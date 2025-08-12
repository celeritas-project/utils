//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/RootIO.cc
//---------------------------------------------------------------------------//
#include "RootIO.hh"

#include <G4LogicalVolume.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4Threading.hh>
#include <G4VSensitiveDetector.hh>
#include <TROOT.h>
#include <corecel/Assert.hh>
#include <corecel/io/Logger.hh>

#include "JsonReader.hh"

//---------------------------------------------------------------------------//
/*!
 * Return the static thread local singleton instance.
 */
RootIO* RootIO::Instance()
{
    static G4ThreadLocal RootIO instance;
    return &instance;
}

//---------------------------------------------------------------------------//
/*!
 * Construct thread-local ROOT I/O.
 */
RootIO::RootIO()
{
    ROOT::EnableThreadSafety();

    CELER_VALIDATE(G4Threading::IsWorkerThread(),
                   << "Must be constructed on worker thread");

    auto filename = JsonReader::Instance().at("mctruth").get<std::string>();
    CELER_VALIDATE(!filename.empty(), << "ROOT filename must be non-empty");

    // Append thread id to filename
    std::string thread_filename
        = filename.substr(0, filename.find_last_of("."));
    thread_filename += "-" + std::to_string(G4Threading::G4GetThreadId())
                       + ".root";

    CELER_LOG_LOCAL(status) << "Open file " << thread_filename;
    file_ = TFile::Open(thread_filename.c_str(), "recreate");

    // TODO:  tree_ = new TTree("steps", "steps", this->SplitLevel(), file_);

    // Map physical volumes to be scored
    auto const& physvol_store = *G4PhysicalVolumeStore::GetInstance();
    CELER_LOG_LOCAL(status) << "Mapping sensitive detectors for I/O";
    for (auto const& physvol : physvol_store)
    {
        CELER_ASSERT(physvol);
        auto const* logvol = physvol->GetLogicalVolume();
        CELER_ASSERT(logvol);
        auto const* sd = logvol->GetSensitiveDetector();
        if (sd)
        {
            auto name = sd->GetName();
            auto pid = physvol->GetInstanceID();
            auto copy_num = physvol->GetCopyNo();
            std::string sd_name = name + "_" + std::to_string(copy_num);
            hist_store_.InsertSensDet(pid, copy_num, sd_name);

            CELER_LOG_LOCAL(status)
                << "Mapped " << name << " with copy number " << copy_num
                << " as sensitive detector";
        }
    }
}

//---------------------------------------------------------------------------//
/*!
 * Destruct by writing histograms and TTree data to thread-local ROOT file, and
 * write/close file.
 */
RootIO::~RootIO()
{
#define RIO_HIST_WRITE(MEMBER) hist.MEMBER.Write()

    static std::mutex mutex_log;
    std::lock_guard<std::mutex> scoped_lock{mutex_log};
    CELER_LOG_LOCAL(status)
        << "Writing '" << file_->GetName() << "' file to disk";

    // TODO: tree_->Write();

    // Write histograms to disk
    for (auto& [ids, hist] : hist_store_.Map())
    {
        std::string dir_name = "histograms/" + hist.sd_name;
        auto hist_sd_dir = file_->mkdir(dir_name.c_str());
        hist_sd_dir->cd();

        RIO_HIST_WRITE(energy);
        RIO_HIST_WRITE(time);
    }
    file_->Close();

#undef RIO_HIST_WRITE
}
