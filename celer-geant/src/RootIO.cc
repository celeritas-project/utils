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
#include <TTree.h>
#include <corecel/Assert.hh>
#include <corecel/io/Logger.hh>

#include "JsonReader.hh"

//---------------------------------------------------------------------------//
/*!
 * Return the static thread local singleton instance.
 */
RootIO* RootIO::Instance()
{
    static thread_local RootIO instance;
    return &instance;
}

//---------------------------------------------------------------------------//
/*!
 * Construct thread-local ROOT I/O.
 */
RootIO::RootIO()
{
    CELER_VALIDATE(G4Threading::IsWorkerThread(),
                   << "Must be constructed on worker thread");

    ROOT::EnableThreadSafety();
    auto json = JsonReader::Instance();
    auto filename = json.at("root_output").get<std::string>();
    CELER_VALIDATE(!filename.empty(), << "ROOT filename must be non-empty");

    // Append thread ID to filename
    std::string thread_filename
        = filename.substr(0, filename.find_last_of("."));
    thread_filename += "-" + std::to_string(G4Threading::G4GetThreadId())
                       + ".root";

    CELER_LOG_LOCAL(status) << "Open file " << thread_filename;
    file_ = TFile::Open(thread_filename.c_str(), "recreate");
    CELER_VALIDATE(!file_->IsZombie(),
                   << "ROOT file \"" << thread_filename << "\" is zombie");

    // Map physical volumes to be scored
    auto const& physvol_store = *G4PhysicalVolumeStore::GetInstance();
    CELER_LOG_LOCAL(status) << "Mapping sensitive detectors for I/O";
    for (auto const& physvol : physvol_store)
    {
        CELER_ASSERT(physvol);
        auto const* logvol = physvol->GetLogicalVolume();
        CELER_ASSERT(logvol);
        auto const* sd = logvol->GetSensitiveDetector();
        if (!sd)
        {
            // Skip non-sensitive logical volumes
            continue;
        }

        auto const name = sd->GetName();
        auto const pid = physvol->GetInstanceID();
        auto const copy_num = physvol->GetCopyNo();
        std::string sd_name = name + "_" + std::to_string(copy_num);
        hist_store_.InsertSensDet(pid, copy_num, sd_name);

        CELER_LOG_LOCAL(info)
            << "Mapped " << name << " with instance ID " << pid
            << " and copy number " << copy_num << " as sensitive detector";
    }
}

//---------------------------------------------------------------------------//
/*!
 * Write Celeritas Output Registry diagnostics as a string during
 * \c RunAction::EndOfRunAction:: .
 *
 * \note Since this is on a worker thread the diagnostics has no record of the
 * total simulation runtime.
 */
void RootIO::StoreDiagnostics(std::string diagnostics)
{
    char const* name = "diagnostics";
    TTree tree(name, name, this->SplitLevel(), nullptr);
    tree.Branch(name, &diagnostics);
    tree.Fill();
    tree.Write();
}

//---------------------------------------------------------------------------//
/*!
 * Finalize I/O by writing data to thread-local ROOT file and close it.
 */
void RootIO::Finalize()
{
#define RIO_HIST_WRITE(MEMBER) hist.MEMBER.Write();

    for (auto& [ids, hist] : hist_store_.Map())
    {
        std::string dir_name = "histograms/" + hist.sd_name;

        auto hist_sd_dir = file_->mkdir(dir_name.c_str());
        hist_sd_dir->cd();

        RIO_HIST_WRITE(energy);
        RIO_HIST_WRITE(time);
    }
    CELER_LOG_LOCAL(info) << "Wrote Geant4 ROOT output to \""
                          << file_->GetName() << "\"";

    file_->Close();
}
