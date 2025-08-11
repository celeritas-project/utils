//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/RootIO.cc
//---------------------------------------------------------------------------//
#include "RootIO.hh"

#include <G4Threading.hh>
#include <TROOT.h>
#include <corecel/Assert.hh>
#include <corecel/io/Logger.hh>

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

    std::string filename = "test";  // \todo: load from input
    CELER_VALIDATE(!filename.empty(), << "ROOT filename must be non-empty");
    std::string thread_filename
        = filename + "-" + std::to_string(G4Threading::G4GetThreadId())
          + ".root";

    CELER_LOG_LOCAL(status) << "Open file " << thread_filename;
    file_ = TFile::Open(thread_filename.c_str(), "recreate");
    tree_ = new TTree("steps", "steps", this->SplitLevel(), this->file_);

    hist_.Initialize();
}

//---------------------------------------------------------------------------//
/*!
 * Destruct by writing histograms and TTree data to thread-local ROOT file, and
 * write/close.
 */
RootIO::~RootIO()
{
    {
        static std::mutex mutex_log;
        std::lock_guard<std::mutex> scoped_lock{mutex_log};
        CELER_LOG_LOCAL(status)
            << "Writing '" << file_->GetName() << "' file to disk";
    }
    tree_->Write();
    auto* hist_dir = file_->mkdir("histograms");
    hist_dir->cd();
    hist_.energy.Write();
    file_->Close();
}
