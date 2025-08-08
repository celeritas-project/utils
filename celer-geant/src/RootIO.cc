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
 * Construct thread-local Root IO.
 */
RootIO::RootIO()
{
    ROOT::EnableThreadSafety();

    CELER_VALIDATE(G4Threading::IsWorkerThread(),
                   << "Must be constructed on worker thread");

    std::string filename = "test";
    CELER_VALIDATE(!filename.empty(), << "ROOT filename must be non-empty");
    std::string thread_filename
        = filename + "-" + std::to_string(G4Threading::G4GetThreadId())
          + ".root";
    file_ = TFile::Open(thread_filename.c_str(), "recreate");
    CELER_LOG_LOCAL(status) << "Open file " << thread_filename;
    tree_ = new TTree("steps", "steps", this->SplitLevel(), this->file_);

    // Initialize histograms
    hist_.energy = TH1D("energy", "energy", 100, 0, 100);
}

//---------------------------------------------------------------------------//
void RootIO::Finalize()
{
    CELER_LOG_LOCAL(status)
        << "Writing '" << file_->GetName() << "' file to disk";
    tree_->Write();
    auto* hist_dir = file_->mkdir("histograms");
    hist_dir->cd();
    hist_.energy.Write();
    file_->Close();
}
