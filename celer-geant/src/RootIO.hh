//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/RootIO.hh
//---------------------------------------------------------------------------//
#pragma once

#include <memory>
#include <mutex>
#include <string>
#include <G4ThreadLocalSingleton.hh>
#include <TFile.h>
#include <TTree.h>
#include <celeritas/ext/RootFileManager.hh>

#include "HistogramStore.hh"

//---------------------------------------------------------------------------//
/*!
 * Thread-local ROOT I/O manager singleton.
 */
class RootIO
{
    friend class G4ThreadLocalSingleton<RootIO>;

  public:
    //!@{
    //! \name type aliases
    using UPRootFile = std::unique_ptr<TFile>;
    using UPRootTree = std::unique_ptr<TTree>;
    //!@}

    // Return a thread-local singleton instance
    static RootIO* Instance();

    // Get reference to thread-local TFile
    TFile& File() { return *file_; }

    // Get reference to thread-local TTree
    // TTree& Tree() { return *tree_; }

    // Get reference to thread-local Histogram data
    HistogramStore& Histograms() { return hist_store_; }

  private:
    //// DATA ////

    TFile* file_;
    // TTree* tree_;

    // Histogram data
    HistogramStore hist_store_;

  private:
    //// HELPER FUNCTIONS ////

    // Construct with filename on worker thread
    RootIO();

    // Destruct by writing histograms/TTree and file to disk
    ~RootIO();

    //! ROOT TTree split level
    static constexpr short int SplitLevel() { return 99; }
};
