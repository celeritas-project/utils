//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/RootIO.hh
//---------------------------------------------------------------------------//
#pragma once

#include <G4ThreadLocalSingleton.hh>
#include <TFile.h>

#include "HistogramStore.hh"

//---------------------------------------------------------------------------//
/*!
 * Thread-local ROOT I/O manager singleton.
 */
class RootIO
{
    friend class G4ThreadLocalSingleton<RootIO>;

  public:
    // Return a thread-local singleton instance
    static RootIO* Instance();

    // Get reference to thread-local TFile
    TFile& File() { return *file_; }

    // Get reference to thread-local Histogram data
    HistogramStore& Histograms() { return hist_store_; }

  private:
    //// DATA ////

    TFile* file_;
    HistogramStore hist_store_;

    //// HELPER FUNCTIONS ////

    // Construct with filename on worker thread
    RootIO();

    // Destruct by writing histograms/TTree and file to disk
    ~RootIO();
};
