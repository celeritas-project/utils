//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/RootIO.hh
//---------------------------------------------------------------------------//
#pragma once

#include <TFile.h>

#include "HistogramStore.hh"

//---------------------------------------------------------------------------//
/*!
 * Thread-local ROOT I/O manager singleton.
 */
class RootIO
{
  public:
    //! Return a thread-local singleton instance
    static RootIO* Instance();

    //! Get reference to thread-local Histogram data
    HistogramStore& Histograms() { return hist_store_; }

    //! Store OutputRegistry diagnostics
    void StoreDiagnostics(std::string diagnostics);

    //! Write data to ROOT file and close it
    void Finalize();

  private:
    //// DATA ////

    TFile* file_;
    HistogramStore hist_store_;

    //// HELPER FUNCTIONS ////

    // Construct with JSON input filename on worker thread
    RootIO();

    // ROOT TTree split level
    static constexpr short int SplitLevel() { return 99; }
};
