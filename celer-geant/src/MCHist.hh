//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file src/MCHist.hh
//---------------------------------------------------------------------------//
#pragma once

#include <TH1D.h>

//---------------------------------------------------------------------------//
/*!
 * Store and initialize histogram for faster MC truth output results.
 */
struct MCHist
{
    TH1D energy;
    TH1D time;

    // Initialize histograms with input data
    // \todo: fetch input parameters
    void Initialize()
    {
        this->energy = TH1D("energy", "energy", 100, 0, 100);
        this->time = TH1D("time", "time", 100, 0, 100e-8);
    }
};
