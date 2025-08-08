//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file src/MCHist.hh
//---------------------------------------------------------------------------//
#pragma once

#include <TH1D.h>

//---------------------------------------------------------------------------//
struct MCHist
{
    TH1D energy;
    TH1D time;
};
