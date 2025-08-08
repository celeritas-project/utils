//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file src/MCTruth.hh
//---------------------------------------------------------------------------//
#pragma once

#include <array>

//---------------------------------------------------------------------------//
struct Step
{
    // Process?
    std::array<double, 3> dir[2];
    std::array<double, 3> pos[2];
    double global_time[2];
};

//---------------------------------------------------------------------------//
struct Hit
{
    // Volume?
    double energy_deposition;
    double global_time;
};
