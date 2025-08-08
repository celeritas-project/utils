//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/MakeCelerOptions.cc
//---------------------------------------------------------------------------//
#include "MakeCelerOptions.hh"

#include <accel/AlongStepFactory.hh>
#include <accel/SetupOptions.hh>

#include "EventAction.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Build options to set up Celeritas.
 */
celeritas::SetupOptions MakeCelerOptions()
{
    celeritas::SetupOptions opts;
    opts.max_num_tracks = 1024 * 16;
    opts.initializer_capacity = 1024 * 128 * 4;
    opts.ignore_processes = {"CoulombScat"};

    // Set along-step factory with zero field
    opts.make_along_step = celeritas::UniformAlongStepFactory();

    return opts;
}
//---------------------------------------------------------------------------//
}  // namespace celeritas
