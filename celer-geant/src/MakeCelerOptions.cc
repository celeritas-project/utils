//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/MakeCelerOptions.cc
//---------------------------------------------------------------------------//
#include "MakeCelerOptions.hh"

#include <G4Electron.hh>
#include <G4Gamma.hh>
#include <G4MuonMinus.hh>
#include <G4MuonPlus.hh>
#include <G4Neutron.hh>
#include <G4Positron.hh>
#include <accel/AlongStepFactory.hh>
#include <accel/SetupOptions.hh>
#include <accel/TrackingManagerConstructor.hh>

#include "EventAction.hh"

//---------------------------------------------------------------------------//
/*!
 * Celeritas runtime options.
 */
celeritas::SetupOptions MakeCelerOptions()
{
    celeritas::SetupOptions opts;
    opts.max_num_tracks = 1024 * 16;
    opts.initializer_capacity = opts.max_num_tracks * 8;
    opts.ignore_processes = {"CoulombScat"};
    opts.offload_particles
        = {G4MuonMinus::Definition(), G4MuonPlus::Definition()};
    opts.sd.ignore_zero_deposition = false;

    // Set along-step factory with zero field
    opts.make_along_step = celeritas::UniformAlongStepFactory();

    return opts;
}
