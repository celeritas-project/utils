//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/MakeCelerOptions.hh
//---------------------------------------------------------------------------//
#pragma once

#include <G4Electron.hh>
#include <G4Gamma.hh>
#include <G4MuonMinus.hh>
#include <G4MuonPlus.hh>
#include <G4Neutron.hh>
#include <G4Positron.hh>
#include <accel/AlongStepFactory.hh>
#include <accel/SetupOptions.hh>
#include <accel/TrackingManagerConstructor.hh>
#include <celeritas/phys/PDGNumber.hh>

#include "JsonReader.hh"

//---------------------------------------------------------------------------/
/*!
 * Helper map for loading vector of \c G4ParticleDefinition from PDGs in the
 * JSON input.
 */
celeritas::SetupOptions::VecG4PD from_json()
{
    using celeritas::PDGNumber;

    celeritas::SetupOptions::VecG4PD result;
    auto& json = JsonReader::Instance();
    auto const input = json.at("offload_particles").get<std::vector<int>>();

    static std::unordered_map<PDGNumber, G4ParticleDefinition*> supported = {
        {celeritas::pdg::gamma(), G4Gamma::Definition()},
        {celeritas::pdg::electron(), G4Electron::Definition()},
        {celeritas::pdg::positron(), G4Positron::Definition()},
        {celeritas::pdg::mu_minus(), G4MuonMinus::Definition()},
        {celeritas::pdg::mu_plus(), G4MuonPlus::Definition()},
    };

    for (auto pdg : input)
    {
        auto it = supported.find(PDGNumber{pdg});
        CELER_VALIDATE(it != supported.end(),
                       << "PDG '" << pdg << "' not available");
        result.push_back(it->second);
    }
    return result;
}

//---------------------------------------------------------------------------/
/*!
 * Celeritas runtime options.
 */
celeritas::SetupOptions MakeCelerOptions()
{
    celeritas::SetupOptions opts;
    opts.max_num_tracks = 1024 * 16;
    opts.initializer_capacity = opts.max_num_tracks * 8;
    opts.offload_particles = from_json();
    opts.sd.ignore_zero_deposition = false;

    // Set along-step factory with zero field
    opts.make_along_step = celeritas::UniformAlongStepFactory();

    return opts;
}
