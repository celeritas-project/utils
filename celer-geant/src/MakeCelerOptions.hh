//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/MakeCelerOptions.hh
//---------------------------------------------------------------------------//
#pragma once

#include <unordered_map>
#include <vector>
#include <G4Electron.hh>
#include <G4Gamma.hh>
#include <G4MuonMinus.hh>
#include <G4MuonPlus.hh>
#include <G4Neutron.hh>
#include <G4OpticalPhoton.hh>
#include <G4Positron.hh>
#include <accel/AlongStepFactory.hh>
#include <accel/SetupOptions.hh>
#include <accel/TrackingManagerConstructor.hh>
#include <celeritas/phys/PDGNumber.hh>
#include <corecel/Assert.hh>
#include <corecel/io/Logger.hh>

#include "JsonReader.hh"
#include "RootOpticalIO.hh"

namespace detail
{
//---------------------------------------------------------------------------/
/*!
 * Static list of valid PDGs for Celeritas offload.
 */
using PDG = int;
using VecPDG = std::vector<PDG>;

static VecPDG const supported_celeritas_pdgs{
    celeritas::pdg::gamma().get(),
    celeritas::pdg::electron().get(),
    celeritas::pdg::positron().get(),
    celeritas::pdg::mu_minus().get(),
    celeritas::pdg::mu_plus().get(),
    -22, /* Optical photon */
};

//! Helper function to verify if PDG is in the list of particles
static bool is_valid_celeritas_pdg(VecPDG const& valid_pdgs, PDG pdg)
{
    return std::any_of(supported_celeritas_pdgs.begin(),
                       supported_celeritas_pdgs.end(),
                       [pdg](PDG this_pdg) { return this_pdg == pdg; });
}

//! Return list of PDGs used given the JSON input
static VecPDG offloaded_pdgs_from_json()
{
    JsonReader::Validate(JsonReader::Instance(), "celeritas");
    auto const& json = JsonReader::Instance().at("celeritas");
    VecPDG result = json.contains("offload_particles")
                        ? json.at("offload_particles").get<VecPDG>()
                        : supported_celeritas_pdgs;
    return result;
}
//---------------------------------------------------------------------------/
}  // namespace detail

//---------------------------------------------------------------------------/
/*!
 * Load vector of \c G4ParticleDefinition from list of PDGs.
 */
inline celeritas::SetupOptions::VecG4PD initialize_pdgs_from_json()
{
    using celeritas::PDGNumber;
    static std::unordered_map<PDGNumber, G4ParticleDefinition*> supported = {
        {celeritas::pdg::gamma(), G4Gamma::Definition()},
        {celeritas::pdg::electron(), G4Electron::Definition()},
        {celeritas::pdg::positron(), G4Positron::Definition()},
        {celeritas::pdg::mu_minus(), G4MuonMinus::Definition()},
        {celeritas::pdg::mu_plus(), G4MuonPlus::Definition()},
        {PDGNumber{-22}, G4OpticalPhoton::Definition()},
    };

    celeritas::SetupOptions::VecG4PD result;
    auto const input = detail::offloaded_pdgs_from_json();
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
inline celeritas::SetupOptions MakeCelerOptions()
{
    using PDG = int;
    using VecPDG = std::vector<PDG>;

    JsonReader::Validate(JsonReader::Instance(), "celeritas");
    auto const& json = JsonReader::Instance().at("celeritas");

    celeritas::SetupOptions opts;
    JsonReader::Validate(json, "max_num_tracks");
    opts.max_num_tracks = json.at("max_num_tracks").get<size_t>();

    JsonReader::Validate(json, "initializer_capacity");
    opts.initializer_capacity = json.at("initializer_capacity").get<size_t>();

    if (json.contains("offload_particles"))
    {
        opts.offload_particles = initialize_pdgs_from_json();
    }
    else
    {
        CELER_LOG(info)
            << "Celeritas' \"offload_particles\" option not present. "
               "Using default list.";
    }

    opts.sd.ignore_zero_deposition = false;

    // Set along-step factory with zero field
    opts.make_along_step = celeritas::UniformAlongStepFactory();

    // Optical physics setup
    opts.optical = [] {
        celeritas::OpticalSetupOptions opt;
        opt.capacity.primaries = std::pow(2, 19);
        opt.capacity.tracks = std::pow(2, 17);
        opt.capacity.generators = std::pow(2, 16);
        opt.generator = celeritas::inp::OpticalDirectGenerator{};
        opt.detectors.callback = celeritas::optical_hits_callback;
        return opt;
    }();

    return opts;
}
