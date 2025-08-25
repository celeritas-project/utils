//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/celer-geant.cc
//! \brief Celeritas-Geant4 offloading application
//---------------------------------------------------------------------------//
#include <iostream>
#include <memory>
#include <G4Electron.hh>
#include <G4Positron.hh>
#include <G4RunManagerFactory.hh>
#include <G4Threading.hh>
#include <G4UImanager.hh>
#include <accel/TrackingManagerConstructor.hh>
#include <accel/TrackingManagerIntegration.hh>
#include <celeritas/ext/EmPhysicsList.hh>
#include <corecel/io/Logger.hh>

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "JsonReader.hh"
#include "MakeCelerOptions.hh"
#include "RootIO.hh"

//---------------------------------------------------------------------------//
/*!
 * Validate minimal set of JSON input keys.
 */
void validate_input()
{
#define VALIDATE_MAIN_KEY(KEY) \
    CELER_VALIDATE(j.contains(#KEY), << "\"" << #KEY << "\" key missing");

#define VALIDATE_HIST_KEY(HIST)        \
    CELER_VALIDATE(jh.contains(#HIST), \
                   << "Histogram \"" << #HIST << "\" missing");

#define VALIDATE_HIST_DEF_KEY(HIST, KEY)                                    \
    CELER_VALIDATE(jh.at(#HIST).contains(#KEY),                             \
                   << "Missing \"" << #KEY << "\" in histogram \"" << #HIST \
                   << "\"");

#define VALIDATE_HIST_DEF(HIST)            \
    VALIDATE_HIST_KEY(HIST)                \
    VALIDATE_HIST_DEF_KEY(HIST, num_bins); \
    VALIDATE_HIST_DEF_KEY(HIST, min);      \
    VALIDATE_HIST_DEF_KEY(HIST, max);

    // Validate all main JSON input terms
    auto const& j = JsonReader::Instance();
    VALIDATE_MAIN_KEY(geometry);
    VALIDATE_MAIN_KEY(root_output);
    VALIDATE_MAIN_KEY(num_threads);
    VALIDATE_MAIN_KEY(histograms);

    // Validate histogram information
    auto const& jh = j.at("histograms");
    VALIDATE_HIST_DEF(energy);
    VALIDATE_HIST_DEF(time);

#undef VALIDATE_MAIN_KEY
#undef VALIDATE_HIST_KEY
#undef VALIDATE_HIST_DEF_KEY
#undef VALIDATE_HIST_DEF
}

//---------------------------------------------------------------------------//
/*!
 * Run a Celeritas-Geant4 execution run for physics validation.
 *
 * See README for details.
 */
int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        // Print help message
        std::cout << "Usage: " << argv[0] << " input.json" << std::endl;
        return EXIT_FAILURE;
    }

    // Load and verify input file
    JsonReader::Construct(argv[1]);
    validate_input();

    auto const& json = JsonReader::Instance();
    auto const num_threads = json.at("num_threads").get<size_t>();
    CELER_VALIDATE(num_threads > 0, << "Number of threads must be positive");

    std::unique_ptr<G4RunManager> run_manager;
    run_manager.reset(
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default));

    run_manager->SetNumberOfThreads(num_threads);

    // Initialize Celeritas
    auto& tmi = celeritas::TrackingManagerIntegration::Instance();
    tmi.SetOptions(MakeCelerOptions());

    // Initialize physics with Celeritas offload
    using PhysicsOptions = celeritas::GeantPhysicsOptions;
    auto phys_opts = PhysicsOptions{PhysicsOptions::deactivated()};
    phys_opts.muon = celeritas::GeantMuonPhysicsOptions{};

    auto physics = std::make_unique<celeritas::EmPhysicsList>(phys_opts);
    physics->RegisterPhysics(new celeritas::TrackingManagerConstructor(&tmi));
    run_manager->SetUserInitialization(physics.release());

    // Initialize geometry and actions
    run_manager->SetUserInitialization(
        new DetectorConstruction(json.at("geometry").get<std::string>()));
    run_manager->SetUserInitialization(new ActionInitialization());

    // Run events
    auto const num_events
        = json.at("particle_gun").at("num_events").get<size_t>();
    CELER_VALIDATE(num_events, << "Number of events must be positive");

    run_manager->Initialize();
    run_manager->BeamOn(num_events);

    return EXIT_SUCCESS;
}
