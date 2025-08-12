//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/celer-geant.cc
//! \brief Celeritas-Geant4 offloading application
//---------------------------------------------------------------------------//
#include <iostream>
#include <memory>
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

//---------------------------------------------------------------------------//
/*!
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

    // Load input file
    JsonReader::Construct(argv[1]);
    auto const& json = JsonReader::Instance();
    auto const num_threads = json.at("num_threads").get<size_t>();

    CELER_VALIDATE(num_threads > 0, << "Number of threads must be positive");

    auto const num_cores = G4Threading::G4GetNumberOfCores();
    if (2 * num_cores < num_threads)
    {
        CELER_LOG(warning)
            << "Number of threads (" << num_threads
            << ") is larger than the number of available cores (" << num_cores
            << "), assuming hyperthreading";
    }

    std::unique_ptr<G4RunManager> run_manager;
    run_manager.reset(
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default));

    run_manager->SetNumberOfThreads(num_threads);

    // Initialize Celeritas
    auto& tmi = celeritas::TrackingManagerIntegration::Instance();
    tmi.SetOptions(celeritas::MakeCelerOptions());

    // Initialize physics with Celeritas offload
    celeritas::GeantPhysicsOptions phys_opts;
    auto physics = std::make_unique<celeritas::EmPhysicsList>(phys_opts);
    physics->RegisterPhysics(new celeritas::TrackingManagerConstructor(&tmi));
    run_manager->SetUserInitialization(physics.release());

    // Initialize geometry and actions
    run_manager->SetUserInitialization(
        new DetectorConstruction(json.at("geometry").get<std::string>()));
    run_manager->SetUserInitialization(new ActionInitialization());

    // Run events
    run_manager->Initialize();
    run_manager->BeamOn(20);

    return EXIT_SUCCESS;
}
