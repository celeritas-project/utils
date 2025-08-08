//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/celer-geant.cc
//! \brief Celeritas-Geant4 offloading application
//---------------------------------------------------------------------------//
#include <iostream>
#include <memory>
#include <FTFP_BERT.hh>
#include <G4RunManagerFactory.hh>
#include <G4Threading.hh>
#include <G4UImanager.hh>
#include <accel/TrackingManagerConstructor.hh>
#include <accel/TrackingManagerIntegration.hh>

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "MakeCelerOptions.hh"

//---------------------------------------------------------------------------//
/*!
 * See README for details.
 */
int main(int argc, char* argv[])
{
    if (argc > 2)
    {
        // Print help message
        std::cout << "Usage: " << argv[0] << " input.json" << std::endl;
        return EXIT_FAILURE;
    }

    std::string gdml_input
        = "/Users/4s2/celeritas-project/celeritas/app/data/simple-cms.gdml";

    std::unique_ptr<G4RunManager> run_manager;
    run_manager.reset(
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default));
    run_manager->SetNumberOfThreads(5);

    // Initialize Celeritas
    auto& tmi = celeritas::TrackingManagerIntegration::Instance();

    // Initialize physics with celeritas offload
    auto* physics_list = new FTFP_BERT{/* verbosity = */ 0};  // todo: FIXME
    physics_list->RegisterPhysics(
        new celeritas::TrackingManagerConstructor(&tmi));
    run_manager->SetUserInitialization(physics_list);
    tmi.SetOptions(celeritas::MakeCelerOptions());

    // Initialize geometry and actions
    run_manager->SetUserInitialization(new DetectorConstruction(gdml_input));
    run_manager->SetUserInitialization(new ActionInitialization());

    // Run four events
    run_manager->Initialize();
    run_manager->BeamOn(20);

    return EXIT_SUCCESS;
}
