//----------------------------------*-C++-*----------------------------------//
// Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file main.cc
//---------------------------------------------------------------------------//
#include <G4EmStandardPhysics.hh>
#include <G4RunManager.hh>
#include <G4VModularPhysicsList.hh>

#include "src/ActionInitialization.hh"
#include "src/OpticsDetector.hh"
#include "src/PrimaryGenerator.hh"

//---------------------------------------------------------------------------//
/*!
 * Geant4 optical physics testbed application.
 * See README for details.
 */
int main(int argc, char* argv[])
{
    if (argc != 1)
    {
        // Print help message
        std::cout << "Usage: " << argv[0] << std::endl;
        return EXIT_FAILURE;
    }

    // Construct run manager
    G4RunManager run_manager;
    run_manager.SetVerboseLevel(0);

    // TODO: Set up correct physics list
    std::unique_ptr<G4VModularPhysicsList> physics
        = std::make_unique<G4VModularPhysicsList>();
    physics->RegisterPhysics(new G4EmStandardPhysics(/* verbosity = */ 0));

    // Initialize physics, geometry, and actions
    run_manager.SetUserInitialization(physics.release());
    run_manager.SetUserInitialization(new OpticsDetector());
    run_manager.SetUserInitialization(new ActionInitalization());

    // TODO: maybe set up visualization?

    // Initialize and run
    run_manager.Initialize();
    run_manager.BeamOn(1);

    return EXIT_SUCCESS;
}
