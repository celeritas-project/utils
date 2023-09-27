//----------------------------------*-C++-*----------------------------------//
// Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file main.cc
//---------------------------------------------------------------------------//
#include <G4EmStandardPhysics.hh>
#include <G4GDMLParser.hh>
#include <G4RunManager.hh>
#include <G4UIExecutive.hh>
#include <G4UImanager.hh>
#include <G4VModularPhysicsList.hh>
#include <G4VisExecutive.hh>

#include "src/ActionInitialization.hh"
#include "src/OpticalDetector.hh"
#include "src/OpticalPhysics.hh"
#include "src/PrimaryGenerator.hh"

//---------------------------------------------------------------------------//
//! HELPER function for exporting the geomtry as a GDML file.
void export_gdml()
{
    G4GDMLParser parser;
    parser.SetEnergyCutsExport(true);
    parser.SetSDExport(true);
    parser.SetOverlapCheck(true);
    parser.SetOutputFileOverwrite(true);
    parser.Write("geo-optics.gdml",
                 G4TransportationManager::GetTransportationManager()
                     ->GetNavigatorForTracking()
                     ->GetWorldVolume()
                     ->GetLogicalVolume(),
                 false);  // bool appends ptr address to name
}

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

    // Initialize physics, geometry, and actions
    run_manager.SetUserInitialization(new OpticalPhysics());
    run_manager.SetUserInitialization(new OpticalDetector());
    run_manager.SetUserInitialization(new ActionInitalization());
    run_manager.Initialize();

    // Visualization
    auto qt_interface = new G4UIExecutive(1, argv);
    auto vis_manager = new G4VisExecutive();
    vis_manager->Initialize();

    auto ui_manager = G4UImanager::GetUIpointer();
    ui_manager->ApplyCommand("/control/execute vis.mac");

    run_manager.BeamOn(1);
    qt_interface->SessionStart();

    if (false)
    {
        // Generate gdml for testing
        export_gdml();
    }

    return EXIT_SUCCESS;
}
