//----------------------------------*-C++-*----------------------------------//
// Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file main.cc
//---------------------------------------------------------------------------//
#include <G4EmStandardPhysics.hh>
#include <G4GDMLParser.hh>
#include <G4MaterialTable.hh>
#include <G4OpticalParameters.hh>
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
//! Test optical properties export.
void export_properties()
{
    auto const& g4material_table = G4Material::GetMaterialTable();
    auto const num_materials = g4material_table->size();
    std::cout << "materials size: " << num_materials << std::endl;

    for (auto i = 0; i < num_materials; i++)
    {
        auto const* g4material = g4material_table->at(i);

        auto const* mat_prop_table = g4material->GetMaterialPropertiesTable();

        auto mat_properties = mat_prop_table->GetProperties();
        auto mat_prop_names = mat_prop_table->GetMaterialPropertyNames();
        std::cout << "property names size: " << mat_prop_names.size()
                  << std::endl;

        for (auto const& name : mat_prop_names)
        {
            std::cout << "property: " << name << std::endl;
        }

        for (auto const& prop : mat_properties)
        {
            // tbd
        }
    }
}

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
    run_manager.RunInitialization();
    // run_manager.BeamOn(1);

    G4OpticalParameters::Instance()->Dump();

    // export_properties();

    bool const vis_interface = false;
    if (vis_interface)
    {
        // Open visualization
        auto qt_interface = new G4UIExecutive(argc, argv);
        auto vis_manager = new G4VisExecutive();
        vis_manager->Initialize();

        auto ui_manager = G4UImanager::GetUIpointer();
        ui_manager->ApplyCommand("/control/execute vis.mac");
        qt_interface->SessionStart();
    }

    if (false)
    {
        // Generate gdml for testing
        export_gdml();
    }

    return EXIT_SUCCESS;
}
