//----------------------------------*-C++-*----------------------------------//
// Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file src/OpticalDetector.hh
//---------------------------------------------------------------------------//
#include "OpticalDetector.hh"

#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4VPhysicalVolume.hh>

//---------------------------------------------------------------------------//
/*!
 * Construct and set up materials.
 */
OpticalDetector::OpticalDetector() : G4VUserDetectorConstruction()
{
    auto const nist = G4NistManager::Instance();
    geo_mat_.world = nist->FindOrBuildMaterial("G4_Galactic");
    geo_mat_.box = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
}

//---------------------------------------------------------------------------//
/*!
 * Construct geometry: a vacuum world box with a smaller box with a material
 * with optical properties.
 */
G4VPhysicalVolume* OpticalDetector::Construct()
{
    // World
    auto world_def = new G4Box(
        "world_def", geo_sizes_.world, geo_sizes_.world, geo_sizes_.world);

    auto world_lv = new G4LogicalVolume(world_def,  // Solid definition
                                        geo_mat_.world,  // Material definition
                                        "world_lv");  // LV name

    auto world_pv = new G4PVPlacement(0,  // Rotation matrix
                                      G4ThreeVector(),  // Position
                                      world_lv,  // Current LV
                                      "world_pv",  // Name
                                      nullptr,  // Mother LV
                                      false,  // Bool operation
                                      0,  // Copy number
                                      false);  // Overlap check

    // Box with optical properties
    auto box_def
        = new G4Box("box_def", geo_sizes_.box, geo_sizes_.box, geo_sizes_.box);

    auto box_lv = new G4LogicalVolume(box_def, geo_mat_.box, "box_lv");

    auto box_pv = new G4PVPlacement(
        0, G4ThreeVector(), box_lv, "box_pv", world_lv, false, 0, false);

    // TODO: Set up optical properties

    return world_pv;
}
