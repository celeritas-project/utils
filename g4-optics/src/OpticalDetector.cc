//----------------------------------*-C++-*----------------------------------//
// Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file src/OpticalDetector.hh
//---------------------------------------------------------------------------//
#include "OpticalDetector.hh"

#include <G4Box.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalVolume.hh>
#include <G4MaterialPropertyVector.hh>
#include <G4NistManager.hh>
#include <G4OpticalSurface.hh>
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
    //// Optical properties ////

    // For reflection/refraction: need surface model, type, finish
    auto optical_surface = new G4OpticalSurface("optical_surface");
    optical_surface->SetModel(G4OpticalSurfaceModel::unified);
    optical_surface->SetType(G4SurfaceType::dielectric_metal);
    optical_surface->SetFinish(G4OpticalSurfaceFinish::polished);
    optical_surface->SetSigmaAlpha(CLHEP::pi);

    // Set up material properties for reflectivity
    std::vector<double> energies, values;
    energies.push_back(2 * eV);
    energies.push_back(3 * eV);
    values.push_back(1);
    values.push_back(1);
    auto mat_property = std::make_unique<G4MaterialPropertyVector>(
        energies, values, /* spline = */ false);

    // Set up material properties table
    auto surface_properties = new G4MaterialPropertiesTable();
    surface_properties->AddProperty("REFLECTIVITY",
                                    mat_property.release(),
                                    /* create new key = */ true);
    optical_surface->SetMaterialPropertiesTable(surface_properties);

    //// Add properties to materials ////

    geo_mat_.box->SetMaterialPropertiesTable(surface_properties);

    //// Construct geometry /////

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

    new G4LogicalBorderSurface(
        "box_interface", world_pv, box_pv, optical_surface);

    return world_pv;
}
