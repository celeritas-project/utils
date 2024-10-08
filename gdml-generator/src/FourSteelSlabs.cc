//----------------------------------*-C++-*----------------------------------//
// Copyright 2024 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file FourSteelSlabs.cc
//---------------------------------------------------------------------------//
#include "FourSteelSlabs.hh"

#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4SystemOfUnits.hh>

#include "core/SensitiveDetector.hh"

//---------------------------------------------------------------------------//
/*!
 * Construct empty.
 */
FourSteelSlabs::FourSteelSlabs() {}

//---------------------------------------------------------------------------//
/*!
 * Mandatory Construct function.
 */
G4VPhysicalVolume* FourSteelSlabs::Construct()
{
    return this->create_geometry();
}

//---------------------------------------------------------------------------//
// PRIVATE
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * Programmatic geometry definition: World volume with 4 steel slabs.
 *
 * \note G4Replicas are avoided.
 */
G4VPhysicalVolume* FourSteelSlabs::create_geometry()
{
    // Geometry materials
    auto const nist_manager = G4NistManager::Instance();
    auto const world_material
        = nist_manager->FindOrBuildMaterial("G4_Galactic");
    auto const slab_material
        = nist_manager->FindOrBuildMaterial("G4_STAINLESS-STEEL");

    // World definition
    double const world_size = 1000. * cm;
    double const half_world = world_size * 0.5;

    auto const world_solid
        = new G4Box("world_box", half_world, half_world, half_world);

    auto const world_log_vol
        = new G4LogicalVolume(world_solid, world_material, "world_lv");

    auto const world_phys_vol = new G4PVPlacement(nullptr,
                                                  G4ThreeVector(),
                                                  world_log_vol,
                                                  "world_pv",
                                                  nullptr,
                                                  false,
                                                  0,
                                                  true);

    // Slabs definition
    auto const slabs_xy = 0.01 * world_size;
    auto const slabs_z = 0.2 * slabs_xy;

    // Slab 0
    auto const slab_solid = new G4Box("box", slabs_xy, slabs_xy, slabs_z);

    auto const slab_log_vol
        = new G4LogicalVolume(slab_solid, slab_material, "box");

    new G4PVPlacement(
        0, G4ThreeVector(), slab_log_vol, "box", world_log_vol, false, 0, true);

    // Slab 1
    auto const slab_replica_solid
        = new G4Box("boxReplica", slabs_xy, slabs_xy, slabs_z);

    auto const slab_replica_log_vol
        = new G4LogicalVolume(slab_replica_solid, slab_material, "boxReplica");

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, 3 * slabs_z),
                      slab_replica_log_vol,
                      "box",
                      world_log_vol,
                      false,
                      0,
                      true);

    // Slab 2
    auto const slab_replica_2_solid
        = new G4Box("boxReplica2", slabs_xy, slabs_xy, slabs_z);

    auto const slab_replica_2_log_vol = new G4LogicalVolume(
        slab_replica_2_solid, slab_material, "boxReplica");

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, 6 * slabs_z),
                      slab_replica_2_log_vol,
                      "box",
                      world_log_vol,
                      false,
                      0,
                      true);

    // Slab 3
    auto const slab_replica_3_solid
        = new G4Box("boxReplica3", slabs_xy, slabs_xy, slabs_z);

    auto const slab_replica_3_log_vol = new G4LogicalVolume(
        slab_replica_3_solid, slab_material, "boxReplica");

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, 9 * slabs_z),
                      slab_replica_3_log_vol,
                      "box",
                      world_log_vol,
                      false,
                      0,
                      true);

    return world_phys_vol;
}
