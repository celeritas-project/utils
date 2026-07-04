//----------------------------------*-C++-*----------------------------------//
// Copyright 2024-2025 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file NotionalDUNE.cc
//---------------------------------------------------------------------------//
#include "NotionalDUNE.hh"

#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4NistManager.hh>
#include <G4Orb.hh>
#include <G4PVPlacement.hh>
#include <G4SDManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>

#include "core/SensitiveDetector.hh"

//---------------------------------------------------------------------------//
/*!
 * Constructor.
 */
NotionalDUNE::NotionalDUNE(int num_spheres_per_axis, int num_levels)
    : num_spheres_per_axis_(num_spheres_per_axis)
    , num_levels_(num_levels)
{
}

//---------------------------------------------------------------------------//
/*!
 * Build a notional DUNE geometry.
 */
G4VPhysicalVolume* NotionalDUNE::Construct()
{
    // Materials
    auto nist = G4NistManager::Instance();
    // Bulk medium of the TPC: liquid argon
    auto lar_mat = nist->FindOrBuildMaterial("G4_lAr");
    // Anode material. Copper is used as a representative stand-in
    auto anode_mat = nist->FindOrBuildMaterial("G4_Cu");
    // Everything outside the active volume (world + concentric shells): vacuum
    auto world_mat = nist->FindOrBuildMaterial("G4_Galactic");
    world_mat->SetName("vacuum");

    // World: large enough to enclose the outermost concentric box, with one
    // extra layer of vacuum margin so nothing is tangent to the world boundary
    double const outer_edge = box_size_ + 2.0 * num_levels_ * level_thickness_;
    double const world_edge = outer_edge + 2.0 * level_thickness_;
    G4Box* world_box = new G4Box("world_box",
                                 0.5 * world_edge * cm,
                                 0.5 * world_edge * cm,
                                 0.5 * world_edge * cm);
    auto const world_lv = new G4LogicalVolume(world_box, world_mat, "world");
    auto const world_pv = new G4PVPlacement(
        nullptr, G4ThreeVector(), world_lv, "world_pv", nullptr, false, 0, false);

    // Concentric vacuum boxes, built from the outermost (level M) inward so
    // that each shell is placed inside the next-larger one. 
    auto mother_lv = world_lv;
    for (int level = num_levels_; level >= 1; --level)
    {
        double const edge = box_size_ + 2.0 * level * level_thickness_;
        G4Box* shell_box = new G4Box(
            "shell_box", 0.5 * edge * cm, 0.5 * edge * cm, 0.5 * edge * cm);
        auto const shell_lv
            = new G4LogicalVolume(shell_box, world_mat, "shell_lv");
        new G4PVPlacement(nullptr,
                          G4ThreeVector(),
                          shell_lv,
                          "shell_pv",
                          mother_lv,
                          false,
                          level,
                          false);
        mother_lv = shell_lv;
    }

    // Central liquid-argon box that holds the grid of anode spheres.
    // Placed in the innermost shell, or directly in the world when
    // num_levels_ == 0.
    G4Box* inner_box = new G4Box("inner_box",
                                 0.5 * box_size_ * cm,
                                 0.5 * box_size_ * cm,
                                 0.5 * box_size_ * cm);
    auto const inner_lv = new G4LogicalVolume(inner_box, lar_mat, "inner_lv");
    new G4PVPlacement(
        nullptr, G4ThreeVector(), inner_lv, "inner_pv", mother_lv, false, 0, false);

    // Grid of N^3 "anode" spheres, equally spaced and symmetric about the 
    // origin.
    double const pitch = box_size_ / num_spheres_per_axis_;
    double const offset = -0.5 * box_size_ + 0.5 * pitch;

    G4Orb* sphere = new G4Orb("sphere", sphere_radius_ * cm);
    auto const sphere_lv = new G4LogicalVolume(sphere, anode_mat, "sphere_lv");

    int copy_no = 0;
    for (int i = 0; i < num_spheres_per_axis_; ++i)
    {
        for (int j = 0; j < num_spheres_per_axis_; ++j)
        {
            for (int k = 0; k < num_spheres_per_axis_; ++k)
            {
                G4ThreeVector pos((offset + i * pitch) * cm,
                                  (offset + j * pitch) * cm,
                                  (offset + k * pitch) * cm);
                new G4PVPlacement(nullptr,
                                  pos,
                                  sphere_lv,
                                  "sphere_pv",
                                  inner_lv,
                                  false,
                                  copy_no++,
                                  false);
            }
        }
    }

    return world_pv;
}

//---------------------------------------------------------------------------//
/*!
 * Set the bulk liquid-argon as the sensitive detector.
 */
void NotionalDUNE::ConstructSDandField()
{
    auto lar_sd = new SensitiveDetector("lar_sd");
    G4SDManager::GetSDMpointer()->AddNewDetector(lar_sd);
    G4VUserDetectorConstruction::SetSensitiveDetector("inner_lv", lar_sd);
}
