//----------------------------------*-C++-*----------------------------------//
// Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file src/OpticsDetector.hh
//---------------------------------------------------------------------------//
#pragma once

#include <G4Material.hh>
#include <G4SystemOfUnits.hh>
#include <G4VUserDetectorConstruction.hh>

//---------------------------------------------------------------------------//
/*!
 * Construct detector geometry.
 */
class OpticsDetector final : public G4VUserDetectorConstruction
{
  public:
    // Construct and set up materials
    OpticsDetector();

    // Define programmatic geometry
    G4VPhysicalVolume* Construct() final;

  private:
    struct GeometryDef
    {
        double const world{10 * m};
        double const box{2 * m};
    } geo_sizes_;

    struct GeometryMaterial
    {
        G4Material* world;
        G4Material* box;
    } geo_mat_;
};
