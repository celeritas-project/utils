//----------------------------------*-C++-*----------------------------------//
// Copyright 2022-2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file BoxDetector.hh
//! \brief Create the detector geometry.
//---------------------------------------------------------------------------//
#pragma once

#include <string>
#include <G4VUserDetectorConstruction.hh>
#include <G4Material.hh>
#include <G4VPhysicalVolume.hh>

//---------------------------------------------------------------------------//
/*!
 * Construct a programmatic detector geometry.
 */
class BoxDetector : public G4VUserDetectorConstruction
{
  public:
    // Construct
    BoxDetector();

    // Construct geometry
    G4VPhysicalVolume* Construct() final;
    // Set up sensitive detectors and magnetic field
    void ConstructSDandField() final;

  private:
    // "Infinite" box: Pb cube with 500 m side
    G4VPhysicalVolume* create_box();
    // Set volume as sensitive detector
    void set_sd();
};
