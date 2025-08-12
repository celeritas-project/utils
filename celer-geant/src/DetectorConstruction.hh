//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/DetectorConstruction.hh
//---------------------------------------------------------------------------//
#pragma once

#include <string>
#include <G4GDMLParser.hh>
#include <G4VUserDetectorConstruction.hh>

//---------------------------------------------------------------------------//
/*!
 * Construct detector geometry.
 */
class DetectorConstruction final : public G4VUserDetectorConstruction
{
  public:
    // Construct with GDML path
    DetectorConstruction(std::string gdml_path);

    // Load GDML geometry
    G4VPhysicalVolume* Construct() final;

    // Sensitive detectors are the only Celeritas interface with Geant4
    void ConstructSDandField() final;

  private:
    //// DATA ////
    G4GDMLParser parser_;

    //// HELPER FUNCTIONS ////
    void set_sd();
};
