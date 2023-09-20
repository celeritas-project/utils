//----------------------------------*-C++-*----------------------------------//
// Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file SegmentedSimpleCmsDetector.hh
//! \brief Create the detector geometry.
//---------------------------------------------------------------------------//
#pragma once

#include <string>
#include <G4Material.hh>
#include <G4SystemOfUnits.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VUserDetectorConstruction.hh>

//---------------------------------------------------------------------------//
/*!
 * Construct a programmatic detector geometry.
 */
class SegmentedSimpleCmsDetector : public G4VUserDetectorConstruction
{
  public:
    enum class MaterialType
    {
        simple,
        composite
    };

    struct SegmentDefinition
    {
        std::size_t num_theta;
        std::size_t num_r;
        std::size_t num_z;
    };

    // Construct with geometry type
    SegmentedSimpleCmsDetector(MaterialType type, SegmentDefinition segments);

    // Construct geometry
    G4VPhysicalVolume* Construct() final;
    // Set up sensitive detectors and magnetic field
    void ConstructSDandField() final;

  private:
    // Material selection
    struct MaterialList
    {
        G4Material* world;
        G4Material* vacuum_tube;
        G4Material* si_tracker;
        G4Material* em_calorimeter;
        G4Material* had_calorimeter;
        G4Material* sc_solenoid;
        G4Material* muon_chambers;
    };

#if 0
    // Volume sizes
    struct VolumeSize
    {
        double half_z_len = 7 * m;
        double vacuum_tube_r = 30 * cm;
        double si_tracker_r = 125 * cm - vacuum_tube_r;
        double em_calorimeter_r = 175 * cm - si_tracker_r;
        double had_calorimeter_r = 275 * cm - em_calorimeter_r;
        double sc_solenoid_r = 375 * cm - em_calorimeter_r;
        double mu_chambers_r = 700 * cm - em_calorimeter_r;
    } const static volume_size_;

    struct VolumeRadiusLimit
    {
        double vacuum_tube[2] = {0, volume_size_.vacuum_tube_r};
        double s_tracker[2]
            = {vacuum_tube[1], vacuum_tube[1] + volume_size_.si_tracker_r};
    } const static volume_r_lim_;
#endif

    // Select simple/composite geometry
    MaterialType geometry_type_;

    // Number of segments
    SegmentDefinition num_segments_;

  private:
    // Segmented simple CMS
    G4VPhysicalVolume* segmented_simple_cms();
    // Set sensitive detectors
    void set_sd();
    // Build material list based on MaterialType
    MaterialList build_materials();
};
