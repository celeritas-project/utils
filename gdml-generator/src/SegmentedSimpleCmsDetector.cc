//----------------------------------*-C++-*----------------------------------//
// Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file SegmentedSimpleCmsDetector.cc
//---------------------------------------------------------------------------//
#include "SegmentedSimpleCmsDetector.hh"

#include <G4Box.hh>
#include <G4Colour.hh>
#include <G4LogicalVolume.hh>
#include <G4NistManager.hh>
#include <G4PVDivision.hh>
#include <G4PVPlacement.hh>
#include <G4PVReplica.hh>
#include <G4SDManager.hh>
#include <G4VisAttributes.hh>
#include <stdlib.h>

#include "core/SensitiveDetector.hh"

//---------------------------------------------------------------------------//
/*!
 * Construct with geometry type enum and number of segments.
 */
SegmentedSimpleCmsDetector::SegmentedSimpleCmsDetector(
    MaterialType type, SegmentDefinition segments)
    : geometry_type_(type), num_segments_(segments)
{
    if (num_segments_.num_r < 1 || num_segments_.num_theta < 1
        || num_segments_.num_z < 1)
    {
        std::cout << "Number of segments must be at least 1 in every axis"
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    materials_ = this->build_materials();
}

//---------------------------------------------------------------------------//
/*!
 * Mandatory Construct function.
 */
G4VPhysicalVolume* SegmentedSimpleCmsDetector::Construct()
{
    return this->segmented_simple_cms();
}

//---------------------------------------------------------------------------//
/*!
 * Set sensitive detectors and (TODO) magnetic field.
 */
void SegmentedSimpleCmsDetector::ConstructSDandField()
{
    this->set_sd();

    // TODO: Add magnetic field
}

//---------------------------------------------------------------------------//
// PRIVATE
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * Define list of materials.
 */
SegmentedSimpleCmsDetector::MaterialList
SegmentedSimpleCmsDetector::build_materials()
{
    MaterialList materials;
    G4NistManager* nist = G4NistManager::Instance();

    switch (geometry_type_)
    {
        case MaterialType::simple:
            // Load materials
            materials.world = nist->FindOrBuildMaterial("G4_Galactic");
            materials.vacuum_tube = nist->FindOrBuildMaterial("G4_Galactic");
            materials.si_tracker = nist->FindOrBuildMaterial("G4_Si");
            materials.em_calorimeter = nist->FindOrBuildMaterial("G4_Pb");
            materials.had_calorimeter = nist->FindOrBuildMaterial("G4_C");
            materials.sc_solenoid = nist->FindOrBuildMaterial("G4_Ti");
            materials.muon_chambers = nist->FindOrBuildMaterial("G4_Fe");

            // Update names
            materials.world->SetName("vacuum");
            materials.vacuum_tube->SetName("vacuum");
            materials.si_tracker->SetName("Si");
            materials.em_calorimeter->SetName("Pb");
            materials.had_calorimeter->SetName("C");
            materials.sc_solenoid->SetName("Ti");
            materials.muon_chambers->SetName("Fe");
            break;

        case MaterialType::composite:
            // Load materials
            materials.world = nist->FindOrBuildMaterial("G4_Galactic");
            materials.vacuum_tube = nist->FindOrBuildMaterial("G4_Galactic");
            materials.si_tracker
                = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
            materials.em_calorimeter
                = nist->FindOrBuildMaterial("G4_LEAD_OXIDE");
            materials.had_calorimeter = nist->FindOrBuildMaterial("G4_C");
            materials.sc_solenoid = nist->FindOrBuildMaterial("G4_Ti");
            materials.muon_chambers = nist->FindOrBuildMaterial("G4_Fe");

            // Update names
            materials.world->SetName("vacuum");
            materials.vacuum_tube->SetName("vacuum");
            materials.si_tracker->SetName("SiO2");
            materials.em_calorimeter->SetName("Pb3O4");
            materials.had_calorimeter->SetName("C");
            materials.sc_solenoid->SetName("Ti");
            materials.muon_chambers->SetName("Fe");
            break;
        default:
            break;
    }

    return materials;
}

//---------------------------------------------------------------------------//
/*!
 * Programmatic geometry definition: Segmented simple CMS mock up.
 *
 * This is a set of single-element concentric cylinders that can be split
 * multiple times in each direction, to generate a large numbers of volumes.
 *
 * - The World volume is a box, and its dimensions are expressed in cartesian
 * coordinates [x. y, z].
 * - **All** other volumes are concentric cylinders, and their dimensions are
 * expressed as [inner radius, outer radius, length]
 *
 * | Volume                       | Composition      | Dimensions [cm]    |
 * | ---------------------------- | ---------------- | ------------------ |
 * | world                        | H                | [1000, 1000, 2000] |
 * | vacuum tube                  | H                | [0, 30, 1400]      |
 * | silicon tracker              | Si or SiO2       | [30, 125, 1400]    |
 * | electromagnetic calorimeter  | Pb or Pb3O4      | [125, 175, 1400]   |
 * | hadron calorimeter           | C                | [175, 275, 1400]   |
 * | superconducting solenoid     | Ti               | [275, 375, 1400]   |
 * | muon chambers                | Fe               | [375, 700, 1400]   |
 *
 * Table shows the full cylinder values. Each of those cylinders can be
 * segmented in every direction. The sum of all segments will add up to the
 * above values.
 *
 * E.g.: If the silicon tracker is segmented in 2 in the radial axis, it will
 * become si_tracker_1 (r = [30, 77.5]) and si_tracker_2 (r = [77.5, 125]).
 */
G4VPhysicalVolume* SegmentedSimpleCmsDetector::segmented_simple_cms()
{
    // Size of World volume
    double const world_size = 20 * m;

    struct VolumeGap
    {
        double overlap{0};
        double millimeter{1 * mm};
        double tolerance{1e-9 * mm};
    } const volume_gaps;

    struct CylinderRadius
    {
        double vacuum_tube{30 * cm};
        double si_tracker{125 * cm};
        double em_calo{175 * cm};
        double had_calo{275 * cm};
        double sc_solenoid{375 * cm};
        double muon_chambers{700 * cm};
    } const radius;

    // World volume
    G4Box* world_def
        = new G4Box("world_def", world_size / 2, world_size / 2, world_size);

    auto const world_lv
        = new G4LogicalVolume(world_def, materials_.world, "world");

    auto const world_pv = new G4PVPlacement(0,  // Rotation matrix
                                            G4ThreeVector(),  // Position
                                            world_lv,  // Current LV
                                            "world_pv",  // Name
                                            nullptr,  // Mother LV
                                            false,  // Bool operation
                                            0,  // Copy number
                                            false);  // Overlap check

    //// Set up mother tubes needed for the material segments

    // Vacuum tube
    G4Tubs* vacuum_tube_def = new G4Tubs("vacuum_tube_def",
                                         0,  // Inner radius
                                         radius.vacuum_tube,  // Outer radius
                                         half_length_,  // Half-length z
                                         0 * deg,  // Start angle
                                         360 * deg);  // Spanning angle

    auto const vacuum_tube_lv = new G4LogicalVolume(
        vacuum_tube_def, materials_.vacuum_tube, "vacuum_tube_lv");

    new G4PVPlacement(0,
                      G4ThreeVector(),  // Spans -z/2 to +z/2
                      vacuum_tube_lv,
                      "vacuum_tube_pv",
                      world_lv,
                      false,
                      0,
                      false);

    // Si tracker
    G4Tubs* si_tracker_def = new G4Tubs("si_tracker_def",
                                        radius.vacuum_tube,
                                        radius.si_tracker,
                                        half_length_,
                                        0 * deg,
                                        360 * deg);

    auto const si_tracker_lv
        = new G4LogicalVolume(si_tracker_def, materials_.world, "si_tracker");

    auto si_tracker_pv = new G4PVPlacement(0,
                                           G4ThreeVector(),
                                           si_tracker_lv,
                                           "si_tracker_pv",
                                           world_lv,
                                           false,
                                           0,
                                           false);
    // EM calorimeter
    G4Tubs* em_calorimeter_def = new G4Tubs("em_calorimeter_def",
                                            radius.si_tracker,
                                            radius.em_calo,
                                            half_length_,
                                            0 * deg,
                                            360 * deg);

    auto const em_calorimeter_lv = new G4LogicalVolume(
        em_calorimeter_def, materials_.world, "em_calorimeter_lv");

    new G4PVPlacement(0,
                      G4ThreeVector(),
                      em_calorimeter_lv,
                      "em_calorimeter_pv",
                      world_lv,
                      false,
                      0,
                      false);

    // Hadron calorimeter
    G4Tubs* had_calorimeter_def = new G4Tubs("had_calorimeter_def",
                                             radius.em_calo,
                                             radius.had_calo,
                                             half_length_,
                                             0 * deg,
                                             360 * deg);

    auto const had_calorimeter_lv = new G4LogicalVolume(
        had_calorimeter_def, materials_.world, "had_calorimeter_lv");

    new G4PVPlacement(0,
                      G4ThreeVector(),
                      had_calorimeter_lv,
                      "had_calorimeter_pv",
                      world_lv,
                      false,
                      0,
                      false);

    // Superconducting solenoid
    G4Tubs* sc_solenoid_def = new G4Tubs("sc_solenoid_def",
                                         radius.had_calo,
                                         radius.sc_solenoid,
                                         half_length_,
                                         0 * deg,
                                         360 * deg);

    auto const sc_solenoid_lv = new G4LogicalVolume(
        sc_solenoid_def, materials_.world, "sc_solenoid_lv");

    new G4PVPlacement(0,
                      G4ThreeVector(),
                      sc_solenoid_lv,
                      "sc_solenoid_pv",
                      world_lv,
                      false,
                      0,
                      false);

    // Muon chambers
    G4Tubs* muon_chambers_def = new G4Tubs("muon_chambers_def",
                                           radius.sc_solenoid,
                                           radius.muon_chambers,
                                           half_length_,
                                           0 * deg,
                                           360 * deg);

    auto const muon_chambers_lv = new G4LogicalVolume(
        muon_chambers_def, materials_.world, "muon_chambers_lv");

    new G4PVPlacement(0,
                      G4ThreeVector(),
                      muon_chambers_lv,
                      "muon_chambers_pv",
                      world_lv,
                      false,
                      0,
                      false);

    // Add segmented elements
    this->create_segments("si_tracker",
                          radius.vacuum_tube,
                          radius.si_tracker,
                          si_tracker_def,
                          si_tracker_lv,
                          materials_.si_tracker);

    this->create_segments("em_calorimeter",
                          radius.si_tracker,
                          radius.em_calo,
                          em_calorimeter_def,
                          em_calorimeter_lv,
                          materials_.em_calorimeter);

    this->create_segments("had_calorimeter",
                          radius.em_calo,
                          radius.had_calo,
                          had_calorimeter_def,
                          had_calorimeter_lv,
                          materials_.had_calorimeter);

    this->create_segments("sc_solenoid",
                          radius.had_calo,
                          radius.sc_solenoid,
                          sc_solenoid_def,
                          sc_solenoid_lv,
                          materials_.sc_solenoid);

    this->create_segments("muon_chambers",
                          radius.sc_solenoid,
                          radius.muon_chambers,
                          muon_chambers_def,
                          muon_chambers_lv,
                          materials_.muon_chambers);

    return world_pv;
}

//---------------------------------------------------------------------------//
/*!
 * TODO
 */
void SegmentedSimpleCmsDetector::set_sd()
{
    // TODO
}

//---------------------------------------------------------------------------//
/*!
 * Generate segments in r, theta, and z for a given cylinder.
 */
void SegmentedSimpleCmsDetector::create_segments(
    std::string name,
    double inner_r,
    double outer_r,
    G4Tubs* full_culinder_def,
    G4LogicalVolume* full_cylinder_lv,
    G4Material* cyl_material)
{
    std::string name_segment = name + "_segment";
    std::string name_theta = name_segment + "_theta";
    std::string name_r = name_segment + "_r";
    std::string name_z = name_segment + "_z";

    // Theta
    double const segment_theta = 2 * CLHEP::pi / num_segments_.num_theta;

    std::string name_theta_def = name_theta + "_def";
    G4Tubs* segment_theta_def = new G4Tubs(
        name_theta_def, inner_r, outer_r, half_length_, 0, segment_theta);

    std::string name_theta_lv = name_theta + "_lv";
    auto const segment_theta_lv = new G4LogicalVolume(
        full_culinder_def, materials_.world, name_theta_lv);

    // R
    double const segment_r = (outer_r - inner_r) / num_segments_.num_r;

    std::string name_r_def = name_r + "_def";
    G4Tubs* segment_r_def = new G4Tubs(
        name_r_def, inner_r, inner_r + segment_r, half_length_, 0, segment_theta);

    std::string name_r_lv = name_r + "_lv";
    auto const segment_r_lv
        = new G4LogicalVolume(segment_r_def, materials_.world, name_r_lv);

    std::string name_r_pv = name_r + "_pv";
    new G4PVReplica(name_r_pv,
                    segment_r_lv,
                    segment_theta_lv,
                    EAxis::kRho,
                    num_segments_.num_r,
                    0);

    std::string name_theta_pv = name_theta + "_pv";
    new G4PVReplica(name_theta_pv,
                    segment_theta_lv,
                    full_cylinder_lv,
                    EAxis::kPhi,
                    num_segments_.num_theta,
                    segment_theta);

    // Z
    double const segment_z = 2 * half_length_ / num_segments_.num_z;

    std::string name_z_def = name_z + "_def";
    G4Tubs* si_tracker_segment_z_def = new G4Tubs(name_z_def,
                                                  inner_r,
                                                  inner_r + segment_r,
                                                  segment_z / 2,
                                                  0,
                                                  segment_theta);

    std::string name_z_lv = name_z + "_lv";
    auto const si_tracker_segment_z_lv = new G4LogicalVolume(
        si_tracker_segment_z_def, cyl_material, name_z_lv);

    std::string name_z_pv = name_z + "_pv";
    new G4PVDivision(name_z_pv,
                     si_tracker_segment_z_lv,
                     segment_r_lv,
                     EAxis::kZAxis,
                     num_segments_.num_z,
                     0);
}
