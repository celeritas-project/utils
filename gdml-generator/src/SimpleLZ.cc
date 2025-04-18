//----------------------------------*-C++-*----------------------------------//
// Copyright 2024-2025 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file SimpleLZ.cc
//---------------------------------------------------------------------------//
#include "SimpleLZ.hh"

#include <cmath>
#include <G4LogicalVolume.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4SDManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4Tubs.hh>
#include <G4VisAttributes.hh>

#include "core/SensitiveDetector.hh"

namespace
{
// X and Y locations of the LZ PMTs
double pmt_x[253]
    = {0,           92,          46,          -46,         -92,
       -46,         46,          138,         0,           -138,
       -138,        0,           138,         184,         92,
       -92,         -184,        -92,         92,          230,
       184,         46,          -46,         -184,        -230,
       -230,        -184,        -46,         46,          184,
       230,         276,         138,         -138,        -276,
       -138,        138,         276,         0,           -276,
       -276,        0,           276,         322,         230,
       92,          -92,         -230,        -322,        -322,
       -230,        -92,         92,          230,         322,
       368,         184,         -184,        -368,        -184,
       184,         366.828253,  320.988013,  45.84024,    -45.84024,
       -320.988013, -366.828253, -366.828253, -320.988013, -45.84024,
       45.84024,    320.988013,  366.828253,  408.490432,  270.559222,
       137.931209,  -137.931209, -270.559222, -408.490432, -408.490432,
       -270.559222, -137.931209, 137.931209,  270.559222,  408.490432,
       473.546359,  236.830639,  -236.71572,  -473.546359, -236.830639,
       236.71572,   412.663978,  0.069725,    -412.594253, -412.663978,
       -0.069725,   412.594253,  456.929688,  369.274335,  87.655353,
       -87.655353,  -369.274335, -456.929688, -456.929688, -369.274335,
       -87.655353,  87.655353,   369.274335,  456.929688,  494.720323,
       320.643074,  174.077248,  -174.077248, -320.643074, -494.720323,
       -494.720323, -320.643074, -174.077248, 174.077248,  320.643074,
       494.720323,  504.218985,  460.583403,  43.635582,   -43.635582,
       -460.583403, -504.218985, -504.218985, -460.583403, -43.635582,
       43.635582,   460.583403,  504.218985,  552.229751,  405.491978,
       146.737773,  -146.737773, -405.491978, -552.229751, -552.229751,
       -405.491978, -146.737773, 146.737773,  405.491978,  552.229751,
       570.783016,  322.752069,  248.030947,  -248.030947, -322.752069,
       -570.783016, -570.783016, -322.752069, -248.030947, 248.030947,
       322.752069,  570.783016,  590.793211,  507.875591,  82.91762,
       -82.91762,   -507.875591, -590.793211, -590.793211, -507.875591,
       -82.91762,   82.91762,    507.875591,  590.793211,  646.325679,
       397.420482,  248.905197,  -248.905197, -397.420482, -646.325679,
       -646.325679, -397.420482, -248.905197, 248.905197,  397.420482,
       646.325679,  635.947962,  466.425365,  169.522598,  -169.522598,
       -466.425365, -635.947962, -635.947962, -466.425365, -169.522598,
       169.522598,  466.425365,  635.947962,  659.264523,  570.940173,
       329.648599,  0.001489,    -329.615923, -570.938683, -659.264523,
       -570.940173, -329.648599, -0.001489,   329.615923,  570.938683,
       731.430591,  718.915611,  694.099786,  657.407722,  609.467223,
       551.098576,  483.300481,  407.23298,   324.197611,  235.61513,
       143.00121,   47.940499,   -47.940499,  -143.00121,  -235.61513,
       -324.197611, -407.23298,  -483.300481, -551.098576, -609.467223,
       -657.407722, -694.099786, -718.915611, -731.430591, -731.430591,
       -718.915611, -694.099786, -657.407722, -609.467223, -551.098576,
       -483.300481, -407.23298,  -324.197611, -235.61513,  -143.00121,
       -47.940499,  47.940499,   143.00121,   235.61513,   324.197611,
       407.23298,   483.300481,  551.098576,  609.467223,  657.407722,
       694.099786,  718.915611,  731.430591};

double pmt_y[253]
    = {0,           0,           -79.674337,  -79.674337,  0,
       79.674337,   79.674337,   -79.674337,  -159.348674, -79.674337,
       79.674337,   159.348674,  79.674337,   0,           -159.348674,
       -159.348674, 0,           159.348674,  159.348674,  -79.674337,
       -159.348674, -239.023011, -239.023011, -159.348674, -79.674337,
       79.674337,   159.348674,  239.023011,  239.023011,  159.348674,
       79.674337,   0,           -239.023011, -239.023011, 0,
       239.023011,  239.023011,  -159.348674, -318.697349, -159.348674,
       159.348674,  318.697349,  159.348674,  -79.674337,  -239.023011,
       -318.697349, -318.697349, -239.023011, -79.674337,  79.674337,
       239.023011,  318.697349,  318.697349,  239.023011,  79.674337,
       0,           -318.697349, -318.697349, 0,           318.697349,
       318.697349,  -158.856641, -238.254266, -397.110906, -397.110906,
       -238.254266, -158.856641, 158.856641,  238.254266,  397.110906,
       397.110906,  238.254266,  158.856641,  -76.572819,  -315.476682,
       -392.0495,   -392.0495,   -315.476682, -76.572819,  76.572819,
       315.476682,  392.0495,    392.0495,    315.476682,  76.572819,
       0.066348,    -410.070002, -410.136351, -0.066348,   410.070002,
       410.136351,  -238.171147, -476.463062, -238.291915, 238.171147,
       476.463062,  238.291915,  -162.592795, -314.41632,  -477.009115,
       -477.009115, -314.41632,  -162.592795, 162.592795,  314.41632,
       477.009115,  477.009115,  314.41632,   162.592795,  -84.619819,
       -386.130458, -470.750277, -470.750277, -386.130458, -84.619819,
       84.619819,   386.130458,  470.750277,  470.750277,  386.130458,
       84.619819,   -240.724937, -316.303982, -557.028918, -557.028918,
       -316.303982, -240.724937, 240.724937,  316.303982,  557.028918,
       557.028918,  316.303982,  240.724937,  -149.39181,  -403.549089,
       -552.940898, -552.940898, -403.549089, -149.39181,  149.39181,
       403.549089,  552.940898,  552.940898,  403.549089,  149.39181,
       -43.14026,   -472.742462, -515.882722, -515.882722, -472.742462,
       -43.14026,   43.14026,    472.742462,  515.882722,  515.882722,
       472.742462,  43.14026,    -245.349599, -388.967129, -634.316729,
       -634.316729, -388.967129, -245.349599, 245.349599,  388.967129,
       634.316729,  634.316729,  388.967129,  245.349599,  -85.74534,
       -516.861787, -602.607127, -602.607127, -516.861787, -85.74534,
       85.74534,    516.861787,  602.607127,  602.607127,  516.861787,
       85.74534,    -171.416892, -465.038645, -636.455537, -636.455537,
       -465.038645, -171.416892, 171.416892,  465.038645,  636.455537,
       636.455537,  465.038645,  171.416892,  0.018865,    -329.630743,
       -570.930392, -659.264065, -570.949257, -329.633322, -0.018865,
       329.630743,  570.930392,  659.264065,  570.949257,  329.633322,
       -47.940493,  -143.001204, -235.615124, -324.197605, -407.232986,
       -483.300486, -551.098581, -609.467226, -657.407719, -694.099784,
       -718.91561,  -731.43059,  -731.43059,  -718.91561,  -694.099784,
       -657.407719, -609.467226, -551.098581, -483.300486, -407.232986,
       -324.197605, -235.615124, -143.001204, -47.940493,  47.940493,
       143.001204,  235.615124,  324.197605,  407.232986,  483.300486,
       551.098581,  609.467226,  657.407719,  694.099784,  718.91561,
       731.43059,   731.43059,   718.91561,   694.099784,  657.407719,
       609.467226,  551.098581,  483.300486,  407.232986,  324.197605,
       235.615124,  143.001204,  47.940493};
}  // namespace

//---------------------------------------------------------------------------//
/*!
 * Constructor
 */
SimpleLZ::SimpleLZ(int sqrt_num_pmts) : sqrt_num_pmts_(sqrt_num_pmts) {}

//---------------------------------------------------------------------------//
/*!
 * Build a simple LZ geometry, or a notional geometry with a set number of
 * PMTs.
 *
 * This methods operates in two modes. If sqrt_num_pmts_ is zero (the default
 * value) a simplified model of the LZ top PMT detector array is created with
 * 253 PMTs placed in their real locations. If sqrt_num_pmts_ is positive, then
 * a square-pitch sqrt_num_pnts_ by sqrt_num_pmts_ array of detectors created.
 * The PMTs are spaced 1 mm apart, and sized such that the array is inscribed
 * within the world cyclinder. These notional geometries are useful for
 * performance studies.
 */
G4VPhysicalVolume* SimpleLZ::Construct()
{
    // Materials
    auto nist = G4NistManager::Instance();
    auto pmt_mat = nist->FindOrBuildMaterial("G4_Au");
    auto world_mat = nist->FindOrBuildMaterial("G4_Galactic");
    world_mat->SetName("vacuum");

    // World Cylinder
    G4Tubs* world_cylinder = new G4Tubs(
        "world_cylinder", 0, 780.288 * mm, 123 * mm, 0 * deg, 360 * deg);
    auto const world_lv
        = new G4LogicalVolume(world_cylinder, world_mat, "world");
    auto const world_pv = new G4PVPlacement(
        nullptr, G4ThreeVector(), world_lv, "world_pv", nullptr, false, 0, false);

    if (sqrt_num_pmts_ == 0)
    {
        // Simple LZ model
        G4Tubs* pmt
            = new G4Tubs("pmt", 0, 39 * mm, 123 * mm, 0 * deg, 360 * deg);
        auto const pmt_lv = new G4LogicalVolume(pmt, pmt_mat, "pmt_lv");

        for (auto i = 0; i < 253; ++i)
        {
            new G4PVPlacement(nullptr,
                              G4ThreeVector({pmt_x[i] * mm, pmt_y[i] * mm, 0}),
                              pmt_lv,
                              "pmt_pv",
                              world_lv,
                              false,
                              i,
                              false);
        }
    }
    else
    {
        // Notional model
        double a = 780.288 * std::sqrt(2) * mm;
        double pmt_pitch = a / (sqrt_num_pmts_);
        double pmt_radius = (pmt_pitch - 1 * mm) / 2;
        double offset = a / 2 - pmt_pitch / 2;

        G4Tubs* pmt
            = new G4Tubs("pmt", 0, pmt_radius, 123 * mm, 0 * deg, 360 * deg);
        auto const pmt_lv = new G4LogicalVolume(pmt, pmt_mat, "pmt_lv");

        for (auto i = 0; i < sqrt_num_pmts_; ++i)
        {
            for (auto j = 0; j < sqrt_num_pmts_; ++j)
            {
                new G4PVPlacement(nullptr,
                                  G4ThreeVector({pmt_pitch * i - offset,
                                                 pmt_pitch * j - offset,
                                                 0 * mm}),
                                  pmt_lv,
                                  "pmt_pv",
                                  world_lv,
                                  false,
                                  i,
                                  false);
            }
        }
    }

    return world_pv;
}

//---------------------------------------------------------------------------//
/*!
 * Set sensitive detectors.
 */
void SimpleLZ::ConstructSDandField()
{
    auto pmt_sd = new SensitiveDetector("pmt_sd");
    G4SDManager::GetSDMpointer()->AddNewDetector(pmt_sd);
    G4VUserDetectorConstruction::SetSensitiveDetector("pmt_lv", pmt_sd);
}

//---------------------------------------------------------------------------//
// PRIVATE
//---------------------------------------------------------------------------//
