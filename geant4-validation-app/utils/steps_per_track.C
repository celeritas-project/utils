//----------------------------------*-C++-*----------------------------------//
// Copyright 2022-2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file steps_per_track.C
//! \brief Plot Geant4 vs. Celeritas steps per track and their relative error.
//---------------------------------------------------------------------------//
#include <iostream>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TText.h>

//---------------------------------------------------------------------------//
/*!
 * Enums for safely selecting particles and MC codes.
 */
enum MC
{
    g4,
    cel,
    size
};

enum PID
{
    e_plus,
    e_minus,
    photon,
    size
};

//---------------------------------------------------------------------------//
/*!
 * Geant4 vs. Celeritas steps per track for electrons, positrons, gammas,
 * and their relative errors.
 *
 * Make sure you have celeritas standard plotting options loaded into ROOT.
 * Instructions in celeritas-docs/utils/rootlogon.C
 *
 * Usage:
 * $ root steps_per_track.C
 *
 * The G4 bin data in this plot is extracted using
 * \c g4-validation-app/utils/diagnostics.C .
 *
 * \note
 * The data starts at bin[1], with bin[0] being reserved for ROOT's underflow:
 * - bin[0] = underflow, manually set to 0.
 * - bin[1] = number of tracks with 1 step.
 */
void steps_per_track()
{
    // Create histograms pointers for G4 and Celeritas for each particle type
    TH1D* h_steps[MC::size][PID::size];
    static int const n_bins = 180;

    for (int i = 0; i < MC::size; i++)
    {
        for (int j = 0; j < PID::size; j++)
        {
            // Initialize pointers
            h_steps[i][j] = new TH1D("", "", n_bins, 0, n_bins);
        }
    }

    // Histrograms for the relative error
    TH1D* h_error[PID::size];

    for (int j = 0; j < PID::size; j++)
    {
        // Initialize pointers
        h_error[j] = new TH1D("", "", n_bins, 0, n_bins);
    }

    // >>> Geant4 histogram data
    int const g4_positron_steps[200] = {
        0,     3265,  53562, 25525, 19295, 16622, 14570, 13094, 12452, 12371,
        12279, 12774, 12992, 13287, 13104, 13175, 12549, 11891, 11420, 10746,
        10026, 9485,  8848,  8399,  8011,  7243,  6893,  6489,  6071,  5578,
        5277,  4877,  4546,  4230,  3886,  3679,  3408,  3217,  2862,  2639,
        2561,  2305,  2213,  1966,  1917,  1766,  1689,  1528,  1302,  1277,
        1249,  1103,  1058,  983,   910,   820,   791,   744,   698,   569,
        621,   516,   499,   471,   416,   418,   382,   340,   347,   299,
        316,   277,   234,   256,   210,   203,   169,   167,   173,   148,
        123,   111,   121,   107,   88,    95,    79,    67,    70,    68,
        57,    69,    52,    36,    33,    31,    27,    27,    41,    27,
        25,    28,    26,    20,    8,     21,    14,    18,    16,    12,
        18,    6,     7,     5,     12,    6,     6,     7,     7,     8,
        3,     4,     4,     5,     1,     0,     3,     0,     1,     2,
        3,     0,     1,     1,     1,     3,     0,     0,     3,     3,
        0,     1,     1,     1,     2,     1,     0,     1,     0,     0,
        0,     0,     2,     0,     0,     0,     0,     1,     0,     0,
        0,     0,     0,     0,     0,     0,     0,     0,     0,     1,
        0,     0,     0,     0,     0,     0,     0,     1,     0,     0,
        0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
        0,     0,     0,     0,     0,     0,     0,     0,     0,     0};

    int const g4_electron_steps[200]
        = {0,     6077590, 510999, 215683, 140735, 103357, 80746, 68622, 61707,
           57290, 53406,   51653,  50312,  48193,  46529,  44082, 41471, 38682,
           35739, 33044,   30414,  27800,  25273,  23387,  21585, 19964, 18394,
           16822, 15257,   14330,  13107,  11972,  11130,  10546, 9376,  8834,
           8146,  7548,    7115,   6455,   6022,   5658,   5105,  4748,  4474,
           4106,  3878,    3583,   3302,   3087,   2973,   2713,  2418,  2354,
           2097,  1978,    1916,   1715,   1699,   1534,   1464,  1375,  1277,
           1158,  1073,    971,    955,    911,    796,    836,   776,   668,
           671,   670,     604,    567,    500,    528,    441,   475,   406,
           370,   377,     357,    299,    307,    278,    280,   275,   233,
           205,   203,     189,    183,    206,    150,    152,   149,   141,
           122,   120,     114,    91,     93,     88,     81,    80,    83,
           80,    58,      60,     58,     54,     48,     55,    54,    54,
           45,    26,      41,     33,     37,     34,     20,    30,    15,
           28,    17,      18,     19,     15,     17,     9,     16,    18,
           11,    11,      6,      7,      7,      12,     8,     5,     8,
           7,     4,       3,      6,      1,      2,      7,     3,     8,
           1,     6,       2,      5,      3,      1,      1,     1,     4,
           0,     0,       1,      1,      2,      0,      0,     0,     0,
           1,     0,       0,      0,      0,      0,      0,     1,     0,
           1,     0,       0,      0,      0,      0,      0,     0,     0,
           0,     0,       0,      0,      0,      0,      0,     0,     0,
           0,     0};

    int const g4_gamma_steps[200]
        = {0,      1080030, 400098, 448583, 269737, 244710, 183679, 178335,
           135128, 134065,  102571, 104546, 80652,  83727,  64373,  68094,
           52001,  56530,   43407,  47597,  36452,  40158,  31039,  34420,
           26581,  29837,   22875,  26008,  19893,  22865,  17610,  20035,
           15588,  17465,   13500,  15776,  12309,  13961,  10981,  12326,
           9659,   11055,   8629,   9839,   7803,   8889,   7114,   7967,
           6285,   7258,    5727,   6506,   5052,   5730,   4704,   5269,
           4214,   4877,    3971,   4346,   3506,   3881,   3249,   3578,
           2838,   3102,    2654,   2815,   2381,   2639,   2132,   2400,
           1966,   2181,    1782,   1864,   1605,   1728,   1398,   1536,
           1254,   1406,    1123,   1275,   1033,   1142,   873,    958,
           773,    838,     606,    702,    544,    546,    439,    430,
           335,    325,     242,    229,    181,    160,    116,    66,
           39,     23,      5,      1,      8,      1,      1,      0,
           0,      0,       0,      0,      1,      0,      0,      0,
           0,      0,       0,      0,      0,      0,      1,      0,
           0,      0,       0,      0,      0,      0,      0,      0,
           0,      0,       0,      0,      0,      0,      0,      0,
           0,      0,       0,      0,      0,      0,      0,      0,
           0,      0,       0,      0,      0,      0,      0,      0,
           0,      0,       0,      0,      0,      0,      0,      0,
           0,      0,       0,      0,      0,      0,      0,      0,
           0,      0,       0,      0,      0,      0,      0,      0,
           0,      0,       0,      0,      0,      0,      0,      0,
           0,      0,       0,      0,      0,      0,      0,      0};

    // >>> Celeritas histogram data
    int const cel_positron_steps[200] = {
        0,     24611, 42042, 22408, 18300, 15636, 14071, 13030, 12458, 12309,
        12444, 12869, 13166, 13600, 13677, 13539, 13199, 12203, 11526, 10936,
        10189, 9571,  8970,  8395,  7730,  7243,  6852,  6187,  5783,  5382,
        4996,  4618,  4220,  3963,  3621,  3353,  3152,  2929,  2687,  2522,
        2340,  2156,  2106,  1804,  1758,  1676,  1432,  1363,  1216,  1165,
        1064,  986,   966,   913,   784,   738,   710,   658,   614,   554,
        471,   500,   443,   425,   373,   374,   348,   283,   252,   241,
        223,   201,   210,   195,   173,   170,   147,   134,   137,   125,
        104,   99,    118,   79,    70,    80,    68,    80,    60,    54,
        35,    38,    44,    39,    32,    19,    32,    24,    25,    22,
        23,    22,    19,    21,    16,    13,    14,    12,    11,    9,
        10,    9,     6,     11,    4,     7,     7,     7,     5,     3,
        5,     10,    1,     3,     5,     0,     3,     2,     1,     3,
        1,     1,     1,     1,     1,     0,     0,     1,     2,     0,
        0,     0,     0,     1,     0,     1,     1,     1,     0,     0,
        0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
        1,     0,     0,     0,     0,     0,     0,     0,     0,     0,
        0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
        0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
        0,     0,     0,     0,     0,     0,     0,     0,     0,     0};

    int const cel_electron_steps[200]
        = {0,     6031570, 503932, 214939, 143562, 105225, 82600, 69006, 62122,
           57233, 54677,   52449,  51104,  49746,  47104,  45166, 41967, 38864,
           35759, 33129,   30574,  28308,  25924,  23268,  21565, 19845, 18151,
           16701, 15493,   14133,  12759,  11826,  10900,  9988,  9193,  8498,
           7885,  7215,    6644,   6131,   5709,   5441,   4900,  4636,  4294,
           3903,  3594,    3400,   3078,   2889,   2695,   2446,  2343,  2077,
           2026,  1907,    1739,   1624,   1606,   1507,   1320,  1330,  1155,
           1107,  1003,    990,    888,    838,    848,    731,   754,   654,
           635,   605,     496,    534,    485,    486,    428,   390,   400,
           368,   335,     314,    291,    286,    298,    246,   242,   237,
           232,   196,     174,    163,    167,    152,    159,   142,   146,
           113,   107,     105,    104,    98,     120,    78,    74,    67,
           61,    65,      71,     60,     52,     50,     44,    44,    33,
           33,    29,      42,     23,     35,     17,     25,    30,    29,
           23,    18,      17,     15,     12,     21,     19,    16,    4,
           6,     8,       4,      6,      5,      11,     8,     4,     7,
           6,     7,       2,      6,      7,      3,      4,     0,     3,
           2,     3,       0,      2,      1,      2,      0,     0,     3,
           0,     0,       1,      0,      1,      0,      1,     0,     1,
           3,     1,       1,      0,      0,      0,      0,     0,     0,
           0,     0,       1,      1,      0,      1,      0,     0,     0,
           0,     0,       0,      0,      0,      1,      0,     0,     0,
           0,     0};

    int const cel_gamma_steps[200]
        = {0,      1072950, 395342, 451120, 267313, 245985, 182717, 177819,
           133743, 133610,  102587, 104585, 79857,  83240,  63767,  67742,
           52235,  56503,   42934,  47432,  36498,  40721,  31097,  34978,
           26502,  29668,   22837,  26128,  20168,  23090,  17677,  20336,
           15552,  17855,   13792,  15792,  12349,  14158,  11183,  12780,
           9980,   11277,   8854,   10027,  7923,   9201,   7200,   8143,
           6494,   7499,    5923,   6676,   5237,   6105,   4803,   5599,
           4504,   5023,    3908,   4591,   3800,   4149,   3341,   3799,
           3037,   3316,    2750,   3053,   2411,   2778,   2205,   2464,
           2061,   2217,    1886,   2021,   1769,   1973,   1631,   1702,
           1447,   1529,    1245,   1380,   1077,   1204,   853,    970,
           824,    907,     673,    716,    569,    602,    473,    435,
           339,    325,     243,    231,    176,    183,    111,    63,
           43,     19,      8,      3,      0,      2,      0,      1,
           0,      1,       0,      0,      1,      0,      0,      1,
           0,      0,       0,      0,      0,      0,      0,      0,
           0,      0,       0,      0,      0,      0,      0,      0,
           0,      0,       0,      0,      0,      0,      0,      0,
           0,      0,       0,      0,      0,      0,      0,      0,
           0,      0,       0,      0,      0,      0,      0,      0,
           0,      0,       0,      0,      0,      0,      0,      0,
           0,      0,       0,      0,      0,      0,      0,      0,
           0,      0,       0,      0,      0,      0,      0,      0,
           0,      0,       0,      0,      0,      0,      0,      0,
           0,      0,       0,      0,      0,      0,      0,      0};

    // >>> Fill histogram data
    auto& h_steps_g4 = h_steps[MC::g4];
    auto& h_steps_cel = h_steps[MC::cel];

    for (int i = 0; i < h_steps_g4[photon]->GetNbinsX(); i++)
    {
        // Geant4
        h_steps_g4[e_plus]->SetBinContent(i, g4_positron_steps[i]);
        h_steps_g4[e_minus]->SetBinContent(i, g4_electron_steps[i]);
        h_steps_g4[photon]->SetBinContent(i, g4_gamma_steps[i]);

        // Celeritas
        h_steps_cel[e_plus]->SetBinContent(i, cel_positron_steps[i]);
        h_steps_cel[e_minus]->SetBinContent(i, cel_electron_steps[i]);
        h_steps_cel[photon]->SetBinContent(i, cel_gamma_steps[i]);
    }

    // Calculate relative error
    for (int i = 0; i < PID::size; i++)
    {
        h_error[i] = (TH1D*)h_steps_cel[i]->Clone();
        h_error[i]->Add(h_steps_g4[i], -1);
        h_error[i]->Divide(h_steps_cel[i]);
    }

    // Create clones to draw lines over the shaded error bar regions
    TH1D* h_steps_clones[PID::size];
    for (int j = 0; j < PID::size; j++)
    {
        h_steps_clones[j] = (TH1D*)h_steps_cel[j]->Clone();
    }

    h_steps_clones[e_plus]->SetLineColor(kGreen + 2);
    h_steps_clones[e_minus]->SetLineColor(kBlue);
    h_steps_clones[photon]->SetLineColor(kViolet);

    // Create canvas
    auto canvas = new TCanvas("", "", 750, 600);
    canvas->Divide(1, 2);

    // Create top pad
    TPad* pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.11);
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();

    // Set histogram attributes
    for (int i = 0; i < MC::size; i++)
    {
        h_steps[i][e_plus]->SetLineColor(kGreen + 2);
        h_steps[i][e_minus]->SetLineColor(kBlue);
        h_steps[i][photon]->SetLineColor(kViolet);

        h_steps[i][e_plus]->SetFillColorAlpha(kGreen + 2, 0.3);
        h_steps[i][e_minus]->SetFillColorAlpha(kBlue, 0.3);
        h_steps[i][photon]->SetFillColorAlpha(kViolet, 0.3);

        h_steps[i][e_plus]->SetMarkerColor(kGreen + 2);
        h_steps[i][e_minus]->SetMarkerColor(kBlue);
        h_steps[i][photon]->SetMarkerColor(kViolet);
    }

    for (int j = 0; j < PID::size; j++)
    {
        h_steps[MC::g4][j]->SetMarkerSize(1.1);
        h_steps[MC::cel][j]->SetMarkerSize(0);
    }

    // Make X marker a little bit bigger
    h_steps[MC::g4][e_minus]->SetMarkerSize(1.3);

    h_steps[MC::g4][e_plus]->SetMarkerStyle(53);
    h_steps[MC::g4][e_minus]->SetMarkerStyle(52);
    h_steps[MC::g4][photon]->SetMarkerStyle(55);

    for (int j = 0; j < PID::size; j++)
    {
        h_steps[MC::g4][j]->SetLineWidth(2);
        h_steps[MC::cel][j]->SetLineWidth(2);

        h_steps_clones[j]->SetLineWidth(2);
    }

    h_steps_g4[e_minus]->GetXaxis()->SetLabelOffset(99);
    h_steps_g4[e_minus]->GetYaxis()->SetRangeUser(0.5, 1E7);
    h_steps_g4[e_minus]->GetYaxis()->SetNdivisions(-5);

    // Draw plots
    h_steps_g4[e_minus]->Draw("P");
    h_steps_g4[e_plus]->Draw("P sames");
    h_steps_g4[photon]->Draw("P sames");

    h_steps_cel[e_minus]->Draw("E2 sames");
    h_steps_cel[e_plus]->Draw("E2 sames");
    h_steps_cel[photon]->Draw("E2 sames");

    for (int j = 0; j < PID::size; j++)
    {
        h_steps_clones[j]->Draw("sames");
    }

    // Draw legends
    auto legend_g4 = new TLegend(0.57, 0.65, 0.70, 0.86);
    legend_g4->SetHeader("Geant4");
    legend_g4->AddEntry(h_steps_g4[e_minus], "e^{-}", "p");
    legend_g4->AddEntry(h_steps_g4[e_plus], "e^{+}", "p");
    legend_g4->AddEntry(h_steps_g4[photon], "#gamma", "p");
    legend_g4->SetMargin(0.65);
    legend_g4->SetLineColor(kGray);
    // legend_g4->Draw();

    auto legend_cel = new TLegend(0.72, 0.65, 0.86, 0.86);
    legend_cel->SetHeader("Celeritas*");
    legend_cel->AddEntry(h_steps_cel[e_minus], "e^{-}", "lf");
    legend_cel->AddEntry(h_steps_cel[e_plus], "e^{+}", "lf");
    legend_cel->AddEntry(h_steps_cel[photon], "#gamma", "lf");
    legend_cel->SetMargin(0.65);
    legend_cel->SetLineColor(kGray);
    // legend_cel->Draw();

    auto legend_g42 = new TLegend(0.36, 0.73, 0.57, 0.86);
    legend_g42->SetHeader(" Geant4 v11.0.3");
    legend_g42->SetNColumns(3);
    legend_g42->AddEntry(h_steps_g4[e_minus], "e^{-}", "p");
    legend_g42->AddEntry(h_steps_g4[e_plus], "e^{+}", "p");
    legend_g42->AddEntry(h_steps_g4[photon], "#gamma^{ }", "p");
    legend_g42->SetMargin(0.37);
    legend_g42->SetLineColor(kGray);
    legend_g42->Draw();

    auto legend_cel2 = new TLegend(0.58, 0.73, 0.865, 0.86);
    legend_cel2->SetHeader("Celeritas (c0a251de4)");
    legend_cel2->SetNColumns(3);
    legend_cel2->AddEntry(h_steps_cel[e_minus], "e^{-}", "lf");
    legend_cel2->AddEntry(h_steps_cel[e_plus], "e^{+}", "lf");
    legend_cel2->AddEntry(h_steps_cel[photon], "#gamma^{ }", "lf");
    legend_cel2->SetMargin(0.65);
    legend_cel2->SetLineColor(kGray);
    legend_cel2->Draw();

    // Draw title
    auto title_text
        = new TText(0.11, 0.92, "Steps per track per particle type");
    title_text->SetNDC();
    title_text->SetTextColor(kGray);
    title_text->Draw();

    // Redraw axis, so the tick marks stay in front of the histogram lines
    pad1->RedrawAxis();

    // Move back to top canvas
    canvas->cd();

    // Create bottom pad
    TPad* pad2 = new TPad("pad2", "", 0.0, 0.0, 1.0, 0.3);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.33);
    pad2->SetLeftMargin(0.11);
    pad2->Draw();
    pad2->cd();

    h_error[e_minus]->GetXaxis()->SetTitle("Steps per track");
    h_error[e_minus]->GetXaxis()->CenterTitle();
    h_error[e_minus]->GetXaxis()->SetTitleSize(0.15);
    h_error[e_minus]->GetXaxis()->SetLabelSize(0.1155);
    h_error[e_minus]->GetXaxis()->SetLabelOffset(0.04);
    h_error[e_minus]->GetXaxis()->SetTitleOffset(1.2);

    h_error[e_minus]->GetYaxis()->SetTitle("Rel. Diff.");
    h_error[e_minus]->GetYaxis()->CenterTitle();
    h_error[e_minus]->GetYaxis()->SetTitleSize(0.13);
    h_error[e_minus]->GetYaxis()->SetTitleOffset(0.37);
    h_error[e_minus]->GetYaxis()->SetLabelSize(0.117);
    h_error[e_minus]->GetYaxis()->SetLabelOffset(0.011);
    h_error[e_minus]->GetYaxis()->SetRangeUser(-2.5, 1.5);

    h_error[e_minus]->SetLineColor(kBlue);
    h_error[e_plus]->SetLineColor(kGreen + 2);
    h_error[photon]->SetLineColor(kViolet);

    for (int j = 0; j < PID::size; j++)
    {
        h_error[j]->SetLineWidth(1);
    }

    h_error[e_minus]->Draw();
    h_error[e_plus]->Draw("sames");
    h_error[photon]->Draw("sames");

    // Draw line at 0
    auto line = new TLine(0, 0, 125, 0);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->SetLineStyle(7);
    line->Draw();

    // Refresh axis
    pad2->RedrawAxis();

    // canvas->SaveAs("steps-per-track.pdf");
}
