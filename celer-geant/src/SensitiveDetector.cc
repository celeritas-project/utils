//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/SensitiveDetector.cc
//---------------------------------------------------------------------------//
#include "SensitiveDetector.hh"

#include <G4SystemOfUnits.hh>
#include <corecel/Assert.hh>
#include <corecel/io/Logger.hh>

#include "RootIO.hh"

//---------------------------------------------------------------------------//
/*!
 * Construct with sensitive detector name.
 */
SensitiveDetector::SensitiveDetector(std::string sd_name)
    : G4VSensitiveDetector(sd_name)
{
    CELER_VALIDATE(!sd_name.empty(),
                   << "must provide a valid sensitive detector name");
}

//---------------------------------------------------------------------------//
/*!
 * Callback interface between Geant4 and Celeritas.
 */
G4bool SensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    CELER_EXPECT(step);
    auto pre = step->GetPreStepPoint();
    CELER_ASSERT(pre);
    auto pre_th = pre->GetTouchableHandle();
    CELER_ASSERT(pre_th);
    auto phys_vol = pre_th->GetVolume();
    CELER_ASSERT(phys_vol);

    auto rio = RootIO::Instance();
    CELER_ASSERT(rio);
    auto& hists = rio->Histograms().Find(phys_vol->GetInstanceID(),
                                         phys_vol->GetCopyNo());

#define SD_1D_FILL(MEMBER, VALUE) hists.MEMBER.Fill(VALUE);
#define SD_2D_FILL(MEMBER, X, Y) hists.MEMBER.Fill(X, Y);
#define SD_1D_FILL_WEIGHT(MEMBER, VALUE, W)         \
    {                                               \
        auto& h = hists.MEMBER;                     \
        auto const i = h.FindBin(VALUE);            \
        h.SetBinContent(i, h.GetBinContent(i) + W); \
    }

    auto const& pos = pre->GetPosition() / cm;

    SD_1D_FILL_WEIGHT(energy_dep, pos.x(), step->GetTotalEnergyDeposit());
    SD_1D_FILL(step_len, step->GetStepLength() / cm);
    SD_2D_FILL(pos_yz, pos.y() / cm, pos.z() / cm);
    SD_1D_FILL(time, pre->GetGlobalTime());

    return true;

#undef SD_1D_FILL
#undef SD_2D_FILL
#undef SD_1D_FILL_WEIGHT
}
