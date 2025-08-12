//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/SensitiveDetector.cc
//---------------------------------------------------------------------------//
#include "SensitiveDetector.hh"

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

    auto& hists = RootIO::Instance()->Histograms().Find(
        phys_vol->GetInstanceID(), phys_vol->GetCopyNo());

    hists.energy.Fill(step->GetTotalEnergyDeposit());
    hists.time.Fill(pre->GetGlobalTime());

    return true;
}
