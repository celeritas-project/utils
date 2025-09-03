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

    auto& json = JsonReader::Instance().at("celeritas");
    if (json.contains("offload_particles"))
    {
        valid_pdgs_ = json.at("offload_particles").get<std::vector<PDG>>();
    }
}

//---------------------------------------------------------------------------//
/*!
 * Callback interface between Geant4 and Celeritas.
 */
G4bool SensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    CELER_EXPECT(step);
    auto* track = step->GetTrack();
    CELER_ASSERT(track);
    auto* pd = track->GetParticleDefinition();
    CELER_ASSERT(pd);

    if (!this->is_pdg_valid(pd->GetPDGEncoding()))
    {
        // Do not score particles that aren't in the offload list
        return false;
    }

    auto* pre = step->GetPreStepPoint();
    CELER_ASSERT(pre);
    auto& pre_th = pre->GetTouchableHandle();
    CELER_ASSERT(pre_th);
    auto* phys_vol = pre_th->GetVolume();
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
    auto const& len = step->GetStepLength() / cm;

    SD_1D_FILL_WEIGHT(energy_dep, pos.x(), step->GetTotalEnergyDeposit());
    SD_1D_FILL(step_len, step->GetStepLength());
    SD_2D_FILL(pos_yz, pos.y(), pos.z());
    SD_1D_FILL(time, pre->GetGlobalTime());

    return true;

#undef SD_1D_FILL
#undef SD_2D_FILL
#undef SD_1D_FILL_WEIGHT
}

//---------------------------------------------------------------------------//
/*!
 * Only process PDGs that are listed the \c SetupOptions::offload_particles .
 *
 * If the list is empty, it defaults to the Celeritas basic EM list.
 */
bool SensitiveDetector::is_pdg_valid(PDG id) const
{
    return std::any_of(this->valid_pdgs_.begin(),
                       this->valid_pdgs_.end(),
                       [&id](PDG this_pdg) { return id == this_pdg; });
};
