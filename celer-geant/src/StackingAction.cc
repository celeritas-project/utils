//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/StackingAction.cc
//---------------------------------------------------------------------------//
#include "StackingAction.hh"

#include <algorithm>
#include <G4ClassificationOfNewTrack.hh>
#include <G4Track.hh>
#include <corecel/Assert.hh>

#include "JsonReader.hh"
#include "MakeCelerOptions.hh"

//---------------------------------------------------------------------------//
/*!
 * Construct with list of valid PDGs from the JSON input.
 */
StackingAction::StackingAction()
    : G4UserStackingAction()
    , offloaded_pdgs_(detail::offloaded_pdgs_from_json())
{
}

//---------------------------------------------------------------------------//
/*!
 * Assign \c fKill to all non-offloaded particles.
 */
G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(G4Track* const track)
{
    using detail::is_valid_celeritas_pdg;

    auto* pd = track->GetParticleDefinition();
    CELER_ASSERT(pd);
    return is_valid_celeritas_pdg(offloaded_pdgs_, pd->GetPDGEncoding())
               ? fUrgent
               : fKill;
}
