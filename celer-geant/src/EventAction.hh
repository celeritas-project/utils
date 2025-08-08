//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/EventAction.hh
//---------------------------------------------------------------------------//
#pragma once

#include <memory>
#include <G4UserEventAction.hh>

//---------------------------------------------------------------------------//
/*!
 * Print step statistics at the end of every event.
 */
class EventAction final : public G4UserEventAction
{
  public:
    EventAction() = default;

    void BeginOfEventAction(G4Event const*) final;
    void EndOfEventAction(G4Event const* event) final;
};
