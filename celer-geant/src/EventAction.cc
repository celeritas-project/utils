//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/EventAction.cc
//---------------------------------------------------------------------------//
#include "EventAction.hh"

#include <G4Event.hh>
#include <G4SystemOfUnits.hh>
#include <accel/TrackingManagerIntegration.hh>
#include <corecel/Assert.hh>
#include <corecel/io/Logger.hh>

#include "RootIO.hh"

//---------------------------------------------------------------------------//
/*!
 * At the end of each event, copy statistics from the local Celeritas state.
 */
void EventAction::BeginOfEventAction(G4Event const* event)
{
    CELER_LOG_LOCAL(status) << "Begin event " << event->GetEventID();
}

//---------------------------------------------------------------------------//
/*!
 * At the end of each event, copy statistics from the local Celeritas state.
 */
void EventAction::EndOfEventAction(G4Event const* event)
{
    auto& state = celeritas::TrackingManagerIntegration::Instance().GetState();
}
