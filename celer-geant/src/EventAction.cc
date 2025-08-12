//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/EventAction.cc
//---------------------------------------------------------------------------//
#include "EventAction.hh"

#include <G4Event.hh>
#include <corecel/io/Logger.hh>

//---------------------------------------------------------------------------//
/*!
 * Thread-local begin event action.
 */
void EventAction::BeginOfEventAction(G4Event const* event)
{
    CELER_LOG_LOCAL(status) << "Begin event " << event->GetEventID();
}
