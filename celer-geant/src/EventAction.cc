//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/EventAction.cc
//---------------------------------------------------------------------------//
#include "EventAction.hh"

#include <G4Event.hh>
#include <corecel/io/Logger.hh>

#include "JsonReader.hh"

//---------------------------------------------------------------------------//
/*!
 * Construct thread-local event action and set up event logging.
 */
EventAction::EventAction() : G4UserEventAction()
{
    auto json = JsonReader::Instance();
    log_progress_ = (json.contains("log_progress"))
                        ? json.at("log_progress").get<size_t>()
                        : 1;
}

//---------------------------------------------------------------------------//
/*!
 * Thread-local begin of event action.
 */
void EventAction::BeginOfEventAction(G4Event const* event)
{
    if (auto const id = event->GetEventID(); id % log_progress_ == 0)
    {
        CELER_LOG_LOCAL(status) << "Begin event " << id;
    }
}
