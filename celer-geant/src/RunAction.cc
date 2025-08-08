//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/RunAction.cc
//---------------------------------------------------------------------------//
#include "RunAction.hh"

#include <G4Threading.hh>
#include <accel/TrackingManagerIntegration.hh>
#include <corecel/io/Logger.hh>

#include "RootIO.hh"

//---------------------------------------------------------------------------//
/*!
 * Construct empty.
 */
RunAction::RunAction() : G4UserRunAction() {}

//---------------------------------------------------------------------------//
/*!
 * Initialize master and worker threads in Celeritas.
 */
void RunAction::BeginOfRunAction(G4Run const* run)
{
    CELER_LOG_LOCAL(status) << "Begin of run action";
    celeritas::TrackingManagerIntegration::Instance().BeginOfRunAction(run);
}

//---------------------------------------------------------------------------//
/*!
 * Clear local data and return Celeritas to an invalid state.
 */
void RunAction::EndOfRunAction(G4Run const* run)
{
    if (G4Threading::IsWorkerThread())
    {
        RootIO::Instance()->Finalize();
    }
    celeritas::TrackingManagerIntegration::Instance().EndOfRunAction(run);
}
