//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/RunAction.cc
//---------------------------------------------------------------------------//
#include "RunAction.hh"

#include <G4Threading.hh>
#include <accel/TrackingManagerIntegration.hh>
#include <celeritas/global/CoreParams.hh>
#include <corecel/io/Logger.hh>
#include <corecel/io/OutputRegistry.hh>

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

    if (G4Threading::IsWorkerThread())
    {
        // Construct thread-local ROOT I/O
        // Initialization at RunAction ensures geometry/SD data is available
        RootIO::Instance();
    }
}

//---------------------------------------------------------------------------//
/*!
 * Clear local data and return Celeritas to an invalid state.
 */
void RunAction::EndOfRunAction(G4Run const* run)
{
    auto& tmi = celeritas::TrackingManagerIntegration::Instance();
    if (G4Threading::IsWorkerThread())
    {
        // Write diagnostics to ROOT file and close it
        std::ostringstream diagnostics;
        tmi.GetParams().output_reg()->output(&diagnostics);
        auto rio = RootIO::Instance();
        rio->StoreDiagnostics(diagnostics.str());
        rio->Finalize();
    }
    tmi.EndOfRunAction(run);
}
