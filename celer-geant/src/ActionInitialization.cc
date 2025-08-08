//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/ActionInitialization.cc
//---------------------------------------------------------------------------//
#include "ActionInitialization.hh"

#include <corecel/io/Logger.hh>

#include "EventAction.hh"
#include "G4Threading.hh"
#include "PrimaryGeneratorAction.hh"
#include "RootIO.hh"
#include "RunAction.hh"

//---------------------------------------------------------------------------//
/*!
 * Construct empty.
 */
ActionInitialization::ActionInitialization() : G4VUserActionInitialization() {}

//---------------------------------------------------------------------------//
/*!
 * Set up Celeritas offload on master thread and initialize it via the
 * \c G4UserRunAction .
 */
void ActionInitialization::BuildForMaster() const
{
    // RunAction is responsible for initializing Celeritas
    this->SetUserAction(new RunAction());
}

//---------------------------------------------------------------------------//
/*!
 * Set up all worker thread user actions and Celeritas offload interface.
 */
void ActionInitialization::Build() const
{
    this->SetUserAction(new RunAction());
    this->SetUserAction(new PrimaryGeneratorAction());
    this->SetUserAction(new EventAction());

    // Construct thread-local ROOT I/O
    RootIO::Instance();
}
