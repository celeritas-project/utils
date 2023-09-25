//----------------------------------*-C++-*----------------------------------//
// Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file src/ActionInitialization.cc
//---------------------------------------------------------------------------//
#include "ActionInitialization.hh"

#include "PrimaryGenerator.hh"

//---------------------------------------------------------------------------//
/*!
 * Construct empty.
 */
ActionInitalization::ActionInitalization() : G4VUserActionInitialization() {}

//---------------------------------------------------------------------------//
/*!
 * Set up user actions.
 */
void ActionInitalization::Build() const
{
    this->SetUserAction(new PrimaryGenerator());
}
