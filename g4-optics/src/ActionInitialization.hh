//----------------------------------*-C++-*----------------------------------//
// Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file src/ActionInitialization.hh
//---------------------------------------------------------------------------//
#pragma once

#include <G4VUserActionInitialization.hh>

//---------------------------------------------------------------------------//
/*!
 * Generate primaries.
 */
class ActionInitalization final : public G4VUserActionInitialization
{
  public:
    // Construct empty
    ActionInitalization();

    // Set up user actions
    void Build() const final;
};
