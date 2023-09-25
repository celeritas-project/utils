//----------------------------------*-C++-*----------------------------------//
// Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file src/PrimaryGenerator.hh
//---------------------------------------------------------------------------//
#pragma once

#include <G4ParticleGun.hh>
#include <G4VUserPrimaryGeneratorAction.hh>

//---------------------------------------------------------------------------//
/*!
 * Generate primaries.
 */
class PrimaryGenerator final : public G4VUserPrimaryGeneratorAction
{
  public:
    // Construct empty
    PrimaryGenerator();

    // Place primaries in the event simulation
    void GeneratePrimaries(G4Event* event) final;
};
