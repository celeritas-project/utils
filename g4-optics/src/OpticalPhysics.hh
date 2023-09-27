//----------------------------------*-C++-*----------------------------------//
// Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file src/OpticalPhysics.hh
//---------------------------------------------------------------------------//
#pragma once

#include <G4VUserPhysicsList.hh>

class OpticalPhysics : public G4VUserPhysicsList
{
  public:
    // Construct empty
    OpticalPhysics();

    // Initialize minimal list of particles needed for optical physics
    void ConstructParticle() override;

    // Initialize optical physics models
    void ConstructProcess() override;
};
