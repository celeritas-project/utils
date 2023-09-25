//----------------------------------*-C++-*----------------------------------//
// Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file src/PrimaryGenerator.cc
//---------------------------------------------------------------------------//
#include "PrimaryGenerator.hh"

#include <G4ParticleTable.hh>
#include <G4SystemOfUnits.hh>

#include "OpticsDetector.hh"

//---------------------------------------------------------------------------//
/*!
 * Construct empty.
 */
PrimaryGenerator::PrimaryGenerator() : G4VUserPrimaryGeneratorAction() {}

//---------------------------------------------------------------------------//
/*!
 * Generate primary at each new event.
 */
void PrimaryGenerator::GeneratePrimaries(G4Event* event)
{
    G4ParticleGun gun;
    gun.SetParticleDefinition(
        G4ParticleTable::GetParticleTable()->FindParticle(11));  // e-
    gun.SetParticleEnergy(1 * GeV);
    gun.SetParticlePosition(G4ThreeVector(-3 * m, 0, 0));  // Before box
    gun.SetParticleMomentumDirection(G4ThreeVector(1, 0, 0));  // +x

    gun.GeneratePrimaryVertex(event);
}
