//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/PrimaryGeneratorAction.cc
//---------------------------------------------------------------------------//
#include "PrimaryGeneratorAction.hh"

#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4SystemOfUnits.hh>
#include <corecel/Assert.hh>

#include "DetectorConstruction.hh"

//---------------------------------------------------------------------------//
/*!
 * Generate primaries.
 */
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    // CELER_EXPECT(event);
    // CELER_EXPECT(event->GetEventID() > 0);

    G4ParticleGun particle_gun;
    particle_gun.SetParticleDefinition(
        G4ParticleTable::GetParticleTable()->FindParticle(11));
    particle_gun.SetParticleEnergy(10 * MeV);
    particle_gun.SetParticlePosition(G4ThreeVector());  // Origin
    particle_gun.SetParticleMomentumDirection(G4ThreeVector(1, 0, 0));  // +x
    particle_gun.GeneratePrimaryVertex(event);
}
