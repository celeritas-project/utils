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

#include "DetectorConstruction.hh"

//---------------------------------------------------------------------------//
/*!
 * Construct empty.
 */
PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction()
{
}

//---------------------------------------------------------------------------//
/*!
 * Generate a simple primary.
 *
 * TODO: Replace by HepMC3 reader.
 */
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    G4ParticleGun particle_gun;
    particle_gun.SetParticleDefinition(
        G4ParticleTable::GetParticleTable()->FindParticle(11));
    particle_gun.SetParticleEnergy(5 * GeV);
    particle_gun.SetParticlePosition(G4ThreeVector());  // Origin
    particle_gun.SetParticleMomentumDirection(G4ThreeVector(1, 0, 0));  // +x
    particle_gun.GeneratePrimaryVertex(event);
}
