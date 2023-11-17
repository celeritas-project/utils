//----------------------------------*-C++-*----------------------------------//
// Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file src/OpticalPhysics.cc
//---------------------------------------------------------------------------//
#include "OpticalPhysics.hh"

#include <memory>
#include <G4Cerenkov.hh>
#include <G4Electron.hh>
#include <G4Gamma.hh>
#include <G4OpAbsorption.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4OpticalParameters.hh>
#include <G4OpticalPhoton.hh>
#include <G4Positron.hh>
#include <G4ProcessManager.hh>

//---------------------------------------------------------------------------//
/*!
 * Construct empty.
 */
OpticalPhysics::OpticalPhysics() : G4VUserPhysicsList() {}

//---------------------------------------------------------------------------//
/*!
 * Construct particles.
 *
 * TODO: list particles and caveats
 */
void OpticalPhysics::ConstructParticle()
{
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();
    G4Gamma::GammaDefinition();
    G4OpticalPhoton::OpticalPhotonDefinition();
}

//---------------------------------------------------------------------------//
/*!
 * Construct physics.
 *
 * TODO: describe processes
 */
void OpticalPhysics::ConstructProcess()
{
    // Add mandatory transportation
    this->AddTransportation();

    // Add optical processes/models
    auto gamma_proc_mgr = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
    auto params = G4OpticalParameters::Instance();

    // Attenuate photon based on material properties; Used to calculate mfp
    if (params->GetProcessActivation("OpAbsorption"))
    {
        auto absorption_proc = new G4OpAbsorption();
        gamma_proc_mgr->AddDiscreteProcess(absorption_proc);
    }

    // Update velocity, polarization, and direction when crossing boundaries
    if (params->GetProcessActivation("OpBoundary"))
    {
        auto boundary_proc = new G4OpBoundaryProcess();
        gamma_proc_mgr->AddDiscreteProcess(boundary_proc);
        gamma_proc_mgr->SetProcessOrdering(
            boundary_proc, G4ProcessVectorDoItIndex::idxPostStep);
    }

    // Initialize Cherenkov radiation process
    auto e_proc_mgr = G4Electron::Electron()->GetProcessManager();
    auto cherenkov = new G4Cerenkov();
    e_proc_mgr->AddProcess(cherenkov);
    e_proc_mgr->SetProcessOrdering(cherenkov,
                                   G4ProcessVectorDoItIndex::idxPostStep);
}
