//----------------------------------*-C++-*----------------------------------//
// Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file src/OpticalPhysics.cc
//---------------------------------------------------------------------------//
#include "OpticalPhysics.hh"

#include <memory>
#include <G4OpAbsorption.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4OpticalParameters.hh>
#include <G4OpticalPhoton.hh>
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
    auto mgr = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
    auto params = G4OpticalParameters::Instance();

    // Attenuate photon based on material properties; Used to calculate mfp
    if (params->GetProcessActivation("OpAbsorption"))
    {
        auto absorption_proc = std::make_unique<G4OpAbsorption>();
        mgr->AddDiscreteProcess(absorption_proc.release());
    }

    // Update velocity, polarization, and direction when crossing boundaries
    if (params->GetProcessActivation("OpBoundary"))
    {
        auto boundary_proc = std::make_unique<G4OpBoundaryProcess>();
        mgr->AddDiscreteProcess(boundary_proc.release());
    }
}
