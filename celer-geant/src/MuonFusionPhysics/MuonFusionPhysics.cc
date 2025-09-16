// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Authors: Ara Knaian (ara@nklabs.com)
//          Sridhar Tripathy
//          Kevin Lynch (krlynch@fnal.gov)
//
// Class Description:
//
// Muonic atom physics list
//

#include "MuonFusionPhysics.hh"

#include <G4Deuteron.hh>
#include <G4Electron.hh>
#include <G4GenericIon.hh>
#include <G4GenericMuonicAtom.hh>
#include <G4IonTable.hh>
#include <G4MuonMinus.hh>
#include <G4MuonPlus.hh>
#include <G4MuonicAtomDecay.hh>
#include <G4Neutron.hh>
#include <G4PionMinus.hh>
#include <G4PionPlus.hh>
#include <G4Positron.hh>
#include <G4ProcessManager.hh>
#include <G4ProcessTable.hh>
#include <G4Proton.hh>
#include <G4StepLimiter.hh>
#include <G4Triton.hh>
#include <G4hMultipleScattering.hh>
#include <G4ionIonisation.hh>

#include "MuonCatalyzedDDFusion.hh"
#include "MuonCatalyzedDTFusion.hh"
#include "MuonCatalyzedTTFusion.hh"
#include "MuonMinusAtomicCapture.hh"
#include "MuonStripping.hh"
#include "MuonicAtomSpinFlip.hh"
#include "MuonicAtomTransfer.hh"

MuonFusionPhysics::MuonFusionPhysics(G4String const& name)
    : G4VPhysicsConstructor(name)
{
}

MuonFusionPhysics::~MuonFusionPhysics() {}

void MuonFusionPhysics::ConstructParticle()
{
    G4GenericIon::GenericIonDefinition();
    G4GenericMuonicAtom::GenericMuonicAtomDefinition();
}

void MuonFusionPhysics::ConstructProcess()
{
    // ********** Muon minus processes

    G4ProcessManager* muonMinusProcessManager
        = G4MuonMinus::MuonMinus()->GetProcessManager();

    // Muon Minus Atomic Capture at Rest
    muonMinusProcessManager->AddRestProcess(new MuonMinusAtomicCapture());

    // Remove standard muMinusCaptureAtRest process (replaced by muonic atom
    // process)
    muonMinusProcessManager->RemoveProcess(
        G4ProcessTable::GetProcessTable()->FindProcess("muMinusCaptureAtRest",
                                                       "mu-"));

    // Limit step size of mu- so we can accurately track its trajectory in the
    // magnetic field
    muonMinusProcessManager->AddDiscreteProcess(new G4StepLimiter());

    // ********** Muonic Atom Processes

    auto* generic_muonic_atom = G4GenericMuonicAtom::GenericMuonicAtom();
    G4ProcessManager* genericMuonicAtomProcessManager
        = generic_muonic_atom->GetProcessManager();

    // Muonic Atom Decay
    // genericMuonicAtomProcessManager->AddRestProcess(new
    // G4MuonicAtomDecay());

    // Muon Catalyzed Fusion
    genericMuonicAtomProcessManager->AddRestProcess(
        new MuonCatalyzedDDFusion());
    genericMuonicAtomProcessManager->AddRestProcess(
        new MuonCatalyzedDTFusion());
    genericMuonicAtomProcessManager->AddRestProcess(
        new MuonCatalyzedTTFusion());

    // Spin flip for muonic deuterium from spin=3/2 to spin=1/2
    // genericMuonicAtomProcessManager->AddRestProcess(new
    // MuonicAtomSpinFlip());

    // Muon Transfer (TODO check and add)
    // genericMuonicAtomProcessManager->AddRestProcess(new
    // MuonicAtomTransfer());

    // Muon Stripping from MuHe3 and MuHe4
    // genericMuonicAtomProcessManager->AddDiscreteProcess(new
    // MuonStripping());

    // Add basic electromagnetic physics processes for muonic atoms
    // TODO make sure multiple scattering for muonic atom is using the correct
    // (reduced by one) charge
    // TODO add hadronic processes for muonic atom (generic ion has hElastic
    // and hInelastic in addition to these)
    // TODO do we need to modify muonic atom decay to apply to moving muonic
    // atom, or does it already?
    // TODO migrate to adding all processes using the physics list helper

    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

    // ph->RegisterProcess(fMuonStripping, generic_muonic_atom);
    ph->RegisterProcess(new G4ionIonisation(), generic_muonic_atom);
    ph->RegisterProcess(new G4hMultipleScattering(), generic_muonic_atom);
    // generic_muonic_atom->GetProcessManager()->DumpInfo();
}
