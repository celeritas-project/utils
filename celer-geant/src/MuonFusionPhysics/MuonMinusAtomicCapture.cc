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

#include "MuonMinusAtomicCapture.hh"

#include <G4EmCaptureCascade.hh>
#include <G4HadProjectile.hh>
#include <G4HadSecondary.hh>
#include <G4HadronicInteraction.hh>
#include <G4HadronicProcessStore.hh>
#include <G4HadronicProcessType.hh>
#include <G4IonTable.hh>
#include <G4MuonMinus.hh>
#include <G4MuonMinusBoundDecay.hh>
#include <G4ParticleDefinition.hh>
#include <G4RandomDirection.hh>
#include <G4SystemOfUnits.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MuonMinusAtomicCapture::MuonMinusAtomicCapture(G4String const& name)
    : G4VRestProcess(name, fHadronic)
    ,
    //  : G4HadronicProcess(name, fHadronAtRest),// name, process type
    fElementSelector(new G4ElementSelector())
    , fEmCascade(new G4EmCaptureCascade())
    ,  // Owned by InteractionRegistry
    theTotalResult(new G4ParticleChange())
{
    // Modify G4VProcess flags to emulate G4VRest instead of G4VDiscrete
    //  enableAtRestDoIt = true;
    //  enablePostStepDoIt = false;
    SetProcessSubType(fMuAtomicCapture);
    G4HadronicProcessStore::Instance()->RegisterExtraProcess(this);

    if (verboseLevel > 0)
        G4cout << "MuonMinusAtomicCapture (with Q1S) is created." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MuonMinusAtomicCapture::~MuonMinusAtomicCapture()
{
    G4HadronicProcessStore::Instance()->DeRegisterExtraProcess(this);
    delete theTotalResult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool MuonMinusAtomicCapture::IsApplicable(G4ParticleDefinition const& p)
{
    return (&p == G4MuonMinus::MuonMinus());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MuonMinusAtomicCapture::PreparePhysicsTable(G4ParticleDefinition const& p)
{
    G4HadronicProcessStore::Instance()->RegisterParticleForExtraProcess(this,
                                                                        &p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MuonMinusAtomicCapture::BuildPhysicsTable(G4ParticleDefinition const& p)
{
    G4HadronicProcessStore::Instance()->PrintInfo(&p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double MuonMinusAtomicCapture::AtRestGetPhysicalInteractionLength(
    G4Track const&, G4ForceCondition* condition)
{
    *condition = NotForced;
    return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange*
MuonMinusAtomicCapture::AtRestDoIt(G4Track const& track, G4Step const&)
{
    // if primary is not Alive then do nothing (how?)
    theTotalResult->Initialize(track);

    G4Nucleus* nucleus = &targetNucleus;
    // the call below actually sets the nucleus params;
    // G4Nucleus targetNucleus; is a member of G4HadronicProcess
    // G4Element* elm =
    fElementSelector->SelectZandA(track, nucleus);

    thePro.Initialise(track);  // thePro was G4HadProjectile from
                               // G4HadronicProcess

    // save track time and start capture from zero time
    thePro.SetGlobalTime(0.0);
    G4double time0 = track.GetGlobalTime();

    // Do the electromagnetic cascade in the nuclear field.
    // EM cascade should keep G4HadFinalState object,
    // because it will not be deleted at the end of this method
    //
    result = fEmCascade->ApplyYourself(thePro, *nucleus);
    G4double ebound = result->GetLocalEnergyDeposit();  // may need to carry
                                                        // this over; review
    G4double edep = 0.0;
    G4int nSecondaries = result->GetNumberOfSecondaries();
    thePro.SetBoundEnergy(ebound);

    G4int Z = nucleus->GetZ_asInt();
    G4int A = nucleus->GetA_asInt();
    if (verboseLevel > 0)
        G4cout << "MuonMinusAtomicCapture: Z = " << Z << " A = " << A << G4endl;

    // Muon transfer from deutrium to tritium during deexcitation cascade
    // In reality this effect would apply to other atoms and isotopes as well
    // If the G4element selector selected capture onto deuterium or tritium,
    // this code checks the relative fraction of each present in the mixture
    // and adjusts the fraction captured onto deiterium or tritium to account
    // for the muon transfer from deuterim to tritium during the deexcitation
    // casade.  TODO this should later be extended, perhaps using transportable
    // muonic atoms, to give correct results in the presence of impurities,
    // in the presence of protium, for other elements besides hydrogen, etc.

    // If the muon was initially captured onto deuterium or tritium
    if ((Z == 1) && ((A == 2) || (A == 3)))
    {
        if (verboseLevel > 0)
            G4cout << "MuonMinusAtomicCapture: capture onto deuterium or "
                      "tritium detected."
                   << G4endl;

        // Get material information
        G4Material* mat = track.GetMaterial();
        G4ElementVector const* elementVector
            = mat->GetElementVector();  // std vector of g4 element
        G4double const* fractionVector = mat->GetFractionVector();

        // Get relative mole fraction of D and T in mixture
        G4double deuteriumMoleFraction = 0.0;
        G4double tritiumMoleFraction = 0.0;
        G4double deuteriumMassFraction = 0.0;
        G4double tritiumMassFraction = 0.0;
        for (unsigned int i = 0; i < elementVector->size(); i++)
        {
            G4IsotopeVector const* isotopeVector
                = ((*elementVector)[i])->GetIsotopeVector();
            G4double const* relativeAbundance
                = ((*elementVector)[i])->GetRelativeAbundanceVector();
            G4double const atomicMass = ((*elementVector)[i])->GetA();
            for (unsigned int j = 0; j < isotopeVector->size(); j++)
            {
                G4int thisZ = ((*isotopeVector)[j])->GetZ();
                G4int thisN = ((*isotopeVector)[j])->GetN();
                G4double thisA = ((*isotopeVector)[j])->GetA();
                if ((thisZ == 1) && (thisN == 2))
                {
                    deuteriumMassFraction += fractionVector[i] * thisA
                                             * relativeAbundance[j]
                                             / atomicMass;
                }
                if ((thisZ == 1) && (thisN == 3))
                {
                    tritiumMassFraction += fractionVector[i] * thisA
                                           * relativeAbundance[j] / atomicMass;
                }
            }
        }
        if (tritiumMassFraction == 0)
        {
            deuteriumMoleFraction = 1.0;
        }
        else
        {
            G4double deuteriumMassRatio = deuteriumMassFraction
                                          / tritiumMassFraction;
            G4double deuteriumMoleRatio
                = 1.5 * deuteriumMassRatio;  // 1.5 = atomic weight of tritium
                                             // divided by atomic weight of
                                             // deuterium
            deuteriumMoleFraction = deuteriumMoleRatio
                                    / (deuteriumMoleRatio + 1.0);
        }
        tritiumMoleFraction = 1 - deuteriumMoleFraction;

        if (verboseLevel > 0)
        {
            G4cout << "deuteriumMassFraction: " << deuteriumMassFraction
                   << G4endl;
            G4cout << "tritiumMassFraction: " << tritiumMassFraction << G4endl;
            G4cout << "deuteriumMoleFraction: " << deuteriumMoleFraction
                   << G4endl;
            G4cout << "tritiumMoleFraction: " << tritiumMoleFraction << G4endl;
        }

        // Using Q1S formula, compute relative probability of having muonic
        // deuterium vs. muonic tritium at end of deexcitation cascade Using
        // data-fit Q1S formula from V.R. Bom 2005 q1s = 1/(1 + 7.2*C_t) p_dmu
        // = C_d*q1s p_tmu = 1 - p_dmu

        G4double q1s = 1 / (1 + 2.9 * tritiumMoleFraction);
        G4double deuteriumProbability = deuteriumMoleFraction * q1s;
        G4double randomNumber = G4UniformRand();

        if (verboseLevel > 0)
        {
            G4cout << "q1s: " << q1s << G4endl;
            G4cout << "deuteriumProbability: " << deuteriumProbability
                   << G4endl;
            G4cout << "randomNumber: " << randomNumber << G4endl;
        }

        // Select deuterium or tritium capture
        if (randomNumber <= deuteriumProbability)
            A = 2;
        else
            A = 3;
    }

    if (verboseLevel > 0)
        G4cout << "MuonMinusAtomicCapture, after Q1s correction, Z = " << Z
               << " A = " << A << G4endl;

    // creating the muonic atom
    ++nSecondaries;

    G4IonTable* itp = G4IonTable::GetIonTable();
    G4ParticleDefinition* muonicAtom = itp->GetMuonicAtom(Z, A);
    G4DynamicParticle* dp
        = new G4DynamicParticle(muonicAtom, G4RandomDirection(), 0.);

    // If this is a muonic hydrogen atom, set the initial spin
    // Yamashita, T., Kino, Y., Okutsu, K. et al. Roles of resonant muonic
    // molecule in new kinetics model and muon catalyzed fusion in compressed
    // gas. Sci Rep 12, 6393 (2022). https://doi.org/10.1038/s41598-022-09487-0
    if ((Z == 1) && ((A == 1) || (A == 2) || (A == 3)))
    {
        G4double upFraction;
        G4int spinUp;
        G4int spinDown;

        if (A == 2)
        {
            upFraction = 2. / 3.;
            spinUp = 3;  // spin 3/2
            spinDown = 1;  // spin 1/2
        }
        else
        {
            upFraction = 3. / 4.;
            spinUp = 2;  // spin 1
            spinDown = 0;  // spin 0
        }

        if (G4UniformRand() < upFraction)
        {
            dp->SetSpin(spinUp);
        }
        else
        {
            dp->SetSpin(spinDown);
        }
    }

    G4HadSecondary hadSec(dp);
    result->AddSecondary(hadSec);

    // Fill results
    //
    theTotalResult->ProposeTrackStatus(fStopAndKill);
    theTotalResult->ProposeLocalEnergyDeposit(edep);
    theTotalResult->SetNumberOfSecondaries(nSecondaries);
    G4double w = track.GetWeight();
    theTotalResult->ProposeWeight(w);

#ifdef G4VERBOSE
    if (GetVerboseLevel() > 1)
    {
        G4cout << __func__ << " nSecondaries " << nSecondaries << G4endl;
    }
#endif

    for (G4int i = 0; i < nSecondaries; ++i)
    {
        G4HadSecondary* sec = result->GetSecondary(i);

        // add track global time to the reaction time
        G4double time = sec->GetTime();
        if (time < 0.0)
        {
            time = 0.0;
        }
        time += time0;

#ifdef G4VERBOSE
        if (GetVerboseLevel() > 1)
        {
            G4cout << __func__ << " " << i << " Resulting secondary "
                   << sec->GetParticle()->GetPDGcode() << " "
                   << sec->GetParticle()->GetDefinition()->GetParticleName()
                   << G4endl;
        }
#endif

        // create secondary track
        G4Track* t = new G4Track(sec->GetParticle(), time, track.GetPosition());
        t->SetWeight(w * sec->GetWeight());

        t->SetTouchableHandle(track.GetTouchableHandle());
        theTotalResult->AddSecondary(t);
    }
    result->Clear();

    // fixme: needs to be done at the MuonicAtom level
    // if (epReportLevel != 0) { // G4HadronicProcess::
    //   CheckEnergyMomentumConservation(track, *nucleus);
    // }
    return theTotalResult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MuonMinusAtomicCapture::ProcessDescription(std::ostream& outFile) const
{
    outFile << "Stopping of mu- using default element selector, EM cascade"
            << "G4MuonicAtom is created\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
