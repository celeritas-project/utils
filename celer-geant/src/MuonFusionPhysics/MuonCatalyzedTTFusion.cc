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
// Muon-catalyzed tritium-tritium fusion process
//

#include "MuonCatalyzedTTFusion.hh"

#include <G4EmCaptureCascade.hh>
#include <G4GenericMuonicAtom.hh>
#include <G4HadDecayGenerator.hh>
#include <G4HadProjectile.hh>
#include <G4HadSecondary.hh>
#include <G4HadronicInteraction.hh>
#include <G4HadronicProcessStore.hh>
#include <G4HadronicProcessType.hh>
#include <G4IonTable.hh>
#include <G4LinInterpolation.hh>
#include <G4MuonMinus.hh>
#include <G4MuonMinusBoundDecay.hh>
#include <G4ParticleDefinition.hh>
#include <G4RandomDirection.hh>
#include <G4SystemOfUnits.hh>
#include <corecel/io/Logger.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MuonCatalyzedTTFusion::MuonCatalyzedTTFusion(G4String const& name)
    : VMuonCatalyzedFusionProcess(name, fHadronic)
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
        G4cout << "MuonCatalyzedTTFusion is created." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MuonCatalyzedTTFusion::~MuonCatalyzedTTFusion()
{
    G4HadronicProcessStore::Instance()->DeRegisterExtraProcess(this);
    delete theTotalResult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool MuonCatalyzedTTFusion::IsApplicable(G4ParticleDefinition const& p)
{
    return (&p == G4GenericMuonicAtom::GenericMuonicAtom());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MuonCatalyzedTTFusion::PreparePhysicsTable(G4ParticleDefinition const& p)
{
    G4HadronicProcessStore::Instance()->RegisterParticleForExtraProcess(this,
                                                                        &p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MuonCatalyzedTTFusion::BuildPhysicsTable(G4ParticleDefinition const&)
{
    //  G4HadronicProcessStore::Instance()->PrintInfo(&p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double MuonCatalyzedTTFusion::AtRestGetPhysicalInteractionLength(
    G4Track const& track, G4ForceCondition* condition)
{
    *condition = NotForced;

    // check if this is the beginning of tracking
    if (theNumberOfInteractionLengthLeft < 0.)
    {
        ResetNumberOfInteractionLengthLeft();
    }

    // Check particle
    // Must be a MuH3 for this process to be applicable
    if (!particleMatches(&track, 1, 3))
    {
        if (verboseLevel > 1)
            G4cout << "MuonCatalyzedTTFusion GPIL: Muonic ion is not tritium, "
                      "returning infinity.";
        return std::numeric_limits<G4double>::infinity();
    }

    G4double meanCycleTime = GetMeanCycleTime(track);
    G4double interactionLength = theNumberOfInteractionLengthLeft
                                 * meanCycleTime;

    if (verboseLevel > 1)
        G4cout << "MuonCatalyzedTTFusion GPIL: Interaction Length: "
               << interactionLength / microsecond << " microsecond." << G4endl;

    return (interactionLength);
}

G4double MuonCatalyzedTTFusion::GetMeanCycleTime(G4Track const& track)
{
    // update density information
    updateMaterialInfo(&track);

    // Process does not happen unless tritium is present
    if (tritiumPhi < std::numeric_limits<G4double>::epsilon())
        return (std::numeric_limits<G4double>::infinity());

    // By default, process does not happen
    G4double meanCycleTime = std::numeric_limits<G4double>::infinity();

    // The cycling rate used is the mean of the two theory vales given in Table
    // 2 of Bogdanoda, doi:10.1134/s1063776109020034 This has reasonable
    // agreement with experimental data - total range for all theory and dat is
    // from 1.8/microsecond to 2.96/microsecond
    // TODO determine if temperature depenedence is expected in this cycling
    // rate - I think it is not because it is a non-resonant process
    G4double lambda_c = tritiumPhi * 2.8 / microsecond;

    if (lambda_c > std::numeric_limits<G4double>::epsilon())  // if rate is
                                                              // zero, keep
                                                              // cycle time at
                                                              // infinity so
                                                              // process does
                                                              // not trigger
        meanCycleTime = 1 / lambda_c;

    if (verboseLevel > 0)
    {
        G4cout << "MuonCatalyzedTTFusion GetMeanCycleTime: Temperature: "
               << temperature << " Tritium Phi: " << tritiumPhi << G4endl;
        G4cout << "MuonCatalyzedTTFusion GetMeanCycleTime: Mean Cycle Time: "
               << meanCycleTime / microsecond << " microsecond" << G4endl;
        G4cout << "MuonCatalyzedTTFusion GetMeanCycleTime: Global time: "
               << track.GetGlobalTime() / microsecond << " microsecond "
               << G4endl;
    }
    return (meanCycleTime);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange*
MuonCatalyzedTTFusion::AtRestDoIt(G4Track const& track, G4Step const& step)
{
    theTotalResult->Initialize(track);

    CELER_LOG_LOCAL(info) << "TT mucf: "
                          << track.GetParticleDefinition()->GetParticleName();
    auto const post_step = step.GetPostStepPoint();
    CELER_LOG_LOCAL(info) << "DT mucf: " << post_step->GetPosition().x() / cm
                          << " cm and " << post_step->GetKineticEnergy()
                          << " MeV";

    G4double finalGlobalTime = track.GetGlobalTime();
    G4double finalLocalTime = track.GetLocalTime();

    G4double meanCycleTime = GetMeanCycleTime(track);
    G4double interactionTime = theNumberOfInteractionLengthLeft * meanCycleTime;
    finalGlobalTime += interactionTime;
    finalLocalTime += interactionTime;

    // Main reference is doi:10.1134/s1063776109020034
    // T+T fusion reaction
    // MuH3 + t = mu- + alpha + neutron + neutron + 11.3 MeV   (non-sticking
    // branch) MuHe3 + t = MuAlpha + neutron + neutron + 11.3 MeV Assume the
    // muon is emited with the energy distibution of the dt reaction
    // TODO In reality, the neutrons and alpha are correlated (see Fig 9) but I
    // don't think the precise correlation is known, so for simplicity this
    // process emits them using a pure phase space generator.  This gets the
    // neutron energy spectrum slightly wrong.

    // Calculate kinetic energy of resultant muon by interpolating on CDF of
    // distribution function
    G4double muonKineticEnergy
        = muonEnergyInterpolator->CubicSplineInterpolation(G4UniformRand());

    // Physics parameters
    // From doi:10.1134/s1063776109020034
    G4double initialStickingFraction = 0.14;  // TODO not sure if 14% should be
                                              // the initial or final sticking
                                              // fraction, could refer to
                                              // Bogdanova reference to see
                                              // calculation
    // CODATA values
    G4double tritiumMass = 3.0155007134 * amu;
    G4double neutronMass = 1.00866491595 * amu;
    G4double alphaMass = 4.001506179127 * amu;
    G4double muonMass = 0.11344426775770586 * amu;

    // Check if alpha sticking has occured
    G4bool sticking = (G4UniformRand() <= initialStickingFraction);

    // Run the phase space generator to get product energies and momenta
    // (I did a manual check with the mass values here and this does indeed
    // result in the expected 11.3 MeV yield per Bogdonova) NOTE that muon
    // binding energy to both the tritium ion and alpha particle are ignored in
    // computation of the product energies and momenta

    G4double initialMass;
    std::vector<G4double> masses;
    if (sticking)
    {
        // sticking case
        initialMass = 2 * tritiumMass + muonMass;
        masses = {neutronMass, neutronMass, alphaMass + muonMass};
        // G4cout << "sticking case" << G4endl;
    }
    else
    {
        // no sticking case - muon handled seperately
        initialMass = 2 * tritiumMass;
        masses = {neutronMass, neutronMass, alphaMass};
        // G4cout << "no sticking case" << G4endl;
    }
    std::vector<G4LorentzVector> finalState;
    G4HadDecayGenerator generator;
    generator.Generate(initialMass, masses, finalState);

    // G4cout << "initial mass: " << initialMass/amu << " amu " << G4endl;
    // G4double totalMass = 0.0;
    // for (G4int i=0; i<3; i++) {
    //   G4cout << i << " mass in: " << masses[i]/amu << " amu" << G4endl;
    //   totalMass += masses[i];
    //  }
    // G4cout << "Mass difference: " << (initialMass - totalMass)/amu << "
    // energy: " << (initialMass - totalMass)*c_squared/MeV << " MeV" <<G4endl;

    // Create the muon and alpha particle, or muonic alpha particle, and add to
    // the list of secondaries

    // G4cout << "Phase space generator done." << G4endl;
    // for (G4int i=0; i<3; i++) {
    //   G4cout << "i " << i << " px:" << finalState[i].px() << " py:" <<
    //   finalState[i].py() << " pz:" << finalState[i].pz() << " mass: " <<
    //   finalState[i].m()/amu << " KE: " <<
    //   (finalState[i].e()-finalState[i].m())*c_light*c_light/MeV << " MeV "
    //   << G4endl;
    // }

    if (sticking)
    {
        // case with sticking
        if (verboseLevel > 1)
            G4cout << "MuonCatalyzedTTFusion: Alpha sticking occured."
                   << G4endl;
        theTotalResult->SetNumberOfSecondaries(3);
        G4IonTable* itp = G4IonTable::GetIonTable();
        G4ParticleDefinition* muonicAlphaDefinition = itp->GetMuonicAtom(2, 4);
        G4DynamicParticle* muonicAlpha = new G4DynamicParticle(
            muonicAlphaDefinition,
            finalState[2].vect().unit(),
            (finalState[2].e() - finalState[2].m()) * c_light * c_light);
        muonicAlpha->SetCharge(+1 * eplus);
        theTotalResult->AddSecondary(muonicAlpha, finalGlobalTime, true);
    }
    else
    {
        // case without sticking - muon goes in a random direction
        if (verboseLevel > 1)
            G4cout << "MuonCatalyzedFusion: Alpha sticking did not occur."
                   << G4endl;
        theTotalResult->SetNumberOfSecondaries(4);  // TODO test
        G4DynamicParticle* muon = new G4DynamicParticle(
            G4MuonMinus::MuonMinus(), G4RandomDirection(), muonKineticEnergy);
        theTotalResult->AddSecondary(muon, finalGlobalTime, true);
        // Alpha particle with equal and opposite momentum to neutron plus muon
        G4DynamicParticle* alpha = new G4DynamicParticle(
            G4Alpha::Alpha(),
            finalState[2].vect().unit(),
            (finalState[2].e() - finalState[2].m()) * c_light * c_light);
        theTotalResult->AddSecondary(alpha, finalGlobalTime, true);
    }

    // Create the neutrons and add to the list of secondaries
    G4DynamicParticle* neutron0 = new G4DynamicParticle(
        G4Neutron::Neutron(),
        finalState[0].vect().unit(),
        (finalState[0].e() - finalState[0].m()) * c_light * c_light);
    theTotalResult->AddSecondary(neutron0, finalGlobalTime, true);
    G4DynamicParticle* neutron1 = new G4DynamicParticle(
        G4Neutron::Neutron(),
        finalState[1].vect().unit(),
        (finalState[1].e() - finalState[1].m()) * c_light * c_light);
    theTotalResult->AddSecondary(neutron1, finalGlobalTime, true);

    // Kill primary particle (the MuH3)
    theTotalResult->ProposeTrackStatus(fStopAndKill);
    theTotalResult->ProposeLocalTime(finalLocalTime);

    ClearNumberOfInteractionLengthLeft();
    return theTotalResult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MuonCatalyzedTTFusion::ProcessDescription(std::ostream& outFile) const
{
    outFile << "Model of muon-catalyzed TT fusion process.";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
