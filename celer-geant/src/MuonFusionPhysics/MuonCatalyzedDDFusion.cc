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
// Muon-catalyzed deuterium-deuterium fusion process
//

#include "MuonCatalyzedDDFusion.hh"

#include <G4EmCaptureCascade.hh>
#include <G4GenericMuonicAtom.hh>
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

#include "VMuonCatalyzedFusionProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MuonCatalyzedDDFusion::MuonCatalyzedDDFusion(G4String const& name)
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

    // DD fusion rates from
    // Balin, D.V., Ganzha, V.A., Kozlov, S.M. et al. High precision study of
    // muon catalyzed fusion in D2 and HD gas. Phys. Part. Nuclei 42, 185–214
    // (2011). https://doi.org/10.1134/S106377961102002X

    G4int static const numLambdaDDOneHalf = 21;

    G4double temperatureDDmuFromSpinOneHalf[numLambdaDDOneHalf]
        = {0.0,
           36.43103676856282 * kelvin,
           49.84753686046977 * kelvin,
           62.30013430590361 * kelvin,
           72.8249264583911 * kelvin,
           94.80254726485413 * kelvin,
           118.65929788178391 * kelvin,
           134.8637502794038 * kelvin,
           174.9275151614521 * kelvin,
           212.16105692405534 * kelvin,
           260.90506935936025 * kelvin,
           325.03565799325077 * kelvin,
           380.6089687138337 * kelvin,
           434.29674495520254 * kelvin,
           486.0836155123639 * kelvin,
           564.7393523983244 * kelvin,
           668.3489596576321 * kelvin,
           764.2819309563899 * kelvin,
           860.2110594538995 * kelvin,
           935.0316806939586 * kelvin,
           998.3392693938174 * kelvin};

    G4double lamdaDDmuFromSpinOneHalf[numLambdaDDOneHalf]
        = {0.0,
           0.005464463375221662e6 / second,
           0.020292872692477815e6 / second,
           0.04715981762065269e6 / second,
           0.098103833770665e6 / second,
           0.2721222062034343e6 / second,
           0.5362926959236205e6 / second,
           0.7584940318094278e6 / second,
           1.2388537948623028e6 / second,
           1.6110396634730852e6 / second,
           1.9710284809214524e6 / second,
           2.2225696043387786e6 / second,
           2.3059026708109145e6 / second,
           2.3141137763784285e6 / second,
           2.283275296360035e6 / second,
           2.19792475923249e6 / second,
           2.05207892345204e6 / second,
           1.921378528091838e6 / second,
           1.7996962265613146e6 / second,
           1.7144033314524956e6 / second,
           1.6473195500592115e6 / second};

    G4int static const numLambdaDDThreeHalves = 28;

    G4double temperatureDDmuFromSpinThreeHalves[numLambdaDDThreeHalves]
        = {0.0 * kelvin,
           5.4317995646105715 * kelvin,
           8.023769006655044 * kelvin,
           15.563345056011997 * kelvin,
           21.23403876489806 * kelvin,
           24.19235525928977 * kelvin,
           30.970416194589205 * kelvin,
           36.77176514592074 * kelvin,
           40.6299375993124 * kelvin,
           48.3052926261127 * kelvin,
           54.996250066448454 * kelvin,
           65.48389514020138 * kelvin,
           83.57516295078466 * kelvin,
           104.55813870078731 * kelvin,
           127.50862868996978 * kelvin,
           141.88390769335217 * kelvin,
           157.2358986807023 * kelvin,
           176.44862398894298 * kelvin,
           208.18055529758965 * kelvin,
           248.59273415949275 * kelvin,
           289.00491302139596 * kelvin,
           347.6774429488119 * kelvin,
           437.0821353934164 * kelvin,
           532.1831402218961 * kelvin,
           626.2766906564209 * kelvin,
           748.1620201862349 * kelvin,
           868.0888020131156 * kelvin,
           992.7877025236317 * kelvin};

    G4double lamdaDDmuFromSpinThreeHalves[numLambdaDDThreeHalves]
        = {0.0e6 / second,
           3.002924371750031 / second,
           3.670220083632165 / second,
           3.976719989803767 / second,
           4.169019528475477 / second,
           3.9765902952616323 / second,
           3.820175795347777 / second,
           3.705860143810432 / second,
           3.6516939388136374 / second,
           3.6395545296699483 / second,
           3.6875501565621276 / second,
           3.8256690797323616 / second,
           4.119986345246231 / second,
           4.378188003927343 / second,
           4.519125621813277 / second,
           4.533939620625851 / second,
           4.506654771061913 / second,
           4.419191653948062 / second,
           4.202281855381299 / second,
           3.8650011112100278 / second,
           3.5277203670387562 / second,
           3.0879607598755188 / second,
           2.527498765500099 / second,
           2.0992156842652023 / second,
           1.785176202044242 / second,
           1.500779127953113 / second,
           1.312603875721246 / second,
           1.1754591026674808 / second};

    lambdaDDmuFromSpinOneHalfInterpolator = new G4DataInterpolation(
        temperatureDDmuFromSpinOneHalf,
        lamdaDDmuFromSpinOneHalf,
        numLambdaDDOneHalf,
        (lamdaDDmuFromSpinOneHalf[1] - lamdaDDmuFromSpinOneHalf[0])
            / (temperatureDDmuFromSpinOneHalf[1]
               - temperatureDDmuFromSpinOneHalf[0]),
        (lamdaDDmuFromSpinOneHalf[numLambdaDDOneHalf - 1]
         - lamdaDDmuFromSpinOneHalf[numLambdaDDOneHalf - 2])
            / (temperatureDDmuFromSpinOneHalf[numLambdaDDOneHalf - 1]
               - temperatureDDmuFromSpinOneHalf[numLambdaDDOneHalf - 2]));

    lambdaDDmuFromSpinThreeHalvesInterpolator = new G4DataInterpolation(
        temperatureDDmuFromSpinThreeHalves,
        lamdaDDmuFromSpinThreeHalves,
        numLambdaDDThreeHalves,
        (lamdaDDmuFromSpinThreeHalves[1] - lamdaDDmuFromSpinThreeHalves[0])
            / (temperatureDDmuFromSpinThreeHalves[1]
               - temperatureDDmuFromSpinThreeHalves[0]),
        (lamdaDDmuFromSpinThreeHalves[numLambdaDDOneHalf - 1]
         - lamdaDDmuFromSpinThreeHalves[numLambdaDDThreeHalves - 2])
            / (temperatureDDmuFromSpinThreeHalves[numLambdaDDThreeHalves - 1]
               - temperatureDDmuFromSpinThreeHalves[numLambdaDDThreeHalves - 2]));

    if (verboseLevel > 0)
        G4cout << "MuonCatalyzedDDFusion is created." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MuonCatalyzedDDFusion::~MuonCatalyzedDDFusion()
{
    G4HadronicProcessStore::Instance()->DeRegisterExtraProcess(this);
    delete theTotalResult;

    delete lambdaDDmuFromSpinOneHalfInterpolator;
    delete lambdaDDmuFromSpinThreeHalvesInterpolator;
}

// TODO should these three functions be moved to the virtual base class?

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool MuonCatalyzedDDFusion::IsApplicable(G4ParticleDefinition const& p)
{
    return (&p == G4GenericMuonicAtom::GenericMuonicAtom());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MuonCatalyzedDDFusion::PreparePhysicsTable(G4ParticleDefinition const& p)
{
    G4HadronicProcessStore::Instance()->RegisterParticleForExtraProcess(this,
                                                                        &p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MuonCatalyzedDDFusion::BuildPhysicsTable(G4ParticleDefinition const&)
{
    //  G4HadronicProcessStore::Instance()->PrintInfo(&p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double MuonCatalyzedDDFusion::AtRestGetPhysicalInteractionLength(
    G4Track const& track, G4ForceCondition* condition)
{
    *condition = NotForced;

    // check if this is the beginning of tracking
    if (theNumberOfInteractionLengthLeft < 0.)
    {
        ResetNumberOfInteractionLengthLeft();
    }

    // Check particle
    // Must be a MuH2 for this process to be applicable

    if (!particleMatches(&track, 1, 2))
    {
        if (verboseLevel > 1)
            G4cout << "MuonCatalyzedDDFusion GPIL: Muonic ion is not "
                      "deuterium, returning infinity.";
        return std::numeric_limits<G4double>::infinity();
    }

    G4double meanCycleTime = GetMeanCycleTime(track);
    G4double interactionLength = theNumberOfInteractionLengthLeft
                                 * meanCycleTime;

    if (verboseLevel > 1)
        G4cout << "MuonCatalyzedDDFusion GPIL: Interaction Length: "
               << interactionLength / microsecond << " microsecond." << G4endl;

    return (interactionLength);
}

// TODO the checking we are doing in these functions is overkill because GPIL
// alrady checks the muonic atom type - remove in all 3 to simplify

G4double MuonCatalyzedDDFusion::GetMeanCycleTime(G4Track const& track)
{
    G4double lambda_c;

    // update density information
    updateMaterialInfo(&track);

    // Process does not happen unless deuterium is present
    if (deuteriumPhi < std::numeric_limits<G4double>::epsilon())
        return (std::numeric_limits<G4double>::infinity());

    // By default, process does not happen
    G4double meanCycleTime = std::numeric_limits<G4double>::infinity();

    // Pure D2 fusion
    // TODO Assumes linear density dependence TODO fix
    // TODO: add proper allowance for back-decay rates and fusion times now
    // that you are using Faifman rates
    // TODO: should this by deuteriumPhi or DDPhi?  Can DD fusion happen to a
    // DT molecule or HD molecule???????????

    G4int spin = (G4int)round(2.0 * track.GetDynamicParticle()->GetSpin());
    // G4cout << "MuCF spin = " << spin << G4endl;
    if (spin == 3)
    {
        lambda_c = deuteriumPhi
                   * lambdaDDmuFromSpinThreeHalvesInterpolator
                         ->CubicSplineInterpolation(temperature);
        // G4cout << "MuCF from spin 3/2 phi = " << phi << " lambda_c = " <<
        // lambda_c/(1/microsecond ) << "/us temperature: " <<
        // temperature/kelvin << " kelvin"  <<  G4endl;
    }
    else
    {
        lambda_c
            = deuteriumPhi
              * lambdaDDmuFromSpinOneHalfInterpolator->CubicSplineInterpolation(
                  temperature);
        // G4cout << "MuCF from spin 1/2 phi = " << phi << " lambda_c = " <<
        // lambda_c/(1/microsecond ) << "/us temperature: " <<
        // temperature/kelvin << " kelvin"  <<  G4endl;
    }

    if (lambda_c > std::numeric_limits<G4double>::epsilon())  // if rate is
                                                              // zero, keep
                                                              // cycle time at
                                                              // infinity so
                                                              // process does
                                                              // not trigger
        meanCycleTime = 1 / lambda_c;

    if (verboseLevel > 0)
    {
        G4cout << "MuonCatalyzedDDFusion GetMeanCycleTime: Temperature: "
               << temperature << " Phi: " << deuteriumPhi << G4endl;
        G4cout << "MuonCatalyzedDDFusion GetMeanCycleTime: Mean Cycle Time: "
               << meanCycleTime / microsecond << " microsecond" << G4endl;
        G4cout << "MuonCatalyzedDDFusion GetMeanCycleTime: Global time: "
               << track.GetGlobalTime() / microsecond << " microsecond "
               << G4endl;
    }

    return (meanCycleTime);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange*
MuonCatalyzedDDFusion::AtRestDoIt(G4Track const& track, G4Step const&)
{
    theTotalResult->Initialize(track);

    G4double finalGlobalTime = track.GetGlobalTime();
    G4double finalLocalTime = track.GetLocalTime();

    G4double meanCycleTime = GetMeanCycleTime(track);
    G4double interactionTime = theNumberOfInteractionLengthLeft * meanCycleTime;
    finalGlobalTime += interactionTime;
    finalLocalTime += interactionTime;

    // Calculate kinetic energy of resultant muon by interpolating on CDF of
    // distribution function Note: this uses the dt distribution.  TODO use
    // proper distribution for DD also
    G4double muonKineticEnergy
        = muonEnergyInterpolator->CubicSplineInterpolation(G4UniformRand());

    // Temperature dependent branching ratio between
    // n/He-3 channel and t/p channel
    G4double branchingRatio;
    if (temperature / kelvin < 50)
        branchingRatio = 1.0;
    else if (temperature / kelvin < 100)
        branchingRatio = 1.0 + (0.44 / 50) * (temperature / kelvin - 50);
    else
        branchingRatio = 1.44;

    if (verboseLevel > 0)
    {
        G4cout << "Branching ratio: " << branchingRatio << G4endl;
        G4cout << "Probability: " << (branchingRatio / (branchingRatio + 1))
               << G4endl;
    }

    G4double randomNumber = G4UniformRand();
    if (randomNumber < (branchingRatio / (branchingRatio + 1)))
    {
        // n/He-3 channel
        // 3.3 MeV total energy
        // He-3 mass is 3 times the mass of the neutron
        // so neutron KE is 3/4 and He-3 KE is 1/4
        // Determine if there was sticking
        // TODO: The number used here is a final sticking fraction
        // TODO: Properly, this should be the initial sticking fraction
        randomNumber = G4UniformRand();
        G4double initialStickingFraction = 0.122;
        if (randomNumber <= initialStickingFraction)
        {
            // sticking case
            if (verboseLevel > 1)
                G4cout << "MuonCatalyzedFusion: DD He-3 channel, sticking "
                          "occured."
                       << G4endl;
            theTotalResult->SetNumberOfSecondaries(2);
            // neutron
            G4double neutronKineticEnergy = 3.3 * (0.75) * MeV;
            G4DynamicParticle* neutron = new G4DynamicParticle(
                G4Neutron::Neutron(), G4RandomDirection(), neutronKineticEnergy);
            theTotalResult->AddSecondary(neutron, finalGlobalTime, true);
            G4ThreeVector neutronMomentum = neutron->GetMomentum();
            // muonic He-3
            G4IonTable* itp = G4IonTable::GetIonTable();
            G4ParticleDefinition* muonicHeliumDefinition
                = itp->GetMuonicAtom(2, 3);  // MuHe3
            G4DynamicParticle* muonicHelium = new G4DynamicParticle(
                muonicHeliumDefinition, -neutronMomentum);
            muonicHelium->SetCharge(+1 * eplus);
            theTotalResult->AddSecondary(muonicHelium, finalGlobalTime, true);
        }
        else
        {
            // no sticking case
            if (verboseLevel > 1)
                G4cout << "MuonCatalyzedFusion: DD He-3 channel." << G4endl;
            theTotalResult->SetNumberOfSecondaries(3);
            // neutron
            G4double neutronKineticEnergy = 3.3 * (0.75) * MeV;
            G4DynamicParticle* neutron = new G4DynamicParticle(
                G4Neutron::Neutron(), G4RandomDirection(), neutronKineticEnergy);
            theTotalResult->AddSecondary(neutron, finalGlobalTime, true);
            G4ThreeVector neutronMomentum = neutron->GetMomentum();
            // muon
            G4DynamicParticle* muon
                = new G4DynamicParticle(G4MuonMinus::MuonMinus(),
                                        G4RandomDirection(),
                                        muonKineticEnergy);
            theTotalResult->AddSecondary(muon, finalGlobalTime, true);
            G4ThreeVector muonMomentum = muon->GetMomentum();
            // He-3
            G4IonTable* itp = G4IonTable::GetIonTable();
            G4ParticleDefinition* heliumThreeDefinition
                = itp->GetIon(2, 3, 0);  // He3
            G4DynamicParticle* heliumThree = new G4DynamicParticle(
                heliumThreeDefinition, -neutronMomentum - muonMomentum);
            heliumThree->SetCharge(+2 * eplus);
            theTotalResult->AddSecondary(heliumThree, finalGlobalTime, true);
        }
    }
    else
    {
        // t/p channel
        // TODO: Does not take into account sticking to tritium
        // Sticking to tritium would typically cause a subsequent DT fusion,
        // and then release a free muon so sticking to muon would typically
        // increase the yield somewhat
        //
        // 4.03 MeV is the total energy output
        // triton mass is 3X proton mass
        // from mass conservation proton has 3/4 of KE, and triton has 1/4 of
        // KE proton KE is therefore 3.02 MeV
        if (verboseLevel > 1)
            G4cout << "MuonCatalyzedFusion: DD Proton/Triton channel."
                   << G4endl;
        theTotalResult->SetNumberOfSecondaries(3);
        // proton
        G4double protonKineticEnergy = 4.03 * (0.75) * MeV;
        G4DynamicParticle* proton = new G4DynamicParticle(
            G4Proton::Proton(), G4RandomDirection(), protonKineticEnergy);
        theTotalResult->AddSecondary(proton, finalGlobalTime, true);
        G4ThreeVector protonMomentum = proton->GetMomentum();
        // muon
        G4DynamicParticle* muon = new G4DynamicParticle(
            G4MuonMinus::MuonMinus(), G4RandomDirection(), muonKineticEnergy);
        theTotalResult->AddSecondary(muon, finalGlobalTime, true);
        G4ThreeVector muonMomentum = muon->GetMomentum();
        // triton
        G4DynamicParticle* triton = new G4DynamicParticle(
            G4Triton::Triton(), -(protonMomentum + muonMomentum));
        theTotalResult->AddSecondary(triton, finalGlobalTime, true);
    }

    // Kill primary particle (the MuH2)
    //
    theTotalResult->ProposeTrackStatus(fStopAndKill);
    theTotalResult->ProposeLocalTime(finalLocalTime);

    ClearNumberOfInteractionLengthLeft();
    return theTotalResult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MuonCatalyzedDDFusion::ProcessDescription(std::ostream& outFile) const
{
    outFile << "Model of muon catalyzed D-D fusion.\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
