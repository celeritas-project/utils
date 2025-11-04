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
// Muon-catalyzed deuterium-tritium fusion process
//

// TODO We can use a G4MaterialPropertiesTable to encode the fraction of
// DD/DT/T2/HD/HT/TT in the gas also could encode ortho/para fraction that way
// if no material properties table, assume something default (maybe equilibrium
// state)

// TODO add nontrivial density dependence, for example from this paper:
// doi:10.1016/s0375-9474(99)00168-2

// TODO plot results versus experimental data (as was done in 2022 paper)

#include "MuonCatalyzedDTFusion.hh"

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
#include <corecel/io/Logger.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MuonCatalyzedDTFusion::MuonCatalyzedDTFusion(G4String const& name)
    : VMuonCatalyzedFusionProcess(name, fHadronic)
    ,
    //  : G4HadronicProcess(name, fHadronAtRest),// name, process type
    fElementSelector(new G4ElementSelector())
    , fEmCascade(new G4EmCaptureCascade())
    ,  // Owned by InteractionRegistry
    theTotalResult(new G4ParticleChange())
{
    verboseLevel = 2;
    // Modify G4VProcess flags to emulate G4VRest instead of G4VDiscrete
    //  enableAtRestDoIt = true;
    //  enablePostStepDoIt = false;
    SetProcessSubType(fMuAtomicCapture);
    G4HadronicProcessStore::Instance()->RegisterExtraProcess(this);

    // Theoretical Molecular formation rates versus temperature
    // from Faifman, M.P., Strizh, T.A., Armour, E.A.G. et al.
    // Quadrupole corrections to matrix elements of transitions in resonant
    // reactions of muonic molecule formation. Hyperfine Interact 101, 179–189
    // (1996). https://doi.org/10.1007/BF02227621

    // F=0, DD
    G4int static const numDD0 = 15;

    G4double lambdaDD0temp[numDD0] = {9.719 * kelvin,
                                      61.224 * kelvin,
                                      113.677 * kelvin,
                                      149.956 * kelvin,
                                      203.175 * kelvin,
                                      301.385 * kelvin,
                                      407.915 * kelvin,
                                      491.843 * kelvin,
                                      609.109 * kelvin,
                                      731.772 * kelvin,
                                      859.935 * kelvin,
                                      1041.085 * kelvin,
                                      1205.668 * kelvin,
                                      1370.401 * kelvin,
                                      1487.748 * kelvin};

    G4double lambdaDD0[numDD0] = {0.073e8 / second,
                                  1.614e8 / second,
                                  2.294e8 / second,
                                  2.404e8 / second,
                                  2.386e8 / second,
                                  2.195e8 / second,
                                  2.077e8 / second,
                                  2.143e8 / second,
                                  2.447e8 / second,
                                  2.934e8 / second,
                                  3.514e8 / second,
                                  4.286e8 / second,
                                  4.847e8 / second,
                                  5.271e8 / second,
                                  5.502e8 / second};

    lambdaDD0Interpolator = new G4DataInterpolation(
        lambdaDD0temp,
        lambdaDD0,
        numDD0,
        (lambdaDD0[1] - lambdaDD0[0]) / (lambdaDD0temp[1] - lambdaDD0temp[0]),
        (lambdaDD0[numDD0 - 1] - lambdaDD0[numDD0 - 2])
            / (lambdaDD0temp[numDD0 - 1] - lambdaDD0temp[numDD0 - 2]));

    // F=0, DT
    G4int static const numDT0 = 17;

    G4double lambdaDT0temp[numDT0] = {2.800 * kelvin,
                                      40.599 * kelvin,
                                      90.979 * kelvin,
                                      155.357 * kelvin,
                                      194.476 * kelvin,
                                      233.474 * kelvin,
                                      277.910 * kelvin,
                                      319.385 * kelvin,
                                      392.354 * kelvin,
                                      462.421 * kelvin,
                                      554.546 * kelvin,
                                      669.030 * kelvin,
                                      831.184 * kelvin,
                                      991.254 * kelvin,
                                      1175.567 * kelvin,
                                      1349.074 * kelvin,
                                      1500.362 * kelvin};

    G4double lambdaDT0[numDT0] = {0.000e8 / second,
                                  0.001e8 / second,
                                  0.020e8 / second,
                                  0.039e8 / second,
                                  0.113e8 / second,
                                  0.296e8 / second,
                                  0.627e8 / second,
                                  1.104e8 / second,
                                  2.224e8 / second,
                                  3.435e8 / second,
                                  4.958e8 / second,
                                  6.518e8 / second,
                                  8.015e8 / second,
                                  8.860e8 / second,
                                  9.303e8 / second,
                                  9.388e8 / second,
                                  9.307e8 / second};

    lambdaDT0Interpolator = new G4DataInterpolation(
        lambdaDT0temp,
        lambdaDT0,
        numDT0,
        (lambdaDT0[1] - lambdaDT0[0]) / (lambdaDT0temp[1] - lambdaDT0temp[0]),
        (lambdaDT0[numDT0 - 1] - lambdaDT0[numDT0 - 2])
            / (lambdaDT0temp[numDT0 - 1] - lambdaDT0temp[numDT0 - 2]));

    // F=0, HD
    G4int static const numHD0 = 20;

    G4double lambdaHD0temp[numHD0]
        = {1.400 * kelvin,    65.789 * kelvin,   141.358 * kelvin,
           205.626 * kelvin,  237.603 * kelvin,  273.559 * kelvin,
           308.014 * kelvin,  344.996 * kelvin,  409.061 * kelvin,
           474.303 * kelvin,  555.983 * kelvin,  628.326 * kelvin,
           714.709 * kelvin,  813.915 * kelvin,  926.103 * kelvin,
           1058.264 * kelvin, 1183.899 * kelvin, 1262.237 * kelvin,
           1368.676 * kelvin, 1500.476 * kelvin};

    G4double lambdaHD0[numHD0]
        = {0.000e8 / second,  0.010e8 / second,  0.039e8 / second,
           0.159e8 / second,  0.361e8 / second,  0.765e8 / second,
           1.260e8 / second,  2.003e8 / second,  3.581e8 / second,
           5.361e8 / second,  7.470e8 / second,  9.158e8 / second,
           10.810e8 / second, 12.260e8 / second, 13.361e8 / second,
           14.124e8 / second, 14.456e8 / second, 14.512e8 / second,
           14.476e8 / second, 14.295e8 / second};

    lambdaHD0Interpolator = new G4DataInterpolation(
        lambdaHD0temp,
        lambdaHD0,
        numHD0,
        (lambdaHD0[1] - lambdaHD0[0]) / (lambdaHD0temp[1] - lambdaHD0temp[0]),
        (lambdaHD0[numHD0 - 1] - lambdaHD0[numHD0 - 2])
            / (lambdaHD0temp[numHD0 - 1] - lambdaHD0temp[numHD0 - 2]));

    // F=1, DD
    G4int static const numDD1 = 15;

    G4double lambdaDD1temp[numDD1] = {0.0 * kelvin,
                                      40.97238 * kelvin,
                                      86.18093 * kelvin,
                                      154.00251 * kelvin,
                                      216.16937 * kelvin,
                                      309.43132 * kelvin,
                                      391.40911 * kelvin,
                                      484.73130 * kelvin,
                                      592.21748 * kelvin,
                                      733.66109 * kelvin,
                                      900.53147 * kelvin,
                                      1070.18645 * kelvin,
                                      1251.11005 * kelvin,
                                      1383.94917 * kelvin,
                                      1501.23102 * kelvin};

    G4double lambdaDD1[numDD1] = {0.0e8 / second,
                                  0.02590e8 / second,
                                  0.03788e8 / second,
                                  0.11774e8 / second,
                                  0.17032e8 / second,
                                  0.33172e8 / second,
                                  0.61735e8 / second,
                                  1.20516e8 / second,
                                  2.05375e8 / second,
                                  3.27240e8 / second,
                                  4.47630e8 / second,
                                  5.39123e8 / second,
                                  6.07188e8 / second,
                                  6.38302e8 / second,
                                  6.57099e8 / second};

    lambdaDD1Interpolator = new G4DataInterpolation(
        lambdaDD1temp,
        lambdaDD1,
        numDD1,
        (lambdaDD1[1] - lambdaDD1[0]) / (lambdaDD1temp[1] - lambdaDD1temp[0]),
        (lambdaDD1[numDD1 - 1] - lambdaDD1[numDD1 - 2])
            / (lambdaDD1temp[numDD1 - 1] - lambdaDD1temp[numDD1 - 2]));

    // F=1, DT
    G4int static const numDT1 = 17;

    G4double lambdaDT1temp[numDT1] = {
        0.0 * kelvin,
        42.3405 * kelvin,
        90.3261 * kelvin,
        149.6118 * kelvin,
        193.4004 * kelvin,
        257.0266 * kelvin,
        312.1948 * kelvin,
        377.2730 * kelvin,
        445.1485 * kelvin,
        527.1292 * kelvin,
        630.2549 * kelvin,
        731.9188 * kelvin,
        825.0831 * kelvin,
        935.1662 * kelvin,
        1079.0999 * kelvin,
        1268.1724 * kelvin,
        1502.3866 * kelvin,
    };

    G4double lambdaDT1[numDT1] = {0.0,
                                  0.0254e8 / second,
                                  0.0505e8 / second,
                                  0.2126e8 / second,
                                  0.7603e8 / second,
                                  2.4341e8 / second,
                                  4.2458e8 / second,
                                  6.4969e8 / second,
                                  8.3905e8 / second,
                                  10.1734e8 / second,
                                  11.6117e8 / second,
                                  12.3352e8 / second,
                                  12.6055e8 / second,
                                  12.6414e8 / second,
                                  12.3869e8 / second,
                                  11.8140e8 / second,
                                  10.9640e8 / second};

    lambdaDT1Interpolator = new G4DataInterpolation(
        lambdaDT1temp,
        lambdaDT1,
        numDT1,
        (lambdaDT1[1] - lambdaDT1[0]) / (lambdaDT1temp[1] - lambdaDT1temp[0]),
        (lambdaDT1[numDT1 - 1] - lambdaDT1[numDT1 - 2])
            / (lambdaDT1temp[numDT1 - 1] - lambdaDT0temp[numDT1 - 2]));

    // F=0, HD
    G4int static const numHD1 = 21;

    G4double lambdaHD1temp[numHD1]
        = {0.0 * kelvin,        57.92486 * kelvin,   108.78618 * kelvin,
           151.19847 * kelvin,  182.34991 * kelvin,  197.94409 * kelvin,
           227.75862 * kelvin,  274.69662 * kelvin,  320.26078 * kelvin,
           380.06082 * kelvin,  442.66685 * kelvin,  515.12296 * kelvin,
           595.97369 * kelvin,  690.86987 * kelvin,  795.54230 * kelvin,
           883.16703 * kelvin,  976.39400 * kelvin,  1089.35418 * kelvin,
           1261.56836 * kelvin, 1405.53035 * kelvin, 1505.73924 * kelvin};

    G4double lambdaHD1[numHD1] = {0.0,
                                  0.02523e8 / second,
                                  0.05075e8 / second,
                                  0.26916e8 / second,
                                  0.77688e8 / second,
                                  1.16141e8 / second,
                                  2.20563e8 / second,
                                  4.45963e8 / second,
                                  6.98879e8 / second,
                                  10.28768e8 / second,
                                  13.44890e8 / second,
                                  16.33464e8 / second,
                                  18.64233e8 / second,
                                  20.37175e8 / second,
                                  21.30299e8 / second,
                                  21.56089e8 / second,
                                  21.47469e8 / second,
                                  21.07136e8 / second,
                                  20.10173e8 / second,
                                  19.14696e8 / second,
                                  18.48278e8 / second};

    lambdaHD1Interpolator = new G4DataInterpolation(
        lambdaHD1temp,
        lambdaHD1,
        numHD1,
        (lambdaHD1[1] - lambdaHD1[0]) / (lambdaHD1temp[1] - lambdaHD1temp[0]),
        (lambdaHD1[numHD1 - 1] - lambdaHD1[numHD1 - 2])
            / (lambdaHD1temp[numHD1 - 1] - lambdaHD1temp[numHD1 - 2]));

    if (verboseLevel > 0)
        G4cout << "MuonCatalyzedDTFusion is created." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MuonCatalyzedDTFusion::~MuonCatalyzedDTFusion()
{
    G4HadronicProcessStore::Instance()->DeRegisterExtraProcess(this);
    delete theTotalResult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool MuonCatalyzedDTFusion::IsApplicable(G4ParticleDefinition const& p)
{
    return (&p == G4GenericMuonicAtom::GenericMuonicAtom());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MuonCatalyzedDTFusion::PreparePhysicsTable(G4ParticleDefinition const& p)
{
    G4HadronicProcessStore::Instance()->RegisterParticleForExtraProcess(this,
                                                                        &p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MuonCatalyzedDTFusion::BuildPhysicsTable(G4ParticleDefinition const&)
{
    //  G4HadronicProcessStore::Instance()->PrintInfo(&p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double MuonCatalyzedDTFusion::AtRestGetPhysicalInteractionLength(
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
            G4cout << "MuonCatalyzedDTFusion GPIL: Muonic ion is not "
                      "deuterium or tritium, returning infinity.";
        return std::numeric_limits<G4double>::infinity();
    }

    G4double meanCycleTime = GetMeanCycleTime(track);
    G4double interactionLength = theNumberOfInteractionLengthLeft
                                 * meanCycleTime;

    if (verboseLevel > 1)
        G4cout << "MuonCatalyzedDTFusion GPIL: Interaction Length: "
               << interactionLength / microsecond << " microsecond." << G4endl;

    return (interactionLength);
}

G4double MuonCatalyzedDTFusion::GetMeanCycleTime(G4Track const& track)
{
    // Update density information
    updateMaterialInfo(&track);

    // Process does not happen unless deuterium and tritum are present
    if ((deuteriumPhi < std::numeric_limits<G4double>::epsilon())
        || (tritiumPhi < std::numeric_limits<G4double>::epsilon()))
        return (std::numeric_limits<G4double>::infinity());

    // By default, nothing happens
    G4double meanCycleTime = std::numeric_limits<G4double>::infinity();

    equilibrateHydrogens(&track);  // TODO eventually allow non-equilbrium
                                   // mixtures
    G4DynamicParticle const* p = track.GetDynamicParticle();
    G4int spin = (G4int)round(p->GetSpin());

    G4double lambda_c;
    if (spin == 1)
    {
        lambda_c
            = (DDPhi
               * lambdaDD1Interpolator->CubicSplineInterpolation(temperature))
              + (DTPhi
                 * lambdaDT1Interpolator->CubicSplineInterpolation(temperature))
              + (HDPhi
                 * lambdaHD1Interpolator->CubicSplineInterpolation(temperature));
        // G4cout << "Temperature: " << temperature/kelvin << " kelvin " <<
        // G4endl; G4cout << "Spin 1 DDPhi:" << DDPhi << " DTPhi: " << DTPhi <<
        // " HDPhi: " << HDPhi << G4endl; G4cout << "DD1: " <<
        // lambdaDD1Interpolator->CubicSplineInterpolation(temperature)/(1/microsecond)
        // <<
        //           "DT1: " <<
        //           lambdaDT1Interpolator->CubicSplineInterpolation(temperature)/(1/microsecond)
        //           << "HD1: " <<
        //           lambdaHD1Interpolator->CubicSplineInterpolation(temperature)/(1/microsecond)
        //           << G4endl;
        // G4cout << "lambda_c " << lambda_c / (1/microsecond) << G4endl;
    }
    else  // spin == 0
    {
        lambda_c
            = (DDPhi
               * lambdaDD0Interpolator->CubicSplineInterpolation(temperature))
              + (DTPhi
                 * lambdaDT0Interpolator->CubicSplineInterpolation(temperature))
              + (HDPhi
                 * lambdaHD0Interpolator->CubicSplineInterpolation(temperature));
        //    G4cout << "Temperature: " << temperature/kelvin << " kelvin " <<
        //    G4endl;

        //          G4cout << "Spin 0 DDPhi:" << DDPhi << " DTPhi: " << DTPhi
        //          << " HDPhi: " << HDPhi << G4endl;
        //   G4cout << "DD0: " <<
        //   lambdaDD0Interpolator->CubicSplineInterpolation(temperature)/(1/microsecond)
        //   <<
        //           "DT0: " <<
        //           lambdaDT0Interpolator->CubicSplineInterpolation(temperature)/(1/microsecond)
        //           <<
        //         "HD0: " <<
        //         lambdaHD0Interpolator->CubicSplineInterpolation(temperature)/(1/microsecond)
        //         << G4endl;
        //   G4cout << "lambda_c " << lambda_c / (1/microsecond) << G4endl;
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
        G4cout << "MuonCatalyzedDTFusion GetMeanCycleTime: spin: " << spin
               << G4endl;
        G4cout << "MuonCatalyzedDTFusion GetMeanCycleTime: Temperature: "
               << temperature << " DDPhi: " << DDPhi << " DTPhi: " << DTPhi
               << " HDPhi: " << HDPhi << G4endl;
        G4cout << "MuonCatalyzedDTFusion GetMeanCycleTime: Mean Cycle Time: "
               << meanCycleTime / microsecond << " microsecond" << G4endl;
        G4cout << "MuonCatalyzedDTFusion GetMeanCycleTime: Global time: "
               << track.GetGlobalTime() / microsecond << " microsecond "
               << G4endl;
    }

    return (meanCycleTime);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange*
MuonCatalyzedDTFusion::AtRestDoIt(G4Track const& track, G4Step const& step)
{
    theTotalResult->Initialize(track);

    // CELER_LOG_LOCAL(info) << "DT mucf: "
    //                       <<
    //                       track.GetParticleDefinition()->GetParticleName();
    // auto const post_step = step.GetPostStepPoint();
    // CELER_LOG_LOCAL(info) << "DT mucf: " << post_step->GetPosition().x() /
    // cm
    //                       << " cm and " << post_step->GetKineticEnergy()
    //                       << " MeV";

    G4double finalGlobalTime = track.GetGlobalTime();
    G4double finalLocalTime = track.GetLocalTime();

    G4double meanCycleTime = GetMeanCycleTime(track);
    G4double interactionTime = theNumberOfInteractionLengthLeft * meanCycleTime;
    finalGlobalTime += interactionTime;
    finalLocalTime += interactionTime;

    // Calculate kinetic energy of resultant muon by interpolating on CDF of
    // distribution function
    G4double muonKineticEnergy
        = muonEnergyInterpolator->CubicSplineInterpolation(G4UniformRand());

    // Physics parameters
    // (initial sticking from https://arxiv.org/pdf/2112.08399.pdf)
    G4double initialStickingFraction = 0.00857;
    G4double neutronKineticEnergy = 14.1 * MeV;

    // In both cases, there is a 14.1 MeV Neutron, random direction
    G4DynamicParticle* neutron = new G4DynamicParticle(
        G4Neutron::Neutron(), G4RandomDirection(), neutronKineticEnergy);
    G4ThreeVector neutronMomentum = neutron->GetMomentum();

    // Do random sampling to determine if there is alpha sticking
    G4double randomNumber = G4UniformRand();
    if (randomNumber <= initialStickingFraction)
    {
        if (verboseLevel > 1)
            G4cout << "MuonCatalyzedFusion: Alpha sticking occured." << G4endl;
        // case with alpha sticking
        // Muonic alpha with equal and opposite momentum to neutron
        theTotalResult->SetNumberOfSecondaries(2);
        G4IonTable* itp = G4IonTable::GetIonTable();
        G4ParticleDefinition* muonicAlphaDefinition = itp->GetMuonicAtom(2, 4);
        G4DynamicParticle* muonicAlpha
            = new G4DynamicParticle(muonicAlphaDefinition, -neutronMomentum);
        muonicAlpha->SetCharge(+1 * eplus);
        theTotalResult->AddSecondary(muonicAlpha, finalGlobalTime, true);
    }
    else
    {
        // case without alpha sticking
        // Muon in random direction
        if (verboseLevel > 1)
            G4cout << "MuonCatalyzedFusion: Alpha sticking did not occur."
                   << G4endl;
        theTotalResult->SetNumberOfSecondaries(3);
        G4DynamicParticle* muon = new G4DynamicParticle(
            G4MuonMinus::MuonMinus(), G4RandomDirection(), muonKineticEnergy);
        theTotalResult->AddSecondary(muon, finalGlobalTime, true);
        G4ThreeVector muonMomentum = muon->GetMomentum();
        // Alpha particle with equal and opposite momentum to neutron plus muon
        G4DynamicParticle* alpha = new G4DynamicParticle(
            G4Alpha::Alpha(), -(neutronMomentum + muonMomentum));
        theTotalResult->AddSecondary(alpha, finalGlobalTime, true);
    }
    // Add the neutron to the list of secondaries
    theTotalResult->AddSecondary(neutron, finalGlobalTime, true);

    // Kill primary particle
    //
    theTotalResult->ProposeTrackStatus(fStopAndKill);
    theTotalResult->ProposeLocalTime(finalLocalTime);

    ClearNumberOfInteractionLengthLeft();
    return theTotalResult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MuonCatalyzedDTFusion::ProcessDescription(std::ostream& outFile) const
{
    outFile << "Model of muon catalyzed DT fusion.";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
