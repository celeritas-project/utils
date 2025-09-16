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
// Authors: Sridhar Tripathy
//          Ara Knaian (ara@nklabs.com)
//          Kevin Lynch (krlynch@fnal.gov)
//          Andrzej Adamczak (andrzej.adamczak@ifj.edu.pl)
//
// Class Description:
//
// Muon Transfer between Muonic Atoms

#include "MuonicAtomTransfer.hh"

#include <iostream>
#include <G4EmCaptureCascade.hh>
#include <G4GenericMuonicAtom.hh>
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

#include "VMuonCatalyzedFusionProcess.hh"
using namespace std;

MuonicAtomTransfer::MuonicAtomTransfer(G4String const& name)
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

    // The exchange rates were computed by A. Adamczak
    // using the method described in "A. Adamczak, PRA 74, 042718 (2006)"

    G4int static const numHTemp = 24;
    minHTemp = 5. * kelvin;
    maxHTemp = 1500. * kelvin;

    G4double HTemp[numHTemp]
        = {5.0 * kelvin,    7.0 * kelvin,    10.0 * kelvin,   15.0 * kelvin,
           20.0 * kelvin,   30.0 * kelvin,   50.0 * kelvin,   70.0 * kelvin,
           100.0 * kelvin,  150.0 * kelvin,  200.0 * kelvin,  300.0 * kelvin,
           400.0 * kelvin,  500.0 * kelvin,  600.0 * kelvin,  700.0 * kelvin,
           800.0 * kelvin,  900.0 * kelvin,  1000.0 * kelvin, 1100.0 * kelvin,
           1200.0 * kelvin, 1300.0 * kelvin, 1400.0 * kelvin, 1500.0 * kelvin};

    {
        G4double table[numHTemp]
            = {2.050700e+10 / second, 2.038600e+10 / second,
               2.022200e+10 / second, 1.999100e+10 / second,
               1.979600e+10 / second, 1.948400e+10 / second,
               1.904100e+10 / second, 1.873300e+10 / second,
               1.840600e+10 / second, 1.805100e+10 / second,
               1.781800e+10 / second, 1.752100e+10 / second,
               1.733500e+10 / second, 1.720400e+10 / second,
               1.710400e+10 / second, 1.702500e+10 / second,
               1.696000e+10 / second, 1.690500e+10 / second,
               1.685800e+10 / second, 1.681800e+10 / second,
               1.678200e+10 / second, 1.675100e+10 / second,
               1.672300e+10 / second, 1.669900e+10 / second};
        addExchangeRate(
            PROTIUM, DEUTERIUM, DEUTERIUM, DEUTERIUM, HTemp, table, numHTemp);
    }

    {
        G4double table[numHTemp]
            = {1.016100e+10 / second, 1.010600e+10 / second,
               1.003200e+10 / second, 9.927600e+09 / second,
               9.839400e+09 / second, 9.697100e+09 / second,
               9.493100e+09 / second, 9.349400e+09 / second,
               9.195500e+09 / second, 9.025900e+09 / second,
               8.912800e+09 / second, 8.766900e+09 / second,
               8.673600e+09 / second, 8.606900e+09 / second,
               8.555900e+09 / second, 8.515100e+09 / second,
               8.481400e+09 / second, 8.453100e+09 / second,
               8.428800e+09 / second, 8.407800e+09 / second,
               8.389500e+09 / second, 8.373500e+09 / second,
               8.359300e+09 / second, 8.346800e+09 / second};
        addExchangeRate(
            PROTIUM, PROTIUM, DEUTERIUM, DEUTERIUM, HTemp, table, numHTemp);
    }

    {
        G4double table[numHTemp]
            = {4.357500e+09 / second, 4.332600e+09 / second,
               4.298900e+09 / second, 4.250700e+09 / second,
               4.210500e+09 / second, 4.146600e+09 / second,
               4.057400e+09 / second, 3.996100e+09 / second,
               3.931900e+09 / second, 3.863300e+09 / second,
               3.819200e+09 / second, 3.765300e+09 / second,
               3.733800e+09 / second, 3.713400e+09 / second,
               3.699700e+09 / second, 3.690200e+09 / second,
               3.683700e+09 / second, 3.679400e+09 / second,
               3.676800e+09 / second, 3.675600e+09 / second,
               3.675500e+09 / second, 3.676300e+09 / second,
               3.677800e+09 / second, 3.680000e+09 / second};
        addExchangeRate(
            PROTIUM, PROTIUM, TRITIUM, TRITIUM, HTemp, table, numHTemp);
    }

    {
        G4double table[numHTemp]
            = {8.822800e+09 / second, 8.766500e+09 / second,
               8.690500e+09 / second, 8.582400e+09 / second,
               8.492200e+09 / second, 8.348900e+09 / second,
               8.149400e+09 / second, 8.013400e+09 / second,
               7.872400e+09 / second, 7.723700e+09 / second,
               7.629900e+09 / second, 7.517800e+09 / second,
               7.453700e+09 / second, 7.413100e+09 / second,
               7.386100e+09 / second, 7.367500e+09 / second,
               7.354900e+09 / second, 7.346400e+09 / second,
               7.341200e+09 / second, 7.338400e+09 / second,
               7.337800e+09 / second, 7.338800e+09 / second,
               7.341300e+09 / second, 7.345000e+09 / second};
        addExchangeRate(
            PROTIUM, TRITIUM, TRITIUM, TRITIUM, HTemp, table, numHTemp);
    }

    {
        G4double table[numHTemp] = {
            0.15034E+09 / second, 0.14850E+09 / second, 0.14617E+09 / second,
            0.14312E+09 / second, 0.14077E+09 / second, 0.13740E+09 / second,
            0.13347E+09 / second, 0.13137E+09 / second, 0.12984E+09 / second,
            0.12931E+09 / second, 0.12995E+09 / second, 0.13271E+09 / second,
            0.13636E+09 / second, 0.14042E+09 / second, 0.14470E+09 / second,
            0.14913E+09 / second, 0.15366E+09 / second, 0.15828E+09 / second,
            0.16296E+09 / second, 0.16770E+09 / second, 0.17251E+09 / second,
            0.17736E+09 / second, 0.18227E+09 / second, 0.18722E+09 / second};
        addExchangeRate(
            DEUTERIUM, DEUTERIUM, TRITIUM, TRITIUM, HTemp, table, numHTemp);
    }

    {
        G4double table[numHTemp] = {
            0.30520E+09 / second, 0.30117E+09 / second, 0.29607E+09 / second,
            0.28941E+09 / second, 0.28433E+09 / second, 0.27703E+09 / second,
            0.26846E+09 / second, 0.26385E+09 / second, 0.26044E+09 / second,
            0.25907E+09 / second, 0.26013E+09 / second, 0.26531E+09 / second,
            0.27226E+09 / second, 0.27998E+09 / second, 0.28812E+09 / second,
            0.29650E+09 / second, 0.30508E+09 / second, 0.31380E+09 / second,
            0.32265E+09 / second, 0.33161E+09 / second, 0.34067E+09 / second,
            0.34984E+09 / second, 0.35909E+09 / second, 0.36843E+09 / second};
        addExchangeRate(
            DEUTERIUM, TRITIUM, TRITIUM, TRITIUM, HTemp, table, numHTemp);
    }

    {
        G4double table[numHTemp]
            = {1.031500e+10 / second, 1.024900e+10 / second,
               1.016200e+10 / second, 1.003800e+10 / second,
               9.934300e+09 / second, 9.769500e+09 / second,
               9.538500e+09 / second, 9.379600e+09 / second,
               9.212700e+09 / second, 9.032900e+09 / second,
               8.915400e+09 / second, 8.766800e+09 / second,
               8.673500e+09 / second, 8.607700e+09 / second,
               8.558000e+09 / second, 8.518400e+09 / second,
               8.486000e+09 / second, 8.458600e+09 / second,
               8.435200e+09 / second, 8.414900e+09 / second,
               8.397000e+09 / second, 8.381300e+09 / second,
               8.367300e+09 / second, 8.354800e+09 / second};
        addExchangeRate(
            PROTIUM, DEUTERIUM, TRITIUM, DEUTERIUM, HTemp, table, numHTemp);
    }

    {
        G4double table[numHTemp]
            = {4.390700e+09 / second, 4.364000e+09 / second,
               4.327900e+09 / second, 4.276300e+09 / second,
               4.233200e+09 / second, 4.164500e+09 / second,
               4.068500e+09 / second, 4.002700e+09 / second,
               3.934100e+09 / second, 3.861400e+09 / second,
               3.815300e+09 / second, 3.759700e+09 / second,
               3.727800e+09 / second, 3.707500e+09 / second,
               3.694000e+09 / second, 3.684700e+09 / second,
               3.678400e+09 / second, 3.674200e+09 / second,
               3.671700e+09 / second, 3.670500e+09 / second,
               3.670300e+09 / second, 3.671000e+09 / second,
               3.672400e+09 / second, 3.674500e+09 / second};
        addExchangeRate(
            PROTIUM, DEUTERIUM, TRITIUM, TRITIUM, HTemp, table, numHTemp);
    }

    // Data from Table 3
    // A. V. Kravtsov; A. I. Mikhailov (2000). Temperature dependence of the
    // formation rates of hydrogen-helium mesic molecules in collisions of slow
    // hydrogen atoms with helium. , 90(1), 45–49. doi:10.1134/1.559092

    G4int static const numHeTemp = 16;
    minHeTemp = 15 * kelvin;
    maxHeTemp = 500 * kelvin;

    G4double heTemp[numHeTemp] = {
        15 * kelvin,
        20 * kelvin,
        25 * kelvin,
        30 * kelvin,
        35 * kelvin,
        40 * kelvin,
        50 * kelvin,
        100 * kelvin,
        150 * kelvin,
        200 * kelvin,
        250 * kelvin,
        300 * kelvin,
        350 * kelvin,
        400 * kelvin,
        450 * kelvin,
        500 * kelvin,
    };

    G4double lambdapHe3[numHeTemp] = {0.52e8 / second,
                                      0.52e8 / second,
                                      0.51e8 / second,
                                      0.51e8 / second,
                                      0.51e8 / second,
                                      0.50e8 / second,
                                      0.50e8 / second,
                                      0.47e8 / second,
                                      0.46e8 / second,
                                      0.45e8 / second,
                                      0.44e8 / second,
                                      0.43e8 / second,
                                      0.42e8 / second,
                                      0.42e8 / second,
                                      0.41e8 / second,
                                      0.40e8 / second};

    lambdapHe3Interpolator = new G4DataInterpolation(
        heTemp,
        lambdapHe3,
        numHeTemp,
        (lambdapHe3[1] - lambdapHe3[0]) / (heTemp[1] - heTemp[0]),
        (lambdapHe3[numHeTemp - 1] - lambdapHe3[numHeTemp - 2])
            / (heTemp[numHeTemp - 1] - heTemp[numHeTemp - 2]));

    G4double lambdapHe4[numHeTemp] = {0.33e8 / second,
                                      0.33e8 / second,
                                      0.32e8 / second,
                                      0.32e8 / second,
                                      0.32e8 / second,
                                      0.31e8 / second,
                                      0.31e8 / second,
                                      0.29e8 / second,
                                      0.28e8 / second,
                                      0.27e8 / second,
                                      0.27e8 / second,
                                      0.26e8 / second,
                                      0.25e8 / second,
                                      0.25e8 / second,
                                      0.25e8 / second,
                                      0.24e8 / second};

    lambdapHe4Interpolator = new G4DataInterpolation(
        heTemp,
        lambdapHe4,
        numHeTemp,
        (lambdapHe4[1] - lambdapHe4[0]) / (heTemp[1] - heTemp[0]),
        (lambdapHe4[numHeTemp - 1] - lambdapHe4[numHeTemp - 2])
            / (heTemp[numHeTemp - 1] - heTemp[numHeTemp - 2]));

    G4double lambdadHe3[numHeTemp] = {2.40e8 / second,
                                      2.34e8 / second,
                                      2.28e8 / second,
                                      2.24e8 / second,
                                      2.20e8 / second,
                                      2.16e8 / second,
                                      2.09e8 / second,
                                      1.85e8 / second,
                                      1.70e8 / second,
                                      1.60e8 / second,
                                      1.52e8 / second,
                                      1.45e8 / second,
                                      1.39e8 / second,
                                      1.34e8 / second,
                                      1.30e8 / second,
                                      1.26e8 / second};

    lambdadHe3Interpolator = new G4DataInterpolation(
        heTemp,
        lambdadHe3,
        numHeTemp,
        (lambdadHe3[1] - lambdadHe3[0]) / (heTemp[1] - heTemp[0]),
        (lambdadHe3[numHeTemp - 1] - lambdadHe3[numHeTemp - 2])
            / (heTemp[numHeTemp - 1] - heTemp[numHeTemp - 2]));

    G4double lambdadHe4[numHeTemp] = {12.7e8 / second,
                                      11.8e8 / second,
                                      11.0e8 / second,
                                      10.4e8 / second,
                                      9.8e8 / second,
                                      9.3e8 / second,
                                      8.5e8 / second,
                                      6.2e8 / second,
                                      5.1e8 / second,
                                      4.4e8 / second,
                                      3.8e8 / second,
                                      3.5e8 / second,
                                      3.2e8 / second,
                                      2.9e8 / second,
                                      2.7e8 / second,
                                      2.6e8 / second};

    lambdadHe4Interpolator = new G4DataInterpolation(
        heTemp,
        lambdadHe4,
        numHeTemp,
        (lambdadHe4[1] - lambdadHe4[0]) / (heTemp[1] - heTemp[0]),
        (lambdadHe4[numHeTemp - 1] - lambdadHe4[numHeTemp - 2])
            / (heTemp[numHeTemp - 1] - heTemp[numHeTemp - 2]));

    G4double lambdatHe3[numHeTemp] = {51.2e8 / second,
                                      45.5e8 / second,
                                      41.1e8 / second,
                                      37.6e8 / second,
                                      34.7e8 / second,
                                      32.3e8 / second,
                                      28.6e8 / second,
                                      18.8e8 / second,
                                      14.4e8 / second,
                                      11.9e8 / second,
                                      10.2e8 / second,
                                      9.0e8 / second,
                                      8.1e8 / second,
                                      7.3e8 / second,
                                      6.7e8 / second,
                                      6.2e8 / second};

    lambdatHe3Interpolator = new G4DataInterpolation(
        heTemp,
        lambdatHe3,
        numHeTemp,
        (lambdatHe3[1] - lambdatHe3[0]) / (heTemp[1] - heTemp[0]),
        (lambdatHe3[numHeTemp - 1] - lambdatHe3[numHeTemp - 2])
            / (heTemp[numHeTemp - 1] - heTemp[numHeTemp - 2]));

    G4double lambdatHe4[numHeTemp] = {1.89e8 / second,
                                      1.86e8 / second,
                                      1.84e8 / second,
                                      1.82e8 / second,
                                      1.79e8 / second,
                                      1.77e8 / second,
                                      1.74e8 / second,
                                      1.60e8 / second,
                                      1.51e8 / second,
                                      1.43e8 / second,
                                      1.37e8 / second,
                                      1.31e8 / second,
                                      1.27e8 / second,
                                      1.22e8 / second,
                                      1.19e8 / second,
                                      1.15e8 / second};

    lambdatHe4Interpolator = new G4DataInterpolation(
        heTemp,
        lambdatHe4,
        numHeTemp,
        (lambdatHe4[1] - lambdatHe4[0]) / (heTemp[1] - heTemp[0]),
        (lambdatHe4[numHeTemp - 1] - lambdatHe4[numHeTemp - 2])
            / (heTemp[numHeTemp - 1] - heTemp[numHeTemp - 2]));
}

MuonicAtomTransfer::~MuonicAtomTransfer()
{
    G4HadronicProcessStore::Instance()->DeRegisterExtraProcess(this);
    delete theTotalResult;
}

G4bool MuonicAtomTransfer::IsApplicable(G4ParticleDefinition const& p)
{
    return (&p == G4GenericMuonicAtom::GenericMuonicAtom());
}

void MuonicAtomTransfer::PreparePhysicsTable(G4ParticleDefinition const& p)
{
    G4HadronicProcessStore::Instance()->RegisterParticleForExtraProcess(this,
                                                                        &p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MuonicAtomTransfer::BuildPhysicsTable(G4ParticleDefinition const& p)
{
    G4HadronicProcessStore::Instance()->PrintInfo(&p);
}

G4double MuonicAtomTransfer::AtRestGetPhysicalInteractionLength(
    G4Track const& track, G4ForceCondition* condition)
{
    *condition = NotForced;

    if (theNumberOfInteractionLengthLeft < 0.)
    {
        ResetNumberOfInteractionLengthLeft();
    }

    G4double interactionLength
        = theNumberOfInteractionLengthLeft
          * GetMeanTransferTime(track, false, nullptr, nullptr);
    return (interactionLength);
}

// If selectElement is false, just returns the mean transfer time  (use in get
// physical interaction length) if selectElement is true, does two passes and
// tells you to which element the transfer occured  (use in post-step do it)

// TODO can remove this class from VMuonCatatalyzedFusionProcess now if you
// like (but don't have to)

G4double MuonicAtomTransfer::GetMeanTransferTime(G4Track const& track,
                                                 G4bool selectElement,
                                                 G4int* newZ,
                                                 G4int* newA)
{
    // Default is nothing happens
    G4double meanTransferTime = std::numeric_limits<G4double>::infinity();

    // Check particle
    G4DynamicParticle const* particle = track.GetDynamicParticle();
    G4ParticleDefinition const* particleDefinition = particle->GetDefinition();
    G4int startingA = particleDefinition->GetAtomicMass();
    G4int startingZ = particleDefinition->GetAtomicNumber();

    // Update info on hydrogen isoprotologue concentrations
    updateMaterialInfo(&track);
    equilibrateHydrogens(&track);  // TODO eventually allow non-equilbrium
                                   // mixtures

    // Get material information
    G4Material* mat = track.GetMaterial();
    G4ElementVector const* elementVector = mat->GetElementVector();  // std
                                                                     // vector
                                                                     // of g4
                                                                     // element
    G4double const* fractionVector = mat->GetFractionVector();
    G4double density = mat->GetDensity();
    temperature = mat->GetTemperature();

    G4double lambda = 0.0;
    G4double totalLambda = 0.0;
    G4double thresholdLambda = 0.0;

    // Do two passes if selectElement is true, so we can return a specific
    // target
    G4int passes;
    if (selectElement)
        passes = 2;
    else
        passes = 1;
    for (G4int pass = 1; pass <= passes; pass++)
    {
        // If the second pass, randomly select a specific target isotope in
        // propotion to their relative rates
        if (pass == 2)
        {
            totalLambda = lambda;
            thresholdLambda = G4UniformRand() * totalLambda;
            lambda = 0.0;
        }

        // Handle transfers between hydrogen isotopes seperately,
        // because we use hydrogen isoprotologue concentrations
        // for these to take into account molecular effects

        // If this is a muonic hydrogen
        if (startingZ == 1)
        {
            G4double HTemp = temperature;
            if (HTemp < minHTemp)
                HTemp = minHTemp;
            if (HTemp > maxHTemp)
                HTemp = maxHTemp;

            // Iterate over all hydrogen exchange rates
            for (auto& rate : rateTable)
            {
                if (rate.startIsotope == (startingA - 1))
                {
                    lambda += moleculeDensity[rate.isotope1][rate.isotope2]
                              * rate.interpolator->CubicSplineInterpolation(
                                  HTemp);

                    // If the second pass
                    if (pass == 2)
                    {
                        if (lambda > thresholdLambda)
                        {
                            // If this is the selected transfer, record the
                            // final state and exit
                            *newZ = 1;
                            *newA = rate.endIsotope + 1;
                            return (1 / totalLambda);
                        }
                    }
                }
            }
        }

        // Go through each constituent isotope of the material and calculate
        // their contribution to the muon transfer rate
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

                G4double thisMassFraction = fractionVector[i] * thisA
                                            * relativeAbundance[j] / atomicMass;
                G4double thisPhi
                    = thisMassFraction * density / (thisA / Avogadro)
                      / liquidHydrogenAtomicNumberDensity;  // do the units
                                                            // work here to use
                                                            // thisA?  do you
                                                            // need an explicit
                                                            // amu?
                // G4cout << "pass " << pass << " atomicMass:" <<
                // atomicMass/(gram/mole) << "g/mol relativeAbundance:" <<
                // relativeAbundance[j] << " fractionVector:"
                // <<fractionVector[i] << "thisZ:" << thisZ << " thisA:" <<
                // thisA/(gram/mole) << "g/mol thisN:" << thisN << "
                // thisMassFraction:" << thisMassFraction << " thisPhi:" <<
                // thisPhi << G4endl;

                // We handle hydrogen, helium, and Z>3 seperately
                if (thisZ == 1)
                {
                    // Note: transfer to hydrogen is handled above as a special
                    // case
                }
                else if (thisZ == 2)
                {
                    // ************************ TRANSFER TO HELIUM
                    // *************************************

                    G4double heTemp = temperature;
                    if (heTemp < minHeTemp)
                        heTemp = minHeTemp;
                    if (heTemp > maxHeTemp)
                        heTemp = maxHeTemp;

                    if (thisN == 3)
                    {
                        // Helium-3
                        if ((startingZ == 1) && (startingA == 1))
                            lambda += thisPhi
                                      * lambdapHe3Interpolator
                                            ->CubicSplineInterpolation(heTemp);
                        else if ((startingZ == 1) && (startingA == 2))
                            lambda += thisPhi
                                      * lambdadHe3Interpolator
                                            ->CubicSplineInterpolation(heTemp);
                        else if ((startingZ == 1) && (startingA == 3))
                            lambda += thisPhi
                                      * lambdatHe3Interpolator
                                            ->CubicSplineInterpolation(heTemp);
                        // TODO add in transfer rate from Helium-3 to Helium-4
                        // (is this known?????)
                    }
                    else if (thisN == 4)
                    {
                        // Helium-4
                        if ((startingZ == 1) && (startingA == 1))
                            lambda += thisPhi
                                      * lambdapHe4Interpolator
                                            ->CubicSplineInterpolation(heTemp);
                        else if ((startingZ == 1) && (startingA == 2))
                            lambda += thisPhi
                                      * lambdadHe4Interpolator
                                            ->CubicSplineInterpolation(heTemp);
                        else if ((startingZ == 1) && (startingA == 3))
                            lambda += thisPhi
                                      * lambdatHe4Interpolator
                                            ->CubicSplineInterpolation(heTemp);
                    }
                }
                else
                {
                    // ************************ TRANSFER TO HEAVIER ELEMENTS
                    // **************************** generic transfer to Z > 3
                    // isotope uses very approximate relationship
                    // TODO may want to add special cases for transfer from
                    // hydrogens to specific elements (e.g. nitrogen, oxygen,
                    // carbon, neon, argon) for which there is better data
                    // TODO if you take some data with different materials at
                    // different temperatures at the beam you could validate
                    // this model

                    if (thisZ > startingZ)
                    {
                        lambda += thisPhi * (thisZ - startingZ)
                                  * (temperature / (300 * kelvin))
                                  * (4.0e11 / second / 36.0);
                        // G4cout << "-------- thisZ: " << thisZ << "
                        // startingZ: " << startingZ << " lambda:" <<
                        // lambda/(1/second) << "/second thisPhi" << thisPhi <<
                        // G4endl;
                        //  Approximate relationship from plot in Schellenberg,
                        //  Muon Cat. Fusion 6
                        //  TODO Temp dependence from hydrogen to oxygen looks
                        //  linear
                        //  https://www.sciencedirect.com/science/article/pii/S037596012030534X
                        //  is this justified in general, or is this specific
                        //  to oxygen?
                        //  TODO is there a more correct / more justifed
                        //  expression to use here?
                    }
                }

                // For each element/isotope, check if this is the one that was
                // randomly selected If so, return it for the postStepDoIt
                if (pass == 2)
                {
                    if (lambda > thresholdLambda)
                    {
                        *newZ = thisZ;
                        *newA = thisN;
                        // G4cout << "********* pass 2 returned value " <<
                        // *newZ << " " << *newA << " " << "totalLambda " <<
                        // totalLambda/(1/second) << "/second" << G4endl;
                        return (1 / totalLambda);
                    }
                }
            }
        }
    }

    if (lambda > std::numeric_limits<G4double>::epsilon())  // if rate is zero,
                                                            // keep cycle time
                                                            // at infinity so
                                                            // process does not
                                                            // trigger
        meanTransferTime = 1 / lambda;
    ;

    // G4cout << "at end, returned " << meanTransferTime/microsecond << "
    // microsecond" << G4endl;
    return (meanTransferTime);
}

G4VParticleChange*
MuonicAtomTransfer::AtRestDoIt(G4Track const& track, G4Step const&)
{
    theTotalResult->Initialize(track);

    // G4cout << "((((((((((((At start of atRestDoIt))))))))))))" << G4endl;

    G4double finalGlobalTime = track.GetGlobalTime();
    G4double finalLocalTime = track.GetLocalTime();
    G4int newZ;
    G4int newA;
    G4double interactionTime = theNumberOfInteractionLengthLeft
                               * GetMeanTransferTime(track, true, &newZ, &newA);
    finalGlobalTime += interactionTime;
    finalLocalTime += interactionTime;

    theTotalResult->SetNumberOfSecondaries(1);

    // Create the new muonic atom
    G4IonTable* itp = G4IonTable::GetIonTable();
    G4ParticleDefinition* muonicAtom = itp->GetMuonicAtom(newZ, newA);
    G4DynamicParticle* dp
        = new G4DynamicParticle(muonicAtom, G4RandomDirection(), 0.);

    // TODO is the initial spin random with the same distribution as capture or
    // is it dependent on the previous spin?

    // If hydrogen, set the spin
    // Yamashita, T., Kino, Y., Okutsu, K. et al. Roles of resonant muonic
    // molecule in new kinetics model and muon catalyzed fusion in compressed
    // gas. Sci Rep 12, 6393 (2022). https://doi.org/10.1038/s41598-022-09487-0

    if (newZ == 1)
    {
        G4double upFraction;
        G4int spinUp;
        G4int spinDown;

        if (newA == 2)
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

    // Replace the old particle with the new one
    theTotalResult->AddSecondary(dp, finalGlobalTime, true);
    theTotalResult->ProposeLocalEnergyDeposit(0.0);
    theTotalResult->ProposeTrackStatus(fStopAndKill);
    theTotalResult->ProposeLocalTime(finalLocalTime);
    ClearNumberOfInteractionLengthLeft();

    // G4cout << "New finalGlobalTime: " << finalGlobalTime/nanosecond << " ns
    // " << G4endl;

    return theTotalResult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MuonicAtomTransfer::ProcessDescription(std::ostream& outFile) const
{
    outFile << "Transfer" << endl;
}

// Register a transfer rate table
void MuonicAtomTransfer::addExchangeRate(int startIsotope,
                                         int isotope1,
                                         int isotope2,
                                         int endIsotope,
                                         G4double* tableTemperature,
                                         G4double* table,
                                         int numTemperatures)
{
    exchangeRate rate;
    rate.interpolator = new G4DataInterpolation(
        tableTemperature,
        table,
        numTemperatures,
        (table[1] - table[0]) / (tableTemperature[1] - tableTemperature[0]),
        (table[numTemperatures - 1] - table[numTemperatures - 2])
            / (tableTemperature[numTemperatures - 1]
               - tableTemperature[numTemperatures - 2]));
    rate.startIsotope = startIsotope;
    rate.isotope1 = isotope1;
    rate.isotope2 = isotope2;
    rate.endIsotope = endIsotope;

    rateTable.push_back(rate);
}
