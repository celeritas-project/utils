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
//          Andrzej Adamczak (andrzej.adamczak@ifj.edu.pl)
//
// Class Description:
//
// Muonic atom spin-flip process for MuH1, MuH2, and MuH3
//

#include "MuonicAtomSpinFlip.hh"

#include <G4DataInterpolation.hh>
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
#include <G4PhysicalConstants.hh>
#include <G4RandomDirection.hh>
#include <G4SystemOfUnits.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MuonicAtomSpinFlip::MuonicAtomSpinFlip(G4String const& name)
    : VMuonCatalyzedFusionProcess(name, fHadronic)
    ,
    //  : G4HadronicProcess(name, fHadronAtRest),// name, process type
    fElementSelector(new G4ElementSelector())
    , fEmCascade(new G4EmCaptureCascade())
    ,  // Owned by InteractionRegistry
    theTotalResult(new G4ParticleChange())
{
    // TODO is this important / correct?
    SetProcessSubType(fMuAtomicCapture);
    G4HadronicProcessStore::Instance()->RegisterExtraProcess(this);

    // The spin flip rates were computed by A. Adamczak
    // using the method described in "A. Adamczak, PRA 74, 042718 (2006)"

    // Spin flip only occurs in collision with another molecule containing
    // the same istope as one of its constitutents.  The spin flip rates
    // are indexed by muonic atom isotope, the other isotope on the molecule,
    // and the flip direction (up or down)

    G4int static const numTemperatures = 24;

    G4double tableTemperature[numTemperatures]
        = {5.0 * kelvin,    7.0 * kelvin,    10.0 * kelvin,   15.0 * kelvin,
           20.0 * kelvin,   30.0 * kelvin,   50.0 * kelvin,   70.0 * kelvin,
           100.0 * kelvin,  150.0 * kelvin,  200.0 * kelvin,  300.0 * kelvin,
           400.0 * kelvin,  500.0 * kelvin,  600.0 * kelvin,  700.0 * kelvin,
           800.0 * kelvin,  900.0 * kelvin,  1000.0 * kelvin, 1100.0 * kelvin,
           1200.0 * kelvin, 1300.0 * kelvin, 1400.0 * kelvin, 1500.0 * kelvin};

    {
        G4double table[numTemperatures]
            = {0.000000e+00 / second, 0.000000e+00 / second,
               0.000000e+00 / second, 0.000000e+00 / second,
               0.000000e+00 / second, 0.000000e+00 / second,
               9.926000e+01 / second, 3.754200e+03 / second,
               5.492500e+04 / second, 4.691300e+05 / second,
               1.436400e+06 / second, 4.812900e+06 / second,
               9.430300e+06 / second, 1.463500e+07 / second,
               2.047700e+07 / second, 2.643900e+07 / second,
               3.308700e+07 / second, 3.971900e+07 / second,
               4.657800e+07 / second, 5.365500e+07 / second,
               6.094100e+07 / second, 6.842900e+07 / second,
               7.610900e+07 / second, 8.398000e+07 / second};
        addInterpolator(
            DEUTERIUM, PROTIUM, UP, tableTemperature, table, numTemperatures);
    }

    {
        G4double table[numTemperatures]
            = {2.765100e+07 / second, 2.681600e+07 / second,
               2.574400e+07 / second, 2.456900e+07 / second,
               2.377600e+07 / second, 2.310500e+07 / second,
               2.249400e+07 / second, 2.242300e+07 / second,
               2.254100e+07 / second, 2.319600e+07 / second,
               2.398500e+07 / second, 2.614500e+07 / second,
               2.867800e+07 / second, 3.122100e+07 / second,
               3.433500e+07 / second, 3.729700e+07 / second,
               4.094000e+07 / second, 4.439600e+07 / second,
               4.799100e+07 / second, 5.171600e+07 / second,
               5.556500e+07 / second, 5.952800e+07 / second,
               6.360000e+07 / second, 6.777200e+07 / second};
        addInterpolator(
            DEUTERIUM, PROTIUM, DOWN, tableTemperature, table, numTemperatures);
    }

    {
        G4double table[numTemperatures]
            = {0.000000e+00 / second, 0.000000e+00 / second,
               0.000000e+00 / second, 0.000000e+00 / second,
               0.000000e+00 / second, 0.000000e+00 / second,
               0.000000e+00 / second, 0.000000e+00 / second,
               0.000000e+00 / second, 1.950000e+02 / second,
               2.752200e+04 / second, 2.393500e+06 / second,
               2.084800e+07 / second, 7.599500e+07 / second,
               1.858000e+08 / second, 3.499600e+08 / second,
               5.766300e+08 / second, 8.443300e+08 / second,
               1.148800e+09 / second, 1.480700e+09 / second,
               1.832500e+09 / second, 2.197200e+09 / second,
               2.569800e+09 / second, 2.946000e+09 / second};
        addInterpolator(
            PROTIUM, DEUTERIUM, UP, tableTemperature, table, numTemperatures);
    }

    {
        G4double table[numTemperatures]
            = {1.286600e+10 / second, 1.249500e+10 / second,
               1.201600e+10 / second, 1.148700e+10 / second,
               1.112300e+10 / second, 1.079100e+10 / second,
               1.040600e+10 / second, 1.022300e+10 / second,
               1.002800e+10 / second, 9.858600e+09 / second,
               9.717500e+09 / second, 9.589000e+09 / second,
               9.517600e+09 / second, 9.405000e+09 / second,
               9.396400e+09 / second, 9.319600e+09 / second,
               9.346700e+09 / second, 9.311200e+09 / second,
               9.278200e+09 / second, 9.246500e+09 / second,
               9.215800e+09 / second, 9.185700e+09 / second,
               9.156100e+09 / second, 9.126700e+09 / second};
        addInterpolator(
            PROTIUM, DEUTERIUM, DOWN, tableTemperature, table, numTemperatures);
    }

    {
        G4double table[numTemperatures]
            = {0.000000e+00 / second, 0.000000e+00 / second,
               0.000000e+00 / second, 0.000000e+00 / second,
               0.000000e+00 / second, 0.000000e+00 / second,
               0.000000e+00 / second, 0.000000e+00 / second,
               0.000000e+00 / second, 1.041000e+03 / second,
               9.192100e+04 / second, 6.172900e+06 / second,
               5.088800e+07 / second, 1.856300e+08 / second,
               4.418800e+08 / second, 8.260100e+08 / second,
               1.324500e+09 / second, 1.940000e+09 / second,
               2.665700e+09 / second, 3.427400e+09 / second,
               4.232500e+09 / second, 5.065100e+09 / second,
               5.914000e+09 / second, 6.767900e+09 / second};
        addInterpolator(
            PROTIUM, PROTIUM, UP, tableTemperature, table, numTemperatures);
    }

    {
        G4double table[numTemperatures]
            = {2.847900e+10 / second, 2.761900e+10 / second,
               2.650000e+10 / second, 2.524600e+10 / second,
               2.440700e+10 / second, 2.333500e+10 / second,
               2.164800e+10 / second, 2.048400e+10 / second,
               1.968700e+10 / second, 1.911000e+10 / second,
               1.900200e+10 / second, 1.873400e+10 / second,
               1.862100e+10 / second, 1.861900e+10 / second,
               1.852500e+10 / second, 1.840800e+10 / second,
               1.826100e+10 / second, 1.828100e+10 / second,
               1.840900e+10 / second, 1.837200e+10 / second,
               1.833600e+10 / second, 1.830200e+10 / second,
               1.826600e+10 / second, 1.822900e+10 / second};
        addInterpolator(
            PROTIUM, PROTIUM, DOWN, tableTemperature, table, numTemperatures);
    }

    {
        G4double table[numTemperatures]
            = {0.000000e+00 / second, 0.000000e+00 / second,
               0.000000e+00 / second, 0.000000e+00 / second,
               0.000000e+00 / second, 0.000000e+00 / second,
               0.000000e+00 / second, 0.000000e+00 / second,
               0.000000e+00 / second, 4.508500e+01 / second,
               1.047200e+04 / second, 1.241900e+06 / second,
               1.220500e+07 / second, 4.866900e+07 / second,
               1.222400e+08 / second, 2.441100e+08 / second,
               4.074300e+08 / second, 6.096600e+08 / second,
               8.443000e+08 / second, 1.105100e+09 / second,
               1.385800e+09 / second, 1.681000e+09 / second,
               1.986700e+09 / second, 2.299100e+09 / second};
        addInterpolator(
            PROTIUM, TRITIUM, UP, tableTemperature, table, numTemperatures);
    }

    {
        G4double table[numTemperatures]
            = {1.315400e+10 / second, 1.276900e+10 / second,
               1.227100e+10 / second, 1.171100e+10 / second,
               1.129000e+10 / second, 1.098900e+10 / second,
               1.056100e+10 / second, 1.038400e+10 / second,
               1.015900e+10 / second, 9.982500e+09 / second,
               9.801800e+09 / second, 9.659100e+09 / second,
               9.580700e+09 / second, 9.525100e+09 / second,
               9.429000e+09 / second, 9.429000e+09 / second,
               9.379700e+09 / second, 9.333700e+09 / second,
               9.290400e+09 / second, 9.249100e+09 / second,
               9.209400e+09 / second, 9.170900e+09 / second,
               9.133600e+09 / second, 9.097200e+09 / second};
        addInterpolator(
            PROTIUM, TRITIUM, DOWN, tableTemperature, table, numTemperatures);
    }

    {
        G4double table[numTemperatures]
            = {0.000000e+00 / second, 0.000000e+00 / second,
               0.000000e+00 / second, 0.000000e+00 / second,
               0.000000e+00 / second, 0.000000e+00 / second,
               0.000000e+00 / second, 0.000000e+00 / second,
               0.000000e+00 / second, 0.000000e+00 / second,
               1.806600e+02 / second, 4.355800e+04 / second,
               6.025100e+05 / second, 2.896300e+06 / second,
               8.240900e+06 / second, 1.764000e+07 / second,
               3.118500e+07 / second, 4.870600e+07 / second,
               6.973500e+07 / second, 9.368000e+07 / second,
               1.199700e+08 / second, 1.480500e+08 / second,
               1.795400e+08 / second, 2.076600e+08 / second};
        addInterpolator(
            TRITIUM, PROTIUM, UP, tableTemperature, table, numTemperatures);
    }

    {
        G4double table[numTemperatures]
            = {1.096200e+09 / second, 1.056300e+09 / second,
               1.004800e+09 / second, 9.451900e+08 / second,
               9.005900e+08 / second, 8.612100e+08 / second,
               8.093200e+08 / second, 7.851100e+08 / second,
               7.586400e+08 / second, 7.371500e+08 / second,
               7.194700e+08 / second, 7.051300e+08 / second,
               6.982700e+08 / second, 6.944300e+08 / second,
               6.882100e+08 / second, 6.896500e+08 / second,
               6.875200e+08 / second, 6.857700e+08 / second,
               6.842900e+08 / second, 6.830200e+08 / second,
               6.818700e+08 / second, 6.808300e+08 / second,
               6.798400e+08 / second, 6.788900e+08 / second};
        addInterpolator(
            TRITIUM, PROTIUM, DOWN, tableTemperature, table, numTemperatures);
    }

    {
        G4double table[numTemperatures] = {
            0.53961E+08 / second, 0.52295E+08 / second, 0.50155E+08 / second,
            0.47665E+08 / second, 0.46755E+08 / second, 0.45978E+08 / second,
            0.46476E+08 / second, 0.47078E+08 / second, 0.47384E+08 / second,
            0.48318E+08 / second, 0.49982E+08 / second, 0.53014E+08 / second,
            0.57235E+08 / second, 0.61389E+08 / second, 0.66476E+08 / second,
            0.71389E+08 / second, 0.76476E+08 / second, 0.81719E+08 / second,
            0.87111E+08 / second, 0.92646E+08 / second, 0.98320E+08 / second,
            0.10413E+09 / second, 0.11008E+09 / second, 0.11616E+09 / second};
        addInterpolator(
            DEUTERIUM, DEUTERIUM, DOWN, tableTemperature, table, numTemperatures);
    }

    {
        G4double table[numTemperatures] = {
            0.00000E+00 / second, 0.00000E+00 / second, 0.00000E+00 / second,
            0.00000E+00 / second, 0.00000E+00 / second, 0.00000E+00 / second,
            0.10661E+03 / second, 0.50077E+04 / second, 0.82913E+05 / second,
            0.77414E+06 / second, 0.25105E+07 / second, 0.86580E+07 / second,
            0.17339E+08 / second, 0.27135E+08 / second, 0.38166E+08 / second,
            0.49370E+08 / second, 0.60892E+08 / second, 0.72660E+08 / second,
            0.84645E+08 / second, 0.96838E+08 / second, 0.10924E+09 / second,
            0.12186E+09 / second, 0.13470E+09 / second, 0.14776E+09 / second};
        addInterpolator(
            DEUTERIUM, DEUTERIUM, UP, tableTemperature, table, numTemperatures);
    }

    {
        G4double table[numTemperatures] = {
            0.31194E+08 / second, 0.30180E+08 / second, 0.28826E+08 / second,
            0.27458E+08 / second, 0.26584E+08 / second, 0.25543E+08 / second,
            0.24806E+08 / second, 0.24433E+08 / second, 0.24599E+08 / second,
            0.25179E+08 / second, 0.25795E+08 / second, 0.27744E+08 / second,
            0.30127E+08 / second, 0.32865E+08 / second, 0.35337E+08 / second,
            0.37844E+08 / second, 0.40382E+08 / second, 0.42959E+08 / second,
            0.45583E+08 / second, 0.48262E+08 / second, 0.51002E+08 / second,
            0.53810E+08 / second, 0.56689E+08 / second, 0.59642E+08 / second};
        addInterpolator(
            DEUTERIUM, TRITIUM, DOWN, tableTemperature, table, numTemperatures);
    }

    {
        G4double table[numTemperatures] = {
            0.00000E+00 / second, 0.00000E+00 / second, 0.00000E+00 / second,
            0.00000E+00 / second, 0.00000E+00 / second, 0.00000E+00 / second,
            0.34793E+02 / second, 0.18009E+04 / second, 0.32221E+05 / second,
            0.31624E+06 / second, 0.10353E+07 / second, 0.37715E+07 / second,
            0.77642E+07 / second, 0.12586E+08 / second, 0.17627E+08 / second,
            0.22847E+08 / second, 0.28167E+08 / second, 0.33556E+08 / second,
            0.39006E+08 / second, 0.44524E+08 / second, 0.50114E+08 / second,
            0.55791E+08 / second, 0.61557E+08 / second, 0.67431E+08 / second};
        addInterpolator(
            DEUTERIUM, TRITIUM, UP, tableTemperature, table, numTemperatures);
    }

    {
        G4double table[numTemperatures] = {
            0.11373E+10 / second, 0.10941E+10 / second, 0.10365E+10 / second,
            0.97510E+09 / second, 0.93351E+09 / second, 0.87925E+09 / second,
            0.82610E+09 / second, 0.79157E+09 / second, 0.76728E+09 / second,
            0.74172E+09 / second, 0.72146E+09 / second, 0.70314E+09 / second,
            0.69592E+09 / second, 0.69537E+09 / second, 0.69095E+09 / second,
            0.68744E+09 / second, 0.68446E+09 / second, 0.68176E+09 / second,
            0.67927E+09 / second, 0.67691E+09 / second, 0.67466E+09 / second,
            0.67246E+09 / second, 0.67032E+09 / second, 0.66824E+09 / second};
        addInterpolator(
            TRITIUM, DEUTERIUM, DOWN, tableTemperature, table, numTemperatures);
    }

    {
        G4double table[numTemperatures] = {
            0.00000E+00 / second, 0.00000E+00 / second, 0.00000E+00 / second,
            0.00000E+00 / second, 0.00000E+00 / second, 0.00000E+00 / second,
            0.00000E+00 / second, 0.00000E+00 / second, 0.00000E+00 / second,
            0.00000E+00 / second, 0.12913E+03 / second, 0.35420E+05 / second,
            0.50797E+06 / second, 0.24886E+07 / second, 0.72276E+07 / second,
            0.15429E+08 / second, 0.27783E+08 / second, 0.43691E+08 / second,
            0.62936E+08 / second, 0.85024E+08 / second, 0.10943E+09 / second,
            0.13566E+09 / second, 0.16324E+09 / second, 0.19188E+09 / second};
        addInterpolator(
            TRITIUM, DEUTERIUM, UP, tableTemperature, table, numTemperatures);
    }

    {
        G4double table[numTemperatures] = {
            0.26463E+10 / second, 0.25306E+10 / second, 0.23839E+10 / second,
            0.21772E+10 / second, 0.20188E+10 / second, 0.18191E+10 / second,
            0.16665E+10 / second, 0.16025E+10 / second, 0.15370E+10 / second,
            0.14752E+10 / second, 0.14521E+10 / second, 0.14072E+10 / second,
            0.13981E+10 / second, 0.13827E+10 / second, 0.13702E+10 / second,
            0.13595E+10 / second, 0.13500E+10 / second, 0.13413E+10 / second,
            0.13334E+10 / second, 0.13261E+10 / second, 0.13193E+10 / second,
            0.13130E+10 / second, 0.13071E+10 / second, 0.13016E+10 / second};
        addInterpolator(
            TRITIUM, TRITIUM, DOWN, tableTemperature, table, numTemperatures);
    }

    {
        G4double table[numTemperatures] = {
            0.00000E+00 / second, 0.00000E+00 / second, 0.00000E+00 / second,
            0.00000E+00 / second, 0.00000E+00 / second, 0.00000E+00 / second,
            0.00000E+00 / second, 0.00000E+00 / second, 0.00000E+00 / second,
            0.00000E+00 / second, 0.45020E+02 / second, 0.19680E+05 / second,
            0.37377E+06 / second, 0.21045E+07 / second, 0.70457E+07 / second,
            0.16499E+08 / second, 0.31514E+08 / second, 0.52504E+08 / second,
            0.79345E+08 / second, 0.11170E+09 / second, 0.14892E+09 / second,
            0.19041E+09 / second, 0.23535E+09 / second, 0.28304E+09 / second,
        };
        addInterpolator(
            TRITIUM, TRITIUM, UP, tableTemperature, table, numTemperatures);
    }

    if (verboseLevel > 0)
        G4cout << "MuonicAtomSpinFlip Created" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MuonicAtomSpinFlip::~MuonicAtomSpinFlip()
{
    delete theTotalResult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool MuonicAtomSpinFlip::IsApplicable(G4ParticleDefinition const& p)
{
    return (&p == G4GenericMuonicAtom::GenericMuonicAtom());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double MuonicAtomSpinFlip::AtRestGetPhysicalInteractionLength(
    G4Track const& track, G4ForceCondition* condition)
{
    *condition = NotForced;

    // check if this is the beginning of tracking
    if (theNumberOfInteractionLengthLeft < 0.)
    {
        ResetNumberOfInteractionLengthLeft();
    }

    G4double interactionLength = theNumberOfInteractionLengthLeft
                                 * GetMeanSpinFlipTime(track);
    return (interactionLength);
}

G4double MuonicAtomSpinFlip::GetMeanSpinFlipTime(G4Track const& track)
{
    // By default, process does not occur
    G4double meanSpinFlipTime = std::numeric_limits<G4double>::infinity();

    // Update density information
    updateMaterialInfo(&track);
    equilibrateHydrogens(&track);  // TODO eventually allow non-equilbrium
                                   // mixtures

    // Check the type of particle
    G4DynamicParticle const* p = track.GetDynamicParticle();
    G4ParticleDefinition const* pd = p->GetDefinition();

    // Spin flip rate is proportional to density per
    // https://doi.org/10.1016/0375-9601(89)90366-6

    // If this is hydrogen
    if (pd->GetAtomicNumber() == 1)
    {
        G4int spin = (G4int)round(2.0 * p->GetSpin());
        G4int mass = pd->GetAtomicMass();

        // Get isotope number (number of neutrons) from mass
        G4int muonicAtomIsotope = mass - 1;
        if ((muonicAtomIsotope < PROTIUM) || (muonicAtomIsotope > TRITIUM))
            G4cout << "Error, mass out of range." << G4endl;

        // Set spin flip direction; spin is stored multiplied by two
        G4int flipDirection = UP;
        if (((muonicAtomIsotope == PROTIUM) && (spin == 2))
            || ((muonicAtomIsotope == DEUTERIUM) && (spin == 3))
            || ((muonicAtomIsotope == TRITIUM) && (spin == 2)))
            flipDirection = DOWN;

        // Loop over possible molecules with which the muonic atom can collide
        // and sum rates
        G4double lambda = 0;
        for (int otherIsotope = PROTIUM; otherIsotope <= TRITIUM;
             otherIsotope++)
            lambda += moleculeDensity[muonicAtomIsotope][otherIsotope]
                      * spinInterpolator[muonicAtomIsotope][otherIsotope]
                                        [flipDirection]
                                            ->CubicSplineInterpolation(
                                                temperature);
        meanSpinFlipTime = 1 / lambda;
    }
    return (meanSpinFlipTime);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange*
MuonicAtomSpinFlip::AtRestDoIt(G4Track const& track, G4Step const&)
{
    theTotalResult->Initialize(track);

    // NOTE: It does not look like there is a way to change the spin
    // using G4ParticleChange so this function kills the old muonic
    // atom and creates a new one in its place.

    // G4cout << "SpinFlip at rest do it"  << G4endl;

    G4double finalGlobalTime = track.GetGlobalTime();
    G4double finalLocalTime = track.GetLocalTime();
    G4double interactionTime = theNumberOfInteractionLengthLeft
                               * GetMeanSpinFlipTime(track);
    finalGlobalTime += interactionTime;
    finalLocalTime += interactionTime;

    theTotalResult->SetNumberOfSecondaries(1);

    // Check the type of particle
    G4DynamicParticle const* p = track.GetDynamicParticle();
    G4ParticleDefinition const* pd = p->GetDefinition();

    G4int spin = (G4int)round(2.0 * p->GetSpin());
    G4int mass = pd->GetAtomicMass();
    G4int newSpin = 0;

    // Flip the spin

    // If protium or tritum
    if ((mass == 1) || (mass == 3))
    {
        if (spin == 2)
            newSpin = 0;
        else if (spin == 0)
            newSpin = 2;
    }
    // If deuterium
    else if (mass == 2)
    {
        if (spin == 3)
            newSpin = 1;
        else if (spin == 1)
            newSpin = 3;
    }

    // G4cout << "starting spin was "  << spin << " ending spin will be " <<
    // newSpin << G4endl;

    G4IonTable* itp = G4IonTable::GetIonTable();
    G4ParticleDefinition* muonicAtom = itp->GetMuonicAtom(1, mass);
    G4DynamicParticle* dp
        = new G4DynamicParticle(muonicAtom, G4RandomDirection(), 0.);
    dp->SetSpin(newSpin);

    theTotalResult->AddSecondary(dp, finalGlobalTime, true);
    theTotalResult->ProposeLocalEnergyDeposit(0.0);
    theTotalResult->ProposeTrackStatus(fStopAndKill);
    theTotalResult->ProposeLocalTime(finalLocalTime);

    ClearNumberOfInteractionLengthLeft();
    return theTotalResult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MuonicAtomSpinFlip::ProcessDescription(std::ostream& outFile) const
{
    outFile << "Temperature-dependent spin flip in muonic protium, deuterium, "
               "and tritium";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MuonicAtomSpinFlip::addInterpolator(int muonIsotope,
                                         int otherIsotope,
                                         int spin,
                                         G4double* tableTemperature,
                                         G4double* table,
                                         int numTemperatures)
{
    spinInterpolator[muonIsotope][otherIsotope][spin] = new G4DataInterpolation(
        tableTemperature,
        table,
        numTemperatures,
        (table[1] - table[0]) / (tableTemperature[1] - tableTemperature[0]),
        (table[numTemperatures - 1] - table[numTemperatures - 2])
            / (tableTemperature[numTemperatures - 1]
               - tableTemperature[numTemperatures - 2]));
}
