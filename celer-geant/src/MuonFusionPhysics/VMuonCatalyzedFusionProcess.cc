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
// Virtual muon-catalyzed fusion process base class
// Includes functions needed by all muon catalyzed fusion processes

#include "VMuonCatalyzedFusionProcess.hh"

#include <G4DataInterpolation.hh>
#include <G4MaterialTable.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>
#include <G4VParticleChange.hh>

// --------------------------------------------------------------------
VMuonCatalyzedFusionProcess::VMuonCatalyzedFusionProcess()
    : G4VRestProcess("Virtual Muon Catalyzed Fusion Process Base Class")
{
    G4Exception("VMuonCatalyzedFusionProcess::VMuonCatalyzedFusionProcess()",
                "ProcMan102",
                JustWarning,
                "Default constructor is called");
}

// --------------------------------------------------------------------
VMuonCatalyzedFusionProcess::VMuonCatalyzedFusionProcess(G4String const& aName,
                                                         G4ProcessType aType)
    : G4VRestProcess(aName, aType)
{
    // CDF of muon energy after muon-catalyzed d-t fusion
    // From Comprehensive study of nuclear reactions in muon catalyzed fusion:
    // I. dtμ system M. Kamimura, Y. Kino, T. Yamashita
    // https://doi.org/10.48550/arXiv.2112.08399
    // uses the top "red curve" TODO: understand which is correct curve to use
    // cuts off at 80 keV, which is 99% CDF

    G4int static const numMuonEnergyPoints = 21;

    G4double muonEnergy[numMuonEnergyPoints] = {0. * keV,
                                                0.48850540675768084 * keV,
                                                0.8390389347819425 * keV,
                                                1.2521213482687141 * keV,
                                                1.7153033196164724 * keV,
                                                2.253638712180777 * keV,
                                                2.854653691809707 * keV,
                                                3.606073540073316 * keV,
                                                4.470346052913727 * keV,
                                                5.560291219507215 * keV,
                                                6.700556502915258 * keV,
                                                7.953772477101693 * keV,
                                                9.194596305637525 * keV,
                                                10.849180562221111 * keV,
                                                12.353474314071864 * keV,
                                                14.045888515617822 * keV,
                                                15.650634617544647 * keV,
                                                17.38079707555165 * keV,
                                                19.111008546659452 * keV,
                                                19.976130619913615 * keV,
                                                80.0 * keV};

    G4double muonEnergyCDF[numMuonEnergyPoints] = {0,
                                                   0.04169381107491854,
                                                   0.08664495114006499,
                                                   0.14332247557003264,
                                                   0.20456026058631915,
                                                   0.2723127035830618,
                                                   0.34136807817589576,
                                                   0.41563517915309456,
                                                   0.48990228013029324,
                                                   0.5667752442996744,
                                                   0.6306188925081434,
                                                   0.6866449511400652,
                                                   0.7309446254071662,
                                                   0.7778501628664496,
                                                   0.8104234527687297,
                                                   0.8403908794788275,
                                                   0.8618892508143323,
                                                   0.8814332247557004,
                                                   0.8970684039087949,
                                                   0.903583061889251,
                                                   1.0};

    muonEnergyInterpolator
        = new G4DataInterpolation(muonEnergyCDF,
                                  muonEnergy,
                                  numMuonEnergyPoints,
                                  muonEnergy[1] / muonEnergyCDF[1],
                                  0.0);
}

VMuonCatalyzedFusionProcess::~VMuonCatalyzedFusionProcess()
{
    delete muonEnergyInterpolator;
}

bool VMuonCatalyzedFusionProcess::particleMatches(G4Track const* track,
                                                  G4int Z,
                                                  G4int A)
{
    G4DynamicParticle const* particle = track->GetDynamicParticle();
    G4ParticleDefinition const* particleDefinition
        = particle->GetParticleDefinition();
    return ((Z == particleDefinition->GetAtomicNumber())
            && (A == particleDefinition->GetAtomicMass()));
}

void VMuonCatalyzedFusionProcess::updateMaterialInfo(G4Track const* track)
{
    // Get material information
    G4Material* mat = track->GetMaterial();

    // Update density of each hydrogen isotope
    G4ElementVector const* elementVector = mat->GetElementVector();  // std
                                                                     // vector
                                                                     // of g4
                                                                     // element
    G4double const* fractionVector = mat->GetFractionVector();
    G4double density = mat->GetDensity();
    temperature = mat->GetTemperature();

    // Get relative mass fraction of hydrogen isotopes in mixture
    G4double protiumMassFraction = 0.0;
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
            if (thisZ == 1)
            {
                if (thisN == 1)
                    protiumMassFraction += fractionVector[i] * thisA
                                           * relativeAbundance[j] / atomicMass;
                else if (thisN == 2)
                    deuteriumMassFraction += fractionVector[i] * thisA
                                             * relativeAbundance[j]
                                             / atomicMass;
                else if (thisN == 3)
                    tritiumMassFraction += fractionVector[i] * thisA
                                           * relativeAbundance[j] / atomicMass;
            }
        }
    }

    // Compute phi, the atomic number density per volume relative to liquid
    // hydrogen density
    protiumPhi = protiumMassFraction * density / protiumAtomicMass
                 / liquidHydrogenAtomicNumberDensity;
    deuteriumPhi = deuteriumMassFraction * density / deuteriumAtomicMass
                   / liquidHydrogenAtomicNumberDensity;
    tritiumPhi = tritiumMassFraction * density / tritiumAtomicMass
                 / liquidHydrogenAtomicNumberDensity;
}

// Set densities of hydrogen isoprotologue molecule densities
// to long-time equlibirated densities based on isotope densities
// using the model from
// https://www.osti.gov/servlets/purl/6205719-wlGVDf/Cryogenichydrogendatapertinenttomagneticfusionenergy.pdf
void VMuonCatalyzedFusionProcess::equilibrateHydrogens(G4Track const* track)
{
    // Get material information
    G4Material* mat = track->GetMaterial();

    // Update temperature
    temperature = mat->GetTemperature();

    // Compute overall density
    G4double totalPhi = protiumPhi + deuteriumPhi + tritiumPhi;

    G4double HHf = protiumPhi / totalPhi;
    G4double DDf = deuteriumPhi / totalPhi;
    G4double TTf = tritiumPhi / totalPhi;
    G4double DTf = 0.;
    G4double HDf = 0.;
    G4double HTf = 0.;

    G4double xDiff = std::numeric_limits<G4double>::infinity();

    G4double kDTc = kDT(temperature);
    G4double kHTc = kHT(temperature);
    G4double kHDc = kHD(temperature);
    G4double xLast[6];
    G4double x[6];
    G4int maxIter = 1000;
    G4int iter = 0;

    // Compute pairwise equilibra repeatedly until concentrations are stable
    while ((xDiff > 1e-6) && (iter < maxIter))
    {
        xLast[0] = HHf;
        xLast[1] = DDf;
        xLast[2] = TTf;
        xLast[3] = DTf;
        xLast[4] = HTf;
        xLast[5] = HDf;
        pairwiseEquilibrate(
            DDf + DTf / 2, TTf + DTf / 2, kDTc, &DDf, &DTf, &TTf);
        pairwiseEquilibrate(
            HHf + HTf / 2, TTf + HTf / 2, kHTc, &HHf, &HTf, &TTf);
        pairwiseEquilibrate(
            HHf + HDf / 2, DDf + HDf / 2, kHDc, &HHf, &HDf, &DDf);
        x[0] = HHf;
        x[1] = DDf;
        x[2] = TTf;
        x[3] = DTf;
        x[4] = HTf;
        x[5] = HDf;

        xDiff = 0.;
        for (G4int i = 0; i < 6; i++)
        {
            G4double diff = abs(x[i] - xLast[i]);
            if (diff > xDiff)
                xDiff = diff;
        }
        iter++;
        // G4cout << "Ran an iter " << iter << " xDiff is now " << xDiff <<
        // G4endl;
    }

    if (iter == maxIter)
    {
        G4cout << "WARNING: Reached maxIter calculing hydrogen isoprotologue "
                  "equalibrium; results may not be valid"
               << G4endl;
        G4cout << "xDiff is " << xDiff << G4endl;
    }

    HHPhi = HHf * totalPhi;
    DDPhi = DDf * totalPhi;
    TTPhi = TTf * totalPhi;
    DTPhi = DTf * totalPhi;
    HDPhi = HDf * totalPhi;
    HTPhi = HTf * totalPhi;

    // Fill in array form of molecule densities
    // which is a symmetric matrix
    // Calling code must avoid double counting
    moleculeDensity[PROTIUM][PROTIUM] = HHPhi;
    moleculeDensity[PROTIUM][DEUTERIUM] = HDPhi;
    moleculeDensity[PROTIUM][TRITIUM] = HTPhi;
    moleculeDensity[DEUTERIUM][PROTIUM] = HDPhi;
    moleculeDensity[DEUTERIUM][DEUTERIUM] = DDPhi;
    moleculeDensity[DEUTERIUM][TRITIUM] = DTPhi;
    moleculeDensity[TRITIUM][PROTIUM] = HTPhi;
    moleculeDensity[TRITIUM][DEUTERIUM] = DTPhi;
    moleculeDensity[TRITIUM][TRITIUM] = TTPhi;
}

void VMuonCatalyzedFusionProcess::pairwiseEquilibrate(G4double c_a,
                                                      G4double c_b,
                                                      G4double k_ab,
                                                      G4double* f_a,
                                                      G4double* f_ab,
                                                      G4double* f_b)
{
    G4double eps = 1e-6;
    G4double sigma
        = ((c_a + c_b)
           - sqrt((c_a - c_b) * (c_a - c_b) + 16 * c_a * c_b / (k_ab - eps)))
          / (2 * (1 - 4 / (k_ab - eps)));
    *f_a = c_a - sigma;
    *f_ab = 2 * sigma;
    *f_b = c_b - sigma;
}

G4double VMuonCatalyzedFusionProcess::kDT(G4double T)
{
    G4double R = 8.314 * joule / mole / kelvin;
    G4double c_dt_30 = R * (30 * kelvin) * (log(4) - log(2.09));  // TODO could
                                                                  // optimize
    G4double c_dt_100 = R * (100 * kelvin) * (log(4) - log(3.29));

    if (T < (15 * kelvin))
        return (5.924 * exp(-(168.3 * joule / mole) / (R * T)));
    else if (T < (30 * kelvin))
        return (2.995 * exp(-(89.96 * joule / mole) / (R * T)));
    else if (T < (100 * kelvin))
        return (4.0 * exp(-c_dt_30 / (R * T)));
    else
        return (4.0 * exp(-c_dt_100 / (R * T)));
}

G4double VMuonCatalyzedFusionProcess::kHT(G4double T)
{
    G4double R = 8.314 * joule / mole / kelvin;
    G4double c_ht = R * (30 * kelvin) * (log(4) - log(0.034));

    if (T < (30 * kelvin))
        return (10.22 * exp(-1423. / (R * T)));
    else
        return (4.0 * exp(-c_ht / (R * T)));
}

G4double VMuonCatalyzedFusionProcess::kHD(G4double T)
{
    G4double R = 8.314 * joule / mole / kelvin;
    G4double c_hd = R * (30 * kelvin) * (log(4) - log(0.49));

    if (T < (30 * kelvin))
        return (6.785 * exp(-654.3 / (R * T)));
    else
        return (4.0 * exp(-c_hd / (R * T)));
}
