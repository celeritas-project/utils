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
// Base class for muon-catalyzed fusion processes
//

#ifndef VMuonCatalyzedFusionProcess_hh
#define VMuonCatalyzedFusionProcess_hh 1

#include <G4DataInterpolation.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4VProcess.hh>
#include <G4VRestProcess.hh>
#include <G4ios.hh>
#include <globals.hh>

class VMuonCatalyzedFusionProcess : public G4VRestProcess
{
  public:
    VMuonCatalyzedFusionProcess(G4String const& aName,
                                G4ProcessType aType = fNotDefined);
    VMuonCatalyzedFusionProcess(VMuonCatalyzedFusionProcess&);

    virtual ~VMuonCatalyzedFusionProcess();

    bool particleMatches(G4Track const* track, G4int Z, G4int A);
    void updateMaterialInfo(G4Track const* track);

    // Atomic densities
    G4double tritiumPhi = 0.;
    G4double deuteriumPhi = 0.;
    G4double protiumPhi = 0.;

    // Molecular densities
    G4double HHPhi = 0.;
    G4double DDPhi = 0.;
    G4double TTPhi = 0.;
    G4double DTPhi = 0.;
    G4double HDPhi = 0.;
    G4double HTPhi = 0.;

    // Array form of molecule densities,
    // indexed with isotope identifiers
    G4double moleculeDensity[3][3];

    G4double temperature = 0.;

    G4double const protiumAtomicMass = 1.007825031898 * CLHEP::amu;
    G4double const deuteriumAtomicMass = 2.014101777844 * CLHEP::amu;
    G4double const tritiumAtomicMass = 3.016049281320 * CLHEP::amu;
    G4double const liquidHydrogenAtomicNumberDensity = 4.25e22 * (1 / cm3);

    G4DataInterpolation* muonEnergyInterpolator;

  protected:
    void equilibrateHydrogens(G4Track const* track);

    // Isotope identifiers
    int const PROTIUM = 0;
    int const DEUTERIUM = 1;
    int const TRITIUM = 2;

    // Spin direction identifiers
    int const DOWN = 0;
    int const UP = 1;

  private:
    VMuonCatalyzedFusionProcess();
    // Hidden default constructor

    void pairwiseEquilibrate(G4double c_a,
                             G4double c_b,
                             G4double k_ab,
                             G4double* f_a,
                             G4double* f_ab,
                             G4double* f_b);
    G4double kDT(G4double temperature);
    G4double kHT(G4double temperature);
    G4double kHD(G4double temperature);
};

#endif
