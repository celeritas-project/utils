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
// Tempeature-dependent spin-flip of thermalized MuH1, MuH2, or MuH3

#ifndef MuonicAtomSpinFlip_h
#define MuonicAtomSpinFlip_h 1

#include <G4DataInterpolation.hh>
#include <G4ElementSelector.hh>
#include <G4ForceCondition.hh>
#include <G4HadFinalState.hh>
#include <G4HadronicInteraction.hh>
#include <G4HadronicProcessType.hh>
#include <G4ParticleDefinition.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4VRestProcess.hh>
#include <globals.hh>

#include "VMuonCatalyzedFusionProcess.hh"

class G4HadronicInteraction;

class MuonicAtomSpinFlip : public VMuonCatalyzedFusionProcess

{
  public:
    explicit MuonicAtomSpinFlip(G4String const& name = "MuonicAtomSpinFlip");

    ~MuonicAtomSpinFlip();

    G4bool IsApplicable(G4ParticleDefinition const&);

    virtual G4double
    AtRestGetPhysicalInteractionLength(G4Track const& track,
                                       G4ForceCondition* condition);

    virtual G4VParticleChange* AtRestDoIt(G4Track const&, G4Step const&);

    void ProcessDescription(std::ostream& outFile) const;

  protected:
    // set effective lifetime for at-rest process (default is forced action)
    // FIXME: This should be computed by subprocesses via cross-section
    // analogue
    G4double
    GetMeanLifeTime(G4Track const& /*aTrack*/, G4ForceCondition* /*condition*/)
    {
        return -1.0;
    }

  private:
    // hide assignment operator as private
    MuonicAtomSpinFlip& operator=(MuonicAtomSpinFlip const& right);
    MuonicAtomSpinFlip(MuonicAtomSpinFlip const&);

    G4double GetMeanSpinFlipTime(G4Track const& track);

    G4ElementSelector* fElementSelector;
    G4HadronicInteraction* fEmCascade;
    G4ParticleChange* theTotalResult;

    // Array of interpolators for the spin flip rates
    // format is [muon isotope][other isotope on molecule][spin down = 0, spin
    // up = 1] with the isotopes coded protium = 0, deuterium = 1, tritium = 2
    G4DataInterpolation* spinInterpolator[3][3][2];

    // Function to initialize each interpolator
    void addInterpolator(int muonIsotope,
                         int otherIsotope,
                         int spin,
                         G4double* tableTemperature,
                         G4double* table,
                         int numTemperatures);
};

#endif
