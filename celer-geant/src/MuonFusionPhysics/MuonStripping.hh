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
// Muon stripping process
//

#ifndef MuonStripping_h
#define MuonStripping_h 1

/////////////
// Includes
/////////////

#include <G4DynamicParticle.hh>
#include <G4GenericMuonicAtom.hh>
#include <G4Material.hh>
#include <G4Step.hh>
#include <G4VDiscreteProcess.hh>
#include <Randomize.hh>
#include <globals.hh>
#include <templates.hh>

// Class Description:

// Discrete Process -- Stripping of muons from muonic alpha particles
// Class inherits publicly from G4VDiscreteProcess
//
// TODO Temporarily implemented as a process to act on H4 ions, because
// that is what we are using to simulate muonic alpha particles
// until we get the muonic alpha particle transportable, having correct
// implementation for standard electromagnetic processes, etc.

// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

class EventAction;
class MuonStripping : public G4VDiscreteProcess
{
  public:
    ////////////////////////////////
    // Constructors and Destructor
    ////////////////////////////////

    // TODO: not sure if this should be fElectroMagnetic or fUserDefined

    MuonStripping(G4String const& processName = "MuonStripping",
                  G4ProcessType type = fHadronic);
    ~MuonStripping();

  private:
    // G4OpAbsorption(const G4OpAbsorption &right);

    //////////////
    // Operators
    //////////////

    //  G4OpAbsorption& operator=(const G4OpAbsorption &right);

  public:
    ////////////
    // Methods
    ////////////

    G4bool IsApplicable(G4ParticleDefinition const& aParticleType);

    G4double
    GetMeanFreePath(G4Track const& aTrack, G4double, G4ForceCondition*);

    G4VParticleChange* PostStepDoIt(G4Track const& aTrack, G4Step const& aStep);

    void PrintCrossSectionTable();
    G4double AlphaProtonCMFrameEnergy(G4double v);

  private:
    G4double CalculateCrossSection(G4double velocity);
    G4double* CalculateCrossSectionPerVolumeVector(G4Track const& aTrack);
    EventAction* fEventAction;
};

////////////////////
// Inline methods
////////////////////

inline G4bool MuonStripping::IsApplicable(G4ParticleDefinition const& p)
{
    return (&p == G4GenericMuonicAtom::GenericMuonicAtom());
}

#endif /* MuonStripping */
