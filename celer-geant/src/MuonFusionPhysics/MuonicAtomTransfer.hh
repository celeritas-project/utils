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
// Muonic atom transfer physics process
//

#ifndef MuonicAtomTransfer_h
#define MuonicAtomTransfer_h 1

#include <G4ElementSelector.hh>
#include <G4ForceCondition.hh>
#include <G4GenericMuonicAtom.hh>
#include <G4HadFinalState.hh>
#include <G4HadronicInteraction.hh>
#include <G4HadronicProcessType.hh>
#include <G4ParticleDefinition.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4VPhysicsConstructor.hh>
#include <G4VRestProcess.hh>
#include <globals.hh>

#include "MuonicAtomTransfer.hh"
#include "VMuonCatalyzedFusionProcess.hh"
class G4HadronicInteraction;

typedef struct exchangeRate
{
    G4DataInterpolation* interpolator;
    G4int startIsotope;
    G4int isotope1;
    G4int isotope2;
    G4int endIsotope;
} exchangeRate;

class MuonicAtomTransfer : public VMuonCatalyzedFusionProcess

{
  public:
    explicit MuonicAtomTransfer(G4String const& name = "MuonicAtomTransfer");
    virtual ~MuonicAtomTransfer();

    G4bool IsApplicable(G4ParticleDefinition const&);

    virtual void PreparePhysicsTable(G4ParticleDefinition const&);

    virtual void BuildPhysicsTable(G4ParticleDefinition const&);

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
    MuonicAtomTransfer& operator=(MuonicAtomTransfer const& right);
    MuonicAtomTransfer(MuonicAtomTransfer const&);

    G4double GetMeanTransferTime(G4Track const& track,
                                 G4bool selectElement,
                                 G4int* newZ,
                                 G4int* newA);

    G4ElementSelector* fElementSelector;

    G4HadronicInteraction* fEmCascade;

    G4ParticleChange* theTotalResult;

    G4HadFinalState* result;

    // Hydrogen isotpe transfer to helium isotope
    G4DataInterpolation* lambdapHe3Interpolator;
    G4DataInterpolation* lambdapHe4Interpolator;
    G4DataInterpolation* lambdadHe3Interpolator;
    G4DataInterpolation* lambdadHe4Interpolator;
    G4DataInterpolation* lambdatHe3Interpolator;
    G4DataInterpolation* lambdatHe4Interpolator;
    G4double minHeTemp;
    G4double maxHeTemp;

    // Hydrogen isotope transfer to hydrogen isotope
    G4double minHTemp;
    G4double maxHTemp;
    void addExchangeRate(int startingIsotope,
                         int isotope1,
                         int isotope2,
                         int endingIsotope,
                         G4double* tableTemperature,
                         G4double* table,
                         int numTemperatures);
    std::vector<exchangeRate> rateTable;
};
#endif
