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
// Muon stripping in ground-state MuHe3 and MuHe4 due to collisional
// ionization and transfer processes

#include "MuonStripping.hh"

#include <G4Alpha.hh>
#include <G4DynamicParticle.hh>
#include <G4Element.hh>
#include <G4HadDecayGenerator.hh>
#include <G4HadronicProcess.hh>
#include <G4IonTable.hh>
#include <G4Material.hh>
#include <G4MuonMinus.hh>
#include <G4OpProcessSubType.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4ios.hh>
#include <Randomize.hh>

/////////////////////////
// Class Implementation
/////////////////////////

//////////////
// Operators
//////////////

/////////////////
// Constructors
/////////////////

MuonStripping::MuonStripping(G4String const& processName, G4ProcessType type)
    : G4VDiscreteProcess(processName, type)
{
    if (verboseLevel > 0)
    {
        G4cout << GetProcessName() << " is created " << G4endl;
    }

    SetProcessSubType(fHadronInelastic);
}

////////////////
// Destructors
////////////////

MuonStripping::~MuonStripping() {}

////////////
// Methods
////////////

// PostStepDoIt
// -------------
//
G4VParticleChange*
MuonStripping::PostStepDoIt(G4Track const& aTrack, G4Step const&)
{
    aParticleChange.Initialize(aTrack);

    if (verboseLevel > 0)
        G4cout << G4endl << "-------- MuonStripping PostStepDoIt Begin-----"
               << G4endl;

    // STEPS:
    // 1. Get velocity vector of the MuH4 or MuH3 projectile
    // TODO: now that we are using bona-fide muonic atoms rather than H4, it
    // might make sense to get the energy or momentum directly instead.
    G4double velocity = aTrack.GetVelocity();
    G4double initialLabFrameKE = aTrack.GetKineticEnergy();
    G4ThreeVector direction = aTrack.GetMomentumDirection();

    // 1.5 Figure out what atom in the material the projectile hit
    G4Material* material = aTrack.GetMaterial();
    G4int numberOfElements = material->GetNumberOfElements();
    if (verboseLevel > 0)
        G4cout << "Material is " << material->GetName() << G4endl;
    if (verboseLevel > 0)
        G4cout << "Number of Elements is " << numberOfElements << G4endl;

    G4double* crossSectionPerVolumeVector
        = CalculateCrossSectionPerVolumeVector(aTrack);

    // Compute total cross section
    G4double totalCrossSectionPerVolume = 0.0;
    for (int i = 0; i < numberOfElements; i++)
    {
        totalCrossSectionPerVolume += crossSectionPerVolumeVector[i];
    }
    if (verboseLevel > 0)
        G4cout << "totalCrossSectionPerVolume is "
               << totalCrossSectionPerVolume << G4endl;

    // Detemine which element the interaction is with, weighted by the
    // relative cross sections at the current velocity
    G4double cumulativeProbability = 0.0;
    G4double randomNumber = G4UniformRand();
    if (verboseLevel > 0)
        G4cout << "random number is " << randomNumber << G4endl;
    G4int elementIndex = numberOfElements - 1;
    for (int i = 0; i < numberOfElements; i++)
    {
        cumulativeProbability += crossSectionPerVolumeVector[i]
                                 / totalCrossSectionPerVolume;
        if (verboseLevel > 0)
            G4cout << "cumulative probability is " << cumulativeProbability
                   << G4endl;
        if (randomNumber < cumulativeProbability)
        {
            elementIndex = i;
            if (verboseLevel > 0)
                G4cout << "element " << elementIndex << " is selected. "
                       << G4endl;
            break;
        }
    }
    // Free memory allocated by CalculateCrossSectionVector
    delete[] crossSectionPerVolumeVector;

    // Now we need to randomly select an isotope from the element to interact
    // with, based on its relative abundance in the material
    G4Element const* element = material->GetElement(elementIndex);
    G4int numberOfIsotopes = element->GetNumberOfIsotopes();
    if (verboseLevel > 0)
        G4cout << "Element name is " << element->GetName() << G4endl;
    if (verboseLevel > 0)
        G4cout << "Number of isotopes is " << numberOfIsotopes << G4endl;

    G4double* relativeAbundanceVector = element->GetRelativeAbundanceVector();
    randomNumber = G4UniformRand();
    if (verboseLevel > 0)
        G4cout << "random number is " << randomNumber << G4endl;
    cumulativeProbability = 0.0;
    G4int isotopeIndex;
    for (int i = 0; i < numberOfIsotopes; i++)
    {
        cumulativeProbability += relativeAbundanceVector[i];
        if (verboseLevel > 0)
            G4cout << "cumulative probability is " << cumulativeProbability
                   << G4endl;
        if (randomNumber < cumulativeProbability)
        {
            isotopeIndex = i;
            if (verboseLevel > 0)
                G4cout << "isotope " << isotopeIndex << " is selected. "
                       << G4endl;
            break;
        }
    }

    // Get information about the target atom
    G4int Z = element->GetIsotope(isotopeIndex)->GetZ();  // atomic number
    G4int N = element->GetIsotope(isotopeIndex)->GetN();  // atomic mass
    G4double A = element->GetIsotope(isotopeIndex)->GetA();  // atomic mass per
                                                             // mole

    if (verboseLevel > 0)
        G4cout << "Z = " << Z << G4endl;
    if (verboseLevel > 0)
        G4cout << "N = " << N << G4endl;
    if (verboseLevel > 0)
        G4cout << "A = " << A / (g / mole) << " gram/mole" << G4endl;

    // 2. For this relative velocity, compute the total mass-energy of
    // combination of the muonic alpha and the target atom, taking into account
    // the rest masses of the particles, the center-of-momentum frame kinetic
    // energy of the particles, and the mass defect from the bonding energy of
    // the muon and alpha.  This gives the "initial mass"

    // TODO now that we are using bona-fide muonic atoms rather than H4's, if
    // we end up keeping this phase-space generator approach (unlikely) then
    // this could could just use the mass reported by GEANT4 rather than
    // computing it seperately.

    G4int projectileA = aTrack.GetParticleDefinition()->GetAtomicMass();

    G4double targetMass = A / Avogadro;
    G4double muonMass = 1.883531594e-28 * kg;
    G4double nucleusMass = (G4double)projectileA / 4.0 * 6.644657230e-27 * kg;
    G4double muonToElectronMassRatio = 206.7682826;
    G4double heliumSecondIonizationEnergy = 54.417760 * eV;

    // Compute muonic alpha mass, including missing mass from muon binding
    // energy (important because that binding energy is absorbed from KE in the
    // interaction)

    G4double heliumMuonIonizationEnergy = heliumSecondIonizationEnergy
                                          * muonToElectronMassRatio;
    G4double muonicAtomMass = nucleusMass + muonMass
                              - (heliumMuonIonizationEnergy / c_squared);

    // Compute center-of-momentum kinetic energy between muonic alpha and
    // target atom Would be nice to make this calc relativistic, but not high
    // prioirity, this is only gamma = 1.0009
    G4double muonicAlphaCMVelocity = velocity
                                     / (1 + (muonicAtomMass / targetMass));
    G4double targetCMVelocity = velocity - muonicAlphaCMVelocity;
    G4double kineticEnergyCM
        = 0.5 * muonicAtomMass * muonicAlphaCMVelocity * muonicAlphaCMVelocity
          + 0.5 * targetMass * targetCMVelocity * targetCMVelocity;

    // 3. Use the Raubold-Lynch phase space method to calculate final
    // four-momentum vectors for the particles in the center of momentum frame.
    // You will need to pass in the initial mass, and also the rest masses of
    // the particles generated; the rest mass of the target atom, of an alpha
    // particle, and of a muon.  The phase-space generator is an approximation
    // that assumes no interactions between the generated particles; if they
    // are interacting then you need matrix elements or differential cross
    // sections.  For now, as a first approximation, we assume they don't
    // interact; later, we can address this with better QM modeling of the
    // stripping process or experimental data.  (perhaps reproducing the model
    // in the paper)

    // TODO: Do QM or CTMC modeling to get matrix elements or differential
    // cross sections so you don't need to assume no interaction between
    // generated particles.  This might make a significant difference

    // Compute "initial mass" for phase space generator, sum of CM frame KE and
    // masses of the incoming particles
    G4double initialMass = muonicAtomMass + targetMass
                           + (kineticEnergyCM / c_squared);

    // Create vector of masses of secondaries
    std::vector<G4double> masses = {muonMass, nucleusMass, targetMass};

    // Initialize vector to hold particle final states
    std::vector<G4LorentzVector> finalState;

    // Call the phase space generator to create randomized final states
    // with energy and mass conservation

    G4double massAccumulator = 0.0;
    for (G4double mass : masses)
    {
        massAccumulator += mass;
    }

    if (verboseLevel > 0)
        G4cout << "excess mass in eV: "
               << (initialMass - massAccumulator) * c_light * c_light / eV
               << G4endl;

    // If the excess mass is less than zero, which is nonphysical,
    // do not do the interaction.  This can happen (due to inaccuracy
    // in the cross section at very low energies) if postStopDoIt is
    // triggered for a center-of-mass energy less than the helium
    // muon ionization energy.
    if (initialMass - massAccumulator < 0.0)
    {
        if (verboseLevel > 0)
            G4cout << "C/M KE less than helium muon ionization energy.  Muon "
                      "not stripped.";
        return &aParticleChange;
    }

    G4HadDecayGenerator generator;
    generator.Generate(initialMass, masses, finalState);

    if (verboseLevel > 0)
    {
        G4cout << "after call to generate." << G4endl;
        G4cout << "massesSize: " << masses.size() << G4endl;
        G4cout << "finalStateSize: " << finalState.size() << G4endl;
        G4cout << "finalStates: " << G4endl;
        for (G4LorentzVector f : finalState)
        {
            G4cout << "px: " << f.px() * c_light / (amu * meter / second)
                   << " py: " << f.py() * c_light / (amu * meter / second)
                   << " pz: " << f.pz() * c_light / (amu * meter / second)
                   << " e: " << f.e() / amu << " m: " << f.m() / amu << G4endl;
            G4cout
                << "classical KE: "
                << 0.5 / f.m() * c_light * c_light
                       * (f.px() * f.px() + f.py() * f.py() + f.pz() * f.pz())
                       / eV
                << G4endl;
            G4cout << "vx: " << (f.px() * c_light / f.m()) / (meter / second)
                   << " vy: " << (f.py() * c_light / f.m()) / (meter / second)
                   << " vz: " << (f.pz() * c_light / f.m()) / (meter / second)
                   << G4endl;
            G4double E = sqrt(f.m() * f.m() + f.px() * f.px() + f.py() * f.py()
                              + f.pz() * f.pz())
                         * c_light * c_light;
            G4cout << "Calc energy: " << (E / (c_light * c_light)) / amu
                   << G4endl;

            G4cout << "Lorentz factor: " << f.e() / f.m() << G4endl;
        }
        G4cout << G4endl;
    }

    // 4. Boost the four-vectors of the particles into the lab frame.

    // Calculate boost vector; velocity of the center of momentum frame in the
    // lab frame as a fraction of the speed of light.  (Done below as a
    // classical calcualation which is OK as even at 3.5 MeV this is only
    // gamma=1.0009)

    G4ThreeVector labFrameBoost = direction * targetCMVelocity / c_light;

    if (verboseLevel > 0)
    {
        G4cout
            << "boost vx: "
            << (direction.x() * targetCMVelocity) / (meter / second)
            << " vy: " << direction.y() * targetCMVelocity / (meter / second)
            << " vz: " << (direction.z() * targetCMVelocity) / (meter / second)
            << G4endl;
    }

    // Boost the final states back up in the lab frame
    for (size_t i = 0; i < finalState.size(); i++)
    {
        finalState[i].boost(labFrameBoost);
    }

    if (verboseLevel > 0)
    {
        G4cout << "boosted final states:" << G4endl;
        for (G4LorentzVector f : finalState)
        {
            G4cout << "px: " << f.px() * c_light / (amu * meter / second)
                   << " py: " << f.py() * c_light / (amu * meter / second)
                   << " pz: " << f.pz() * c_light / (amu * meter / second)
                   << " e: " << f.e() / amu << " m: " << f.m() / amu << G4endl;
            G4cout
                << "classical KE: "
                << 0.5 / f.m() * c_light * c_light
                       * (f.px() * f.px() + f.py() * f.py() + f.pz() * f.pz())
                       / eV
                << G4endl;
            G4cout << "vx: " << (f.px() * c_light / f.m()) / (meter / second)
                   << " vy: " << (f.py() * c_light / f.m()) / (meter / second)
                   << " vz: " << (f.pz() * c_light / f.m()) / (meter / second)
                   << G4endl;
            G4double E = sqrt(f.m() * f.m() + f.px() * f.px() + f.py() * f.py()
                              + f.pz() * f.pz())
                         * c_light * c_light;
            G4cout << "Calc energy: " << (E / (c_light * c_light)) / amu
                   << G4endl;
            G4cout << "Lorentz factor: " << f.e() / f.m() << G4endl;
            G4ThreeVector momentumDir = f.vect().unit();
            G4cout << "mdx: " << momentumDir.x() << " mdy: " << momentumDir.y()
                   << " mdz: " << momentumDir.z()
                   << " T: " << (f.e() - f.m()) * c_light * c_light / eV
                   << G4endl;
        }
        G4cout << G4endl;
    }

    G4ThreeVector momentumDir;
    G4double kineticEnergy;

    // 5. Generate the requisite secondary particles (target atom, alpha, and
    // muon) and destroy the primary.
    aParticleChange.SetNumberOfSecondaries(3);

    // The newly freed muon
    momentumDir = finalState[0].vect().unit();
    // kineticEnergy = (finalState[0].e() - finalState[0].m()) * c_light *
    // c_light;
    //  The Raubold-Lynch approch is almost certainly wrong. CTMC shows the
    //  muon velocity in the lab frame is typically between 0 and the muonic
    //  alpha velocity in the lab frame.   Very soon, we should replace
    //  raubold-lynch with a more correct description of the collision physics.
    //  For now, this is a huge hack to just say the muon velocity is random
    //  between 0 and the incoming lab frame velocity.  Of course this is very
    //  bad, breaks energy and momentum conservation. etc.  But it is likely
    //  much more correct than what was here before.
    //  TODO FIX FIX FIX!!!! DO NOT RELEASE!
    kineticEnergy = initialLabFrameKE * 0.1 / 4.1
                    * G4UniformRand();  // scale by mass ratio of muon vs.
                                        // muonic alpha

    G4DynamicParticle* aParticle0 = new G4DynamicParticle(
        G4MuonMinus::MuonMinus(), momentumDir, kineticEnergy);
    aParticleChange.AddSecondary(aParticle0);

    // The helium ion
    momentumDir = finalState[1].vect().unit();
    kineticEnergy = (finalState[1].e() - finalState[1].m()) * c_light * c_light;
    assert((projectileA == 3) || (projectileA == 4));
    G4DynamicParticle* aParticle1;
    if (projectileA == 4)
        aParticle1 = new G4DynamicParticle(
            G4Alpha::Alpha(), momentumDir, kineticEnergy);
    else
        aParticle1
            = new G4DynamicParticle(G4IonTable::GetIonTable()->GetIon(2, 3, 0),
                                    momentumDir,
                                    kineticEnergy);
    aParticle1->SetCharge(+2 * eplus);  // alpha++ or He3++, charge = +2 amu
    aParticleChange.AddSecondary(aParticle1);

    // The (now dynamic) target atom
    G4double excitEnergy = 0.0 * keV;
    G4ParticleDefinition* ion
        = G4IonTable::GetIonTable()->GetIon(Z, N, excitEnergy);
    momentumDir = finalState[2].vect().unit();
    kineticEnergy = (finalState[2].e() - finalState[2].m()) * c_light * c_light;
    G4DynamicParticle* aParticle2
        = new G4DynamicParticle(ion, momentumDir, kineticEnergy);
    aParticle2->SetCharge(0);  // charge = 0
    aParticleChange.AddSecondary(aParticle2);

    //
    // Kill the projectile muonic atom
    //
    aParticleChange.ProposeMomentumDirection(0., 0., 0.);
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);

    // TODO: check all this for conservation of mass-energy, momentum, and look
    // at results of a few runs to make sure it looks reasonable

    if (verboseLevel > 0)
    {
        G4cout << "\n** Muon stripped! **" << G4endl;
        G4cout << "-------- MuonStripping PostStepDoIt End -----" << G4endl
               << G4endl;
    }

    // TODO I don't know if this function is keeping track of time correctly.
    // Does it need to explicatly keep track of time the way that the muon
    // catalyzed fusion function does, and clear the number of interaction
    // lengths left here? Or is that only for at-rest processes?

    return &aParticleChange;
}

// GetMeanFreePath
// ---------------
//
G4double MuonStripping::GetMeanFreePath(G4Track const& aTrack,
                                        G4double,
                                        G4ForceCondition*)
{
    // MuonStripping should only trigger for MuHe3 and MuHe4.
    // This process could be extended to work with
    // higher-Z muonic atoms and with muonic hydrogen with
    // these high (keV-to-MeV) energies, but this
    // is not essential as there is not a known process that
    // would create them naturally in the muon catalyzed fusion
    // chain, and no planned experiment to create them artificially.

    G4ParticleDefinition const* particle = aTrack.GetParticleDefinition();
    G4int Z = particle->GetAtomicNumber();
    G4int A = particle->GetAtomicMass();
    if (!((Z == 2) && ((A == 3) || (A == 4))))
        return (std::numeric_limits<G4double>::infinity());

    // Calculate cross section per volume vector for each element, and sum to
    // get total cross section per volume
    G4Material* material = aTrack.GetMaterial();
    G4int numberOfElements = material->GetNumberOfElements();
    G4double* crossSectionPerVolumeVector
        = CalculateCrossSectionPerVolumeVector(aTrack);
    G4double totalCrossSectionPerVolume = 0.0;
    for (int i = 0; i < numberOfElements; i++)
    {
        totalCrossSectionPerVolume += crossSectionPerVolumeVector[i];
    }
    delete[] crossSectionPerVolumeVector;

    // Compute mean free path -- lambda = 1 / sum(sigma * n);
    G4double meanFreePath = 1.0 / totalCrossSectionPerVolume;

    if (verboseLevel > 1)
        G4cout << "In GetMeanFreePath " << G4endl;
    if (verboseLevel > 1)
        G4cout << "Total Cross Section Per Volume is "
               << totalCrossSectionPerVolume << G4endl;
    if (verboseLevel > 1)
        G4cout << "Mean Free Path is " << meanFreePath / mm << " mm." << G4endl;

    return (meanFreePath);
}

G4double*
MuonStripping::CalculateCrossSectionPerVolumeVector(G4Track const& aTrack)
{
    // Theoretical justification for scaling factors:
    // Scaling and formulary of cross sections for ion-atom impact ionization
    // I. Kaganovich, 2006, New J. Phys. 8 278
    // The scaling factors in this paper roughly match the data in Figure 4 of
    // Chou for stripping in gas or foil, even though the paper is about fully
    // stripped ions

    // We scale the velocity down by the ratio of the sqrt(Z_proton + 1) /
    // sqrt(Z + 1) = sqrt(2 / Z+1) And then scale up the cross section by
    // Z*(atoms/volume) which is numerically equal to the electron number
    // density

    if (verboseLevel > 1)
        G4cout << "Entering calculate cross section vector." << G4endl;

    // Get velocity of particle
    G4double velocity = aTrack.GetVelocity();

    if (verboseLevel > 1)
        G4cout << "Velocity is " << velocity / (m / s) << " m/s" << G4endl;

    // Get element vector and number density of each element
    G4Material* material = aTrack.GetMaterial();
    G4int numberOfElements = material->GetNumberOfElements();
    G4double const* atomsPerVolumeVector = material->GetVecNbOfAtomsPerVolume();

    if (verboseLevel > 1)
        G4cout << "Material is " << material->GetName() << G4endl;
    if (verboseLevel > 1)
        G4cout << "Number of Elements is " << numberOfElements << G4endl;

    // Allocate a new array to return, with cross section per element
    // (Must be deleted by caller!)
    G4double* crossSectionPerVolumeVector = new G4double[numberOfElements];

    // For each element in the material
    for (int i = 0; i < numberOfElements; i++)
    {
        if (verboseLevel > 1)
            G4cout << "Element " << i << G4endl;

        // Get atomic number of this element
        G4double Z = material->GetElement(i)->GetZ();
        if (verboseLevel > 1)
            G4cout << "Z is " << Z << G4endl;

        // Calculate the cross section per atom for this element
        // (using the sqrt(2/Z+1) velocity scaling and the p-(He-mu) cross
        // section )
        G4double muonStrippingCrossSectionPerAtom
            = Z * CalculateCrossSection(velocity * sqrt(2.0 / (Z + 1.0)));
        if (verboseLevel > 1)
            G4cout << "muonStrippingCrossSectionPerAtom is "
                   << muonStrippingCrossSectionPerAtom << G4endl;

        // Add the cross section for this element to the total cross section
        crossSectionPerVolumeVector[i] = muonStrippingCrossSectionPerAtom
                                         * atomsPerVolumeVector[i];
        if (verboseLevel > 1)
            G4cout << "crossSectionPerVolumeVector[i] is "
                   << crossSectionPerVolumeVector[i] << G4endl;
    }

    return (crossSectionPerVolumeVector);
}

// CalculateCrossSection
// ---------------------
//
G4double MuonStripping::CalculateCrossSection(G4double velocity)
{
    // Physical constants and unit conversions
    static G4double const unitLengthAU = 5.2917721092e-11 * meter;
    static G4double const unitVelocityAU = 2.1876912633e6 * meter / second;
    static double const electronMuonMassRatioSquared
        = pow(1.0 / 206.7682826, 2.0);

    // Uses the cross section from
    // Muon reactiviation in muon-catalyzed d-t fusion from accurate p-He+
    // stripping and excitation cross sections C.D. Stodden, H.J. Monkhurst, K.
    // Szalewicz, T.G. Winter Physical Review A, Volume 4, Number 3, February
    // 1, 1990, p. 1281

    // TODO: Currently this uses only the stripping, not the excitation cross
    // sections
    // TODO: Model the whole excitation / deexcitation chain
    // TODO: Charge transfer is also possible
    // TODO: confirm that this model in geant4 gives the same results for
    // reactivation probability calculated in paper

    // Convert particle velocity to A.U.
    double velocityAU = velocity / unitVelocityAU;

    // Use piecewise power law fit from paper for electron stripping cross
    // section
    double A1, A2, A3, B1, B2;

    // Data from Table IV, values in A.U. (atomic units)
    if (velocityAU < 4.0)
    {
        A1 = 4.4;
        A2 = 0.303;
        A3 = 3.88;
        B1 = 5.3;
        B2 = 4.75;
    }
    else
    {
        A1 = 0.24;
        A2 = 341.2;
        A3 = -2.4;
        B1 = 7.28;
        B2 = 1.16;
    }

    // Equation 16
    double electronStrippingCrossSectionAU
        = A1 * pow(velocityAU, 8.0)
          / ((A2 + pow(velocityAU, B1)) * (A3 + pow(velocityAU, B2)));

    // Convert to muon stripping cross section using scaling law used in paper
    // Equation 26
    double muonStrippingCrossSectionAU = electronStrippingCrossSectionAU
                                         * electronMuonMassRatioSquared;
    G4double muonStrippingCrossSection = muonStrippingCrossSectionAU
                                         * unitLengthAU * unitLengthAU;

    return (muonStrippingCrossSection);
}

// PrintCrossSectionTable
// ----------------------
//
void MuonStripping::PrintCrossSectionTable()
{
    double v;
    G4double velocity;
    G4double crossSection;
    G4double cmenergy;
    G4double electronCrossSection;

    G4cout << "-----------------------------------" << G4endl;
    G4cout << "Muon Stripping Cross Section Table:" << G4endl;
    G4cout << "Alpha/Proton Ecm (keV), velocity (m/s) , "
              "MuonStrippingCrossSection (m^2), ElectronStrippingCrossSection "
              "(1e-18 cm^2)"
           << G4endl;
    for (v = 0.0; v <= 130.0e5; v += 1.0e5)
    {
        velocity = v * meter / second;
        crossSection = CalculateCrossSection(velocity);
        cmenergy = AlphaProtonCMFrameEnergy(velocity);
        electronCrossSection = crossSection * 206.768282 * 206.768282;
        G4cout << cmenergy / (keV) << " , " << (velocity / (meter / second))
               << " , " << (crossSection / meter2) << " , "
               << (electronCrossSection * 1e18 / cm2) << G4endl;
    }
    G4cout << "End Muon Stripping Cross Section Table" << G4endl;
    G4cout << "--------------------------------------" << G4endl;
}

// For plotting the cross section table to compare with plot in paper
// nonrelativistic calculation (for plotting only)
G4double MuonStripping::AlphaProtonCMFrameEnergy(G4double v)
{
    G4double proton_mass = amu;
    G4double alpha_mass = amu * 4.0;

    G4double v_alpha = v / (1 + (alpha_mass / proton_mass));
    G4double v_proton = v - v_alpha;
    G4double Ecm = (0.5 * alpha_mass * v_alpha * v_alpha)
                   + (0.5 * proton_mass * v_proton * v_proton);

    return (Ecm);
}
