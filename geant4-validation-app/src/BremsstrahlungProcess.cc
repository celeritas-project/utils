//----------------------------------*-C++-*----------------------------------//
// Copyright 2021-2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file BremsstrahlungProcess.cc
//---------------------------------------------------------------------------//
#include "BremsstrahlungProcess.hh"

#include <G4Electron.hh>
#include <G4EmParameters.hh>
#include <G4Gamma.hh>
#include <G4ParticleDefinition.hh>
#include <G4Positron.hh>
#include <G4SeltzerBergerModel.hh>
#include <G4SystemOfUnits.hh>
#include <G4eBremsstrahlungRelModel.hh>

//---------------------------------------------------------------------------//
/*!
 * Construct with model selection.
 */
BremsstrahlungProcess::BremsstrahlungProcess(ModelSelection selection)
    : is_initialized_(false)
    , G4VEnergyLossProcess("eBrem")
    , model_selection_(selection)
{
    G4VEnergyLossProcess::SetProcessSubType(
        G4EmProcessSubType::fBremsstrahlung);
    G4VEnergyLossProcess::SetSecondaryParticle(G4Gamma::Gamma());
    G4VEnergyLossProcess::SetIonisation(false);
}

//---------------------------------------------------------------------------//
/*!
 * Empty destructor.
 */
BremsstrahlungProcess::~BremsstrahlungProcess() {}

//---------------------------------------------------------------------------//
/*!
 * Define applicability based on particle definition.
 */
bool BremsstrahlungProcess::IsApplicable(G4ParticleDefinition const& particle)
{
    return (&particle == G4Electron::Electron()
            || &particle == G4Positron::Positron());
}

//---------------------------------------------------------------------------//
/*!
 * Print documentation in html format.
 */
void BremsstrahlungProcess::ProcessDescription(std::ostream& output) const
{
    output << "  Bremsstrahlung";
    G4VEnergyLossProcess::ProcessDescription(output);
}

//---------------------------------------------------------------------------//
// PROTECTED
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * Initialise process by constructing models based on \c ModelSelection .
 */
void BremsstrahlungProcess::InitialiseEnergyLossProcess(
    G4ParticleDefinition const*, G4ParticleDefinition const*)
{
    if (is_initialized_)
    {
        // Nothing to do
        return;
    }

    auto const& em_parameters = G4EmParameters::Instance();

    double energy_min = em_parameters->MinKinEnergy();
    double energy_max = em_parameters->MaxKinEnergy();
    double sb_energy_limit = 1 * GeV;
    double energy_limit = std::min(energy_max, sb_energy_limit);
    G4VEmFluctuationModel* fluctuation_model = nullptr;

    std::size_t model_index = 0;

    if (model_selection_ == ModelSelection::seltzer_berger
        || model_selection_ == ModelSelection::all)
    {
        if (!G4VEnergyLossProcess::EmModel(model_index))
        {
            G4VEnergyLossProcess::SetEmModel(new G4SeltzerBergerModel());
        }

        auto em_model = G4VEnergyLossProcess::EmModel(model_index);
        em_model->SetLowEnergyLimit(energy_min);
        em_model->SetHighEnergyLimit(energy_limit);
        em_model->SetSecondaryThreshold(em_parameters->BremsstrahlungTh());
        em_model->SetLPMFlag(false);
        G4VEnergyLossProcess::AddEmModel(1, em_model, fluctuation_model);

        model_index++;
    }

    if (model_selection_ == ModelSelection::relativistic
        || model_selection_ == ModelSelection::all)
    {
        if (energy_max > energy_limit)
        {
            if (!G4VEnergyLossProcess::EmModel(model_index))
            {
                G4VEnergyLossProcess::SetEmModel(
                    new G4eBremsstrahlungRelModel());
            }

            auto em_model = G4VEnergyLossProcess::EmModel(model_index);
            em_model->SetLowEnergyLimit(energy_limit);
            em_model->SetHighEnergyLimit(energy_max);
            em_model->SetSecondaryThreshold(em_parameters->BremsstrahlungTh());
            em_model->SetLPMFlag(em_parameters->LPM());
            G4VEnergyLossProcess::AddEmModel(1, em_model, fluctuation_model);
        }
    }

    is_initialized_ = true;
}

//---------------------------------------------------------------------------//
/*!
 * Print class parameters.
 */
void BremsstrahlungProcess::StreamProcessInfo(std::ostream& output) const
{
    if (EmModel(0))
    {
        auto const& param = G4EmParameters::Instance();
        double const energy_threshold = param->BremsstrahlungTh();

        output << "      LPM flag: " << param->LPM() << " for E > "
               << EmModel(0)->HighEnergyLimit() / GeV << " GeV";

        if (energy_threshold < std::numeric_limits<double>::max())
        {
            output << ",  VertexHighEnergyTh(GeV)= " << energy_threshold / GeV;
        }
        output << std::endl;
    }
}
