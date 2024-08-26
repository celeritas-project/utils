//----------------------------------*-C++-*----------------------------------//
// Copyright 2021-2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file SteppingAction.cc
//---------------------------------------------------------------------------//
#include "SteppingAction.hh"

#include <G4SystemOfUnits.hh>
#include <G4VProcess.hh>

#include "JsonReader.hh"
#include "RootData.hh"
#include "RootIO.hh"

//---------------------------------------------------------------------------//
/*!
 * Construct and set up I/O options.
 */
SteppingAction::SteppingAction()
    : G4UserSteppingAction(), root_io_(RootIO::instance())
{
    auto const& json_sim = JsonReader::instance()->json().at("simulation");
    store_step_ = json_sim.at("step_info").get<bool>();
    store_primary_ = json_sim.at("primary_info").get<bool>();
    store_secondary_ = json_sim.at("secondary_info").get<bool>();
}

//---------------------------------------------------------------------------//
/*!
 * Fetch data at every new step and populate Event object.
 */
void SteppingAction::UserSteppingAction(G4Step const* step)
{
    if (!root_io_)
    {
        return;
    }

    auto const parent_id = step->GetTrack()->GetParentID();

    if (store_primary_ && parent_id == 0)
    {
        store_track_data(step);
    }

    if (store_secondary_ && parent_id != 0)
    {
        store_track_data(step);
    }
}

//---------------------------------------------------------------------------//
// PRIVATE
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * Store track data.
 */
void SteppingAction::store_track_data(G4Step const* step)
{
    root_io_->track_.energy_dep += step->GetTotalEnergyDeposit() / MeV;
    root_io_->track_.number_of_steps++;

    if (store_step_)
    {
        this->store_step_data(step);
    }
}

//---------------------------------------------------------------------------//
/*!
 * Populate step information in RootIO::track_.
 */
void SteppingAction::store_step_data(G4Step const* step)
{
    rootdata::Step this_step;

    auto const& post_step = step->GetPostStepPoint();
    if (post_step->GetStepStatus() == fUndefined)
    {
        // Post step status is undefined; GetProcessDefinedStep() is a nullptr
        // Only caused by geantinos
        this_step.process_id = rootdata::ProcessId::not_mapped;
    }
    else
    {
        // Post step is defined; find its ID
        this_step.process_id = rootdata::to_process_name_id(
            step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());
    }

    this_step.kinetic_energy = post_step->GetKineticEnergy() / MeV;
    this_step.energy_loss = step->GetTotalEnergyDeposit() / MeV;
    this_step.length = step->GetStepLength() / cm;
    this_step.global_time = post_step->GetGlobalTime() / s;
    auto pos = post_step->GetPosition() / cm;
    auto dir = post_step->GetMomentumDirection();
    auto pol = post_step->GetPolarization();
    this_step.position = {pos.x(), pos.y(), pos.z()};
    this_step.direction = {dir.x(), dir.y(), dir.z()};
    this_step.polarization = {pol.x(), pol.y(), pol.z()};

    root_io_->track_.steps.push_back(std::move(this_step));

    root_io_->data_limits_.max_time
        = std::max(root_io_->data_limits_.max_time, this_step.global_time);
    root_io_->data_limits_.max_length
        = std::max(root_io_->data_limits_.max_length, this_step.length);
}
