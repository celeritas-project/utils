//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file EventReader.cc
//---------------------------------------------------------------------------//
#include "EventReader.hh"

#include <cmath>
#include <HepMC3/GenEvent.h>
#include <HepMC3/Reader.h>
#include <HepMC3/ReaderFactory.h>
#include <HepMC3/Setup.h>
#include <celeritas/Constants.hh>
#include <corecel/Assert.hh>
#include <corecel/io/Logger.hh>

//---------------------------------------------------------------------------//
/*
 * Make unit vector from momentum.
 */
CELER_FUNCTION std::array<double, 3>
make_unit_vector(std::array<double, 3> const& v)
{
    using Array3 = std::array<double, 3>;

    auto norm = [](Array3 v) -> double {
        auto dot_product = [](Array3 x, Array3 y) -> double {
            double result{};
            for (int i = 0; i != 3; ++i)
            {
                result = std::fma(x[i], y[i], result);
            }
            return result;
        };

        return std::sqrt(dot_product(v, v));
    };

    Array3 result{v};
    double const scale_factor = 1. / norm(result);
    for (auto& el : result)
    {
        el *= scale_factor;
    }
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Construct from a filename.
 */
EventReader::EventReader(std::string const& filename)
{
    // Fetch total number of events by opening a temporary reader
    num_events_ = [&filename] {
        SPReader temp_reader = open_hepmc3(filename);
        CELER_ASSERT(temp_reader);
        unsigned int result = 0;
#if HEPMC3_VERSION_CODE < 3002000
        HepMC3::GenEvent evt;
        temp_reader->read_event(evt);
#else
        temp_reader->skip(0);
#endif
        CELER_VALIDATE(!temp_reader->failed(),
                       << "event file '" << filename
                       << "' did not contain any events");
        do
        {
            result++;
#if HEPMC3_VERSION_CODE < 3002000
            temp_reader->read_event(evt);
#else
            temp_reader->skip(1);
#endif
        } while (!temp_reader->failed());
        CELER_LOG(debug) << "HepMC3 file has " << result << " events";
        return result;
    }();

    CELER_LOG(info) << "HepMC3 Event Reader: found " << num_events_
                    << " event(s) in file '" << filename << "'";

    // Determine the input file format and construct the appropriate reader
    reader_ = open_hepmc3(filename);
    CELER_ENSURE(reader_);
}

//---------------------------------------------------------------------------//
/*!
 * Read a single event from the event record.
 *
 * \note
 * Units are converted manually from HepMC3 to Celeritas rather than using \c
 * GenEvent::set_units(). This method converts the momentum units of all
 * particles and the length units of all daughter vertices; however, it does
 * *not* convert the length units of the root vertex. This means any particle
 * without a vertex will not have the correct units for position and time. This
 * may also be true for any vertices that do not have a position assigned, as
 * their positions will be inherited from the vertices of their ancestors (or
 * by falling back on the event position).
 */
auto EventReader::operator()() -> result_type
{
    // Parse the next event from the record
    HepMC3::GenEvent evt;
    {
        reader_->read_event(evt);
    }
    // There are no more events
    if (reader_->failed())
    {
        return {};
    }

    auto const event_id = evt.event_number();
    CELER_LOG(debug) << "Reading event " << event_id;

    result_type result;
    for (auto const& par : evt.particles())
    {
        if (par->data().status != 1 || par->end_vertex())
        {
            // Skip particles that aren't leaves on the tree of generated
            // particles
            continue;
        }

        // Get the PDG code as particle id
        Primary primary;

        // Set the registered ID of the particle
        primary.pdg = par->pid();

        // Set the event and track number
        primary.event_id = event_id;

        // Get the position of the vertex
        auto pos = par->production_vertex()->position();
        HepMC3::Units::convert(pos, evt.length_unit(), HepMC3::Units::CM);
        auto to_cm = [](double v) {
            return static_cast<double>(v) * celeritas::units::centimeter;
        };
        primary.position = {to_cm(pos.x()), to_cm(pos.y()), to_cm(pos.z())};

        // Get the lab-frame time [time]
        primary.time = to_cm(pos.t()) / celeritas::constants::c_light;

        // Get the direction of the primary
        auto mom = par->momentum();
        HepMC3::Units::convert(mom, evt.momentum_unit(), HepMC3::Units::MEV);
        primary.direction = make_unit_vector({mom.px(), mom.py(), mom.pz()});

        // Get the kinetic energy of the primary
        primary.energy = mom.e() - mom.m();

        result.push_back(primary);
    }

    CELER_VALIDATE(!result.empty(),
                   << "event " << event_id
                   << " did not contain any primaries suitable for "
                      "simulation");

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Wrapper function for HepMC3::deduce_reader to avoid duplicate symbols.
 *
 * HepMC3 through 3.2.6 has a ReaderFactory.h that includes function
 * *definitions* without \c inline keywords, leading to duplicate symbols.
 * Reusing this function rather than including ReaderFactory multiple times in
 * Celeritas is the easiest way to work around the problem.
 *
 * It also sets the debug level from the environment, prints a status
 * message,and validates the file.
 */
std::shared_ptr<HepMC3::Reader> open_hepmc3(std::string const& filename)
{
    CELER_LOG(info) << "Opening HepMC3 input file at " << filename;

    auto result = HepMC3::deduce_reader(filename);
    CELER_VALIDATE(result, << "failed to deduce event input file type");
    return result;
}
