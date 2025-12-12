//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file JsonEventReader.cc
//---------------------------------------------------------------------------//
#include "JsonEventReader.hh"

#include <array>
#include <corecel/Assert.hh>
#include <corecel/io/Logger.hh>
#include <nlohmann/json.hpp>

//---------------------------------------------------------------------------//
/*!
 * Construct with input filename.
 */
JsonEventReader::JsonEventReader(std::string const& filename)
    : infile_(filename)
{
    CELER_EXPECT(!filename.empty());

    CELER_VALIDATE(infile_,
                   << "failed to open file '" << filename << "' for reading");

    // Count the number of events in the file
    std::string line;
    while (std::getline(infile_, line))
    {
        ++num_events_;
    }

    // Reset stream to beginning
    infile_.clear();
    infile_.seekg(0, std::ios::beg);

    CELER_LOG(info) << "Json Event Reader: found " << num_events_
                    << " event(s) in file '" << filename << "'";
}

//---------------------------------------------------------------------------//
/*!
 * Read single event from the file.
 */
auto JsonEventReader::operator()() -> result_type
{
    CELER_EXPECT(infile_);

    result_type result;

    std::string line;
    if (!std::getline(infile_, line))
    {
        // No more events
        return result;
    }
    auto event = nlohmann::json::parse(line);

    auto event_id = event.at("event_id").get<unsigned int>();
    for (auto const& j : event.at("primaries"))
    {
        Primary p;
        p.event_id = event_id;
        p.pdg = j.at("pdg").get<int>();
        p.energy = j.at("energy")[0].get<double>();  // [value, string_unit]
        p.position = j.at("position").get<std::array<double, 3>>();
        p.direction = j.at("direction").get<std::array<double, 3>>();
        p.time = j.at("time").get<double>();

        result.push_back(std::move(p));
    }

    CELER_ENSURE(!result.empty());
    return result;
}
