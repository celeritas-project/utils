//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file EventReader.hh
//---------------------------------------------------------------------------//
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Primary.hh"

namespace HepMC3
{
class Reader;
}

//---------------------------------------------------------------------------//
/*!
 * Read a HepMC3 event record file and create primary particles.
 */
class EventReader
{
  public:
    //!@{
    //! \name Type aliases
    using result_type = std::vector<Primary>;
    //!@}

  public:
    // Construct from a filename
    EventReader(std::string const& filename);
    ~EventReader() = default;

    // Read a single event from the event record
    result_type operator()();

    //! Get total number of events
    unsigned int num_events() const { return num_events_; }

  private:
    using SPReader = std::shared_ptr<HepMC3::Reader>;

    // HepMC3 event record reader
    SPReader reader_;

    // Total number of events in file
    unsigned int num_events_;
};

//---------------------------------------------------------------------------//
// Wrapper function for HepMC3::deduce_reader to avoid duplicate symbols
std::shared_ptr<HepMC3::Reader> open_hepmc3(std::string const& filename);
