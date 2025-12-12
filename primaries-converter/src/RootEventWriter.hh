//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/io/RootEventWriter.hh
//---------------------------------------------------------------------------//
#pragma once

#include <celeritas/ext/RootFileManager.hh>

#include "EventReader.hh"
#include "Primary.hh"

//---------------------------------------------------------------------------//
/*!
 * Export primary data to ROOT.
 *
 * One TTree entry represents one primary.
 */
class RootEventWriter
{
  public:
    //!@{
    //! \name Type aliases
    using SPRootFileManager = std::shared_ptr<celeritas::RootFileManager>;
    using VecPrimary = std::vector<Primary>;
    //!@}

    // Construct with ROOT output filename
    RootEventWriter(SPRootFileManager root_file_manager);

    //! Prevent copying and moving
    ~RootEventWriter() = default;

    // Export primaries to ROOT
    void operator()(VecPrimary const& primaries);

  private:
    SPRootFileManager tfile_mgr_;
    unsigned int event_id_;  // Contiguous event id
    celeritas::UPRootTreeWritable ttree_;
    Primary primary_;  // Temporary object stored to the ROOT TTree
    bool warned_mismatched_events_{false};

    //// HELPER FUNCTIONS ////

    // Hardcoded TTree name and title
    char const* tree_name() { return "primaries"; }
};
