//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/io/RootEventWriter.cc
//---------------------------------------------------------------------------//
#include "RootEventWriter.hh"

#include <TFile.h>
#include <TTree.h>

#include "corecel/Assert.hh"
#include "corecel/io/Logger.hh"
#include "celeritas/ext/ScopedRootErrorHandler.hh"

//---------------------------------------------------------------------------//
/*!
 * Construct with ROOT output filename.
 */
RootEventWriter::RootEventWriter(SPRootFileManager root_file_manager)
    : tfile_mgr_(std::move(root_file_manager)), event_id_(-1)
{
    CELER_EXPECT(tfile_mgr_);

    celeritas::ScopedRootErrorHandler scoped_root_error;

    CELER_LOG(info) << "Creating event tree '" << this->tree_name() << "' at "
                    << tfile_mgr_->filename();

    ttree_ = tfile_mgr_->make_tree(this->tree_name(), this->tree_name());
    ttree_->Branch("event_id", &primary_.event_id);
    ttree_->Branch("particle", &primary_.pdg);
    ttree_->Branch("energy", &primary_.energy);
    ttree_->Branch("time", &primary_.time);
    ttree_->Branch("pos", &primary_.position);
    ttree_->Branch("dir", &primary_.direction);
    scoped_root_error.throw_if_errors();
}

//---------------------------------------------------------------------------//
/*!
 * Export primaries to ROOT.
 */
void RootEventWriter::operator()(VecPrimary const& primaries)
{
    CELER_EXPECT(!primaries.empty());
    celeritas::ScopedRootErrorHandler scoped_root_error;

    // Increment contiguous event id
    event_id_++;

    for (auto const& p : primaries)
    {
        primary_.event_id = event_id_;
        primary_.pdg = p.pdg;
        primary_.energy = p.energy;
        primary_.time = p.time;
        primary_.position = p.position;
        primary_.direction = p.direction;
        ttree_->Fill();
    }

    scoped_root_error.throw_if_errors();
}
