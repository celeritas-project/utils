//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file src/HistogramStore.cc
//---------------------------------------------------------------------------//
#include "HistogramStore.hh"

#include <corecel/Assert.hh>

//---------------------------------------------------------------------------//
/*!
 * Add physical volume ID and copy number to the map and initialize histograms
 * associated to it.
 */
void HistogramStore::InsertSensDet(PhysVolId pid,
                                   CopyNumber cid,
                                   std::string name)
{
    sensdet_map_.insert({{pid, cid}, SDHistograms::Initialize(name)});
}

//---------------------------------------------------------------------------//
/*!
 * Return histogram data for a given sensitive detector.
 */
SDHistograms& HistogramStore::Find(PhysVolId pv_id, CopyNumber copy_num)
{
    auto iter = sensdet_map_.find(SensDetId{pv_id, copy_num});
    CELER_ASSERT(iter != sensdet_map_.end());
    return iter->second;
}
