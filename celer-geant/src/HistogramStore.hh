//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file src/HistogramStore.hh
//---------------------------------------------------------------------------//
#pragma once

#include <string>
#include <TH1D.h>

//---------------------------------------------------------------------------//
/*!
 * Histogram list.
 */
struct SDHistograms
{
    std::string sd_name;
    TH1D energy;
    TH1D time;

    // Initialize histograms TODO: with input data
    static SDHistograms Initialize(std::string name)
    {
        SDHistograms result;
        result.sd_name = name;
        std::string energy_name = name + "_energy";
        std::string time_name = name + "_time";
        result.energy
            = TH1D(energy_name.c_str(), energy_name.c_str(), 100, 0, 100);
        result.time
            = TH1D(time_name.c_str(), time_name.c_str(), 100, -100, 100);
        return result;
    }
};

//---------------------------------------------------------------------------//
/*!
 * Helper struct used to index physical volumes to a vector index.
 * (std::map<SensDetId, index> map)
 */
struct SensDetId
{
    size_t physvol_id;
    size_t copy_number;
};

inline bool operator<(SensDetId const& lhs, SensDetId const& rhs)
{
    return std::make_tuple(lhs.physvol_id, lhs.copy_number)
           < std::make_tuple(rhs.physvol_id, rhs.copy_number);
}

//---------------------------------------------------------------------------//
/*!
 * Store and initialize histogram for faster MC truth output results.
 */
class HistogramStore
{
  public:
    //!@{
    //! \name Type aliases
    using PhysVolId = size_t;
    using CopyNumber = size_t;
    //!@}

    // Construct empty
    HistogramStore() = default;

    // Insert new block of histograms into memory
    void InsertSensDet(PhysVolId pv_id, CopyNumber copy_num, std::string name);

    // Get histogram data for a given physical volume ID and copy number
    SDHistograms& Find(PhysVolId pv_id, CopyNumber copy_num);

    // Access full SD map
    std::map<SensDetId, SDHistograms>& Map() { return sensdet_map_; }

  private:
    std::map<SensDetId, SDHistograms> sensdet_map_;
};
