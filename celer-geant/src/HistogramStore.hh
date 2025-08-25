//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file src/HistogramStore.hh
//---------------------------------------------------------------------------//
#pragma once

#include <string>
#include <TH1D.h>
#include <corecel/Assert.hh>

#include "JsonReader.hh"

//---------------------------------------------------------------------------//
/*!
 * Histogram storage for sensitive detectors.
 */
struct SDHistograms
{
    std::string sd_name;  //!< Sensitive detector name

    //!@{
    //! Histograms
    TH1D energy;
    TH1D time;
    TH1D pos_x;
    //!@}

    //! Initialize histograms using the SD name and JSON input data
    static SDHistograms Initialize(std::string sd_name)
    {
        struct HistDef
        {
            size_t nbins;
            int min;
            int max;
        };
        auto const& j = JsonReader::Instance().at("histograms");
        auto histogram_def = [&](std::string name) -> HistDef {
            auto jx = j.at(name);
            HistDef h;
            h.nbins = jx.at("num_bins").get<size_t>();
            h.min = jx.at("min").get<int>();
            h.max = jx.at("max").get<int>();
            return h;
        };

#define SDH_INIT_TH1D(HIST)                                                 \
    {                                                                       \
        std::string hname = sd_name + "_" + #HIST;                          \
        auto const hd = histogram_def(#HIST);                               \
        result.HIST                                                         \
            = TH1D(hname.c_str(), hname.c_str(), hd.nbins, hd.min, hd.max); \
    }

        SDHistograms result;
        result.sd_name = sd_name;
        SDH_INIT_TH1D(energy);
        SDH_INIT_TH1D(time);
        SDH_INIT_TH1D(pos_x);
        return result;

#undef SDH_INIT_TH1D
    }
};

//---------------------------------------------------------------------------//
/*!
 * Helper struct for indexing physical volumes to an object.
 * (e.g. std::map<SensDetId, Object> map)
 */
struct SensDetId
{
    size_t physvol_id;
    size_t copy_number;
};

//---------------------------------------------------------------------------//
//! Overload \c operator< for \c map.find(sensdetid) compatibility.
inline bool operator<(SensDetId const& lhs, SensDetId const& rhs)
{
    return std::make_tuple(lhs.physvol_id, lhs.copy_number)
           < std::make_tuple(rhs.physvol_id, rhs.copy_number);
}

//---------------------------------------------------------------------------//
/*!
 * Histogram storage manager class.
 *
 * This class stores a set of histograms for every sensitive detector in the
 * geometry and allows an easy way to access them using the physical volume
 * instance ID and copy number.
 */
class HistogramStore
{
  public:
    //!@{
    //! \name Type aliases
    using PhysVolId = size_t;
    using CopyNumber = size_t;
    //!@}

    //! Construct empty
    HistogramStore() = default;

    //! Map and initialize histograms for a sensitive detector
    void InsertSensDet(PhysVolId pv_id, CopyNumber copy_num, std::string name);

    //! Get histogram data for a given physical volume ID and copy number
    SDHistograms& Find(PhysVolId pv_id, CopyNumber copy_num);

    //! Access full SD map
    std::map<SensDetId, SDHistograms>& Map() { return sensdet_map_; }

  private:
    std::map<SensDetId, SDHistograms> sensdet_map_;
};
