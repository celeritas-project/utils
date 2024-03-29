//----------------------------------*-C++-*----------------------------------//
// Copyright 2022 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file GeometryStore.hh
//! \brief Store detector geometry.
//---------------------------------------------------------------------------//
#pragma once

#include <iostream>
#include <map>
#include <string>
#include <G4PhysicalVolumeStore.hh>

//---------------------------------------------------------------------------//
/*!
 * Store volume information.
 */
struct Volume
{
    int         logical_volume_id;
    int         physical_volume_id;
    int         material_id;
    std::string material_name;
    std::string physical_volume_name;
    std::string logical_volume_name;
    int         copy_num;
    int         num_replicas;
};

//---------------------------------------------------------------------------//
/*!
 * Map detector geometry for comparison purposes.
 */
class GeometryStore
{
  public:
    //!@{
    //! Type aliases
    //!@}

    // Constructor
    GeometryStore();
    // Default destructor
    ~GeometryStore() = default;

    // Get constructed map
    std::vector<Volume> get_volumes() const;

    // Save a text output file with the data loaded in this->ids_volumes_
    void save(const std::string filename);

    // Verify if volume ID numbering scheme has any discontinuity
    bool continuous_volume_ids();

  private:
    // Recursive loop over logical volumes
    void loop_volumes();

  private:
    // Map volume id with volume information
    G4PhysicalVolumeStore* phys_vol_store_;
    std::vector<Volume>    volumes_;
};

//---------------------------------------------------------------------------//
// Free functions
//---------------------------------------------------------------------------//

// Define operator << to print a full table with GeoTestMap data.
std::ostream& operator<<(std::ostream& os, std::vector<Volume> map);
