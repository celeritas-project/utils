//----------------------------------*-C++-*----------------------------------//
// Copyright 2021 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file GeometryStore.hh
//---------------------------------------------------------------------------//
#include "GeometryStore.hh"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>

//---------------------------------------------------------------------------//
/*!
 * Construct.
 */
GeometryStore::GeometryStore()
    : phys_vol_store_(G4PhysicalVolumeStore::GetInstance())
{
    std::cout << phys_vol_store_->size() << std::endl;
    this->loop_volumes();
}

//---------------------------------------------------------------------------//
/*!
 * Get volumes vector.
 */
std::vector<Volume> GeometryStore::get_volumes() const
{
    return volumes_;
}

//---------------------------------------------------------------------------//
/*!
 * Save data stored in this->ids_volumes map as a text output file.
 */
void GeometryStore::save(const std::string filename)
{
    std::ofstream output;
    output.open(filename);
    output << volumes_ << std::endl;
    output.close();
}

//---------------------------------------------------------------------------//
// PRIVATE
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * Recursive loop to store all logical volumes.
 */
void GeometryStore::loop_volumes()
{
    for (const auto& phys_vol : *phys_vol_store_)
    {
        const auto& logical_volume = phys_vol->GetLogicalVolume();

        Volume volume;
        volume.name          = phys_vol->GetName();
        volume.volume_id     = phys_vol->GetInstanceID();
        volume.material_id   = logical_volume->GetMaterial()->GetIndex();
        volume.material_name = logical_volume->GetMaterial()->GetName();
        volume.copy_num      = phys_vol->GetCopyNo();

        EAxis  dummy_a;
        double dummy_b, dummy_c;
        bool   dummy_d;
        phys_vol->GetReplicationData(
            dummy_a, volume.num_replicas, dummy_b, dummy_c, dummy_d);
        volumes_.push_back(volume);
    }
}

//---------------------------------------------------------------------------//
// FREE FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * Define operator << to print a full table with the GeoTestMap data.
 */
std::ostream& operator<<(std::ostream& os, std::vector<Volume> list)
{
    size_t width_ids      = 6;
    size_t width_volume   = 0;
    size_t width_material = 0;

    for (const auto& it : list)
    {
        width_volume   = std::max(width_volume, it.name.size());
        width_material = std::max(width_material, it.material_name.size());
    }

    // Title
    os << std::endl;
    os << "| " << std::left << std::setw(width_ids) << "Vol ID"
       << " | " << std::left << std::setw(width_ids) << "Copy"
       << " | " << std::left << std::setw(width_ids) << "Repl"
       << " | " << std::left << std::setw(width_ids) << "Mat ID"
       << " | " << std::setw(width_material) << "Material"
       << " | " << std::setw(width_volume) << "Volume"
       << " |" << std::endl;

    // Dashed line
    os << "| ";

    for (int i = 0; i < width_ids; i++)
        os << "-";
    os << " | ";
    for (int i = 0; i < width_ids; i++)
        os << "-";
    os << " | ";
    for (int i = 0; i < width_ids; i++)
        os << "-";
    os << " | ";
    for (int i = 0; i < width_ids; i++)
        os << "-";
    os << " | ";
    for (int i = 0; i < width_material; i++)
        os << "-";
    os << " | ";
    for (int i = 0; i < width_volume; i++)
        os << "-";

    os << " | ";
    os << std::endl;

    // Table content
    for (const auto& key : list)
    {
        os << "| " << std::left << std::setw(width_ids) << key.volume_id
           << " | " << std::left << std::setw(width_ids) << key.copy_num
           << " | " << std::left << std::setw(width_ids) << key.num_replicas
           << " | " << std::left << std::setw(width_ids) << key.material_id
           << " | " << std::setw(width_material) << key.material_name << " | "
           << std::setw(width_volume) << key.name << " |" << std::endl;
    }

    return os;
}
