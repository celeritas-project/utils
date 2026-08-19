//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/RootOpticalIO.hh
//---------------------------------------------------------------------------//
#pragma once

#include <G4TransportationManager.hh>
#include <G4VPhysicalVolume.hh>
#include <accel/SetupOptions.hh>
#include <celeritas/optical/DetectorData.hh>
#include <corecel/Assert.hh>
#include <corecel/Types.hh>
#include <corecel/math/Quantity.hh>
#include <geocel/g4/Convert.hh>

#include "RootIO.hh"

namespace
{
//---------------------------------------------------------------------------//
/*!
 * Convert MeV to nm using the relation
 * \f[
 *  \lambda_{\text{nm}} = \frac{hc}{E_{\text{MeV}}} .
 * \f]
 */
inline double mev_to_nm(double photon_energy_mev)
{
    static double const hplanck_clight = 1239.8 * 1e-6;  // [MeV * nm]
    return (photon_energy_mev > 0) ? (hplanck_clight / photon_energy_mev * 1e6)
                                   : -1.0;
}

//---------------------------------------------------------------------------//
/*!
 * Get the Geant4 physical volume from a given (x, y, z) point in the world
 * coordinate system via the G4TransportationManager.
 */
inline G4VPhysicalVolume* physvol_from_coordinate(G4ThreeVector const pos)
{
    auto* tm = G4TransportationManager::GetTransportationManager();
    CELER_ASSERT(tm);
    auto* nav = tm->GetNavigatorForTracking();
    CELER_ASSERT(nav);
    auto* result = nav->LocateGlobalPointAndSetup(pos);
    CELER_ENSURE(result);
    return result;
}

//---------------------------------------------------------------------------//
}  // namespace

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Callback function for Celeritas optical detector hits.
 */
inline void optical_hits_callback(
    celeritas::Span<celeritas::optical::DetectorHit const> hits)
{
    using celeritas::native_to_geant;
    using celeritas::real_type;
    using celeritas::units::MevEnergy;

#define OHC_1D_FILL(MEMBER, VALUE) data.MEMBER.Fill(VALUE);
#define OHC_2D_FILL(MEMBER, X, Y) data.MEMBER.Fill(X, Y);
#define OHC_1D_FILL_WEIGHT(MEMBER, VALUE, WEIGHT)        \
    {                                                    \
        auto& h = data.MEMBER;                           \
        auto const i = h.FindBin(VALUE);                 \
        h.SetBinContent(i, h.GetBinContent(i) + WEIGHT); \
    }

    auto rio = RootIO::Instance();
    for (auto const& hit : hits)
    {
        // Locate volume from hit
        auto* phys_vol = physvol_from_coordinate(
            native_to_geant<lengthunits::ClhepLength>(hit.position));
        auto& data = rio->Data().Find(phys_vol->GetInstanceID(),
                                      phys_vol->GetCopyNo());

        // Fill histograms
        OHC_1D_FILL_WEIGHT(energy_dep_x, hit.position[0], hit.energy.value());
        OHC_1D_FILL_WEIGHT(energy_dep_y, hit.position[1], hit.energy.value());
        OHC_1D_FILL_WEIGHT(energy_dep_z, hit.position[2], hit.energy.value());
        OHC_2D_FILL(pos_xy, hit.position[0], hit.position[1]);
        OHC_1D_FILL(time, hit.time);
    }

#undef OHC_1D_FILL
#undef OHC_2D_FILL
#undef OHC_1D_FILL_WEIGHT
}

//---------------------------------------------------------------------------//
}  // namespace celeritas
