//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file Primary.hh
//---------------------------------------------------------------------------//
#pragma once

#include <array>

//---------------------------------------------------------------------------//
/*!
 * Source particle definition. Equivalent to Celeritas \c Primary , but without
 * specific types such as OpaqueIds.
 */
struct Primary
{
    unsigned int event_id{};
    int pdg{};
    double energy{};
    std::array<double, 3> position{0, 0, 0};
    std::array<double, 3> direction{0, 0, 0};
    double time{};
};
