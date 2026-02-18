//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/Stopwatch.cc
//---------------------------------------------------------------------------//
#include "Stopwatch.hh"

//---------------------------------------------------------------------------//
/*!
 * Construct empty.
 */
Stopwatch::Stopwatch() {}

//---------------------------------------------------------------------------//
/*!
 * Start stopwatch.
 */
void Stopwatch::start()
{
    cpu_start_ = std::clock();  // [ms]
    wall_start_ = std::chrono::high_resolution_clock::now();
}

//---------------------------------------------------------------------------//
/*!
 * Stop stopwatch.
 */
void Stopwatch::stop()
{
    cpu_stop_ = std::clock();  // [ms]
    wall_stop_ = std::chrono::high_resolution_clock::now();
}

//---------------------------------------------------------------------------//
/*!
 * Return CPU execution time duration in [s].
 */
double Stopwatch::cpu() const
{
    return (cpu_stop_ - cpu_start_) / CLOCKS_PER_SEC;
}

//---------------------------------------------------------------------------//
/*!
 * Return wall execution time duration in [s].
 */
double Stopwatch::wall() const
{
    return static_cast<double>(
               std::chrono::duration_cast<std::chrono::microseconds>(
                   wall_stop_ - wall_start_)
                   .count())
           / CLOCKS_PER_SEC;
}
