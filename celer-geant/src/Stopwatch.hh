//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celer-geant/src/Stopwatch.hh
//---------------------------------------------------------------------------//
#pragma once

#include <chrono>
#include <sys/time.h>

//---------------------------------------------------------------------------//
/*!
 * Stopwatch helper class.
 *
 * Measure wall and CPU times during execution.
 * \code
 * Stopwatch stopwatch;
 * stopwatch.start();
 * // Do stuff
 * stopwatch.stop();
 * auto cpu_time = stopwatch.duration_cpu();
 * auto wall_time = stopwatch.duration_wall();
 * \endcode
 */
class Stopwatch
{
  public:
    // Construct empty
    Stopwatch();
    ~Stopwatch() = default;

    // Start timer
    void start();

    // Stop timer
    void stop();

    // Return CPU execution time duration in [s]
    double cpu() const;

    // Return wall execution time duration in [s]
    double wall() const;

  private:
    using WallTime = std::chrono::high_resolution_clock::time_point;

    double cpu_start_;
    double cpu_stop_;
    WallTime wall_start_;
    WallTime wall_stop_;
};
