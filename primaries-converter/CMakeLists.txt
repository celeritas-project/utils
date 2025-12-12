#------------------------------- -*- cmake -*- -------------------------------#
# Copyright Celeritas contributors: see top-level COPYRIGHT file for details
# SPDX-License-Identifier: (Apache-2.0 OR MIT)
#-----------------------------------------------------------------------------#
# HepMC3/Jsonl to ROOT converter
#-----------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.18)
project(celer-geant-validation VERSION 0.0.1 LANGUAGES CXX)
cmake_policy(VERSION 3.12...3.30)

if(DEFINED CMAKE_BUILD_TYPE AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
  message(STATUS "Set default CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
endif()

#-----------------------------------------------------------------------------#
# Find packages
find_package(Celeritas 0.7 REQUIRED)
find_package(ROOT REQUIRED)
find_package(HepMC3 REQUIRED)

if(NOT Celeritas_FOUND)
  message(FATAL_ERROR
    "No Celeritas installation found.\n"
    "Export Celeritas_ROOT=/path/to/celeritas/install or use CMAKE_INSTALL_PREFIX")
endif()

if(NOT CELERITAS_USE_Geant4)
  message(SEND_ERROR "This example requires Geant4 support "
    "to be enabled in Celeritas")
endif()

#-----------------------------------------------------------------------------#
# Add executable
add_executable(convert
  primaries-converter.cc
  src/EventReader.cc
  src/JsonEventReader.cc
  src/RootEventWriter.cc
)

target_include_directories(convert PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)

celeritas_target_link_libraries(convert
  Celeritas::celeritas
  ROOT::Tree
  HepMC3::HepMC3
  nlohmann_json::nlohmann_json
)
