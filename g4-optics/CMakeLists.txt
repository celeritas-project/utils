#---------------------------------*-CMake-*----------------------------------#
# Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
# See the Celeritas top-level COPYRIGHT file for details.
# SPDX-License-Identifier: (Apache-2.0 OR MIT)
#----------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9)
project(g4-optics)

cmake_policy(VERSION 3.12...3.16)
set(CMAKE_CXX_STANDARD 17)

#----------------------------------------------------------------------------#
# Find packages
#----------------------------------------------------------------------------#
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------#
# Add executable
#----------------------------------------------------------------------------#
set(sources
  src/ActionInitialization.cc
  src/OpticalPhysics.cc
  src/OpticalDetector.cc
  src/PrimaryGenerator.cc
)

add_executable(main main.cc ${sources})

target_include_directories(main PRIVATE
 ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(main PRIVATE ${Geant4_LIBRARIES})
