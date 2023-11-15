#---------------------------------*-CMake-*----------------------------------#
# Copyright 2023 UT-Battelle, LLC, and other Celeritas developers.
# See the Celeritas top-level COPYRIGHT file for details.
# SPDX-License-Identifier: (Apache-2.0 OR MIT)
#----------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.9)
project(g4-optics)

cmake_policy(VERSION 3.12...3.16)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE "Debug")

#----------------------------------------------------------------------------#
# Find packages
#----------------------------------------------------------------------------#
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------#
# Copy vis.mac to build directory
#----------------------------------------------------------------------------#
configure_file(${PROJECT_SOURCE_DIR}/vis.mac ${PROJECT_BINARY_DIR} COPYONLY)

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
