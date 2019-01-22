#
# Copyright (C) 2015   Michael Lahnert
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

project (NUMSIM)

cmake_minimum_required (VERSION 2.8)

# Force to set a build type
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
    "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Set compiler
set (CMAKE_CC_COMPILER "/usr/bin/gcc")
set (CMAKE_CXX_COMPILER "/usr/bin/g++")

# Set some good compiler flags
set (GOOD_FLAGS "-Wall -Wextra -pedantic -std=c++11")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GOOD_FLAGS}")

# Use debug visualisation?
option (DEBUG_VISU
  "Enable debug visualization" ON)

# Use additional VTK visualisation?
option (VTK_VISU
  "Enable VTK visualization" ON)

# Use debug print for grid u?
option (USE_DEBUG_PRINT_U
  "Enable debug print for grid u" OFF)
# Use debug print for grid v?
option (USE_DEBUG_PRINT_V
  "Enable debug print for grid v" OFF)
# Use debug print for grid p?
option (USE_DEBUG_PRINT_P
  "Enable debug print for grid p" OFF)

# Use debug print for field of cell types?
option (USE_DEBUG_PRINT_TYPES
  "Enable debug print for field of cell types" OFF)

# Print in main.cpp after each single step?
option (USE_STEP_BY_STEP
  "Enable print/visualisation in each step" OFF)
  
# Find sdl2
INCLUDE(FindPkgConfig)
PKG_SEARCH_MODULE(SDL2 REQUIRED sdl2)
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${SDL2_INCLUDE_DIRS})

# Find mpi
find_package (MPI REQUIRED)
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${MPI_INCLUDE_PATH})

# Using glob is discouraged in cmake.
set (HEADERS
  # ${CMAKE_CURRENT_SOURCE_DIR}/src/communicator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/compute.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/geometry.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/grid.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/iterator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/parameter.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/solver.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/typedef.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/zeitgeist.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/argvparser.hpp
  )

set (SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
  # ${CMAKE_CURRENT_SOURCE_DIR}/src/communicator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/compute.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/geometry.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/grid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/iterator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/parameter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/solver.cpp
  )

# Add header and source of debug visualisation if it shall be built
# Additionally set a flag
if (DEBUG_VISU)
  set (HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/visu.hpp
    )

  set (SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/visu.cpp
    )

  add_definitions(-DUSE_DEBUG_VISU)
endif (DEBUG_VISU)

# Add header and source of VTK visualisation if it shall be built
# Aditionally set a flag
if (VTK_VISU)
  set (HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/vtk.hpp
    )

  set (SOURCES
    ${SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/vtk.cpp
    )

  add_definitions(-DUSE_VTK)
endif (VTK_VISU)

if (USE_DEBUG_PRINT_U)
  add_definitions(-DUSE_DEBUG_PRINT_U)
endif (USE_DEBUG_PRINT_U)

if (USE_DEBUG_PRINT_V)
  add_definitions(-DUSE_DEBUG_PRINT_V)
endif (USE_DEBUG_PRINT_V)

if (USE_DEBUG_PRINT_P)
  add_definitions(-DUSE_DEBUG_PRINT_P)
endif (USE_DEBUG_PRINT_P)

if (USE_DEBUG_PRINT_TYPES)
  add_definitions(-DUSE_DEBUG_PRINT_TYPES)
endif (USE_DEBUG_PRINT_TYPES)

if (USE_STEP_BY_STEP)
  add_definitions(-DUSE_STEP_BY_STEP)
endif (USE_STEP_BY_STEP)

# Build executable and link libs against it
add_executable (NumSim ${SOURCES})
target_link_libraries (NumSim ${SDL2_LIBRARIES} ${MPI_LIBRARIES})

if (MPI_COMPILE_FLAGS)
  set_target_properties(NumSim PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if (MPI_LINK_FLAGS)
  set_target_properties (NumSim PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()