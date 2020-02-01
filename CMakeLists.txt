#############################################################################################################################################
 # Saras
 # 
 # Copyright (C) 2019, Mahendra K. Verma
 #
 # All rights reserved.
 # 
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted provided that the following conditions are met:
 #     1. Redistributions of source code must retain the above copyright
 #        notice, this list of conditions and the following disclaimer.
 #     2. Redistributions in binary form must reproduce the above copyright
 #        notice, this list of conditions and the following disclaimer in the
 #        documentation and/or other materials provided with the distribution.
 #     3. Neither the name of the copyright holder nor the
 #        names of its contributors may be used to endorse or promote products
 #        derived from this software without specific prior written permission.
 # 
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 # ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 # WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 # DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 # ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 # (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 # LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 # ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 # (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 # SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 #
 ############################################################################################################################################
 ##
 ##! \file CMakeLists.txt
 #
 #   \brief Root level CMakeLists file.
 #
 #   \author Roshan Samuel
 #   \date Nov 2019
 #   \copyright New BSD License
 #
 ############################################################################################################################################
 ##

# Cmake setup
cmake_minimum_required (VERSION 2.8.3)
project (saras)

get_filename_component (INSTALL_DIR ${CMAKE_BINARY_DIR} PATH)
get_filename_component (PARENT_DIR ${INSTALL_DIR} PATH)

set (HEADER_DIR ${PARENT_DIR}/lib)

include_directories (${HEADER_DIR}/parallel
                     ${HEADER_DIR}/boundary
                     ${HEADER_DIR}/initial
                     ${HEADER_DIR}/poisson
                     ${HEADER_DIR}/force
                     ${HEADER_DIR}/field
                     ${HEADER_DIR}/grid
                     ${HEADER_DIR}/io)

message (STATUS "The parent directory is ${PARENT_DIR}")
message (STATUS "Install folder is ${INSTALL_DIR}")
message (STATUS "All headers are in ${HEADER_DIR}")

# Search for OpenMP package to enable hybrid parallelization and set the correct OpenMP compiler flag
find_package (OpenMP)
message (STATUS "Compiler flag for OpenMP is ${OpenMP_C_FLAGS}")

# Check if user defined flags are set
if (NOT DEFINED PLANAR)
    message (STATUS "PLANAR flag not set. Compiling for default 3D run")
    set (PLANAR OFF)
endif ()

# Add compiler flag for 2D/3D runs as requested by user while running cmake
if (PLANAR)
    message (STATUS "Compiling Saras for 2D simulations")
    add_definitions(-DPLANAR)
else ()
    message (STATUS "Compiling Saras for 3D simulations")
endif ()

# Add compiler flag for test runs as requested by user while running cmake
if (TEST_RUN)
    message (STATUS "Compiling Saras for running unit tests")
    add_definitions(-DTEST_RUN)
endif ()

# Add compiler flag for timing runs for scaling if requested by user
if (TIME_RUN)
    message (STATUS "Compiling Saras with additional timing calls")
    add_definitions(-DTIME_RUN)
endif ()

# Add compiler flag for using double precision calculations
if (REAL_DOUBLE)
    message (STATUS "Compiling Saras to solve with double precision calculations")
    add_definitions(-DREAL_DOUBLE)
endif ()

# Add compiler flag for using single precision calculations
if (REAL_SINGLE)
    message (STATUS "Compiling Saras to solve with single precision calculations")
    add_definitions(-DREAL_SINGLE)
endif ()

# Set compiler flags for normal and debug runs
set (CMAKE_CXX_FLAGS "-Wall ${OpenMP_C_FLAGS} -O3")
set (CMAKE_CXX_FLAGS_DEBUG "-Wall -g ${OpenMP_C_FLAGS} -DBZ_DEBUG")

add_subdirectory (lib)
add_subdirectory (src)

# Remove the PLANAR variable from cache to force user to manually set the PLANAR flag each time a 2D run is needed
unset (PLANAR CACHE)

# Remove the TEST_RUN variable from cache to force user to manually set the TEST_RUN flag each time unit-tests must be performed
unset (TEST_RUN CACHE)

# Remove the TIME_RUN variable from cache to force user to manually set the TIME_RUN flag each time scaling run must be performed
unset (TIME_RUN CACHE)

# Remove the REAL_DOUBLE variable from cache to force user to manually set the precision each time the solver is compiled
unset (REAL_DOUBLE CACHE)

# Remove the REAL_SINGLE variable from cache to force user to manually set the precision each time the solver is compiled
unset (REAL_SINGLE CACHE)
