#!/bin/bash

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
 ##! \file testPoisson.sh
 #
 #   \brief Shell script to automatically compile and run tests on the Poisson library of SARAS
 #
 #   \author Roshan Samuel
 #   \date Mar 2020
 #   \copyright New BSD License
 #
 ############################################################################################################################################
 ##

# Test of Poisson library with Dirichlet BC
PROC=4

# If build directory doesn't exist, create it
if [ ! -d build ]; then
    mkdir build
fi

# Switch to build directory
cd build

# Run cmake with necessary flags for 2D Poisson test
#CC=mpicc CXX=mpicxx cmake ../../ -DPLANAR=ON -DTEST_POISSON=ON

# Run cmake with necessary flags for 3D Poisson test
CC=mpicc CXX=mpicxx cmake ../../ -DTEST_POISSON=ON

# Compile
make -j8

# Remove pre-existing executatbles
rm -f ../../tests/mgTest/saras

# Move the executable to the directory where the test will be performed
mv ../../saras ../../tests/mgTest/

# Switch to mgTest directory
cd ../../tests/mgTest/

# Run the test case
mpirun -np $PROC ./saras
