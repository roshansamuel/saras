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
 #   \brief CMakeLists file where the executable is linked for the Saras solver.
 #
 #   \author Roshan Samuel
 #   \date Nov 2019
 #   \copyright New BSD License
 #
 ############################################################################################################################################
 ##

include_directories (
    ${PARENT_DIR}/src/solvers/hydro
    ${PARENT_DIR}/src/solvers/scalar
)

set (EXECUTABLE_OUTPUT_PATH ${PARENT_DIR})

set (HYDRO_SOURCE hydro/hydro.cc hydro/hydro_d2.cc hydro/hydro_d3.cc)

set (SCALAR_SOURCE scalar/scalar.cc scalar/scalar_d2.cc scalar/scalar_d3.cc)

set (SOURCES main.cc ${HYDRO_SOURCE} ${SCALAR_SOURCE})

add_executable (saras ${SOURCES})

target_link_libraries(saras field grid parser probes initial reader writer tseries boundary parallel poisson force yaml-cpp hdf5 debug /usr/local/lib/libblitz.a)
