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
 #   \brief CMakeLists file where the testing run executable is linked for the unit tests of Saras solver.
 #
 #   \author Roshan Samuel
 #   \date Nov 2019
 #   \copyright New BSD License
 #
 ############################################################################################################################################
 ##

include_directories (
    ${PARENT_DIR}/src/solvers/hydro
    ${PARENT_DIR}/src/tests/global
    ${PARENT_DIR}/src/tests/units
)

set (EXECUTABLE_OUTPUT_PATH ${PARENT_DIR})

set (HYDRO_SOURCE ${PARENT_DIR}/src/solvers/hydro/hydro.cc ${PARENT_DIR}/src/solvers/hydro/hydro_d2.cc ${PARENT_DIR}/src/solvers/hydro/hydro_d3.cc)

set (TEST_SOURCE units/field_test.cc units/differ_test.cc units/nlin_test.cc units/poisson_test.cc units/hydro_test.cc)

set (UNIT_SOURCE global/unittest.cc)

set (SOURCES main.cc ${UNIT_SOURCE} ${HYDRO_SOURCE} ${TEST_SOURCE})

add_executable(saras_test ${SOURCES})

target_link_libraries(saras_test field grid parser probes reader writer parallel poisson yaml-cpp hdf5 debug /usr/local/lib/libblitz.a)
