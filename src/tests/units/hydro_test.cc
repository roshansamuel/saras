/********************************************************************************************************************************************
 * Saras
 * 
 * Copyright (C) 2019, Mahendra K. Verma
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     1. Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *     2. Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *     3. Neither the name of the copyright holder nor the
 *        names of its contributors may be used to endorse or promote products
 *        derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ********************************************************************************************************************************************
 */
/*! \file hydro_test.cc
 *
 *  \brief Definitions of functions to test the hydro class, specifically periodic data transfer.
 *  \sa alltests.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include <iostream>
#include "unittest.h"
#include "alltests.h"
#include "hydro.h"

void hydroTest(grid &gridData, parser &inputData, parallel &mpiData) {
    // ERROR IN COMPUTED SOLUTION FROM POISSON SOLVER
    real errorVal, errorTolerance;

    // POISSON SOLVER INSTANCE
    hydro *nseSolver;

    // CREATE NEW INSTANCE OF THE POISSON SOLVER USING THE GRID
#ifndef PLANAR
    nseSolver = new hydro_d3(gridData, inputData, mpiData);
#else
    nseSolver = new hydro_d2(gridData, inputData, mpiData);
#endif

    if (rootRank == 0) {
        std::cout << "\033[35mTesting NSE solver - Hydro\033[0m" << std::endl;
    }

    errorTolerance = 1.0e-10;
    if (rootRank == 0) {
        std::cout << std::setw(60) << std::left << "\033[34mTesting periodic data transfer operations\033[0m";
    }

    errorVal = nseSolver->testPeriodic();
    printResult(errorVal, errorTolerance);
}
