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
/*! \file poisson_test.cc
 *
 *  \brief Definitions of functions to test if the multi-grid solver is solving the Poisson equation correctly.
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
#include "poisson.h"
#include "sfield.h"

static void taylorGreen(sfield &rho, sfield &P_analytic, grid &mesh);

void poissonTest(grid &gridData, parser &inputData) {
    // ERROR IN COMPUTED SOLUTION FROM POISSON SOLVER
    real errorVal, errorTolerance;

    // POISSON SOLVER INSTANCE
    poisson *mgSolver;

    // SCALAR FIELDS TO TEST POISSON SOLVER
    sfield rho(gridData, "rho", true, true, true);
    sfield P_analytic(gridData, "P_anl", true, true, true);
    sfield P_calculat(gridData, "P_cal", true, true, true);

    // CREATE NEW INSTANCE OF THE POISSON SOLVER USING THE GRID
#ifndef PLANAR
    mgSolver = new multigrid_d3(gridData, inputData);
#else
    mgSolver = new multigrid_d2(gridData, inputData);
#endif

    if (rootRank == 0) {
        std::cout << "\033[35mTesting Poisson solver\033[0m" << std::endl;
    }

    errorTolerance = 1.0e-10;
    if (rootRank == 0) {
        std::cout << std::setw(60) << std::left << "\033[34mTesting prolongation operations\033[0m";
    }

    errorVal = mgSolver->testProlong();
    printResult(errorVal, errorTolerance);

    if (rootRank == 0) {
        std::cout << std::setw(60) << std::left << "\033[34mTesting strided data transfer operations\033[0m";
    }

    errorVal = mgSolver->testTransfer();
    printResult(errorVal, errorTolerance);

    taylorGreen(rho, P_analytic, gridData);

    if (inputData.probType == 2) {
        if (rootRank == 0) {
            std::cout << std::setw(60) << std::left << "\033[34mTesting periodic data transfer operations\033[0m";
        }

        errorVal = mgSolver->testPeriodic();
        printResult(errorVal, errorTolerance);
    }
}

static void taylorGreen(sfield &rho, sfield &P_analytic, grid &mesh) {
    // ANALYTIC SOLUTION AND RHS FOR THE POISSON SOLVER
#ifndef PLANAR
    for (int i=rho.F.F.lbound(0); i <= rho.F.F.ubound(0); i++) {
        for (int j=rho.F.F.lbound(1); j <= rho.F.F.ubound(1); j++) {
            for (int k=rho.F.F.lbound(2); k <= rho.F.F.ubound(2); k++) {
                P_analytic.F.F(i, j, k) = sin(1.0*M_PI*mesh.xStaggr(i))*
                                          cos(2.0*M_PI*mesh.yStaggr(j))*
                                          cos(4.0*M_PI*mesh.zStaggr(k));

                rho.F.F(i, j, k) = -21.0*M_PI*M_PI*sin(1.0*M_PI*mesh.xStaggr(i))*
                                                   cos(2.0*M_PI*mesh.yStaggr(j))*
                                                   cos(4.0*M_PI*mesh.zStaggr(k));
            }
        }
    }
#else
    int j = 0;
    for (int i=rho.F.F.lbound(0); i <= rho.F.F.ubound(0); i++) {
        for (int k=rho.F.F.lbound(2); k <= rho.F.F.ubound(2); k++) {
            P_analytic.F.F(i, j, k) = sin(1.0*M_PI*mesh.xStaggr(i))*
                                      cos(4.0*M_PI*mesh.zStaggr(k));

            rho.F.F(i, j, k) = -17.0*M_PI*M_PI*sin(1.0*M_PI*mesh.xStaggr(i))*
                                               cos(4.0*M_PI*mesh.zStaggr(k));
        }
    }
#endif
}
