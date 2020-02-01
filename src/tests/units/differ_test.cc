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
/*! \file differ_test.cc
 *
 *  \brief Definitions of functions to test correctness of finite differencing.
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
#include "sfield.h"

static void taylorGreen(sfield &F, sfield &dF, grid &mesh);

void differTest(grid &gridData) {
    // ERROR IN COMPUTED DIVERGENCE
    real errorTolerance;

    sfield F(gridData, "F", true, true, true);
    sfield dF(gridData, "dF", true, true, true);

    errorTolerance = gridData.dXi*gridData.dXi*4.0*M_PI*M_PI;

    if (rootRank == 0) {
        std::cout << std::setw(60) << std::left << "\033[34mTesting computation of derivatives in sfield\033[0m";
    }

    taylorGreen(F, dF, gridData);

    F.calcDerivatives1();

    testError(dF.F.F, F.dF_dx, 1, errorTolerance);
}

static void taylorGreen(sfield &F, sfield &dF, grid &mesh) {
#ifndef PLANAR
    for (int i=F.F.F.lbound(0); i <= F.F.F.ubound(0); i++) {
        for (int j=F.F.F.lbound(1); j <= F.F.F.ubound(1); j++) {
            for (int k=F.F.F.lbound(2); k <= F.F.F.ubound(2); k++) {
                F.F.F(i, j, k) = sin(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                 cos(2.0*M_PI*mesh.yStaggr(j)/mesh.yLen)*
                                 cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
                dF.F.F(i, j, k) = 2.0*M_PI*cos(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                           cos(2.0*M_PI*mesh.yStaggr(j)/mesh.yLen)*
                                           cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
            }
        }
    }
#else
    int j = 0;
    for (int i=F.F.F.lbound(0); i <= F.F.F.ubound(0); i++) {
        for (int k=F.F.F.lbound(2); k <= F.F.F.ubound(2); k++) {
            F.F.F(i, j, k) = sin(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                             cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
            dF.F.F(i, j, k) = 2.0*M_PI*cos(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                       cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
        }
    }
#endif
}
