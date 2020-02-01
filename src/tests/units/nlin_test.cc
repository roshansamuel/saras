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
/*! \file nlin_test.cc
 *
 *  \brief Definitions of functions to test if non-linear terms are being computed correctly.
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
#include "vfield.h"

static void taylorGreen(vfield &V, vfield &H, grid &mesh);

void nlinTest(grid &gridData) {
    // ERROR IN COMPUTED NON-LINEAR TERMS
    real errorTolerance;

    // VECTOR FIELDS FOR VELOCITY AND NON-LINEAR TERMS
    vfield V(gridData, "V");
    vfield H_analytic(gridData, "H_anl");
    vfield H_calculat(gridData, "H_cal");

    // TOLERANCE IN ERROR OF COMPUTED VALUES
#ifdef PLANAR
    errorTolerance = std::max(gridData.dXi, gridData.dZt);
#else
    errorTolerance = std::max(std::max(gridData.dXi, gridData.dEt), gridData.dZt);
#endif
    errorTolerance *= 2.0*M_PI;

    if (rootRank == 0) {
        std::cout << std::setw(60) << std::left << "\033[34mTesting non-linear term calculation\033[0m";
    }

    // SET VELOCITY FIELD AND ANALYTIC VALUES OF NON-LINEAR TERMS
    taylorGreen(V, H_analytic, gridData);

    // COMPUTE NON-LINEAR TERMS
    V.computeNLin(V, H_calculat);

    testError(H_analytic.Vx.F.F, H_calculat.Vx.F.F, 1, errorTolerance);
}

static void taylorGreen(vfield &V, vfield &H, grid &mesh) {
    // X-Velocity and X-component of non-linear term
#ifndef PLANAR
    for (int i=V.Vx.F.F.lbound(0); i <= V.Vx.F.F.ubound(0); i++) {
        for (int j=V.Vx.F.F.lbound(1); j <= V.Vx.F.F.ubound(1); j++) {
            for (int k=V.Vx.F.F.lbound(2); k <= V.Vx.F.F.ubound(2); k++) {
                V.Vx.F.F(i, j, k) = sin(2.0*M_PI*mesh.xColloc(i)/mesh.xLen)*
                                    cos(2.0*M_PI*mesh.yStaggr(j)/mesh.yLen)*
                                    cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);

                H.Vx.F.F(i, j, k) = M_PI*sin(4.0*M_PI*mesh.xColloc(i)/mesh.xLen)*
                                     pow(cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen), 2.0);
            }
        }
    }
#else
    int j = 0;
    for (int i=V.Vx.F.F.lbound(0); i <= V.Vx.F.F.ubound(0); i++) {
        for (int k=V.Vx.F.F.lbound(2); k <= V.Vx.F.F.ubound(2); k++) {
            V.Vx.F.F(i, j, k) = sin(2.0*M_PI*mesh.xColloc(i)/mesh.xLen)*
                                cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);

            H.Vx.F.F(i, j, k) = M_PI*sin(4.0*M_PI*mesh.xColloc(i)/mesh.xLen);
        }
    }
#endif

    // Y-Velocity and Y-component of non-linear term
#ifndef PLANAR
    for (int i=V.Vy.F.F.lbound(0); i <= V.Vy.F.F.ubound(0); i++) {
        for (int j=V.Vy.F.F.lbound(1); j <= V.Vy.F.F.ubound(1); j++) {
            for (int k=V.Vy.F.F.lbound(2); k <= V.Vy.F.F.ubound(2); k++) {
                V.Vy.F.F(i, j, k) = -cos(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                     sin(2.0*M_PI*mesh.yColloc(j)/mesh.yLen)*
                                     cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);

                H.Vy.F.F(i, j, k) = M_PI*sin(4.0*M_PI*mesh.yColloc(j)/mesh.yLen)*
                                     pow(cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen), 2.0);
            }
        }
    }
#else
    V.Vy.F.F = 0.0;
#endif

    // Z-Velocity and Z-component of non-linear term
#ifndef PLANAR
    V.Vz.F.F = 0.0;
    H.Vz.F.F = 0.0;
#else
    for (int i=V.Vz.F.F.lbound(0); i <= V.Vz.F.F.ubound(0); i++) {
        for (int k=V.Vz.F.F.lbound(2); k <= V.Vz.F.F.ubound(2); k++) {
            V.Vz.F.F(i, j, k) = -cos(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                 sin(2.0*M_PI*mesh.zColloc(k)/mesh.zLen);

            H.Vz.F.F(i, j, k) = M_PI*sin(4.0*M_PI*mesh.zColloc(k)/mesh.zLen);
        }
    }
#endif
}
