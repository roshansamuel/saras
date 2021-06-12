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
/*! \file channelSine.cc
 *
 *  \brief Definitions for functions of class initial
 *  \sa initial.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include <cstdlib>
#include "initial.h"

/**
 ********************************************************************************************************************************************
 * \brief   Function to generate sinusoidal perturbation for channel flow
 *
 *          The sinusoidal function is multiplied with the equation of a parabola in order to satisfy the channel flow BCs.
 *          This is not divergence-free and has performed poorly in tests so far.
 *
 ********************************************************************************************************************************************
 */
void channelSine::initializeField(vfield &uField) {
    real kx = 10.0;

    if (mesh.rankData.rank == 0) std::cout << "Imposing sinusoidal perturbation initial condition for channel flow" << std::endl << std::endl;

#ifdef PLANAR
    // VELOCITY PERTURBATION FOR PERIODIC CHANNEL FLOW
    // X-VELOCITY
    for (int i=uField.Vx.F.lbound(0); i <= uField.Vx.F.ubound(0); i++) {
        for (int k=uField.Vx.F.lbound(2); k <= uField.Vx.F.ubound(2); k++) {
            uField.Vx.F(i, 0, k) = 4.0*(mesh.zStaggr(k) - mesh.zStaggr(k)*mesh.zStaggr(k));
        }
    }

    // Z-VELOCITY
    for (int i=uField.Vz.F.lbound(0); i <= uField.Vz.F.ubound(0); i++) {
        for (int k=uField.Vz.F.lbound(2); k <= uField.Vz.F.ubound(2); k++) {
            uField.Vz.F(i, 0, k) = 0.1*(mesh.zStaggr(k) - mesh.zStaggr(k)*mesh.zStaggr(k))*sin(2.0*M_PI*kx*mesh.xColloc(i)/mesh.xLen);
        }
    }

#else
    // VELOCITY PERTURBATION FOR PERIODIC CHANNEL FLOW
    // X-VELOCITY
    for (int i=uField.Vx.F.lbound(0); i <= uField.Vx.F.ubound(0); i++) {
        for (int j=uField.Vx.F.lbound(1); j <= uField.Vx.F.ubound(1); j++) {
            for (int k=uField.Vx.F.lbound(2); k <= uField.Vx.F.ubound(2); k++) {
                uField.Vx.F(i, j, k) = 4.0*(mesh.zStaggr(k) - mesh.zStaggr(k)*mesh.zStaggr(k));
            }
        }
    }

    // Y-VELOCITY
    uField.Vy.F = 0.0;

    // Z-VELOCITY
    for (int i=uField.Vz.F.lbound(0); i <= uField.Vz.F.ubound(0); i++) {
        for (int j=uField.Vz.F.lbound(1); j <= uField.Vz.F.ubound(1); j++) {
            for (int k=uField.Vz.F.lbound(2); k <= uField.Vz.F.ubound(2); k++) {
                uField.Vz.F(i, j, k) = 0.1*(mesh.zStaggr(k) - mesh.zStaggr(k)*mesh.zStaggr(k))*sin(2.0*M_PI*kx*mesh.xColloc(i)/mesh.xLen);
            }
        }
    }
#endif
}
