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
/*! \file boundary.cc
 *
 *  \brief Definitions for functions of class boundary
 *  \sa boundary.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "boundary.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the boundary class
 *
 *          The class constructor initializes the mesh and other constant values for imposing the boundary condition.
 *          This constructor is called along with the constructors of any of the available classes of BCs.
 *
 * \param   mesh is a const reference to the global data contained in the grid class.
 * \param   inField is a reference to the field to which the boundary conditions must be applied.
 * \param   bcWall is a const integer which specifies the wall to which the BC must be applied.
 ********************************************************************************************************************************************
 */
boundary::boundary(const grid &mesh, field &inField, const int bcWall):
                          mesh(mesh), dField(inField), wallNum(bcWall) {
    // By default, rankFlag is true. i.e., the BC will be applied on all sub-domains.
    // This works only for Z-direction (in pencil decomposition), or both Z and Y directions (in slab decomposition).
    // This has to be changed appropriately.
    rankFlag = true;

    // By default, shiftVal is 1. i.e., the BC will be applied correctly only on the left wall along a given direction.
    // This has to be changed appropriately for the wall on the other side.
    shiftVal = 1;

    // Update rankFlag for the left wall (along X)
    if (wallNum == 0) rankFlag = mesh.rankData.xRank == 0;
    // Update rankFlag and shiftVal for the right wall (along X)
    if (wallNum == 1) {
        rankFlag = mesh.rankData.xRank == mesh.rankData.npX - 1;
        shiftVal = -1;
    }
#ifndef PLANAR
    // Update rankFlag for the front wall (along Y)
    if (wallNum == 2) rankFlag = mesh.rankData.yRank == 0;
    // Update rankFlag and shiftVal for the front wall (along Y)
    if (wallNum == 3) {
        rankFlag = mesh.rankData.yRank == mesh.rankData.npY - 1;
        shiftVal = -1;
    }
#endif
    // Update shiftVal for the top wall (along Z)
    if (wallNum == 5) {
        shiftVal = -1;
    }

    // Find the dimension along which the BC is being applied (X -> 0, Y -> 1, Z -> 2) using the wallNum
    shiftDim = (int) wallNum/2;
}

/**
 ********************************************************************************************************************************************
 * \brief   Prototype function to impose the boundary conditions on the given field
 *
 *          Based on the values of wallNum, shiftDim, shiftVal and fieldVal, the appropriate BC (neumann/dirichlet/mixed)
 *          is imposed on the field.
 ********************************************************************************************************************************************
 */
void boundary::imposeBC() { };
