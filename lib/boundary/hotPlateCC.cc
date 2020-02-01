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
 * \brief   Constructor of the hotPlateCC class
 *
 *          The constructor initializes the base boundary class using part of the arguments supplied to it.
 *          The radius of the heating plate, denoted by patchRadius, is also set in the initialization list.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   inField is a reference to scalar field to which the boundary conditions must be applied.
 * \param   bcWall is a const integer which specifies the wall to which the BC must be applied.
 * \param   plateRad is the const real value of the radius of the heating plate.
 ********************************************************************************************************************************************
 */
hotPlateCC::hotPlateCC(const grid &mesh, field &inField, const int bcWall, const real plateRad):
                            boundary(mesh, inField, bcWall), patchRadius(plateRad) {
    setXYZ();
    createPatch(patchRadius);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to set the x, y and z arrays by reference to the values in grid class
 *
 *          Depending on whether the variable is collocated or staggered along a given direction,
 *          the function sets the reference to each grid array, x, y and z, to either the
 *          array of locations of face-centers or cell-centers.
 ********************************************************************************************************************************************
 */
void hotPlateCC::setXYZ() {
    dField.xStag ? x.reference(mesh.xStaggr) : x.reference(mesh.xColloc);
    dField.yStag ? y.reference(mesh.yStaggr) : y.reference(mesh.yColloc);
    dField.zStag ? z.reference(mesh.zStaggr) : z.reference(mesh.zColloc);
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to create a heating patch for non-homogeneous boundary conditions
 *
 *          The current implementation is for a conducting circular heating plate surrounded by adiabatic insulated walls.
 *
 ********************************************************************************************************************************************
 */
void hotPlateCC::createPatch(const real patchRadius) {
    real patchCentX, patchCentY, patchCentZ;

    patchCentX = mesh.xLen/2.0;
    patchCentY = mesh.yLen/2.0;
    patchCentZ = 0.0;

    wallData.resize(blitz::TinyVector<int, 3>(dField.fWalls(wallNum).ubound() - dField.fWalls(wallNum).lbound() + 1));
    wallData.reindexSelf(dField.fWalls(wallNum).lbound());
    wallData = 0.0;

    wallMask.resize(blitz::TinyVector<int, 3>(dField.fWalls(wallNum).ubound() - dField.fWalls(wallNum).lbound() + 1));
    wallMask.reindexSelf(dField.fWalls(wallNum).lbound());
    wallMask = true;

    for (int iX = dField.fWalls(wallNum).lbound(0); iX <= dField.fWalls(wallNum).ubound(0); iX++) {
        for (int iY = dField.fWalls(wallNum).lbound(1); iY <= dField.fWalls(wallNum).ubound(1); iY++) {
            for (int iZ = dField.fWalls(wallNum).lbound(2); iZ <= dField.fWalls(wallNum).ubound(2); iZ++) {
                real ptRadius = sqrt(pow((x(iX) - patchCentX), 2) + pow((y(iY) - patchCentY), 2) + pow((z(iZ) - patchCentZ), 2));
                if (ptRadius <= patchRadius) {
                    wallData(iX, iY, iZ) = 1.0;
                    wallMask(iX, iY, iZ) = false;
                }
            }
        }
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to impose Mixed BC of a heating plate on a cell centered variable
 *
 *          For Saras solver, the wall passes through the cell centers of the variables.
 *          Hence the variable is lying on the wall for this case.
 *          Accordingly the value of the variable is directly set to 1.0 on the conducting sections of the wall.
 *          At the rest of the wall, Neumann BC is imposed using a boolean wallMask variable.
 *
 ********************************************************************************************************************************************
 */
void hotPlateCC::imposeBC() {
    // First impose Neumann BC everywhere for adiabatic wall
    dField.F(dField.fWalls(wallNum)) = dField.F(dField.shift(shiftDim, dField.fWalls(wallNum), 1));

    // Now in the area of the circular patch, make all values 0 in order to apply conducting BC
    dField.F(dField.fWalls(wallNum)) = dField.F(dField.fWalls(wallNum))*wallMask;

    // Finally apply conducting BC in the circular patch alone
    dField.F(dField.fWalls(wallNum)) = dField.F(dField.fWalls(wallNum)) +
            wallData*(2.0*wallData - dField.F(dField.shift(shiftDim, dField.fWalls(wallNum), 1)));
}
