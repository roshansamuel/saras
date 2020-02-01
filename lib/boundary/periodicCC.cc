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
/*! \file periodicCC.cc
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
 * \brief   Constructor of the periodicCC class
 *
 *          The constructor simply initializes the base boundary class using all the arguments supplied to it.
 *          Since periodic BC is being implemented, no additional values are necessary for the object.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   inField is a reference to the field to which the boundary conditions must be applied.
 * \param   bcWall is a const integer which specifies the wall to which the BC must be applied.
 ********************************************************************************************************************************************
 */
periodicCC::periodicCC(const grid &mesh, field &inField, const int bcWall):
                            boundary(mesh, inField, bcWall) { }

/**
 ********************************************************************************************************************************************
 * \brief   Function to impose periodic BC on a cell centered variable
 *
 *          For Saras solver, the wall passes through the cell centers of the variables.
 *          Hence the variable is lying on the wall for this case.
 *          This BC is used mainly in the Z-direction, because Saras supports pencil/slab decomposition only.
 *          Along Z, no MPI structures exist at all.
 *          Hence no MPI data transfer takes place, and this BC imposes periodic BC along Z.
 *
 ********************************************************************************************************************************************
 */
inline void periodicCC::imposeBC() {
    // The BC is applied for all ranks and no rankFlag is used
    // NOTE: The second point from the wall of the opposite side is read to impose periodic BC.
    // This corresponds to a bulk that is same as core - and is implemented as Method 3 in setBulkSlice function of field.cc
    if (shiftVal > 0) {
        // If shiftVal = 1, the wall is either left (0), front (2), or bottom (4) wall
        dField.F(dField.fWalls(wallNum)) = dField.F(dField.shift(shiftDim, dField.fWalls(wallNum + 1), -2));
    } else {
        // If shiftVal = -1, the wall is either right (1), back (3), or top (5) wall
        dField.F(dField.fWalls(wallNum)) = dField.F(dField.shift(shiftDim, dField.fWalls(wallNum - 1), 2));
    }
}
