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
/*! \file neumannCC.cc
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
 * \brief   Constructor of the neumannCC class
 *
 *          The constructor initializes the base boundary class using part of the arguments supplied to it.
 *          The value of the derivative of the field at the boundary, denoted by fieldValue, is also set in the initialization list.
 *
 * \param   mesh is a const reference to the global data contained in the grid class.
 * \param   inField is a reference to field to which the boundary conditions must be applied.
 * \param   bcWall is a const integer which specifies the wall to which the BC must be applied.
 * \param   bcValue is the const real value of the derivative of the variable at the boundary.
 ********************************************************************************************************************************************
 */
neumannCC::neumannCC(const grid &mesh, field &inField, const int bcWall, const real bcValue):
                            boundary(mesh, inField, bcWall), fieldValue(bcValue) { }

/**
 ********************************************************************************************************************************************
 * \brief   Function to impose Neumann BC on a cell centered variable
 *
 *          For Saras solver, the wall passes through the cell centers of the variables.
 *          Hence the variable is lying on the wall for this case.
 *          Accordingly the derivative of the variable is set on the wall.
 *
 ********************************************************************************************************************************************
 */
inline void neumannCC::imposeBC() {
    if (rankFlag) {
        // This implementation assumes that the derivative at boundary is 0, and needs update
        dField.F(dField.fWalls(wallNum)) = dField.F(dField.shift(shiftDim, dField.fWalls(wallNum), shiftVal));

        //dField.F(dField.fWalls(4)) = dField.F(dField.fWalls(4)) +
        //        wallData*(2.0*wallData - dField.F(dField.shift(2, dField.fWalls(4), 1)));
    }
}
