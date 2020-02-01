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
/*! \file plainvf.cc
 *
 *  \brief Definitions for functions of class plainvf - plain vector field
 *  \sa plainvf.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "plainvf.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the plainvf class
 *
 *          Three blitz arrays to store the data of the three component scalar fields are initialized.
 *          The name for the plain vector field as given by the user is also assigned.
 *
 * \param   gridData is a const reference to the global data contained in the grid class
 * \param   refV is a const reference to a sample vfield according to whose components the components of plainvf is resized
 ********************************************************************************************************************************************
 */
plainvf::plainvf(const grid &gridData, const vfield &refV): gridData(gridData) {
    Vx.resize(refV.Vx.fSize);
    Vx.reindexSelf(refV.Vx.flBound);

    mpiVxData = new mpidata(Vx, gridData.rankData);
    mpiVxData->createSubarrays(refV.Vx.fSize, refV.Vx.cuBound + 1, gridData.padWidths, refV.Vx.xStag, refV.Vx.yStag);

    Vy.resize(refV.Vy.fSize);
    Vy.reindexSelf(refV.Vy.flBound);

    mpiVyData = new mpidata(Vy, gridData.rankData);
    mpiVyData->createSubarrays(refV.Vy.fSize, refV.Vy.cuBound + 1, gridData.padWidths, refV.Vy.xStag, refV.Vy.yStag);

    Vz.resize(refV.Vz.fSize);
    Vz.reindexSelf(refV.Vz.flBound);

    mpiVzData = new mpidata(Vz, gridData.rankData);
    mpiVzData->createSubarrays(refV.Vz.fSize, refV.Vz.cuBound + 1, gridData.padWidths, refV.Vz.xStag, refV.Vz.yStag);
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given plain vector field
 *
 *          The unary operator += adds a given plain vector field to the entire field stored as plainvf and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another plainvf to be deducted from the member fields
 *
 * \return  A pointer to itself is returned by the plain vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
plainvf& plainvf::operator += (plainvf &a) {
    Vx += a.Vx;
    Vy += a.Vy;
    Vz += a.Vz;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given plain vector field
 *
 *          The unary operator -= subtracts a given plain vector field from the entire field stored as plainvf and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another plainvf to be deducted from the member fields
 *
 * \return  A pointer to itself is returned by the plain vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
plainvf& plainvf::operator -= (plainvf &a) {
    Vx -= a.Vx;
    Vy -= a.Vy;
    Vz -= a.Vz;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given vector field
 *
 *          The unary operator += adds a given vector field to the entire field stored as plainvf and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another vfield to be deducted from the member fields
 *
 * \return  A pointer to itself is returned by the plain vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
plainvf& plainvf::operator += (vfield &a) {
    Vx += a.Vx.F;
    Vy += a.Vy.F;
    Vz += a.Vz.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given vector field
 *
 *          The unary operator -= subtracts a given vector field from the entire field stored as plainvf and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another vfield to be deducted from the member fields
 *
 * \return  A pointer to itself is returned by the plain vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
plainvf& plainvf::operator -= (vfield &a) {
    Vx -= a.Vx.F;
    Vy -= a.Vy.F;
    Vz -= a.Vz.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to multiply a scalar value to the plain vector field
 *
 *          The unary operator *= multiplies a real value to all the fields (Vx, Vy and Vz) stored in plainvf and returns
 *          a pointer to itself.
 *
 * \param   a is a real number to be multiplied to the plain vector field
 *
 * \return  A pointer to itself is returned by the plain vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
plainvf& plainvf::operator *= (real a) {
    Vx *= a;
    Vy *= a;
    Vz *= a;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign another plain vector field to the plain vector field
 *
 *          The operator = assigns all the three scalar sub-fields of a plainvf to all the corresponding fields (Vx, Vy and Vz).
 *
 * \param   a is a plainvf to be assigned to the plain vector field
 ********************************************************************************************************************************************
 */
void plainvf::operator = (plainvf &a) {
    Vx = a.Vx;
    Vy = a.Vy;
    Vz = a.Vz;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign another vector field to the plain vector field
 *
 *          The operator = assigns all the three scalar sub-fields of a plainvf to all the corresponding fields (Vx, Vy and Vz).
 *
 * \param   a is a vfield to be assigned to the plain vector field
 ********************************************************************************************************************************************
 */
void plainvf::operator = (vfield &a) {
    Vx = a.Vx.F;
    Vy = a.Vy.F;
    Vz = a.Vz.F;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a scalar value to the plain vector field
 *
 *          The operator = assigns a real value to all the fields (Vx, Vy and Vz) stored in plainvf.
 *
 * \param   a is a real number to be assigned to the plain vector field
 ********************************************************************************************************************************************
 */
void plainvf::operator = (real a) {
    Vx = a;
    Vy = a;
    Vz = a;
}
