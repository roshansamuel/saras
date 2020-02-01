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
/*! \file sfield.cc
 *
 *  \brief Definitions for functions of class sfield - scalar field
 *  \sa sfield.h
 *  \author Roshan Samuel, Ali Asad
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "plainsf.h"
#include "plainvf.h"
#include "sfield.h"
#include "vfield.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the sfield class
 *
 *          The instance of the field class to store the data of the scalar field is initialized, and the necessary grid
 *          transformation derivatives along each direction are chosen according to the grid staggering.
 *          The arrays to store the output from various operators like derivatives, convective derivatives, etc. are also
 *          allocated.
 *          Finally, an instance of the <B>mpidata</B> class is initialized to store the sub-arrays to be send/received
 *          across the processors during MPI communication.
 *
 * \param   gridData is a const reference to the global data contained in the grid class
 * \param   fieldName is a string value set by the user to name and identify the scalar field
 ********************************************************************************************************************************************
 */
sfield::sfield(const grid &gridData, std::string fieldName):
               gridData(gridData),
               F(gridData, fieldName, true, true, true),
               derS(gridData, F)
{
    this->fieldName = fieldName;

    interTempF.resize(F.fSize);
    interTempF.reindexSelf(F.flBound);

    derivTempF.resize(F.fSize);
    derivTempF.reindexSelf(F.flBound);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the diffusion term
 *
 *          It is assumed that the velocity is specified at face-centers, as required by the \ref sfield#computeNLin
 *          "computeNLin" function of sfield.
 *
 * \param   H is a pointer to a scalar field (sfield) to which the output of the function is to be written
 ********************************************************************************************************************************************
 */
void sfield::computeDiff(plainsf &H) {
    derivTempF = 0.0;
    derS.calcDerivative2xx(derivTempF);
    H.F(F.fCore) += derivTempF(F.fCore);
    
#ifndef PLANAR
    derivTempF = 0.0;
    derS.calcDerivative2yy(derivTempF);
    H.F(F.fCore) += derivTempF(F.fCore);
#endif

    derivTempF = 0.0;
    derS.calcDerivative2zz(derivTempF);
    H.F(F.fCore) += derivTempF(F.fCore);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the convective derivative of the scalar field
 *
 *          The function computes for the operator \f$ (\mathbf{u}.\nabla)f \f$ at the grid nodes of the scalar field f.
 *          To do so, the function needs the vector field (vfield) of velocity. It is assumed that the velocity is always
 *          specified at face-centers, and is interpolated accordingly to the scalar field grid points.
 *
 * \param   V is a const reference to a vector field (vfield) that specifies the convection velocity at each point
 ********************************************************************************************************************************************
 */
void sfield::computeNLin(const vfield &V, plainsf &H) {
    interTempF = 0.0;
    for (unsigned int i=0; i < F.VxIntSlices.size(); i++) {
        interTempF(F.fCore) += V.Vx.F(F.VxIntSlices(i));
    }

    derivTempF = 0.0;
    derS.calcDerivative1_x(derivTempF);
    H.F(F.fCore) -= interTempF(F.fCore)*derivTempF(F.fCore)/F.VxIntSlices.size();

#ifndef PLANAR
    interTempF = 0.0;
    for (unsigned int i=0; i < F.VyIntSlices.size(); i++) {
        interTempF(F.fCore) += V.Vy.F(F.VyIntSlices(i));
    }

    derivTempF = 0.0;
    derS.calcDerivative1_y(derivTempF);
    H.F(F.fCore) -= interTempF(F.fCore)*derivTempF(F.fCore)/F.VyIntSlices.size();
#endif

    interTempF = 0.0;
    for (unsigned int i=0; i < F.VzIntSlices.size(); i++) {
        interTempF(F.fCore) += V.Vz.F(F.VzIntSlices(i));
    }

    derivTempF = 0.0;
    derS.calcDerivative1_z(derivTempF);
    H.F(F.fCore) -= interTempF(F.fCore)*derivTempF(F.fCore)/F.VzIntSlices.size();
}

/**
 ********************************************************************************************************************************************
 * \brief   Operator to compute the gradient of the scalar field
 *
 *          The gradient operator computes the gradient of the cell centered scalar field, and stores it into a face-centered staggered
 *          plain vector field as defined by the tensor operation:
 *          \f$ \nabla f = \frac{\partial f}{\partial x}i + \frac{\partial f}{\partial y}j + \frac{\partial f}{\partial z}k \f$.
 *
 * \param   gradF is a reference to a plain vector field (plainvf) into which the computed gradient must be written.
 * \param   V is a const reference to a vector field (vfield) whose core slices are used to compute gradient, since plainvf doesn't have them
 ********************************************************************************************************************************************
 */
void sfield::gradient(plainvf &gradF, const vfield &V) {
    blitz::firstIndex i;
    blitz::secondIndex j;
    blitz::thirdIndex k;

    blitz::Range xColl, yColl, zColl;

    xColl = blitz::Range(gridData.collocCoreDomain.lbound(0), gridData.collocCoreDomain.ubound(0));
    yColl = blitz::Range(gridData.collocCoreDomain.lbound(1), gridData.collocCoreDomain.ubound(1));
    zColl = blitz::Range(gridData.collocCoreDomain.lbound(2), gridData.collocCoreDomain.ubound(2));

    gradF.Vx(V.Vx.fCore) = gridData.xi_xColloc(xColl)(i)*(F.F(V.Vx.fCRgt) - F.F(V.Vx.fCore))/gridData.dXi;
#ifndef PLANAR
    gradF.Vy(V.Vy.fCore) = gridData.et_yColloc(yColl)(j)*(F.F(V.Vy.fCBak) - F.F(V.Vy.fCore))/gridData.dEt;
#endif
    gradF.Vz(V.Vz.fCore) = gridData.zt_zColloc(zColl)(k)*(F.F(V.Vz.fCTop) - F.F(V.Vz.fCore))/gridData.dZt;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to synchronise data across all processors when performing parallel computations
 *
 *          This function calls the \ref mpidata#syncData "syncData" function of mpidata class to perform perform data-transfer and thus update
 *          the sub-domain boundary pads.
 ********************************************************************************************************************************************
 */
void sfield::syncData() {
    F.syncData();
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given plain scalar field
 *
 *          The unary operator += adds a given plain scalar field to the entire field stored as sfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to a plainsf to be added to the member field
 *
 * \return  A pointer to itself is returned by the scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
sfield& sfield::operator += (plainsf &a) {
    F.F += a.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given plain scalar field
 *
 *          The unary operator -= subtracts a given plain scalar field from the entire field stored as sfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to a plainsf to be deducted from the member field
 *
 * \return  A pointer to itself is returned by the scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
sfield& sfield::operator -= (plainsf &a) {
    F.F -= a.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given scalar field
 *
 *          The unary operator += adds a given scalar field to the entire field stored as sfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another sfield to be added to the member field
 *
 * \return  A pointer to itself is returned by the scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
sfield& sfield::operator += (sfield &a) {
    F.F += a.F.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given scalar field
 *
 *          The unary operator -= subtracts a given scalar field from the entire field stored as sfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another sfield to be deducted from the member field
 *
 * \return  A pointer to itself is returned by the scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
sfield& sfield::operator -= (sfield &a) {
    F.F -= a.F.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to multiply a scalar value to the scalar field
 *
 *          The unary operator *= multiplies a real value to the entire field stored as sfield and returns
 *          a pointer to itself.
 *
 * \param   a is a real number to be multiplied to the scalar field
 *
 * \return  A pointer to itself is returned by the scalar field class to which the operator belongs
 ********************************************************************************************************************************************
 */
sfield& sfield::operator *= (real a) {
    F.F *= a;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a plain scalar field to the scalar field
 *
 *          The operator = copies the contents of the input plain scalar field to itself.
 *
 * \param   a is the plainsf to be assigned to the scalar field
 ********************************************************************************************************************************************
 */
void sfield::operator = (plainsf &a) {
    F.F = a.F;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a scalar field to the scalar field
 *
 *          The operator = copies the contents of the input scalar field to itself.
 *
 * \param   a is the scalar field to be assigned to the scalar field
 ********************************************************************************************************************************************
 */
void sfield::operator = (sfield &a) {
    F.F = a.F.F;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a scalar value to the scalar field
 *
 *          The operator = assigns a real value to all the scalar field.
 *
 * \param   a is a real number to be assigned to the scalar field
 ********************************************************************************************************************************************
 */
void sfield::operator = (real a) {
    F.F = a;
}
