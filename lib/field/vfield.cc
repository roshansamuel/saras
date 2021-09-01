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
/*! \file vfield.cc
 *
 *  \brief Definitions for functions of class vfield - vector field
 *  \sa vfield.h
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
#include <math.h>

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the vfield class
 *
 *          Three instances of field class are initialized.
 *          Each instance corresponds to a component of the vector field.
 *          The fields are initialized with appropriate grid to place the components on the cell faces.
 *          The name of the vector field as given by the user is also assigned.
 *
 * \param   gridData is a const reference to the global data in the grid class
 * \param   fieldName is a string value used to name and identify the vector field
 ********************************************************************************************************************************************
 */
vfield::vfield(const grid &gridData, std::string fieldName):
               gridData(gridData),
               Vx(gridData, "Vx", false, true, true),
               Vy(gridData, "Vy", true, false, true),
               Vz(gridData, "Vz", true, true, false),
               derVx(gridData, Vx), derVy(gridData, Vy), derVz(gridData, Vz)
{
    this->fieldName = fieldName;

    interTempX.resize(Vx.fSize);
    interTempX.reindexSelf(Vx.flBound);

    derivTempX.resize(Vx.fSize);
    derivTempX.reindexSelf(Vx.flBound);

#ifndef PLANAR
    interTempY.resize(Vy.fSize);
    interTempY.reindexSelf(Vy.flBound);

    derivTempY.resize(Vy.fSize);
    derivTempY.reindexSelf(Vy.flBound);
#endif

    interTempZ.resize(Vz.fSize);
    interTempZ.reindexSelf(Vz.flBound);

    derivTempZ.resize(Vz.fSize);
    derivTempZ.reindexSelf(Vz.flBound);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the diffusion term
 *
 *          The diffusion term (grad-squared) is caulculated here.
 *          The second derivatives of each component field are calculated along x, y and z.
 *          These terms are added to the corresponding components of the given plain
 *          vector field (plainvf), which is usually the RHS of the PDE being solved.
 *
 * \param   H is a reference to the plainvf into which the output is written
 ********************************************************************************************************************************************
 */
void vfield::computeDiff(plainvf &H) {
    derivTempX = 0.0;
    derVx.calcDerivative2xx(derivTempX);
    H.Vx(Vx.fCore) += derivTempX(Vx.fCore);

#ifndef PLANAR
    derivTempX = 0.0;
    derVx.calcDerivative2yy(derivTempX);
    H.Vx(Vx.fCore) += derivTempX(Vx.fCore);
#endif

    derivTempX = 0.0;
    derVx.calcDerivative2zz(derivTempX);
    H.Vx(Vx.fCore) += derivTempX(Vx.fCore);

#ifndef PLANAR
    derivTempY = 0.0;
    derVy.calcDerivative2xx(derivTempY);
    H.Vy(Vy.fCore) += derivTempY(Vy.fCore);

    derivTempY = 0.0;
    derVy.calcDerivative2yy(derivTempY);
    H.Vy(Vy.fCore) += derivTempY(Vy.fCore);

    derivTempY = 0.0;
    derVy.calcDerivative2zz(derivTempY);
    H.Vy(Vy.fCore) += derivTempY(Vy.fCore);
#endif

    derivTempZ = 0.0;
    derVz.calcDerivative2xx(derivTempZ);
    H.Vz(Vz.fCore) += derivTempZ(Vz.fCore);

#ifndef PLANAR
    derivTempZ = 0.0;
    derVz.calcDerivative2yy(derivTempZ);
    H.Vz(Vz.fCore) += derivTempZ(Vz.fCore);
#endif

    derivTempZ = 0.0;
    derVz.calcDerivative2zz(derivTempZ);
    H.Vz(Vz.fCore) += derivTempZ(Vz.fCore);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the convective derivative of the vector field
 *
 *          The function calculates \f$ (\mathbf{u}.\nabla)\mathbf{v} \f$ on the vector field, \f$\mathbf{v}\f$.
 *          To do so, the function needs the vector field (vfield) of velocity, \f$\mathbf{u}\f$.
 *          This value is used in three separate calls to the \ref sfield#computeNLin "computeNLin" function of sfield
 *          to compute the derivatives for the three components of the vector field.
 *          It is assumed that the velocity is specified at face-centers, as required by the \ref sfield#computeNLin
 *          "computeNLin" function of sfield.
 *
 * \param   V is a const reference to the vector field (vfield) that specifies the convection velocity
 * \param   H is a reference to the plainvf into which the output is written
 ********************************************************************************************************************************************
 */
void vfield::computeNLin(const vfield &V, plainvf &H) {
    // Compute non-linear term for the Vx component
    interTempX = 0.0;
    for (unsigned int i=0; i < Vx.VxIntSlices.size(); i++) {
        interTempX(Vx.fCore) += V.Vx.F(Vx.VxIntSlices(i));
    }

    derivTempX = 0.0;
    derVx.calcDerivative1_x(derivTempX);
    H.Vx(Vx.fCore) -= interTempX(Vx.fCore)*derivTempX(Vx.fCore)/Vx.VxIntSlices.size();

#ifndef PLANAR
    interTempX = 0.0;
    for (unsigned int i=0; i < Vx.VyIntSlices.size(); i++) {
        interTempX(Vx.fCore) += V.Vy.F(Vx.VyIntSlices(i));
    }

    derivTempX = 0.0;
    derVx.calcDerivative1_y(derivTempX);
    H.Vx(Vx.fCore) -= interTempX(Vx.fCore)*derivTempX(Vx.fCore)/Vx.VyIntSlices.size();
#endif

    interTempX = 0.0;
    for (unsigned int i=0; i < Vx.VzIntSlices.size(); i++) {
        interTempX(Vx.fCore) += V.Vz.F(Vx.VzIntSlices(i));
    }

    derivTempX = 0.0;    
    derVx.calcDerivative1_z(derivTempX);
    H.Vx(Vx.fCore) -= interTempX(Vx.fCore)*derivTempX(Vx.fCore)/Vx.VzIntSlices.size();

// Compute non-linear term for the Vy component
#ifndef PLANAR
    interTempY = 0.0;
    for (unsigned int i=0; i < Vy.VxIntSlices.size(); i++) {
        interTempY(Vy.fCore) += V.Vx.F(Vy.VxIntSlices(i));
    }

    derivTempY = 0.0;
    derVy.calcDerivative1_x(derivTempY);
    H.Vy(Vy.fCore) -= interTempY(Vy.fCore)*derivTempY(Vy.fCore)/Vy.VxIntSlices.size();

    interTempY = 0.0;
    for (unsigned int i=0; i < Vy.VyIntSlices.size(); i++) {
        interTempY(Vy.fCore) += V.Vy.F(Vy.VyIntSlices(i));
    }

    derivTempY = 0.0;
    derVy.calcDerivative1_y(derivTempY);
    H.Vy(Vy.fCore) -= interTempY(Vy.fCore)*derivTempY(Vy.fCore)/Vy.VyIntSlices.size();

    interTempY = 0.0;
    for (unsigned int i=0; i < Vy.VzIntSlices.size(); i++) {
        interTempY(Vy.fCore) += V.Vz.F(Vy.VzIntSlices(i));
    }

    derivTempY = 0.0;
    derVy.calcDerivative1_z(derivTempY);
    H.Vy(Vy.fCore) -= interTempY(Vy.fCore)*derivTempY(Vy.fCore)/Vy.VzIntSlices.size();
#endif

    // Compute non-linear term for the Vz component
    interTempZ = 0.0;
    for (unsigned int i=0; i < Vz.VxIntSlices.size(); i++) {
        interTempZ(Vz.fCore) += V.Vx.F(Vz.VxIntSlices(i));
    }

    derivTempZ = 0.0;
    derVz.calcDerivative1_x(derivTempZ);
    H.Vz(Vz.fCore) -= interTempZ(Vz.fCore)*derivTempZ(Vz.fCore)/Vz.VxIntSlices.size();

#ifndef PLANAR
    interTempZ = 0.0;
    for (unsigned int i=0; i < Vz.VyIntSlices.size(); i++) {
        interTempZ(Vz.fCore) += V.Vy.F(Vz.VyIntSlices(i));
    }

    derivTempZ = 0.0;
    derVz.calcDerivative1_y(derivTempZ);
    H.Vz(Vz.fCore) -= interTempZ(Vz.fCore)*derivTempZ(Vz.fCore)/Vz.VyIntSlices.size();
#endif

    interTempZ = 0.0;
    for (unsigned int i=0; i < Vz.VzIntSlices.size(); i++) {
        interTempZ(Vz.fCore) += V.Vz.F(Vz.VzIntSlices(i));
    }

    derivTempZ = 0.0;
    derVz.calcDerivative1_z(derivTempZ);
    H.Vz(Vz.fCore) -= interTempZ(Vz.fCore)*derivTempZ(Vz.fCore)/Vz.VzIntSlices.size();
}

/**
 ********************************************************************************************************************************************
 * \brief   Operator is used to calculate time step #dt_out using CFL Condition
 *           
 *          When the parameters specify that time-step must be adaptively
 *          calculated using the Courant Number given in the parameters file,
 *          this function will provide the \f$ dt \f$ using the maximum values
 *          of the components of the vfield.
 *
 * \param   dt_out is a reference to the real value into which the calculated value of time-step is written
 *********************************************************************************************************************************************
 */
void vfield::computeTStp(real &dt_out) {
    real Umax, Vmax, Wmax;
    real delx, dely, delz; 

    delx = gridData.dXi;
    Umax = Vx.fieldMax();
#ifdef PLANAR
    dely = 0.0;
    Vmax = 1.0;
#else
    dely = gridData.dEt;
    Vmax = Vy.fieldMax();
#endif
    delz = gridData.dZt;
    Wmax = Vz.fieldMax();

    dt_out = real(gridData.inputParams.courantNumber*(delx/Umax + dely/Vmax + delz/Wmax));
}

/**
 ********************************************************************************************************************************************
 * \brief   Operator to compute the divergence of the vector field
 *
 *          The operator computes the divergence of a face-centered staggered vector field, and stores it into a cell centered
 *          scalar field as defined by the tensor operation:
 *          \f$ \nabla . \mathbf{v} = \frac{\partial \mathbf{v}}{\partial x} +
 *                                    \frac{\partial \mathbf{v}}{\partial y} +
 *                                    \frac{\partial \mathbf{v}}{\partial z} \f$.
 *
 * \param   divV is a reference to the scalar field (sfield) into which the computed divergence is written.
 ********************************************************************************************************************************************
 */
void vfield::divergence(plainsf &divV, const sfield &P) {
    blitz::firstIndex i;
    blitz::secondIndex j;
    blitz::thirdIndex k;

    blitz::Range xStag, yStag, zStag;

    xStag = blitz::Range(P.F.fCore.lbound(0), P.F.fCore.ubound(0));
    yStag = blitz::Range(P.F.fCore.lbound(1), P.F.fCore.ubound(1));
    zStag = blitz::Range(P.F.fCore.lbound(2), P.F.fCore.ubound(2));

    divV = 0.0;

#ifdef PLANAR
    divV.F(P.F.fCore) = gridData.xi_xStaggr(xStag)(i)*(Vx.F(P.F.fCore) - Vx.F(P.F.fCLft))/gridData.dXi + 
                        gridData.zt_zStaggr(zStag)(k)*(Vz.F(P.F.fCore) - Vz.F(P.F.fCBot))/gridData.dZt;
#else
    divV.F(P.F.fCore) = gridData.xi_xStaggr(xStag)(i)*(Vx.F(P.F.fCore) - Vx.F(P.F.fCLft))/gridData.dXi + 
                        gridData.et_yStaggr(yStag)(j)*(Vy.F(P.F.fCore) - Vy.F(P.F.fCFrt))/gridData.dEt + 
                        gridData.zt_zStaggr(zStag)(k)*(Vz.F(P.F.fCore) - Vz.F(P.F.fCBot))/gridData.dZt;
#endif
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to synchronise data across all processors when performing parallel computations
 *
 *          Each of the individual field components have their own subroutine, \ref sfield#syncData "syncData" to send and
 *          receive data across its MPI decomposed sub-domains.
 *          This function calls the \ref sfield#syncData "syncData" function of its components to update the sub-domain boundary pads.
 ********************************************************************************************************************************************
 */
void vfield::syncData() {
    Vx.syncData();
    Vy.syncData();
    Vz.syncData();
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given plain vector field
 *
 *          The unary operator += adds a given plain vector field to the vfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to the plainvf to be added to the member fields
 *
 * \return  A pointer to itself is returned by the vector field object to which the operator belongs
 ********************************************************************************************************************************************
 */
vfield& vfield::operator += (plainvf &a) {
    Vx.F += a.Vx;
    Vy.F += a.Vy;
    Vz.F += a.Vz;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given plain vector field
 *
 *          The unary operator -= subtracts a given plain vector field from the vfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to the plainvf to be subtracted from the member fields
 *
 * \return  A pointer to itself is returned by the vector field object to which the operator belongs
 ********************************************************************************************************************************************
 */
vfield& vfield::operator -= (plainvf &a) {
    Vx.F -= a.Vx;
    Vy.F -= a.Vy;
    Vz.F -= a.Vz;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given vector field
 *
 *          The unary operator += adds a given vector field to the vfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to the vfield to be added to the member fields
 *
 * \return  A pointer to itself is returned by the vector field object to which the operator belongs
 ********************************************************************************************************************************************
 */
vfield& vfield::operator += (vfield &a) {
    Vx.F += a.Vx.F;
    Vy.F += a.Vy.F;
    Vz.F += a.Vz.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given vector field
 *
 *          The unary operator -= subtracts a given vector field from the vfield and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to the vfield to be subtracted from the member fields
 *
 * \return  A pointer to itself is returned by the vector field object to which the operator belongs
 ********************************************************************************************************************************************
 */
vfield& vfield::operator -= (vfield &a) {
    Vx.F -= a.Vx.F;
    Vy.F -= a.Vy.F;
    Vz.F -= a.Vz.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to multiply a scalar value to the vector field
 *
 *          The unary operator *= multiplies a real value to all the fields (Vx, Vy and Vz) stored in vfield and returns
 *          a pointer to itself.
 *
 * \param   a is a real number to be multiplied to the vector field
 *
 * \return  A pointer to itself is returned by the vector field class to which the operator belongs
 ********************************************************************************************************************************************
 */
vfield& vfield::operator *= (real a) {
    Vx.F *= a;
    Vy.F *= a;
    Vz.F *= a;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a plain vector field to the vector field
 *
 *          The operator = assigns all the three blitz arrays of a plain vector field (plainvf)
 *          to the corresponding arrays in the three fields of the vfield.
 *
 * \param   a is a plainvf to be assigned to the vector field
 ********************************************************************************************************************************************
 */
void vfield::operator = (plainvf &a) {
    Vx.F = a.Vx;
    Vy.F = a.Vy;
    Vz.F = a.Vz;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign another vector field to the vector field
 *
 *          The operator = assigns all the three fields of a given vector field (vfield)
 *          to the corresponding fields of the vfield.
 *
 * \param   a is a vfield to be assigned to the vector field
 ********************************************************************************************************************************************
 */
void vfield::operator = (vfield &a) {
    Vx.F = a.Vx.F;
    Vy.F = a.Vy.F;
    Vz.F = a.Vz.F;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a scalar value to the vector field
 *
 *          The operator = assigns a real value to all the fields (Vx, Vy and Vz) stored in vfield.
 *
 * \param   a is a real number to be assigned to the vector field
 ********************************************************************************************************************************************
 */
void vfield::operator = (real a) {
    Vx.F = a;
    Vy.F = a;
    Vz.F = a;
}
