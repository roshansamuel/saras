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
/*! \file plainsf.h
 *
 *  \brief Class declaration of plainsf - plain scalar field
 *
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#ifndef PLAINSF_H
#define PLAINSF_H

#include "plainvf.h"
#include "sfield.h"
#include "vfield.h"
#include "grid.h"

class plainsf {
    private:
        blitz::firstIndex i;
        blitz::secondIndex j;
        blitz::thirdIndex k;

        const grid &gridData;

    public:
        blitz::Array<real, 3> F;

        blitz::Range xColl, yColl, zColl;

        plainsf(const grid &gridData, const sfield &refF);

        mpidata *mpiHandle;

        plainsf& operator += (plainsf &a);
        plainsf& operator -= (plainsf &a);

        plainsf& operator += (sfield &a);
        plainsf& operator -= (sfield &a);

        plainsf& operator *= (real a);

        void operator = (plainsf &a);
        void operator = (sfield &a);
        void operator = (real a);

/**
 ********************************************************************************************************************************************
 * \brief   Operator to compute the gradient of the plain scalar field
 *
 *          The gradient operator computes the gradient of the cell centered scalar field, and stores it into a face-centered staggered
 *          plain vector field as defined by the tensor operation:
 *          \f$ \nabla f = \frac{\partial f}{\partial x}i + \frac{\partial f}{\partial y}j + \frac{\partial f}{\partial z}k \f$.
 *
 * \param   gradF is a reference to a plain vector field (plainvf) into which the computed gradient must be written.
 * \param   V is a const reference to a vector field (vfield) whose core slices are used to compute gradient, since plainvf doesn't have them
 ********************************************************************************************************************************************
 */
        inline void gradient(plainvf &gradF, const vfield &V) {
            gradF.Vx(V.Vx.fCore) = gridData.xi_xColloc(xColl)(i)*(F(V.Vx.fCRgt) - F(V.Vx.fCore))/gridData.dXi;
#ifndef PLANAR
            gradF.Vy(V.Vy.fCore) = gridData.et_yColloc(yColl)(j)*(F(V.Vy.fCBak) - F(V.Vy.fCore))/gridData.dEt;
#endif
            gradF.Vz(V.Vz.fCore) = gridData.zt_zColloc(zColl)(k)*(F(V.Vz.fCTop) - F(V.Vz.fCore))/gridData.dZt;
        }

/**
 ********************************************************************************************************************************************
 * \brief   Function to synchronise data across all processors when performing parallel computations
 *
 *          This function calls the \ref mpidata#syncData "syncData" function of mpidata class to perform perform data-transfer and thus update
 *          the sub-domain boundary pads.
 ********************************************************************************************************************************************
 */
        inline void syncData() {
            mpiHandle->syncData();
        }

/**
 ********************************************************************************************************************************************
 * \brief   Function to extract the maximum value from the plain scalar field
 *
 *          The function uses the in-built blitz function to obtain the maximum value in an array.
 *          While performing parallel computation, the function performs an <B>MPI_Allreduce()</B> to get
 *          the global maximum from the entire computational domain.
 *
 * \return  The real value of the maximum is returned (it is implicitly assumed that only real values are used)
 ********************************************************************************************************************************************
 */
        inline real fxMax() {
            real localMax, globalMax;

            localMax = blitz::max(F);

            MPI_Allreduce(&localMax, &globalMax, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);

            return globalMax;
        }

        ~plainsf() { };
};

/**
 ********************************************************************************************************************************************
 *  \class plainsf plainsf.h "lib/plainsf.h"
 *  \brief Plain scalar field class to store simple scalar fields with no differentiation or interpolation
 *
 *  The class stores scalar fields in the form of a Blitz array
 ********************************************************************************************************************************************
 */

#endif
