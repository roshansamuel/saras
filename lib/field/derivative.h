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
/*! \file derivative.h
 *
 *  \brief Class declaration of derivative
 *
 *  \author Ali Asad, Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#ifndef DERIVATIVE_H
#define DERIVATIVE_H

#include <blitz/array.h>
#include <blitz/array/stencil-et.h>
#include <blitz/array/stencilops.h>
#include <string>

#include "field.h"
#include "grid.h"

class derivative {
    private: 
        const grid &gridData;

        const field &F;

        real invDelx, invDely, invDelz;

        blitz::firstIndex i;
        blitz::secondIndex j;
        blitz::thirdIndex k;    

        blitz::Range fullRange;

        blitz::Array<real, 1> x_Metric, y_Metric, z_Metric;
        blitz::Array<real, 1> xxMetric, yyMetric, zzMetric;
        blitz::Array<real, 1> x2Metric, y2Metric, z2Metric;

        blitz::Array<real, 3> tempMat;

    public:
        derivative(const grid &gridData, const field &F);

        void calcDerivative1_x(blitz::Array<real, 3> outputMat);
        void calcDerivative1_y(blitz::Array<real, 3> outputMat);
        void calcDerivative1_z(blitz::Array<real, 3> outputMat);

        void calcDerivative2xx(blitz::Array<real, 3> outputMat);
        void calcDerivative2yy(blitz::Array<real, 3> outputMat);
        void calcDerivative2zz(blitz::Array<real, 3> outputMat);
};

/**
 ********************************************************************************************************************************************
 *  \class derivative derivative.h "lib/derivative.h"
 *  \brief Derivative class to perform finite difference operations on the data stored in field
 *
 *  It contains functions to perform the finite difference operations with constant grid spacing.
 *  For many classes of this solver, empty destructors are removed. Refer reference [5] of README for more details.
 ********************************************************************************************************************************************
 */

#endif
