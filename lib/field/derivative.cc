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
/*! \file derivative.cc
 *
 *  \brief Definitions for functions of class derivative
 *  \sa derivative.h
 *  \author Ali Asad, Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "derivative.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the derivative class
 *
 *          The empty constructor of the derivative class only serves to assign values to the two const parameters of the differ class,
 *          namely <B>grid</B> and <B>F</B>.
 *
 * \param   gridData is a const reference to the global data contained in the grid class
 * \param   F is a reference to field on which the derivative operations are to be performed
 ********************************************************************************************************************************************
 */

derivative::derivative(const grid &gridData, const field &F): gridData(gridData), F(F) { 
    tempMat.resize(F.fSize);
    tempMat.reindexSelf(F.flBound);

    invDelx = 1.0/gridData.dXi;
    invDely = 1.0/gridData.dEt;
    invDelz = 1.0/gridData.dZt; 

    fullRange = blitz::Range::all();

    if (F.xStag) {
        x_Metric.reference(gridData.xi_xStaggr);
        xxMetric.reference(gridData.xixxStaggr);
        x2Metric.reference(gridData.xix2Staggr);
    } else {
        x_Metric.reference(gridData.xi_xColloc);
        xxMetric.reference(gridData.xixxColloc);
        x2Metric.reference(gridData.xix2Colloc);
    }

    if (F.yStag) {
        y_Metric.reference(gridData.et_yStaggr);
        yyMetric.reference(gridData.etyyStaggr);
        y2Metric.reference(gridData.ety2Staggr);
    } else {
        y_Metric.reference(gridData.et_yColloc);
        yyMetric.reference(gridData.etyyColloc);
        y2Metric.reference(gridData.ety2Colloc);
    }

    if (F.zStag) {
        z_Metric.reference(gridData.zt_zStaggr);
        zzMetric.reference(gridData.ztzzStaggr);
        z2Metric.reference(gridData.ztz2Staggr);
    } else {
        z_Metric.reference(gridData.zt_zColloc);
        zzMetric.reference(gridData.ztzzColloc);
        z2Metric.reference(gridData.ztz2Colloc);
    }

    tempMat = 0.0;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the first derivative of the field with respect to x [second order accurate central difference]
 *
 *          This function must be called using an output array whose shape and size should be same as that of the field
 *          
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative1_x(blitz::Array<real, 3> outputMat) {
    outputMat(blitz::Range(0, F.F.ubound(0) - 1), fullRange, fullRange) = central12n(F.F, 0);
    outputMat *= invDelx;

    outputMat = x_Metric(i)*outputMat(i, j, k);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the first derivative of the field with respect to y [second order accurate central difference]
 *
 *          This function must be called using an output array whose shape and size should be same as that of the field
 *          
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative1_y(blitz::Array<real, 3> outputMat) {
    outputMat(fullRange, blitz::Range(0, F.F.ubound(1) - 1), fullRange) = central12n(F.F, 1);
    outputMat *= invDely;

    outputMat = y_Metric(j)*outputMat(i, j, k);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the first derivative of the field with respect to z [second order accurate central difference]
 *
 *          This function must be called using an output array whose shape and size should be same as that of the field
 *          
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative1_z(blitz::Array<real, 3> outputMat) {
    outputMat(fullRange, fullRange, blitz::Range(0, F.F.ubound(2) - 1)) = central12n(F.F, 2);
    outputMat *= invDelz;

    outputMat = z_Metric(k)*outputMat(i, j, k);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the second derivative of the field with respect to x [second order accurate central difference]
 *
 *          This function must be called using an output array whose shape and size should be same as that of the field
 *          
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative2xx(blitz::Array<real, 3> outputMat) {
    tempMat(blitz::Range(0, F.F.ubound(0) - 1), fullRange, fullRange) = central12n(F.F, 0);
    tempMat *= invDelx;

    outputMat(blitz::Range(0, F.F.ubound(0) - 1), fullRange, fullRange) = central22n(F.F, 0);
    outputMat *= invDelx*invDelx;

    if (gridData.inputParams.iScheme == 1) {
        outputMat = xxMetric(i)*tempMat(i, j, k) + 0.5*x2Metric(i)*outputMat(i, j, k);
    } else {
        outputMat = xxMetric(i)*tempMat(i, j, k) + x2Metric(i)*outputMat(i, j, k);
    }
    
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the second derivatives of the field with respect to y
 *
 *          This function must be called using an output array whose shape and size should be same as that of the field
 *          
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative2yy(blitz::Array<real, 3> outputMat) {
    tempMat(fullRange, blitz::Range(0, F.F.ubound(1) - 1), fullRange) = central12n(F.F, 1);
    tempMat *= invDely;

    outputMat(fullRange, blitz::Range(0, F.F.ubound(1) - 1), fullRange) = central22n(F.F, 1);
    outputMat *= invDely*invDely;
    
    if (gridData.inputParams.iScheme == 1) {
        outputMat = yyMetric(j)*tempMat(i, j, k) + 0.5*y2Metric(j)*outputMat(i, j, k);
    } else {
        outputMat = yyMetric(j)*tempMat(i, j, k) + y2Metric(j)*outputMat(i, j, k);
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the second derivatives of the field with respect to z
 *
 *          This function must be called using an output array whose shape and size should be same as that of the field
 *          
 ********************************************************************************************************************************************
 */
void derivative::calcDerivative2zz(blitz::Array<real, 3> outputMat) {
    tempMat(fullRange, fullRange, blitz::Range(0, F.F.ubound(2) - 1)) = central12n(F.F, 2);
    tempMat *= invDelz;

    outputMat(fullRange, fullRange, blitz::Range(0, F.F.ubound(2) - 1)) = central22n(F.F, 2);
    outputMat *= invDelz*invDelz;
    
    if (gridData.inputParams.iScheme == 1) {
        outputMat = zzMetric(k)*tempMat(i, j, k) + 0.5*z2Metric(k)*outputMat(i, j, k);
    } else {
        outputMat = zzMetric(k)*tempMat(i, j, k) + z2Metric(k)*outputMat(i, j, k);
    }
}

/**
*********************************************************************************************************************************************
*    FUNCTIONS FOR CROSS DERIVATIVES CAN BE DEVELOPED IN THIS CLASS (SEE DOCUMENTATION FOR MORE DETAILS)
*
*********************************************************************************************************************************************
*/
