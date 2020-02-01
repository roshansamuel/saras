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
/*! \file unittest.cc
 *
 *  \brief Definitions for global functions for unit tests
 *  \sa unittest.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "unittest.h"

void testError(blitz::Array<real, 3> A, blitz::Array<real, 3> B, int errorMom, real errorTol) {
    int errorCnt;
    real errorVal;

    errorCnt = 0;
    errorVal = 0.0;

    if (errorMom == 1) {
        for (int i=A.lbound(0)+1; i <= A.ubound(0)-1; i++) {
            for (int j=A.lbound(1)+1; j <= A.ubound(1)-1; j++) {
                for (int k=A.lbound(2)+1; k <= A.ubound(2)-1; k++) {
                    errorCnt += 1;
                    errorVal += fabs(A(i, j, k) - B(i, j, k));
                }
            }
        }

        errorVal /= errorCnt;

    } else if (errorMom == 2) {
        for (int i=A.lbound(0)+1; i <= A.ubound(0)-1; i++) {
            for (int j=A.lbound(1)+1; j <= A.ubound(1)-1; j++) {
                for (int k=A.lbound(2)+1; k <= A.ubound(2)-1; k++) {
                    errorCnt += 1;
                    errorVal += pow((A(i, j, k) - B(i, j, k)), 2.0);
                }
            }
        }

        errorVal /= errorCnt;
        errorVal = sqrt(errorVal);
    }

    printResult(errorVal, errorTol);
}

void printResult(real computedValue, real errorTolerance) {
    if (computedValue > errorTolerance) {
        if (rootRank == 0) {
            std::cout << "\033[31m[FAILED]\033[0m" << std::endl << std::endl;
        }
    } else {
        if (rootRank == 0) {
            std::cout << "\033[32m[PASSED]\033[0m" << std::endl << std::endl;
        }
    }
}
