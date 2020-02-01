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
/*! \file main.cc
 *
 *  \brief Main file for test run of Saras. All test modules are called from here.
 *
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "unittest.h"
#include "alltests.h"

int rootRank;

int main() {
    // INITIALIZE MPI
    MPI_Init(NULL, NULL);

    int numCases = 0;
    blitz::Array<blitz::TinyVector<int, 4>, 1> testParams;

    // ALL PROCESSES READ THE INPUT PARAMETERS
    parser inputData;

    // INITIALIZE PARALLELIZATION DATA
    parallel mpi(inputData);

    rootRank = mpi.rank;

#ifdef PLANAR
    if (rootRank == 0) {
        std::cout << "\n\033[35m" << std::string(16, ' ') << " RUNNING TESTS FOR 2D CASE\033[0m" << std::endl << std::endl;
    }
    if (mpi.nProc > 1) {
        if (rootRank == 0) {
            std::cout << "\033[35m" << std::string(13, ' ') << " PARALLEL TESTS WITH 4 PROCESSORS\033[0m" << std::endl << std::endl;
        }

        numCases = 1;
        testParams.resize(numCases);
        testParams(0) = 7, 0, 4, 3;
    } else {
        if (rootRank == 0) {
            std::cout << "\033[35m" << std::string(14, ' ') << " SERIAL TESTS WITH 1 PROCESSOR\033[0m" << std::endl << std::endl;
        }

        numCases = 2;
        testParams.resize(numCases);
        testParams(0) = 5, 0, 4, 1;
        //testParams(0) = 4, 0, 5, 1;
        testParams(1) = 5, 0, 5, 2;
    }
#else
    if (rootRank == 0) {
        std::cout << "\n\033[35m" << std::string(16, ' ') << " RUNNING TESTS FOR 3D CASE\033[0m" << std::endl << std::endl;
    }
    if (mpi.nProc > 1) {
        if (rootRank == 0) {
            std::cout << "\033[35m" << std::string(13, ' ') << " PARALLEL TESTS WITH 4 PROCESSORS\033[0m" << std::endl << std::endl;
        }

        if (mpi.npX == 4 and mpi.npY == 1) {
            numCases = 1;
            testParams.resize(numCases);

            testParams(0) = 6, 4, 4, 2;
        } else if (mpi.npX == 1 and mpi.npY == 4) {
            numCases = 1;
            testParams.resize(numCases);

            testParams(0) = 4, 6, 4, 2;
        } else if (mpi.npX == 2 and mpi.npY == 2) {
            numCases = 3;
            testParams.resize(numCases);

            testParams(0) = 5, 5, 4, 1;
            testParams(1) = 6, 7, 4, 3;
            testParams(2) = 7, 4, 5, 2;
        }
    } else {
        if (rootRank == 0) {
            std::cout << "\033[35m" << std::string(14, ' ') << " SERIAL TESTS WITH 1 PROCESSOR\033[0m" << std::endl << std::endl;
        }

        numCases = 2;
        testParams.resize(numCases);
        testParams(0) = 4, 4, 5, 1;
        testParams(1) = 4, 5, 5, 2;
    }
#endif

    for (int i=0; i<numCases; i++) {
        inputData.xInd = testParams(i)(0);
        inputData.yInd = testParams(i)(1);
        inputData.zInd = testParams(i)(2);
        inputData.vcDepth = testParams(i)(3);

        // GRID OBJECT
        grid gridData(inputData, mpi);

        if (rootRank == 0) {
            std::cout << "\033[33m" << std::string(4, ' ') << " Test Case " << i+1 << " of " << numCases << ": "
                      << gridData.globalSize(0) << "x" << gridData.globalSize(1) << "x" << gridData.globalSize(2)
                      << " grid with V-Cycle depth: " << inputData.vcDepth << "\033[0m" << std::endl << std::endl;
        }

        differTest(gridData);
        fieldTest(gridData);
        nlinTest(gridData);
        poissonTest(gridData, inputData);
        hydroTest(gridData, inputData, mpi);
    }

    // FINALIZE AND CLEAN-UP
    MPI_Finalize();

    return 0;
}
