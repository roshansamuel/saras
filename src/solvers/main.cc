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
 *  \brief Main file for Saras. The appropriate solver is chosen and initialized here.
 *
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include <iostream>
#include "parallel.h"
#include "scalar.h"
#include "parser.h"
#include "hydro.h"
#include "grid.h"

int main() {
    struct timeval runStart, runEnd;

    // INITIALIZE MPI
    MPI_Init(NULL, NULL);

    // ALL PROCESSES READ THE INPUT PARAMETERS
    parser inputParams;

    // INITIALIZE PARALLELIZATION DATA
    parallel mpi(inputParams);

    // WRITE CONTENTS OF THE INPUT YAML FILE TO THE STANDARD I/O
    if (mpi.rank == 0) {
        inputParams.writeParams();
    }

    // INITIALIZE GRID DATA
    grid gridData(inputParams, mpi);

    gettimeofday(&runStart, NULL);

    if (inputParams.probType <= 4) {
        // SELECTION OF SOLVERS FOR HYDRODYNAMICS SIMULATIONS
        if (mpi.rank == 0) {
            if (inputParams.probType == 1) {
                std::cout << std::endl << "Solving NSE for lid-driven cavity problem" << std::endl;
            } else if (inputParams.probType == 2) {
                std::cout << std::endl << "Solving NSE for decaying flow simulation" << std::endl;
            } else if (inputParams.probType == 3) {
                std::cout << std::endl << "Solving NSE for channel flow problem" << std::endl;
            } else if (inputParams.probType == 4) {
                std::cout << std::endl << "Solving NSE for forced channel flow problem" << std::endl;
            }
            std::cout << std::endl;
        }

        hydro *nseSolver;

        // CREATE NEW INSTANCE OF THE HYDRODYNAMICS SOLVER
#ifdef PLANAR
        nseSolver = new hydro_d2(gridData, inputParams, mpi);
#else
        nseSolver = new hydro_d3(gridData, inputParams, mpi);
#endif

        nseSolver->solvePDE();

        delete nseSolver;

    } else if (inputParams.probType <= 7) {
        // SELECTION OF SOLVERS FOR SCALAR SIMULATIONS
        if (mpi.rank == 0) {
            if (inputParams.probType == 5) {
                std::cout << std::endl << "Solving NSE for heated bottom-plate problem" << std::endl;
            } else if (inputParams.probType == 6) {
                std::cout << std::endl << "Solving NSE for heated top-plate problem" << std::endl;
            } else {
                if (not inputParams.xPer) {
                    std::cout << std::endl << "Solving NSE for heated sidewall problem" << std::endl;
                } else {
                    std::cout << std::endl << "ERROR: X direction cannot be periodic for heated sidewall problem. ABORTING" << std::endl;

                    MPI_Finalize();
                    exit(0);
                }
            }
            std::cout << std::endl;
        }

        scalar *nseSolver;

        // CREATE NEW INSTANCE OF THE SCALAR SOLVER
#ifdef PLANAR
        nseSolver = new scalar_d2(gridData, inputParams, mpi);
#else
        nseSolver = new scalar_d3(gridData, inputParams, mpi);
#endif

        nseSolver->solvePDE();

        delete nseSolver;

    }    
    
    else {
        if (mpi.rank == 0) {
            std::cout << std::endl << "Invalid problem type. ABORTING" << std::endl;
        }

        MPI_Finalize();
        exit(0);
    }

    gettimeofday(&runEnd, NULL);
    real run_time = ((runEnd.tv_sec - runStart.tv_sec)*1000000u + runEnd.tv_usec - runStart.tv_usec)/1.e6;

    if (mpi.rank == 0) {
        std::cout << std::endl << "Simulation completed" << std::endl;
        std::cout << std::endl;
        std::cout << "Time taken by simulation: " << run_time << std::endl;
    }

    // FINALIZE AND CLEAN-UP
    MPI_Finalize();

    return 0;
}
