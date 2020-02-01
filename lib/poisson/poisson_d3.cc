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
/*! \file poisson3.cc
 *
 *  \brief Definitions for functions of class poisson for 3D
 *  \sa poisson.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "poisson.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the multigrid_d3 class derived from the poisson class
 *
 *          The constructor of the derived multigrid_d3 class frst calls the base poisson class with the arguments passed to it.
 *          It then calls a series of functions in sequence to initialize all the necessary parameters and data structures to
 *          store and manipulate the multi-grid data.
 *          Since the multi-grid solver operates on the staggered grid, it first computes the limits of the full and core
 *          staggered grid, as the grid class does the same for the collocated grid.
 *
 *          It then initializes all the Range objects to obtain the correct slices of the full grid at various
 *          levels of the V-cycle.
 *          It also copies the staggered grid derivatives to local arrays with wide pads, and finally generates the MPI datatypes
 *          for data transfer between sub-domain boundaries.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   solParam is a const reference to the user-set parameters contained in the parser class
 ********************************************************************************************************************************************
 */
multigrid_d3::multigrid_d3(const grid &mesh, const parser &solParam): poisson(mesh, solParam) {
    // GET THE localSizeIndex AS IT WILL BE USED TO SET THE FULL AND CORE LIMITS OF THE STAGGERED POINTS
    setLocalSizeIndex();

    // SET THE FULL AND CORE LIMTS SET ABOVE USING THE localSizeIndex VARIABLE SET ABOVE
    setStagBounds();

    // USING THE FULL AND CORE LIMTS SET ABOVE, CREATE ALL Range OBJECTS
    initMeshRanges();

    // SET VALUES OF COEFFICIENTS USED FOR COMPUTING LAPLACIAN
    setCoefficients();

    // COPY THE STAGGERED GRID DERIVATIVES TO LOCAL ARRAYS
    copyStaggrDerivs();

    // RESIZE AND INITIALIZE NECESSARY DATA-STRUCTURES
    initializeArrays();

    // CREATE THE MPI SUB-ARRAYS NECESSARY TO TRANSFER DATA ACROSS SUB-DOMAINS AT ALL MESH LEVELS
    createMGSubArrays();
}

void multigrid_d3::mgSolve(plainsf &inFn, const plainsf &rhs) {
    pressureData = 0.0;
    residualData = 0.0;
    inputRHSData = 0.0;

    // TRANSFER DATA FROM THE INPUT SCALAR FIELD INTO THE DATA-STRUCTURES USED BY poisson
    inputRHSData(stagCore) = rhs.F(stagCore);

    // PERFORM V-CYCLES AS MANY TIMES AS REQUIRED
    for (int i=0; i<inputParams.vcCount; i++) {
        smoothedPres = 0.0;
        iteratorTemp = 0.0;

        vCycle();
    }

    // RETURN CALCULATED PRESSURE DATA
    inFn.F = pressureData(blitz::RectDomain<3>(inFn.F.lbound(), inFn.F.ubound()));
}

void multigrid_d3::vCycle() {
    vLevel = 0;

    // PRE-SMOOTHING
    swap(inputRHSData, residualData);
    smooth(inputParams.preSmooth);
    swap(residualData, inputRHSData);
    // After above 3 lines, pressureData has the pre-smoothed values of pressure, inputRHSData has original RHS data, and residualData = 0.0

    // Compute Laplacian of the pressure field and subtract it from the RHS of Poisson equation to obtain the residual
#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
    for (int iX = xStr; iX <= xEnd; iX += strideValues(vLevel)) {
        for (int iY = yStr; iY <= yEnd; iY += strideValues(vLevel)) {
            for (int iZ = zStr; iZ <= zEnd; iZ += strideValues(vLevel)) {
                residualData(iX, iY, iZ) =  inputRHSData(iX, iY, iZ) -
                                           (xix2(iX) * (pressureData(iX + strideValues(vLevel), iY, iZ) - 2.0*pressureData(iX, iY, iZ) + pressureData(iX - strideValues(vLevel), iY, iZ))/(hx(vLevel)*hx(vLevel)) +
                                            xixx(iX) * (pressureData(iX + strideValues(vLevel), iY, iZ) - pressureData(iX - strideValues(vLevel), iY, iZ))/(2.0*hx(vLevel)) +
                                            ety2(iY) * (pressureData(iX, iY + strideValues(vLevel), iZ) - 2.0*pressureData(iX, iY, iZ) + pressureData(iX, iY - strideValues(vLevel), iZ))/(hy(vLevel)*hy(vLevel)) +
                                            etyy(iY) * (pressureData(iX, iY + strideValues(vLevel), iZ) - pressureData(iX, iY - strideValues(vLevel), iZ))/(2.0*hy(vLevel)) +
                                            ztz2(iZ) * (pressureData(iX, iY, iZ + strideValues(vLevel)) - 2.0*pressureData(iX, iY, iZ) + pressureData(iX, iY, iZ - strideValues(vLevel)))/(hz(vLevel)*hz(vLevel)) +
                                            ztzz(iZ) * (pressureData(iX, iY, iZ + strideValues(vLevel)) - pressureData(iX, iY, iZ - strideValues(vLevel)))/(2.0*hz(vLevel)));
            }
        }
    }

    // Shift pressureData into smoothedPres
    swap(smoothedPres, pressureData);
    // Now pressureData = 0.0 (since smoothedPres was 0.0), and smoothedPres has the pre-smoothed values of pressure

    // RESTRICTION OPERATIONS
    for (int i=0; i<inputParams.vcDepth; i++) {
        vLevel += 1;
    }

    // SOLVE AT COARSEST MESH RESOLUTION
    solve();

    // PROLONGATION OPERATIONS BACK TO FINE MESH
    for (int i=0; i<inputParams.vcDepth; i++) {
        prolong();
        smooth(inputParams.interSmooth[i]);
    }

    pressureData += smoothedPres;

    // POST-SMOOTHING
    swap(inputRHSData, residualData);
    smooth(inputParams.postSmooth);
    swap(residualData, inputRHSData);
}

void multigrid_d3::smooth(const int smoothCount) {
#ifdef TIME_RUN
    struct timeval begin, end;
#endif

    iteratorTemp = 0.0;

    for(int n=0; n<smoothCount; n++) {
#ifdef TIME_RUN
        gettimeofday(&begin, NULL);
#endif

        // IMPOSE BOUNDARY CONDITION
        imposeBC();

#ifdef TIME_RUN
        gettimeofday(&end, NULL);
        smothTimeTran += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;

        gettimeofday(&begin, NULL);
#endif

#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
        for (int iX = xStr; iX <= xEnd; iX += strideValues(vLevel)) {
            for (int iY = yStr; iY <= yEnd; iY += strideValues(vLevel)) {
                for (int iZ = zStr; iZ <= zEnd; iZ += strideValues(vLevel)) {
                    iteratorTemp(iX, iY, iZ) = (hyhz(vLevel) * xix2(iX) * (pressureData(iX + strideValues(vLevel), iY, iZ) + pressureData(iX - strideValues(vLevel), iY, iZ))*2.0 +
                                                hyhz(vLevel) * xixx(iX) * (pressureData(iX + strideValues(vLevel), iY, iZ) - pressureData(iX - strideValues(vLevel), iY, iZ))*hx(vLevel) +
                                                hzhx(vLevel) * ety2(iY) * (pressureData(iX, iY + strideValues(vLevel), iZ) + pressureData(iX, iY - strideValues(vLevel), iZ))*2.0 +
                                                hzhx(vLevel) * etyy(iY) * (pressureData(iX, iY + strideValues(vLevel), iZ) - pressureData(iX, iY - strideValues(vLevel), iZ))*hy(vLevel) +
                                                hxhy(vLevel) * ztz2(iZ) * (pressureData(iX, iY, iZ + strideValues(vLevel)) + pressureData(iX, iY, iZ - strideValues(vLevel)))*2.0 +
                                                hxhy(vLevel) * ztzz(iZ) * (pressureData(iX, iY, iZ + strideValues(vLevel)) - pressureData(iX, iY, iZ - strideValues(vLevel)))*hz(vLevel) -
                                        2.0 * hxhyhz(vLevel) * residualData(iX, iY, iZ))/
                                        (4.0 * (hyhz(vLevel)*xix2(iX) + hzhx(vLevel)*ety2(iY) + hxhy(vLevel)*ztz2(iZ)));
                }
            }
        }

        swap(iteratorTemp, pressureData);

#ifdef TIME_RUN
        gettimeofday(&end, NULL);
        smothTimeComp += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#endif
    }

#ifdef TIME_RUN
    gettimeofday(&begin, NULL);
#endif

    imposeBC();

#ifdef TIME_RUN
    gettimeofday(&end, NULL);
    smothTimeTran += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#endif
}

void multigrid_d3::solve() {
#ifdef TIME_RUN
    struct timeval begin, end;
#endif

    int iterCount = 0;
    real tempValue;
    real localMax, globalMax;

    iteratorTemp = 0.0;

    while (true) {
#ifdef TIME_RUN
        gettimeofday(&begin, NULL);
#endif

        // GAUSS-SEIDEL ITERATIVE SOLVER - FASTEST IN BENCHMARKS
        for (int iX = xStr; iX <= xEnd; iX += strideValues(vLevel)) {
            for (int iY = yStr; iY <= yEnd; iY += strideValues(vLevel)) {
                for (int iZ = zStr; iZ <= zEnd; iZ += strideValues(vLevel)) {
                    iteratorTemp(iX, iY, iZ) = (hyhz(vLevel) * xix2(iX) * (pressureData(iX + strideValues(vLevel), iY, iZ) + iteratorTemp(iX - strideValues(vLevel), iY, iZ))*2.0 +
                                                hyhz(vLevel) * xixx(iX) * (pressureData(iX + strideValues(vLevel), iY, iZ) - iteratorTemp(iX - strideValues(vLevel), iY, iZ))*hx(vLevel) +
                                                hzhx(vLevel) * ety2(iY) * (pressureData(iX, iY + strideValues(vLevel), iZ) + iteratorTemp(iX, iY - strideValues(vLevel), iZ))*2.0 +
                                                hzhx(vLevel) * etyy(iY) * (pressureData(iX, iY + strideValues(vLevel), iZ) - iteratorTemp(iX, iY - strideValues(vLevel), iZ))*hy(vLevel) +
                                                hxhy(vLevel) * ztz2(iZ) * (pressureData(iX, iY, iZ + strideValues(vLevel)) + iteratorTemp(iX, iY, iZ - strideValues(vLevel)))*2.0 +
                                                hxhy(vLevel) * ztzz(iZ) * (pressureData(iX, iY, iZ + strideValues(vLevel)) - iteratorTemp(iX, iY, iZ - strideValues(vLevel)))*hz(vLevel) -
                                        2.0 * hxhyhz(vLevel) * residualData(iX, iY, iZ))/
                                      (4.0 * (hyhz(vLevel)*xix2(iX) + hzhx(vLevel)*ety2(iY) + hxhy(vLevel)*ztz2(iZ)));
                }
            }
        }

        swap(iteratorTemp, pressureData);

#ifdef TIME_RUN
        gettimeofday(&end, NULL);
        solveTimeComp += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;

        gettimeofday(&begin, NULL);
#endif

        // Only the pads within the domain at the sub-domain boundaries are updated here.
        // Boundary conditions are *NOT* applied while solving at the coarsest level.
        // Boundary conditions are applied only while smoothing the solution.
        updatePads();

#ifdef TIME_RUN
        gettimeofday(&end, NULL);
        solveTimeTran += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;

        gettimeofday(&begin, NULL);
#endif

        // Compute the Laplacian of pressure field and subtract the residual. Find the maximum of the absolute value of this difference
        tempValue = 0.0;
        localMax = -1.0e-10;

        // Problem with Koenig lookup is that when using the function abs with blitz arrays, it automatically computes
        // the absolute of the float values without hitch.
        // When replacing with computing absolute of individual array elements in a loop, ADL chooses a version of
        // abs in the STL which **rounds off** the number.
        // In this case, abs has to be replaced with fabs.
        for (int iX = xStr; iX <= xEnd; iX += strideValues(vLevel)) {
            for (int iY = yStr; iY <= yEnd; iY += strideValues(vLevel)) {
                for (int iZ = zStr; iZ <= zEnd; iZ += strideValues(vLevel)) {
                    tempValue = fabs((xix2(iX) * (pressureData(iX + strideValues(vLevel), iY, iZ) - 2.0*pressureData(iX, iY, iZ) + pressureData(iX - strideValues(vLevel), iY, iZ))/(hx(vLevel)*hx(vLevel)) +
                                      xixx(iX) * (pressureData(iX + strideValues(vLevel), iY, iZ) - pressureData(iX - strideValues(vLevel), iY, iZ))/(2.0*hx(vLevel)) +
                                      ety2(iY) * (pressureData(iX, iY + strideValues(vLevel), iZ) - 2.0*pressureData(iX, iY, iZ) + pressureData(iX, iY - strideValues(vLevel), iZ))/(hy(vLevel)*hy(vLevel)) +
                                      etyy(iY) * (pressureData(iX, iY + strideValues(vLevel), iZ) - pressureData(iX, iY - strideValues(vLevel), iZ))/(2.0*hy(vLevel)) +
                                      ztz2(iZ) * (pressureData(iX, iY, iZ + strideValues(vLevel)) - 2.0*pressureData(iX, iY, iZ) + pressureData(iX, iY, iZ - strideValues(vLevel)))/(hz(vLevel)*hz(vLevel)) +
                                      ztzz(iZ) * (pressureData(iX, iY, iZ + strideValues(vLevel)) - pressureData(iX, iY, iZ - strideValues(vLevel)))/(2.0*hz(vLevel))) - residualData(iX, iY, iZ));

                    if (tempValue > localMax) {
                        localMax = tempValue;
                    }
                }
            }
        }

        MPI_Allreduce(&localMax, &globalMax, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);
        if (globalMax < inputParams.tolerance) {
            break;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        iterCount += 1;
        if (iterCount > maxCount) {
            if (mesh.rankData.rank == 0) {
                std::cout << "ERROR: Jacobi iterations for solution at coarsest level not converging. Aborting" << std::endl;
            }
            MPI_Finalize();
            exit(0);
        }

#ifdef TIME_RUN
        gettimeofday(&end, NULL);
        solveTimeComp += ((end.tv_sec - begin.tv_sec)*1000000u + end.tv_usec - begin.tv_usec)/1.e6;
#endif
    }
}

void multigrid_d3::prolong() {
    // Integer values of starting indices, ending indices, and index increments along each direction
    int xSt, xEn, xIn;
    int ySt, yEn, yIn;
    int zSt, zEn, zIn;

    vLevel -= 1;

    // NOTE: Currently interpolating along X first, then Y and finally Z.
    // Test and see if this order is better or the other order, with Z first, then Y and X is better
    // Depending on the order of variables in memory, one of these will give better performance

    // Maybe the below values can be stored in some array instead of recomputing in each prolongation step

    // INTERPOLATE VARIABLE DATA ALONG X-DIRECTION
    xSt = stagCore.lbound(0) + strideValues(vLevel);
    xEn = stagCore.ubound(0) - strideValues(vLevel);
    xIn = strideValues(vLevel+1);

    ySt = stagCore.lbound(1);
    yEn = stagCore.ubound(1);
    yIn = strideValues(vLevel+1);

    zSt = stagCore.lbound(2);
    zEn = stagCore.ubound(2);
    zIn = strideValues(vLevel+1);

    for (int iX = xSt; iX <= xEn; iX += xIn) {
        for (int iY = ySt; iY <= yEn; iY += yIn) {
            for (int iZ = zSt; iZ <= zEn; iZ += zIn) {
                pressureData(iX, iY, iZ) = (pressureData(iX + strideValues(vLevel), iY, iZ) + pressureData(iX - strideValues(vLevel), iY, iZ))/2.0;
                residualData(iX, iY, iZ) = (residualData(iX + strideValues(vLevel), iY, iZ) + residualData(iX - strideValues(vLevel), iY, iZ))/2.0;
            }
        }
    }

    // INTERPOLATE VARIABLE DATA ALONG Y-DIRECTION
    xSt = stagCore.lbound(0);
    xEn = stagCore.ubound(0);
    xIn = strideValues(vLevel);

    ySt = stagCore.lbound(1) + strideValues(vLevel);
    yEn = stagCore.ubound(1) - strideValues(vLevel);
    yIn = strideValues(vLevel+1);

    zSt = stagCore.lbound(2);
    zEn = stagCore.ubound(2);
    zIn = strideValues(vLevel+1);

    for (int iX = xSt; iX <= xEn; iX += xIn) {
        for (int iY = ySt; iY <= yEn; iY += yIn) {
            for (int iZ = zSt; iZ <= zEn; iZ += zIn) {
                pressureData(iX, iY, iZ) = (pressureData(iX, iY + strideValues(vLevel), iZ) + pressureData(iX, iY - strideValues(vLevel), iZ))/2.0;
                residualData(iX, iY, iZ) = (residualData(iX, iY + strideValues(vLevel), iZ) + residualData(iX, iY - strideValues(vLevel), iZ))/2.0;
            }
        }
    }

    // INTERPOLATE VARIABLE DATA ALONG Z-DIRECTION
    xSt = stagCore.lbound(0);
    xEn = stagCore.ubound(0);
    xIn = strideValues(vLevel);

    ySt = stagCore.lbound(1);
    yEn = stagCore.ubound(1);
    yIn = strideValues(vLevel);

    zSt = stagCore.lbound(2) + strideValues(vLevel);
    zEn = stagCore.ubound(2) - strideValues(vLevel);
    zIn = strideValues(vLevel+1);

    for (int iX = xSt; iX <= xEn; iX += xIn) {
        for (int iY = ySt; iY <= yEn; iY += yIn) {
            for (int iZ = zSt; iZ <= zEn; iZ += zIn) {
                pressureData(iX, iY, iZ) = (pressureData(iX, iY, iZ + strideValues(vLevel)) + pressureData(iX, iY, iZ - strideValues(vLevel)))/2.0;
                residualData(iX, iY, iZ) = (residualData(iX, iY, iZ + strideValues(vLevel)) + residualData(iX, iY, iZ - strideValues(vLevel)))/2.0;
            }
        }
    }
}

void multigrid_d3::setLocalSizeIndex() {
    localSizeIndex = blitz::TinyVector<int, 3>(mesh.sizeIndex(0) - int(log2(inputParams.npX)),
                                               mesh.sizeIndex(1) - int(log2(inputParams.npY)),
                                               mesh.sizeIndex(2));
}

void multigrid_d3::setStagBounds() {
    blitz::TinyVector<int, 3> loBound, upBound;

    // LOWER BOUND AND UPPER BOUND OF STAGGERED CORE - USED TO CONSTRUCT THE CORE SLICE
    loBound = 0, 0, 0;
    upBound = mgSizeArray(localSizeIndex(0)) - 1, mgSizeArray(localSizeIndex(1)) - 1, mgSizeArray(localSizeIndex(2)) - 1;
    stagCore = blitz::RectDomain<3>(loBound, upBound);

    // LOWER BOUND AND UPPER BOUND OF STAGGERED FULL SUB-DOMAIN - USED TO CONSTRUCT THE FULL SUB-DOMAIN SLICE
    loBound = -strideValues(inputParams.vcDepth), -strideValues(inputParams.vcDepth), -strideValues(inputParams.vcDepth);
    upBound = stagCore.ubound() - loBound;
    stagFull = blitz::RectDomain<3>(loBound, upBound);
}

void multigrid_d3::setCoefficients() {
    hx.resize(inputParams.vcDepth + 1);
    hy.resize(inputParams.vcDepth + 1);
    hz.resize(inputParams.vcDepth + 1);

    hxhy.resize(inputParams.vcDepth + 1);
    hyhz.resize(inputParams.vcDepth + 1);
    hzhx.resize(inputParams.vcDepth + 1);

    hxhyhz.resize(inputParams.vcDepth + 1);

    for(int i=0; i<=inputParams.vcDepth; i++) {
        hx(i) = strideValues(i)*mesh.dXi;
        hy(i) = strideValues(i)*mesh.dEt;
        hz(i) = strideValues(i)*mesh.dZt;

        hxhy(i) = pow(strideValues(i), 4.0)*pow(mesh.dXi, 2.0)*pow(mesh.dEt, 2.0);
        hyhz(i) = pow(strideValues(i), 4.0)*pow(mesh.dEt, 2.0)*pow(mesh.dZt, 2.0);
        hzhx(i) = pow(strideValues(i), 4.0)*pow(mesh.dZt, 2.0)*pow(mesh.dXi, 2.0);

        hxhyhz(i) = pow(strideValues(i), 6.0)*pow(mesh.dXi, 2.0)*pow(mesh.dEt, 2.0)*pow(mesh.dZt, 2.0);
    }
}

void multigrid_d3::copyStaggrDerivs() {
    xixx.resize(stagFull.ubound(0) - stagFull.lbound(0) + 1);
    xixx.reindexSelf(stagFull.lbound(0));
    xixx = 0.0;
    xixx(blitz::Range(0, stagCore.ubound(0), 1)) = mesh.xixxStaggr(blitz::Range(0, stagCore.ubound(0), 1));

    xix2.resize(stagFull.ubound(0) - stagFull.lbound(0) + 1);
    xix2.reindexSelf(stagFull.lbound(0));
    xix2 = 0.0;
    xix2(blitz::Range(0, stagCore.ubound(0), 1)) = mesh.xix2Staggr(blitz::Range(0, stagCore.ubound(0), 1));

    etyy.resize(stagFull.ubound(1) - stagFull.lbound(1) + 1);
    etyy.reindexSelf(stagFull.lbound(1));
    etyy = 0.0;
    etyy(blitz::Range(0, stagCore.ubound(1), 1)) = mesh.etyyStaggr(blitz::Range(0, stagCore.ubound(1), 1));

    ety2.resize(stagFull.ubound(1) - stagFull.lbound(1) + 1);
    ety2.reindexSelf(stagFull.lbound(1));
    ety2 = 0.0;
    ety2(blitz::Range(0, stagCore.ubound(1), 1)) = mesh.ety2Staggr(blitz::Range(0, stagCore.ubound(1), 1));

    ztzz.resize(stagFull.ubound(2) - stagFull.lbound(2) + 1);
    ztzz.reindexSelf(stagFull.lbound(2));
    ztzz = 0.0;
    ztzz(blitz::Range(0, stagCore.ubound(2), 1)) = mesh.ztzzStaggr(blitz::Range(0, stagCore.ubound(2), 1));

    ztz2.resize(stagFull.ubound(2) - stagFull.lbound(2) + 1);
    ztz2.reindexSelf(stagFull.lbound(2));
    ztz2 = 0.0;
    ztz2(blitz::Range(0, stagCore.ubound(2), 1)) = mesh.ztz2Staggr(blitz::Range(0, stagCore.ubound(2), 1));
}

void multigrid_d3::initMeshRanges() {
    xMeshRange.resize(inputParams.vcDepth + 1);
    yMeshRange.resize(inputParams.vcDepth + 1);
    zMeshRange.resize(inputParams.vcDepth + 1);

    // Range OBJECTS WITH STRIDE TO ACCESS DIFFERENT POINTS OF THE SAME ARRAY AT DIFFERENT MULTI-GRID LEVELS
    for(int i=0; i<=inputParams.vcDepth; i++) {
        xMeshRange(i) = blitz::Range(stagCore.lbound(0), stagCore.ubound(0), strideValues(i));
        yMeshRange(i) = blitz::Range(stagCore.lbound(1), stagCore.ubound(1), strideValues(i));
        zMeshRange(i) = blitz::Range(stagCore.lbound(2), stagCore.ubound(2), strideValues(i));
    }

    // SET THE LIMTS FOR ARRAY LOOPS IN solve AND smooth FUNCTIONS, AND A FEW OTHER PLACES
    // WARNING: THESE VARIABLES HAVE SO FAR BEEN IMPLEMENTED ONLY IN solve, smooth AND vCycle.
    // THE TEST FUNCTIONS HAVE NOT YET BEEN UPDATED WITH THESE
    xStr = stagCore.lbound(0);
    yStr = stagCore.lbound(1);
    zStr = stagCore.lbound(2);

    xEnd = stagCore.ubound(0);
    yEnd = stagCore.ubound(1);
    zEnd = stagCore.ubound(2);
}

void multigrid_d3::createMGSubArrays() {
    int ptsCount;
    int numPoints;
    int areaVal, lengthVal;

    blitz::Array<int, 1> blockIndx, blockSize;

    recvStatus.resize(4);
    recvRequest.resize(4);

    xMGArray.resize(inputParams.vcDepth + 1);
    yMGArray.resize(inputParams.vcDepth + 1);

    mgSendLft.resize(inputParams.vcDepth + 1);        mgSendRgt.resize(inputParams.vcDepth + 1);
    mgRecvLft.resize(inputParams.vcDepth + 1);        mgRecvRgt.resize(inputParams.vcDepth + 1);
    mgSendFrn.resize(inputParams.vcDepth + 1);        mgSendBak.resize(inputParams.vcDepth + 1);
    mgRecvFrn.resize(inputParams.vcDepth + 1);        mgRecvBak.resize(inputParams.vcDepth + 1);

    for(int i=0; i<=inputParams.vcDepth; i++) {
        /**
         * For transfer of non-contiguous, yet uniformly spaced blocks of data, the MPI_Type_indexed datatype is being used here.
         * The number of such data blocks is represented by numPoints.
         * Since single points are being transferred rather than small blocks, the blockSize array, which holds the block length
         * of each block of data to be sent, is set to 1.
         * The blockIndx variable holds the starting index of each block (here, each single data point) within the full global
         * array of data.
         */

        /***************************************************************************************************
         * Previously xMGArray and yMGArray were defined only if npX > 1 or npY > 1 respectively.
         * This condition remained as a hidden bug in the code for the long time
         * Because for periodic cases, it was implicitly assumed that periodic data transfer will serve
         * But for a sequential case with npX = 1 and npY = 1, this transfer will not happen
         * Now xMGArray and yMGArray are defined irrespective of npX and npY
        \**************************************************************************************************/
        // CREATE X_MG_ARRAY DATATYPE
        numPoints = mgSizeArray(localSizeIndex(1) - i)*mgSizeArray(localSizeIndex(2) - i);
        blockIndx.resize(numPoints);
        blockSize.resize(numPoints);

        blockSize = 1;
        ptsCount = 0;

        lengthVal = (stagFull.ubound(2) - stagFull.lbound(2) + 1);
        for (int j = 0; j < mgSizeArray(localSizeIndex(1) - i); j++) {
            for (int k = 0; k < mgSizeArray(localSizeIndex(2) - i); k++) {
                blockIndx(ptsCount) = j*lengthVal*strideValues(i) + k*strideValues(i);
                ptsCount += 1;
            }
        }
        MPI_Type_indexed(numPoints, blockSize.data(), blockIndx.data(), MPI_FP_REAL, &xMGArray(i));
        MPI_Type_commit(&xMGArray(i));

        // CREATE Y_MG_ARRAY DATATYPE
        numPoints = mgSizeArray(localSizeIndex(2) - i)*mgSizeArray(localSizeIndex(0) - i);
        blockIndx.resize(numPoints);
        blockSize.resize(numPoints);

        blockSize = 1;
        ptsCount = 0;

        areaVal = (stagFull.ubound(1) - stagFull.lbound(1) + 1)*(stagFull.ubound(2) - stagFull.lbound(2) + 1);
        for (int j = 0; j < mgSizeArray(localSizeIndex(0) - i); j++) {
            for (int k = 0; k < mgSizeArray(localSizeIndex(2) - i); k++) {
                blockIndx(ptsCount) = j*strideValues(i)*areaVal + k*strideValues(i);
                ptsCount += 1;
            }
        }
        MPI_Type_indexed(numPoints, blockSize.data(), blockIndx.data(), MPI_FP_REAL, &yMGArray(i));
        MPI_Type_commit(&yMGArray(i));

        mgSendLft(i) =  strideValues(i), 0, 0;
        mgRecvLft(i) = -strideValues(i), 0, 0;
        mgSendRgt(i) = stagCore.ubound(0) - strideValues(i), 0, 0;
        mgRecvRgt(i) = stagCore.ubound(0) + strideValues(i), 0, 0;

        mgSendFrn(i) = 0,  strideValues(i), 0;
        mgRecvFrn(i) = 0, -strideValues(i), 0;
        mgSendBak(i) = 0, stagCore.ubound(1) - strideValues(i), 0;
        mgRecvBak(i) = 0, stagCore.ubound(1) + strideValues(i), 0;
    }
}

void multigrid_d3::imposeBC() {
    updatePads();

    if (not inputParams.xPer) {
        // NEUMANN BOUNDARY CONDITION ON PRESSURE AT LEFT WALL
        if (mesh.rankData.xRank == 0) {
            pressureData(-strideValues(vLevel), yMeshRange(vLevel), zMeshRange(vLevel)) = pressureData(strideValues(vLevel), yMeshRange(vLevel), zMeshRange(vLevel));
        }

        // NEUMANN BOUNDARY CONDITION ON PRESSURE AT RIGHT WALL
        if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
            pressureData(stagCore.ubound(0) + strideValues(vLevel), yMeshRange(vLevel), zMeshRange(vLevel)) = pressureData(stagCore.ubound(0) - strideValues(vLevel), yMeshRange(vLevel), zMeshRange(vLevel));
        }
    } // PERIODIC BOUNDARY CONDITIONS ARE AUTOMATICALLY IMPOSED BY PERIODIC DATA TRANSFER ACROSS PROCESSORS THROUGH updatePads()

    if (not inputParams.yPer) {
        // NEUMANN BOUNDARY CONDITION ON PRESSURE AT FRONT WALL
        if (mesh.rankData.yRank == 0) {
            pressureData(xMeshRange(vLevel), -strideValues(vLevel), zMeshRange(vLevel)) = pressureData(xMeshRange(vLevel), strideValues(vLevel), zMeshRange(vLevel));
        }

        // NEUMANN BOUNDARY CONDITION ON PRESSURE AT BACK WALL
        if (mesh.rankData.yRank == mesh.rankData.npY - 1) {
            pressureData(xMeshRange(vLevel), stagCore.ubound(1) + strideValues(vLevel), zMeshRange(vLevel)) = pressureData(xMeshRange(vLevel), stagCore.ubound(1) - strideValues(vLevel), zMeshRange(vLevel));
        }
    } // PERIODIC BOUNDARY CONDITIONS ARE AUTOMATICALLY IMPOSED BY PERIODIC DATA TRANSFER ACROSS PROCESSORS THROUGH updatePads()

    if (inputParams.zPer) {
        // PERIODIC BOUNDARY CONDITION ON PRESSURE AT BOTTOM WALL
        pressureData(xMeshRange(vLevel), yMeshRange(vLevel), -strideValues(vLevel)) = pressureData(xMeshRange(vLevel), yMeshRange(vLevel), stagCore.ubound(2) - strideValues(vLevel));

        // PERIODIC BOUNDARY CONDITION ON PRESSURE AT TOP WALL
        pressureData(xMeshRange(vLevel), yMeshRange(vLevel), stagCore.ubound(2) + strideValues(vLevel)) = pressureData(xMeshRange(vLevel), yMeshRange(vLevel), strideValues(vLevel));

    } else {
        // NEUMANN BOUNDARY CONDITION ON PRESSURE AT BOTTOM WALL
        pressureData(xMeshRange(vLevel), yMeshRange(vLevel), -strideValues(vLevel)) = pressureData(xMeshRange(vLevel), yMeshRange(vLevel), strideValues(vLevel));

        // NEUMANN BOUNDARY CONDITION ON PRESSURE AT TOP WALL
        pressureData(xMeshRange(vLevel), yMeshRange(vLevel), stagCore.ubound(2) + strideValues(vLevel)) = pressureData(xMeshRange(vLevel), yMeshRange(vLevel), stagCore.ubound(2) - strideValues(vLevel));
    }
}

void multigrid_d3::updatePads() {
    recvRequest = MPI_REQUEST_NULL;
    MPI_Irecv(&pressureData(mgRecvLft(vLevel)), 1, xMGArray(vLevel), mesh.rankData.nearRanks(0), 1, MPI_COMM_WORLD, &recvRequest(0));
    MPI_Irecv(&pressureData(mgRecvRgt(vLevel)), 1, xMGArray(vLevel), mesh.rankData.nearRanks(1), 2, MPI_COMM_WORLD, &recvRequest(1));
    MPI_Irecv(&pressureData(mgRecvFrn(vLevel)), 1, yMGArray(vLevel), mesh.rankData.nearRanks(2), 3, MPI_COMM_WORLD, &recvRequest(2));
    MPI_Irecv(&pressureData(mgRecvBak(vLevel)), 1, yMGArray(vLevel), mesh.rankData.nearRanks(3), 4, MPI_COMM_WORLD, &recvRequest(3));

    MPI_Send(&pressureData(mgSendLft(vLevel)), 1, xMGArray(vLevel), mesh.rankData.nearRanks(0), 2, MPI_COMM_WORLD);
    MPI_Send(&pressureData(mgSendRgt(vLevel)), 1, xMGArray(vLevel), mesh.rankData.nearRanks(1), 1, MPI_COMM_WORLD);
    MPI_Send(&pressureData(mgSendFrn(vLevel)), 1, yMGArray(vLevel), mesh.rankData.nearRanks(2), 4, MPI_COMM_WORLD);
    MPI_Send(&pressureData(mgSendBak(vLevel)), 1, yMGArray(vLevel), mesh.rankData.nearRanks(3), 3, MPI_COMM_WORLD);

    MPI_Waitall(4, recvRequest.dataFirst(), recvStatus.dataFirst());
}

real multigrid_d3::testProlong() {
    vLevel = 0;

    // Fill the residualData array with correct values expected after prolongation
    residualData = 0.0;
    for (int iX = stagCore.lbound(0); iX <= stagCore.ubound(0); iX += strideValues(vLevel)) {
        for (int iY = stagCore.lbound(1); iY <= stagCore.ubound(1); iY += strideValues(vLevel)) {
            for (int iZ = stagCore.lbound(2); iZ <= stagCore.ubound(2); iZ += strideValues(vLevel)) {
                residualData(iX, iY, iZ) = (mesh.rankData.rank + 1)*1000 + iX*100 + iY*10 + iZ;
            }
        }
    }

    // After going one level down the V-Cycle, populate the pressureData array with values at the corresponding stride
    vLevel += 1;
    pressureData = 0.0;
    for (int iX = stagCore.lbound(0); iX <= stagCore.ubound(0); iX += strideValues(vLevel)) {
        for (int iY = stagCore.lbound(1); iY <= stagCore.ubound(1); iY += strideValues(vLevel)) {
            for (int iZ = stagCore.lbound(2); iZ <= stagCore.ubound(2); iZ += strideValues(vLevel)) {
                pressureData(iX, iY, iZ) = (mesh.rankData.rank + 1)*1000 + iX*100 + iY*10 + iZ;
            }
        }
    }

    // Perform prolongation
    prolong();

    pressureData -= residualData;

    return blitz::max(fabs(pressureData));
}

real multigrid_d3::testTransfer() {
    real maxVal = 0.0;

    vLevel = 0;

    pressureData = 0.0;
    residualData = 0.0;

    MPI_Barrier(MPI_COMM_WORLD);
    for (int iX = stagCore.lbound(0); iX <= stagCore.ubound(0); iX += 1) {
        for (int iY = stagCore.lbound(1); iY <= stagCore.ubound(1); iY += 1) {
            for (int iZ = stagCore.lbound(2); iZ <= stagCore.ubound(2); iZ += 1) {
                pressureData(iX, iY, iZ) = (mesh.rankData.rank + 1)*1000 + iX*100 + iY*10 + iZ;
                residualData(iX, iY, iZ) = pressureData(iX, iY, iZ);
            }
        }
    }

    // EXPECTED VALUES IN THE PAD REGIONS IF DATA TRANSFER HAPPENS WITH NO HITCH
    for (int iX = 0; iX <= inputParams.vcDepth; iX++) {
        for (int iY = stagCore.lbound(1); iY <= stagCore.ubound(1); iY += strideValues(iX)) {
            for (int iZ = stagCore.lbound(2); iZ <= stagCore.ubound(2); iZ += strideValues(iX)) {
                residualData(-strideValues(iX), iY, iZ) = (mesh.rankData.nearRanks(0) + 1)*1000 + (stagCore.ubound(0) - strideValues(iX))*100 + iY*10 + iZ;
                residualData(stagCore.ubound(0) + strideValues(iX), iY, iZ) = (mesh.rankData.nearRanks(1) + 1)*1000 + strideValues(iX)*100 + iY*10 + iZ;
            }
        }
    }

    for (int iY = 0; iY <= inputParams.vcDepth; iY++) {
        for (int iX = stagCore.lbound(0); iX <= stagCore.ubound(0); iX += strideValues(iY)) {
            for (int iZ = stagCore.lbound(2); iZ <= stagCore.ubound(2); iZ += strideValues(iY)) {
                residualData(iX, -strideValues(iY), iZ) = (mesh.rankData.nearRanks(2) + 1)*1000 + iX*100 + (stagCore.ubound(1) - strideValues(iY))*10 + iZ;
                residualData(iX, stagCore.ubound(1) + strideValues(iY), iZ) = (mesh.rankData.nearRanks(3) + 1)*1000 + iX*100 + strideValues(iY)*10 + iZ;
            }
        }
    }

    for (int i=0; i<=inputParams.vcDepth; i++) {
        updatePads();
        vLevel += 1;
    }

    pressureData -= residualData;

    for (int iX = pressureData.lbound(0); iX <= pressureData.ubound(0); iX += 1) {
        for (int iY = pressureData.lbound(1); iY <= pressureData.ubound(1); iY += 1) {
            for (int iZ = stagCore.lbound(2); iZ <= stagCore.ubound(2); iZ += 1) {
                if (abs(pressureData(iX, iY, iZ)) > maxVal) {
                    maxVal = abs(pressureData(iX, iY, iZ));
                }
            }
        }
    }

    return maxVal;
}

real multigrid_d3::testPeriodic() {
    real xCoord = 0.0;
    real yCoord = 0.0;
    real zCoord = 0.0;

    vLevel = 0;

    pressureData = 0.0;
    residualData = 0.0;

    for (int iX = stagCore.lbound(0); iX <= stagCore.ubound(0); iX += strideValues(vLevel)) {
        for (int iY = stagCore.lbound(1); iY <= stagCore.ubound(1); iY += strideValues(vLevel)) {
            for (int iZ = stagCore.lbound(2); iZ <= stagCore.ubound(2); iZ += strideValues(vLevel)) {
                pressureData(iX, iY, iZ) = sin(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                           cos(2.0*M_PI*mesh.yStaggr(iY)/mesh.yLen)*
                                           cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);
                residualData(iX, iY, iZ) = pressureData(iX, iY, iZ);
            }
        }
    }

    // EXPECTED VALUES IN THE PAD REGIONS IF DATA TRANSFER HAPPENS WITH NO HITCH
    for (int iX = 0; iX <= inputParams.vcDepth; iX++) {
        for (int iY = stagCore.lbound(1); iY <= stagCore.ubound(1); iY += strideValues(iX)) {
            for (int iZ = stagCore.lbound(2); iZ <= stagCore.ubound(2); iZ += strideValues(iX)) {
                xCoord = mesh.xStaggr(stagCore.lbound(0)) - (mesh.xStaggr(stagCore.lbound(0) + strideValues(iX)) - mesh.xStaggr(stagCore.lbound(0)));
                residualData(stagCore.lbound(0) - strideValues(iX), iY, iZ) = sin(2.0*M_PI*xCoord/mesh.xLen)*
                                                                              cos(2.0*M_PI*mesh.yStaggr(iY)/mesh.yLen)*
                                                                              cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);

                xCoord = mesh.xStaggr(stagCore.ubound(0)) + (mesh.xStaggr(stagCore.ubound(0)) - mesh.xStaggr(stagCore.ubound(0) - strideValues(iX)));
                residualData(stagCore.ubound(0) + strideValues(iX), iY, iZ) = sin(2.0*M_PI*xCoord/mesh.xLen)*
                                                                              cos(2.0*M_PI*mesh.yStaggr(iY)/mesh.yLen)*
                                                                              cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);
            }
        }
    }

    for (int iY = 0; iY <= inputParams.vcDepth; iY++) {
        for (int iX = stagCore.lbound(0); iX <= stagCore.ubound(0); iX += strideValues(iY)) {
            for (int iZ = stagCore.lbound(2); iZ <= stagCore.ubound(2); iZ += strideValues(iY)) {
                yCoord = mesh.yStaggr(stagCore.lbound(1)) - (mesh.yStaggr(stagCore.lbound(1) + strideValues(iY)) - mesh.yStaggr(stagCore.lbound(1)));
                residualData(iX, stagCore.lbound(1) - strideValues(iY), iZ) = sin(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                                              cos(2.0*M_PI*yCoord/mesh.yLen)*
                                                                              cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);

                yCoord = mesh.yStaggr(stagCore.ubound(1)) + (mesh.yStaggr(stagCore.ubound(1)) - mesh.yStaggr(stagCore.ubound(1) - strideValues(iY)));
                residualData(iX, stagCore.ubound(1) + strideValues(iY), iZ) = sin(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                                              cos(2.0*M_PI*yCoord/mesh.yLen)*
                                                                              cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);
            }
        }
    }

    for (int iZ = 0; iZ <= inputParams.vcDepth; iZ++) {
        for (int iX = stagCore.lbound(0); iX <= stagCore.ubound(0); iX += strideValues(iZ)) {
            for (int iY = stagCore.lbound(1); iY <= stagCore.ubound(1); iY += strideValues(iZ)) {
                zCoord = mesh.zStaggr(stagCore.lbound(2)) - (mesh.zStaggr(stagCore.lbound(2) + strideValues(iZ)) - mesh.zStaggr(stagCore.lbound(2)));
                residualData(iX, iY, stagCore.lbound(2) - strideValues(iZ)) = sin(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                                              cos(2.0*M_PI*mesh.yStaggr(iY)/mesh.yLen)*
                                                                              cos(2.0*M_PI*zCoord/mesh.zLen);

                zCoord = mesh.zStaggr(stagCore.ubound(2)) + (mesh.zStaggr(stagCore.ubound(2)) - mesh.zStaggr(stagCore.ubound(2) - strideValues(iZ)));
                residualData(iX, iY, stagCore.ubound(2) + strideValues(iZ)) = sin(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                                              cos(2.0*M_PI*mesh.yStaggr(iY)/mesh.yLen)*
                                                                              cos(2.0*M_PI*zCoord/mesh.zLen);
            }
        }
    }

    for (int i=0; i<=inputParams.vcDepth; i++) {
        imposeBC();
        vLevel += 1;
    }

    pressureData -= residualData;

    return blitz::max(fabs(pressureData));
}

real multigrid_d3::testSolve() {
    vLevel = 0;

    pressureData = 0.0;
    residualData = 0.0;
    smoothedPres = 0.0;

    // WARNING: THE EXACT SOLUTION USED HERE ASSUMES xLen = yLen = zLen = 1.0
    for (int iX = stagCore.lbound(0); iX <= stagCore.ubound(0); iX += strideValues(vLevel)) {
        for (int iY = stagCore.lbound(1); iY <= stagCore.ubound(1); iY += strideValues(vLevel)) {
            for (int iZ = stagCore.lbound(2); iZ <= stagCore.ubound(2); iZ += strideValues(vLevel)) {
                smoothedPres(iX, iY, iZ) = sin(1.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                           cos(2.0*M_PI*mesh.yStaggr(iY)/mesh.yLen)*
                                           cos(4.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);

                residualData(iX, iY, iZ) = -21.0*M_PI*M_PI*sin(1.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                           cos(2.0*M_PI*mesh.yStaggr(iY)/mesh.yLen)*
                                                           cos(4.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);
            }
        }
    }

    solve();

    pressureData -= smoothedPres;

    return blitz::max(fabs(pressureData));
}
