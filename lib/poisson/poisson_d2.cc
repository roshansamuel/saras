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
/*! \file poisson2.cc
 *
 *  \brief Definitions for functions of class poisson for 2D
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
 * \brief   Constructor of the multigrid_d2 class derived from the poisson class
 *
 *          The constructor of the derived multigrid_d2 class frst calls the base poisson class with the arguments passed to it.
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
multigrid_d2::multigrid_d2(const grid &mesh, const parser &solParam): poisson(mesh, solParam) {
    // GET THE localSizeIndex AS IT WILL BE USED TO SET THE FULL AND CORE LIMITS OF THE STAGGERED POINTS
    setLocalSizeIndex();

    // SET THE FULL AND CORE LIMTS SET ABOVE USING THE localSizeIndex VARIABLE SET ABOVE
    setStagBounds();

    // SET VALUES OF COEFFICIENTS USED FOR COMPUTING LAPLACIAN
    setCoefficients();

    // COPY THE STAGGERED GRID DERIVATIVES TO LOCAL ARRAYS
    copyStaggrDerivs();

    // RESIZE AND INITIALIZE NECESSARY DATA-STRUCTURES
    initializeArrays();

    // CREATE THE MPI SUB-ARRAYS NECESSARY TO TRANSFER DATA ACROSS SUB-DOMAINS AT ALL MESH LEVELS
    createMGSubArrays();

    // INITIALIZE DIRICHLET BCs WHEN TESTING THE POISSON SOLVER
#ifdef TEST_POISSON
    initDirichlet();
#endif
}


void multigrid_d2::computeResidual() {
    tmpDataArray(vLevel) = 0.0;

    // Compute Laplacian of the pressure field and subtract it from the RHS of Poisson equation to obtain the residual
    // This residual is temporarily stored into tmpDataArray, from which it will be coarsened into residualData array.
#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
    for (int i = 0; i <= xEnd(vLevel); ++i) {
        for (int k = 0; k <= zEnd(vLevel); ++k) {
            tmpDataArray(vLevel)(i, 0, k) =  residualData(vLevel)(i, 0, k) -
                         (xix2(vLevel)(i) * (pressureData(vLevel)(i + 1, 0, k) - 2.0*pressureData(vLevel)(i, 0, k) + pressureData(vLevel)(i - 1, 0, k))/(hx(vLevel)*hx(vLevel)) +
                          xixx(vLevel)(i) * (pressureData(vLevel)(i + 1, 0, k) - pressureData(vLevel)(i - 1, 0, k))/(2.0*hx(vLevel)) +
                          ztz2(vLevel)(k) * (pressureData(vLevel)(i, 0, k + 1) - 2.0*pressureData(vLevel)(i, 0, k) + pressureData(vLevel)(i, 0, k - 1))/(hz(vLevel)*hz(vLevel)) +
                          ztzz(vLevel)(k) * (pressureData(vLevel)(i, 0, k + 1) - pressureData(vLevel)(i, 0, k - 1))/(2.0*hz(vLevel)));
        }
    }

    updatePads(tmpDataArray);
}


void multigrid_d2::smooth(const int smoothCount) {
    tmpDataArray(vLevel) = 0.0;

    for(int n=0; n<smoothCount; ++n) {
        imposeBC();

        // WARNING: When using the gauss-seidel smoothing as written below, the edges of interior sub-domains after MPI decomposition will not have the updated values
        // As a result, the serial and parallel results will not match when using gauss-seidel smoothing
        if (inputParams.gsSmooth) {
            // GAUSS-SEIDEL ITERATIVE SMOOTHING
            for (int i = 0; i <= xEnd(vLevel); ++i) {
                for (int k = 0; k <= zEnd(vLevel); ++k) {
                    pressureData(vLevel)(i, 0, k) = (hz2(vLevel) * xix2(vLevel)(i) * (pressureData(vLevel)(i + 1, 0, k) + pressureData(vLevel)(i - 1, 0, k))*2.0 +
                                                     hz2(vLevel) * xixx(vLevel)(i) * (pressureData(vLevel)(i + 1, 0, k) - pressureData(vLevel)(i - 1, 0, k))*hx(vLevel) +
                                                     hx2(vLevel) * ztz2(vLevel)(k) * (pressureData(vLevel)(i, 0, k + 1) + pressureData(vLevel)(i, 0, k - 1))*2.0 +
                                                     hx2(vLevel) * ztzz(vLevel)(k) * (pressureData(vLevel)(i, 0, k + 1) - pressureData(vLevel)(i, 0, k - 1))*hz(vLevel) -
                                              2.0 * hzhx(vLevel) * residualData(vLevel)(i, 0, k))/
                                             (4.0 * (hz2(vLevel) * xix2(vLevel)(i) + hx2(vLevel)*ztz2(vLevel)(k)));
                }
            }
        } else {
            // JACOBI ITERATIVE SMOOTHING
#pragma omp parallel for num_threads(inputParams.nThreads) default(none)
            for (int i = 0; i <= xEnd(vLevel); ++i) {
                for (int k = 0; k <= zEnd(vLevel); ++k) {
                    tmpDataArray(vLevel)(i, 0, k) = (hz2(vLevel) * xix2(vLevel)(i) * (pressureData(vLevel)(i + 1, 0, k) + pressureData(vLevel)(i - 1, 0, k))*2.0 +
                                                     hz2(vLevel) * xixx(vLevel)(i) * (pressureData(vLevel)(i + 1, 0, k) - pressureData(vLevel)(i - 1, 0, k))*hx(vLevel) +
                                                     hx2(vLevel) * ztz2(vLevel)(k) * (pressureData(vLevel)(i, 0, k + 1) + pressureData(vLevel)(i, 0, k - 1))*2.0 +
                                                     hx2(vLevel) * ztzz(vLevel)(k) * (pressureData(vLevel)(i, 0, k + 1) - pressureData(vLevel)(i, 0, k - 1))*hz(vLevel) -
                                              2.0 * hzhx(vLevel) * residualData(vLevel)(i, 0, k))/
                                             (4.0 * (hz2(vLevel) * xix2(vLevel)(i) + hx2(vLevel)*ztz2(vLevel)(k)));
                }
            }

            swap(tmpDataArray, pressureData);
        }
    }

    imposeBC();
}


void multigrid_d2::solve() {
    int iterCount = 0;
    real tempValue, localMax, globalMax;

    while (true) {
        imposeBC();

        // GAUSS-SEIDEL ITERATIVE SOLVER
        for (int i = 0; i <= xEnd(vLevel); ++i) {
            for (int k = 0; k <= zEnd(vLevel); ++k) {
                pressureData(vLevel)(i, 0, k) = (hz2(vLevel) * xix2(vLevel)(i) * (pressureData(vLevel)(i + 1, 0, k) + pressureData(vLevel)(i - 1, 0, k))*2.0 +
                                                 hz2(vLevel) * xixx(vLevel)(i) * (pressureData(vLevel)(i + 1, 0, k) - pressureData(vLevel)(i - 1, 0, k))*hx(vLevel) +
                                                 hx2(vLevel) * ztz2(vLevel)(k) * (pressureData(vLevel)(i, 0, k + 1) + pressureData(vLevel)(i, 0, k - 1))*2.0 +
                                                 hx2(vLevel) * ztzz(vLevel)(k) * (pressureData(vLevel)(i, 0, k + 1) - pressureData(vLevel)(i, 0, k - 1))*hz(vLevel) -
                                          2.0 * hzhx(vLevel) * residualData(vLevel)(i, 0, k))/
                                         (4.0 * (hz2(vLevel) * xix2(vLevel)(i) + hx2(vLevel)*ztz2(vLevel)(k)));
            }
        }

        tempValue = 0.0;
        localMax = -1.0e-10;
        for (int i = 0; i <= xEnd(vLevel); ++i) {
            for (int k = 0; k <= zEnd(vLevel); ++k) {
                tempValue =  fabs(residualData(vLevel)(i, 0, k) -
                            (xix2(vLevel)(i) * (pressureData(vLevel)(i + 1, 0, k) - 2.0*pressureData(vLevel)(i, 0, k) + pressureData(vLevel)(i - 1, 0, k))/(hx(vLevel)*hx(vLevel)) +
                             xixx(vLevel)(i) * (pressureData(vLevel)(i + 1, 0, k) - pressureData(vLevel)(i - 1, 0, k))/(2.0*hx(vLevel)) +
                             ztz2(vLevel)(k) * (pressureData(vLevel)(i, 0, k + 1) - 2.0*pressureData(vLevel)(i, 0, k) + pressureData(vLevel)(i, 0, k - 1))/(hz(vLevel)*hz(vLevel)) +
                             ztzz(vLevel)(k) * (pressureData(vLevel)(i, 0, k + 1) - pressureData(vLevel)(i, 0, k - 1))/(2.0*hz(vLevel))));

                if (tempValue > localMax) {
                    localMax = tempValue;
                }
            }
        }

        MPI_Allreduce(&localMax, &globalMax, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);

        if (globalMax < inputParams.mgTolerance) {
            break;
        }

        iterCount += 1;
        if (iterCount > maxCount) {
            if (inputParams.printResidual) if (mesh.rankData.rank == 0) std::cout << "WARNING: Iterations for solution at coarsest level not converging." << std::endl;
            break;
        }
    }

    imposeBC();
}


void multigrid_d2::coarsen() {
    real facePoints, vertPoints;

    int i2, k2;
    int pLevel;

    pLevel = vLevel;
    vLevel += 1;

    // Full weighted restriction operation
    // The residual computed at previous vLevel is stored in tmpDataArray.
    // This data is read for coarsening and written into residualData array.

    /*
     * According to An Introduction to Multigrid Methods by P. Wesseling, Page 64 (Sec 5.2),
     * Restriction can be performed at the edges and corners using the same stencil as in the bulk,
     * But by assuming that values of the field outside the domain are all 0.
     * Hence no special treatments at the corners and edges are needed.
     */

    for (int i = 0; i <= xEnd(vLevel); ++i) {
        i2 = i*2;
        for (int k = 0; k <= zEnd(vLevel); ++k) {
            k2 = k*2;
            facePoints = (tmpDataArray(pLevel)(i2 + 1, 0, k2) + tmpDataArray(pLevel)(i2 - 1, 0, k2) +
                          tmpDataArray(pLevel)(i2, 0, k2 + 1) + tmpDataArray(pLevel)(i2, 0, k2 - 1))*0.125;
            vertPoints = (tmpDataArray(pLevel)(i2 + 1, 0, k2 + 1) +
                          tmpDataArray(pLevel)(i2 + 1, 0, k2 - 1) +
                          tmpDataArray(pLevel)(i2 - 1, 0, k2 + 1) +
                          tmpDataArray(pLevel)(i2 - 1, 0, k2 - 1))*0.0625;

            residualData(vLevel)(i, 0, k) = facePoints + vertPoints + tmpDataArray(pLevel)(i2, 0, k2)*0.25;
        }
    }
}


void multigrid_d2::prolong() {
    int pLevel;
    int i2, k2;

    pLevel = vLevel;
    vLevel -= 1;

    pressureData(vLevel) = 0.0;

    for (int i = 0; i <= xEnd(vLevel); ++i) {
        i2 = i/2;
        if (isOdd(i)) {
            for (int k = 0; k <= zEnd(vLevel); ++k) {
                k2 = k/2;
                if (isOdd(k)) { // Both i and k are odd
                    pressureData(vLevel)(i, 0, k) = (pressureData(pLevel)(i2, 0, k2)     + pressureData(pLevel)(i2, 0, k2 + 1) +
                                                     pressureData(pLevel)(i2 + 1, 0, k2) + pressureData(pLevel)(i2 + 1, 0, k2 + 1))/4.0;
                } else {        // Here i is odd, but k is even
                    pressureData(vLevel)(i, 0, k) = (pressureData(pLevel)(i2, 0, k2) + pressureData(pLevel)(i2 + 1, 0, k2))/2.0;
                }
            }
        } else {
            for (int k = 0; k <= zEnd(vLevel); ++k) {
                k2 = k/2;
                if (isOdd(k)) { // Here i is even, but k is odd
                    pressureData(vLevel)(i, 0, k) = (pressureData(pLevel)(i2, 0, k2) + pressureData(pLevel)(i2, 0, k2 + 1))/2.0;
                } else {        // Both i and k are even
                    pressureData(vLevel)(i, 0, k) = pressureData(pLevel)(i2, 0, k2);
                }
            }
        }
    }
}


real multigrid_d2::computeError(const int normOrder) {
    real residualVal = 0.0;

    real tempValue = 0.0;
    real numValLoc = 0.0;
    real denValLoc = 0.0;
    int valCountLoc = 0;

    // This function is called at the finest grid level only.
    // Moreover it called only under the TEST_POISSON flag
    // Hence it is not written to be very fast

    // Problem with Koenig lookup is that when using the function abs with blitz arrays, it automatically computes
    // the absolute of the float values without hitch.
    // When replacing with computing absolute of individual array elements in a loop, ADL chooses a version of
    // abs in the STL which **rounds off** the number.
    // In this case, abs has to be replaced with fabs.
    for (int i = 0; i <= xEnd(0); ++i) {
        for (int k = 0; k <= zEnd(0); ++k) {
            tempValue = fabs((xix2(0)(i) * (pressureData(0)(i + 1, 0, k) - 2.0*pressureData(0)(i, 0, k) + pressureData(0)(i - 1, 0, k))/hx2(0) +
                              xixx(0)(i) * (pressureData(0)(i + 1, 0, k) - pressureData(0)(i - 1, 0, k))/(2.0*hx(0)) +
                              ztz2(0)(k) * (pressureData(0)(i, 0, k + 1) - 2.0*pressureData(0)(i, 0, k) + pressureData(0)(i, 0, k - 1))/hz2(0) +
                              ztzz(0)(k) * (pressureData(0)(i, 0, k + 1) - pressureData(0)(i, 0, k - 1))/(2.0*hz(0))) - residualData(0)(i, 0, k));

            switch (normOrder) {
                case 1:
                    if (tempValue > numValLoc) numValLoc = tempValue;
                    break;
                case 2:
                    numValLoc += tempValue*tempValue;
                    denValLoc += residualData(0)(i, 0, k)*residualData(0)(i, 0, k);
                    valCountLoc += 1;
                    break;
            }
        }
    }

    real numValGlo = 0.0;
    real denValGlo = 0.0;
    int valCountGlo = 0;
    switch (normOrder) {
        case 1:
            denValLoc = blitz::max(fabs(residualData(0)));
            MPI_Allreduce(&numValLoc, &numValGlo, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(&denValLoc, &denValGlo, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);
            if (denValGlo) {
                residualVal = numValGlo/denValGlo;
            } else {
                residualVal = numValGlo;
            }
            break;
        case 2:
            MPI_Allreduce(&numValLoc, &numValGlo, 1, MPI_FP_REAL, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&denValLoc, &denValGlo, 1, MPI_FP_REAL, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&valCountLoc, &valCountGlo, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (denValGlo) {
                residualVal = sqrt(numValGlo/valCountGlo)/sqrt(denValGlo/valCountGlo);
            } else {
                residualVal = sqrt(numValGlo/valCountGlo);
            }
            break;
    }

    return residualVal;
}


void multigrid_d2::createMGSubArrays() {
    int count;

    recvStatus.resize(2);
    recvRequest.resize(2);

    xMGArray.resize(inputParams.vcDepth + 1);
    mgSendLft.resize(inputParams.vcDepth + 1);        mgSendRgt.resize(inputParams.vcDepth + 1);
    mgRecvLft.resize(inputParams.vcDepth + 1);        mgRecvRgt.resize(inputParams.vcDepth + 1);

    for(int n=0; n<=inputParams.vcDepth; ++n) {
        // CREATE X_MG_ARRAY DATATYPE
        count = stagFull(n).ubound(2) + 2;

        MPI_Type_contiguous(count, MPI_FP_REAL, &xMGArray(n));
        MPI_Type_commit(&xMGArray(n));

        // SET STARTING INDICES OF MEMORY LOCATIONS FROM WHERE TO READ (SEND) AND WRITE (RECEIVE) DATA
        mgSendLft(n) =  1, 0, -1;
        mgRecvLft(n) = -1, 0, -1;
        mgSendRgt(n) = stagCore(n).ubound(0) - 1, 0, -1;
        mgRecvRgt(n) = stagCore(n).ubound(0) + 1, 0, -1;
    }
}


void multigrid_d2::initDirichlet() {
    real xDist, zDist;

    // Generate the walls as 1D Blitz arrays
    xWall.resize(stagFull(0).ubound(2) - stagFull(0).lbound(2) + 1);
    xWall.reindexSelf(stagFull(0).lbound(2));
    xWall = 0.0;

    zWall.resize(stagFull(0).ubound(0) - stagFull(0).lbound(0) + 1);
    zWall.reindexSelf(stagFull(0).lbound(0));
    zWall = 0.0;

    // Compute values at the walls using the (r^2)/4 formula
    // Along X-direction - Left and Right Walls
    xDist = hx(0) + mesh.inputParams.Lx/2.0;

    for (int k=0; k<=stagCore(0).ubound(2); ++k) {
        zDist = hz(0)*(k - stagCore(0).ubound(2)/2);

        xWall(k) = (xDist*xDist + zDist*zDist)/4.0;
    }

    // Along Z-direction - Top and Bottom Walls
    zDist = hz(0) + mesh.inputParams.Lz/2.0;

    // In parallel runs, the domain is divided into slabs along X-axis
    // Hence some adjustments have to be made to get the right extents
    int halfIndX = stagCore(0).ubound(0)*mesh.rankData.npX/2;
    for (int i=0; i<=stagCore(0).ubound(0); ++i) {
        xDist = hx(0)*(mesh.rankData.xRank*stagCore(0).ubound(0) + i - halfIndX);

        zWall(i) = (xDist*xDist + zDist*zDist)/4.0;
    }
}


void multigrid_d2::imposeBC() {
    updatePads(pressureData);

    if (not inputParams.xPer) {
#ifdef TEST_POISSON
        // DIRICHLET BOUNDARY CONDITION AT LEFT AND RIGHT WALLS
        if (zeroBC) {
            if (mesh.rankData.xRank == 0) {
                pressureData(vLevel)(-1, 0, all) = -pressureData(vLevel)(1, 0, all);
            }

            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
                pressureData(vLevel)(stagCore(vLevel).ubound(0) + 1, 0, all) = -pressureData(vLevel)(stagCore(vLevel).ubound(0) - 1, 0, all);
            }
        } else {
            if (mesh.rankData.xRank == 0) {
                pressureData(vLevel)(-1, 0, all) = xWall(all);
            }

            if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
                pressureData(vLevel)(stagCore(vLevel).ubound(0) + 1, 0, all) = xWall(all);
            }
        }
#else
        // NEUMANN BOUNDARY CONDITION AT LEFT AND RIGHT WALLS
        if (mesh.rankData.xRank == 0) {
            pressureData(vLevel)(-1, 0, all) = pressureData(vLevel)(1, 0, all);
        }

        if (mesh.rankData.xRank == mesh.rankData.npX - 1) {
            pressureData(vLevel)(stagCore(vLevel).ubound(0) + 1, 0, all) = pressureData(vLevel)(stagCore(vLevel).ubound(0) - 1, 0, all);
        }
#endif
    } // PERIODIC BOUNDARY CONDITIONS ARE AUTOMATICALLY IMPOSED BY PERIODIC DATA TRANSFER ACROSS PROCESSORS THROUGH updatePads()

    if (inputParams.zPer) {
        // PERIODIC BOUNDARY CONDITION AT BOTTOM WALL
        pressureData(vLevel)(all, 0, -1) = pressureData(vLevel)(all, 0, stagCore(vLevel).ubound(2) - 1);

        // PERIODIC BOUNDARY CONDITION AT TOP WALL
        pressureData(vLevel)(all, 0, stagCore(vLevel).ubound(2) + 1) = pressureData(vLevel)(all, 0, 1);

    } else {
#ifdef TEST_POISSON
        // DIRICHLET BOUNDARY CONDITION AT BOTTOM AND TOP WALLS
        if (zeroBC) {
            pressureData(vLevel)(all, 0, -1) = -pressureData(vLevel)(all, 0, 1);

            pressureData(vLevel)(all, 0, stagCore(vLevel).ubound(2) + 1) = -pressureData(vLevel)(all, 0, stagCore(vLevel).ubound(2) - 1);
        } else {
            pressureData(vLevel)(all, 0, -1) = zWall(all);

            pressureData(vLevel)(all, 0, stagCore(vLevel).ubound(2) + 1) = zWall(all);
        }
#else
        // NEUMANN BOUNDARY CONDITION AT BOTTOM AND TOP WALLS
        pressureData(vLevel)(all, 0, -1) = pressureData(vLevel)(all, 0, 1);

        pressureData(vLevel)(all, 0, stagCore(vLevel).ubound(2) + 1) = pressureData(vLevel)(all, 0, stagCore(vLevel).ubound(2) - 1);
#endif
    }
}


void multigrid_d2::updatePads(blitz::Array<blitz::Array<real, 3>, 1> &data) {
    recvRequest = MPI_REQUEST_NULL;

    // TRANSFER DATA FROM NEIGHBOURING CELL TO IMPOSE SUB-DOMAIN BOUNDARY CONDITIONS
    MPI_Irecv(&(data(vLevel)(mgRecvLft(vLevel))), 1, xMGArray(vLevel), mesh.rankData.nearRanks(0), 1, MPI_COMM_WORLD, &recvRequest(0));
    MPI_Irecv(&(data(vLevel)(mgRecvRgt(vLevel))), 1, xMGArray(vLevel), mesh.rankData.nearRanks(1), 2, MPI_COMM_WORLD, &recvRequest(1));

    MPI_Send(&(data(vLevel)(mgSendLft(vLevel))), 1, xMGArray(vLevel), mesh.rankData.nearRanks(0), 2, MPI_COMM_WORLD);
    MPI_Send(&(data(vLevel)(mgSendRgt(vLevel))), 1, xMGArray(vLevel), mesh.rankData.nearRanks(1), 1, MPI_COMM_WORLD);

    MPI_Waitall(2, recvRequest.dataFirst(), recvStatus.dataFirst());
}


real multigrid_d2::testProlong() {
    vLevel = 0;

    // Fill the residualData array with correct values expected after prolongation
    residualData(vLevel) = 0.0;
    for (int i = 0; i <= xEnd(vLevel); ++i) {
        for (int k = 0; k <= zEnd(vLevel); ++k) {
            residualData(vLevel)(i, 0, k) = (mesh.rankData.rank + 1)*100 + i*10 + k;
        }
    }

    // After going one level down the V-Cycle, populate the pressureData array with values at the corresponding stride
    vLevel += 1;
    pressureData(vLevel) = 0.0;
    for (int i = 0; i <= xEnd(vLevel); ++i) {
        for (int k = 0; k <= zEnd(vLevel); ++k) {
            pressureData(vLevel)(i, 0, k) = (mesh.rankData.rank + 1)*100 + i*10 + k;
        }
    }

    // Perform prolongation
    prolong();

    pressureData(vLevel) -= residualData(vLevel - 1);

    return blitz::max(fabs(pressureData(vLevel)));
}


real multigrid_d2::testTransfer() {
    real maxVal = 0.0;

    vLevel = 0;

    pressureData(vLevel) = 0.0;
    residualData(vLevel) = 0.0;

    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i <= xEnd(vLevel); ++i) {
        for (int k = 0; k <= zEnd(vLevel); ++k) {
            pressureData(vLevel)(i, 0, k) = (mesh.rankData.rank + 1)*100 + i*10 + k;
            residualData(vLevel)(i, 0, k) = pressureData(vLevel)(i, 0, k);
        }
    }

    // EXPECTED VALUES IN THE PAD REGIONS IF DATA TRANSFER HAPPENS WITH NO HITCH
    for (int n = 0; n <= inputParams.vcDepth; n++) {
        for (int k = 0; k <= zEnd(n); ++k) {
            residualData(n)(-1, 0, k) = (mesh.rankData.nearRanks(0) + 1)*100 + (xEnd(n) - 1)*10 + k;
            residualData(n)(xEnd(n) + 1, 0, k) = (mesh.rankData.nearRanks(1) + 1)*100 + 10 + k;
        }
    }

    for (int n=0; n<=inputParams.vcDepth; n++) {
        updatePads(pressureData);
        vLevel += 1;
    }

    pressureData(vLevel) -= residualData(vLevel);

    for (int i = pressureData(vLevel).lbound(0); i <= pressureData(vLevel).ubound(0); i += 1) {
        for (int k = 0; k <= zEnd(vLevel); k += 1) {
            if (abs(pressureData(vLevel)(i, 0, k)) > maxVal) {
                maxVal = abs(pressureData(vLevel)(i, 0, k));
            }
        }
    }

    return maxVal;
}


real multigrid_d2::testPeriodic() {
    real xCoord = 0.0;
    real zCoord = 0.0;

    pressureData(0) = 0.0;
    residualData(0) = 0.0;

    for (int i = 0; i <= xEnd(0); ++i) {
        for (int k = 0; k <= zEnd(0); ++k) {
            pressureData(0)(i, 0, k) = sin(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                       cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
            residualData(0)(i, 0, k) = pressureData(0)(i, 0, k);
        }
    }

    // EXPECTED VALUES IN THE PAD REGIONS IF DATA TRANSFER HAPPENS WITH NO HITCH
    for (int n = 0; n <= inputParams.vcDepth; n++) {
        for (int k = 0; k <= zEnd(n); ++k) {
            xCoord = mesh.xStaggr(0) - (mesh.xStaggr(strideValues(n)) - mesh.xStaggr(0));
            residualData(n)(-1, 0, k) = sin(2.0*M_PI*xCoord/mesh.xLen)*
                                        cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);

            xCoord = mesh.xStaggr(xEnd(0)) + (mesh.xStaggr(xEnd(0)) - mesh.xStaggr(xEnd(0) - strideValues(n)));
            residualData(n)(xEnd(n) + 1, 0, k) = sin(2.0*M_PI*xCoord/mesh.xLen)*
                                                 cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
        }
    }

    for (int n = 0; n <= inputParams.vcDepth; n++) {
        for (int i = 0; i <= xEnd(n); ++i) {
            zCoord = mesh.zStaggr(0) - (mesh.zStaggr(strideValues(n)) - mesh.zStaggr(0));
            residualData(n)(i, 0, -1) = sin(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                        cos(2.0*M_PI*zCoord/mesh.zLen);

            zCoord = mesh.zStaggr(zEnd(0)) + (mesh.zStaggr(zEnd(0)) - mesh.zStaggr(zEnd(0) - strideValues(n)));
            residualData(n)(i, 0, zEnd(n) + 1) = sin(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                                                 cos(2.0*M_PI*zCoord/mesh.zLen);
        }
    }

    vLevel = 0;
    for (int n=0; n<=inputParams.vcDepth; n++) {
        imposeBC();
        vLevel += 1;
    }

    pressureData(vLevel) -= residualData(vLevel);

    return blitz::max(fabs(pressureData(vLevel)));
}
