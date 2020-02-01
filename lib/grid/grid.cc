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
/*! \file grid.cc
 *
 *  \brief Definitions for functions of class grid
 *  \sa grid.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "grid.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the grid class
 *
 *          The class constructor initializes the mesh for computational problem.
 *          The pad widths, global grid limits in the full domain, local grid limits in the MPI decomposed sub-domains,
 *          grid spacings, domain lengths, etc., along each direction are set.
 *          Appropriate stretching functions are chosen according to user preferences and their corresponding grid
 *          transformation derivatives are also computed and stored.
 *
 * \param   solParam is a const reference to the global data contained in the parser class
 * \param   parallelData is a reference to the global data contained in the parallel class
 ********************************************************************************************************************************************
 */
grid::grid(const parser &solParam, parallel &parallelData): inputParams(solParam),
                                                            rankData(parallelData) {
    /** Depending on the finite-difference scheme chosen for calculating derivatives, set the \ref padWidths along all directions. */
    if (inputParams.dScheme == "CD2") {
        padWidths = 1, 1, 1;
    } else {
        if (rankData.rank == 0) {
            std::cout << "Undefined finite differencing scheme in YAML file. ABORTING" << std::endl;
        }
        MPI_Finalize();
        exit(0);
    }

    // THE ARRAY sizeArray HAS ELEMENTS [1, 3, 5, 9, 17, 33 ..... ] - STAGGERED GRID SIZE
    makeSizeArray();

    sizeIndex = inputParams.xInd, inputParams.yInd, inputParams.zInd;
    globalSize = sizeArray(sizeIndex(0)), sizeArray(sizeIndex(1)), sizeArray(sizeIndex(2));

    xLen = inputParams.Lx;
    yLen = inputParams.Ly;
    zLen = inputParams.Lz;

    thBeta = inputParams.betaX, inputParams.betaY, inputParams.betaZ;

    dXi = 1.0/real(globalSize(0) - 1);
    dEt = 1.0/real(globalSize(1) - 1);
    dZt = 1.0/real(globalSize(2) - 1);

#ifdef PLANAR
    padWidths(1) = 1;
    yLen = 1.0;
    dEt = 1.0;
#endif

    // COMPUTE THE LOCAL ARRAY SIZES, collocCoreSize, START AND END INDICES, subarrayStarts AND subarrayEnds
    computeGlobalLimits();

    // SET THE TinyVector AND RectDomain VARIABLES BASED ON VALUES COMPUTED IN computeGlobalLimits, FOR RESIZING ALL LOCAL GRIDS
    setDomainSizes();

    // RESIZE GRID USING THE VARIABLES CONSTRUCTED ABOVE IN setDomainSizes
    resizeGrid();

    // GENERATE THE GLOBAL TRANSFORMED GRID
    globalXiEtaZeta();

    // SET LOCAL TRANSFORMED GRID AS SLICES FROM THE GLOBAL TRANSFORMED GRID GENERATED ABOVE IN globalXiEtaZeta
    xi = xiGlo(blitz::Range(subarrayStarts(0) - padWidths(0), subarrayEnds(0) + padWidths(0)));
    et = etGlo(blitz::Range(subarrayStarts(1) - padWidths(1), subarrayEnds(1) + padWidths(1)));
    zt = ztGlo(blitz::Range(subarrayStarts(2) - padWidths(2), subarrayEnds(2) + padWidths(2)));

    // CREATE UNIFORM GRID WHICH IS DEFAULT ALONG ALL THREE DIRECTIONS
    createUniformGrid();

    // FLAG TO CHECK FOR GRID ANISOTROPY - FALSE BY DEFAULT UNLESS NON-UNIFORM GRID IS CREATED
    bool gridCheck = false;

    // DEPENDING ON THE USER-SET PARAMETERS, SWITCH TO TAN-HYP ALONG SELECTED DIRECTIONS
    if (inputParams.xGrid == 2) {
        createTanHypGrid(0);
        gridCheck = true;
    }

    if (inputParams.yGrid == 2) {
        createTanHypGrid(1);
        gridCheck = true;
    }

    if (inputParams.zGrid == 2) {
        createTanHypGrid(2);
        gridCheck = true;
    }

    if (gridCheck) checkAnisotropy();

    gatherGlobal();
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to resize and initialize the size array from which the dimensions of the grid will be determined
 *
 *          The size array will generate grid sizes according to \f$ 2^N + 2 \f$ to enable multigrid operations on the grid.
 *          The \ref parser#xInd "xInd", \ref parser#yInd "yInd" and \ref parser#zInd "zInd" parameters set by the users in
 *          parameters.yaml and read by the \ref parser class will be used to locate the grid size within the \ref sizeArray and
 *          generate grid accordingly.
 *
 *          Note that the grid sizes stored in \ref sizeArray correspond to the collocated grid.
 *          For multi-grid operations, the number of grid points necessary is \f$ 2^N + 1 \f$, which is 1 less than the values
 *          generated by this function.
 *          However, multi-grid is applied here to compute pressure correction and pressure is calculated on the staggered grid.
 *          Since there are \f$ N - 1 \f$ staggered grid points for \f$ N \f$ collocated points, the sizes become consistent.
 ********************************************************************************************************************************************
 */
void grid::makeSizeArray() {
    int maxIndex = 15;

    sizeArray.resize(maxIndex);
    for (int i=0; i < maxIndex; i++) {
        sizeArray(i) = int(pow(2, i)) + 1;
    }

    sizeArray(0) = 1;
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to set the extent of local sub-domains in terms of the global index of the full domain
 *
 *          Depending on the number of processor divisions along each direction, the limits of the grid for each local
 *          sub-domain is set based on its \ref parallel#xRank "xRank" and \ref parallel#yRank "yRank".
 *          These limits are used to locate the local sub-domains within the full domain later.
 ********************************************************************************************************************************************
 */
void grid::computeGlobalLimits() {
    int xiSt, etSt, ztSt;
    int xiEn, etEn, ztEn;
    int localNx, localNy, localNz;

    // NUMBER OF STAGGERED POINTS IN EACH SUB-DOMAIN EXCLUDING PAD POINTS
    localNx = (globalSize(0) - 1)/rankData.npX + 1;
#ifndef PLANAR
    localNy = (globalSize(1) - 1)/rankData.npY + 1;
#else
    localNy = 1;
#endif
    localNz = (globalSize(2));

    // SETTING GLOBAL LIMITS
    // ADD ONE EXTRA POINT EACH AT FIRST AND LAST SUB-DOMAINS
    // FIRST SET THE LIMITS TO DEFAULT VALUES - THIS ELIMINATES AN EXTRA 'if' CONDITION
    // THEN SET LIMITS FOR LAST RANK IN EACH DIRECTION FIRST AND *FINALLY* SET LIMITS OF 0TH RANK
    // THIS IS NECESSARY TO AVOID ERRORS WHEN A PROCESSOR IS BOTH FIRST AND LAST RANK
    // THIS HAPPENS WHEN THERE ARE NO DIVISIONS ALONG AN AXIS AS ALONG Z-DIRECTION

    // ALONG XI-DIRECTION
    xiSt = rankData.xRank*(localNx - 1);
    xiEn = xiSt + localNx - 1;

    // ALONG ETA-DIRECTION
    etSt = rankData.yRank*(localNy - 1);
    etEn = etSt + localNy - 1;

    // ALONG ZETA-DIRECTION
    ztSt = 0;
    ztEn = ztSt + localNz - 1;

    staggrCoreSize = localNx, localNy, localNz;
    staggrFullSize = staggrCoreSize + 2*padWidths;

    // SUB-ARRAY STARTS AND ENDS FOR *STAGGERED* GRID
    subarrayStarts = xiSt, etSt, ztSt;
    subarrayEnds = xiEn, etEn, ztEn;
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to set all the TinyVector and RectDomain variables for all future references throughout the solver
 *
 *          The function sets the core and full domain sizes for all the sub-domains after MPI decomposition.
 *          Additionally, the pad widths and starting indices of the sub-domains within the global domain are also set.
 ********************************************************************************************************************************************
 */
void grid::setDomainSizes() {
    blitz::TinyVector<int, 3> loBound, upBound;

    // SIZE OF THE STAGGERED DOMAIN IN THE CORE OF THE LOCAL SUB-DOMAIN
    collocCoreSize = staggrCoreSize - 1;

#ifdef PLANAR
    collocCoreSize(1) += 1;
#endif

    collocFullSize = collocCoreSize + 2*padWidths;

    // LOWER BOUND AND UPPER BOUND OF CORE - USED TO CONSTRUCT THE CORE SLICE OF COLLOCATED POINTS
    loBound = 0, 0, 0;
    upBound = collocCoreSize - 1;
    collocCoreDomain = blitz::RectDomain<3>(loBound, upBound);

    // LOWER BOUND AND UPPER BOUND OF CORE - USED TO CONSTRUCT THE CORE SLICE OF STAGGERED POINTS
    upBound = staggrCoreSize - 1;
    staggrCoreDomain = blitz::RectDomain<3>(loBound, upBound);

    // LOWER BOUND AND UPPER BOUND OF FULL SUB-DOMAIN - USED TO CONSTRUCT THE FULL SUB-DOMAIN SLICE
    loBound = -padWidths;
    upBound = collocCoreSize + padWidths - 1;
    collocFullDomain = blitz::RectDomain<3>(loBound, upBound);

    // LOWER BOUND AND UPPER BOUND OF FULL SUB-DOMAIN - USED TO CONSTRUCT THE FULL SUB-DOMAIN SLICE
    upBound = staggrCoreSize + padWidths - 1;
    staggrFullDomain = blitz::RectDomain<3>(loBound, upBound);
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to resize and initialize the grid
 *
 *          The global collocated grid in transformed plane are resized according the global size of full domain.
 *          Then the local collocated grid in transformed plane is resized according to the limits defined in \ref computeGlobalLimits.
 *          Correspondingly, this function is called after the global limits have been set.
 *          After defining the transformed plane coordinates, both the staggered and collocated grids in physical plane are resized.
 *          Finally, the arrays for the grid derivative terms are also resized and initialized to 1.
 ********************************************************************************************************************************************
 */
void grid::resizeGrid() {
    // ALL ARRAYS MUST BE RESIZED AND REINDEXED TO LET THE NEGATIVE PADS HAVE NEGATIVE INDICES
    // THIS IS DONE IN A SINGLE STEP BY INITIALIZING THE ARRAYS WITH A blitz::Range OBJECT WHICH CONTAINS SIZE AND INDEXING INFORMATION
    blitz::Range xStagRange, xCollRange;
    blitz::Range yStagRange, yCollRange;
    blitz::Range zStagRange, zCollRange;

    // RANGE OF THE SUB-DOMAIN FOR STAGGERED AND COLLOCATED GRIDS: CONSTRUCTED FROM LOWER AND UPPER BOUNDS OF FULL SUB-DOMAIN
    xCollRange = blitz::Range(collocFullDomain.lbound(0), collocFullDomain.ubound(0));
    yCollRange = blitz::Range(collocFullDomain.lbound(1), collocFullDomain.ubound(1));
    zCollRange = blitz::Range(collocFullDomain.lbound(2), collocFullDomain.ubound(2));

    xStagRange = blitz::Range(staggrFullDomain.lbound(0), staggrFullDomain.ubound(0));
    yStagRange = blitz::Range(staggrFullDomain.lbound(1), staggrFullDomain.ubound(1));
    zStagRange = blitz::Range(staggrFullDomain.lbound(2), staggrFullDomain.ubound(2));

    // LOCAL XI, ETA AND ZETA ARRAYS
    xi.resize(xStagRange);
    et.resize(yStagRange);
    zt.resize(zStagRange);

    // COLLOCATED GRID POINTS AND THEIR METRICS
    xColloc.resize(xCollRange);
    yColloc.resize(yCollRange);
    zColloc.resize(zCollRange);

    xi_xColloc.resize(xCollRange);          xixxColloc.resize(xCollRange);          xix2Colloc.resize(xCollRange);
    et_yColloc.resize(yCollRange);          etyyColloc.resize(yCollRange);          ety2Colloc.resize(yCollRange);
    zt_zColloc.resize(zCollRange);          ztzzColloc.resize(zCollRange);          ztz2Colloc.resize(zCollRange);

    // STAGGERED GRID POINTS AND THEIR METRICS
    xStaggr.resize(xStagRange);
    yStaggr.resize(yStagRange);
    zStaggr.resize(zStagRange);

    xi_xStaggr.resize(xStagRange);          xixxStaggr.resize(xStagRange);          xix2Staggr.resize(xStagRange);
    et_yStaggr.resize(yStagRange);          etyyStaggr.resize(yStagRange);          ety2Staggr.resize(yStagRange);
    zt_zStaggr.resize(zStagRange);          ztzzStaggr.resize(zStagRange);          ztz2Staggr.resize(zStagRange);

    xColloc = 1.0;          xStaggr = 1.0;
    yColloc = 1.0;          yStaggr = 1.0;
    zColloc = 1.0;          zStaggr = 1.0;

    // BELOW ARE DEFAULT VALUES FOR A UNIFORM GRID OVER DOMAIN OF LENGTH 1.0
    // THESE VALUES ARE OVERWRITTEN AS PER GRID TYPE
    xi_xColloc = 1.0;          xixxColloc = 0.0;          xix2Colloc = 1.0;
    et_yColloc = 1.0;          etyyColloc = 0.0;          ety2Colloc = 1.0;
    zt_zColloc = 1.0;          ztzzColloc = 0.0;          ztz2Colloc = 1.0;

    xi_xStaggr = 1.0;          xixxStaggr = 0.0;          xix2Staggr = 1.0;
    et_yStaggr = 1.0;          etyyStaggr = 0.0;          ety2Staggr = 1.0;
    zt_zStaggr = 1.0;          ztzzStaggr = 0.0;          ztz2Staggr = 1.0;
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to compute global values of xi, eta and zeta in transformed plane
 *
 *          The function populates the \ref xiGlo, \ref etGlo and \ref ztGlo arrays from which the local values of
 *          \ref xi, \ref et and \ref zt in each sub-domain are obtained.
 *          These local values are obtained from the global grid according to the limits defined in \ref computeGlobalLimits.
 ********************************************************************************************************************************************
 */
void grid::globalXiEtaZeta() {
    xiGlo.resize(globalSize(0) + 2*padWidths(0));          xiGlo.reindexSelf(-padWidths(0));
    etGlo.resize(globalSize(1) + 2*padWidths(1));          etGlo.reindexSelf(-padWidths(1));
    ztGlo.resize(globalSize(2) + 2*padWidths(2));          ztGlo.reindexSelf(-padWidths(2));

    // ALONG XI-DIRECTION
    for (int i=-padWidths(0); i<globalSize(0)+padWidths(0); i++) {
        xiGlo(i) = real(i)*dXi;
    }

    // ALONG ETA-DIRECTION
    for (int i=-padWidths(1); i<globalSize(1)+padWidths(1); i++) {
        etGlo(i) = real(i)*dEt;
    }

    // ALONG ZETA-DIRECTION
    for (int i=-padWidths(2); i<globalSize(2)+padWidths(2); i++) {
        ztGlo(i) = real(i)*dZt;
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to generate grid with uniform stretching
 *
 *          The local collocated grids, \ref xColloc, \ref yColloc, and \ref zColloc are equated to their corresponding
 *          transformed plane coordinates, \ref xi, \ref et, and \ref zt respectively.
 *          The corresponding grid derivative terms, \ref xi_xColloc, \ref xixxColloc, \ref ety2Colloc, etc are left
 *          unchanged from their initial value of 1.0, indicating that the grid is uniform.
 *
 *          Similarly, the staggered grids, \ref xStaggr, \ref yStaggr, and \ref zStaggr are also equated to the mid-point
 *          averaged values of the nodes in their corresponding transformed plane coordinates, \ref xi, \ref et, and \ref zt
 *          respectively.
 *          As before, the grid derivative terms for the staggered points are also left as 1.0.
 ********************************************************************************************************************************************
 */
void grid::createUniformGrid() {
    int i;

    // COLLOCATED X-GRID POINTS FROM UNIFORM XI-GRID POINTS AND THEIR METRICS
    for (i = -padWidths(0); i < staggrCoreSize(0) + padWidths(0); i++) {
        xStaggr(i) = xLen*xi(i);
    }

    // STAGGERED X-GRID POINTS FROM UNIFORM XI-GRID POINTS AND THEIR METRICS
    for (i = -padWidths(0); i < collocCoreSize(0) + padWidths(0); i++) {
        xColloc(i) = xLen*(xi(i) + xi(i + 1))/2.0;
    }

#ifndef PLANAR
    // COLLOCATED Y-GRID POINTS FROM UNIFORM ETA-GRID POINTS AND THEIR METRICS
    for (i = -padWidths(1); i < staggrCoreSize(1) + padWidths(1); i++) {
        yStaggr(i) = yLen*et(i);
    }

    // STAGGERED Y-GRID POINTS FROM UNIFORM ETA-GRID POINTS AND THEIR METRICS
    for (i = -padWidths(1); i < collocCoreSize(1) + padWidths(1); i++) {
        yColloc(i) = yLen*(et(i) + et(i + 1))/2.0;
    }
#endif

    // COLLOCATED Z-GRID POINTS FROM UNIFORM ZETA-GRID POINTS AND THEIR METRICS
    for (i = -padWidths(2); i < staggrCoreSize(2) + padWidths(2); i++) {
        zStaggr(i) = zLen*zt(i);
    }

    // STAGGERED Z-GRID POINTS FROM UNIFORM ZETA-GRID POINTS AND THEIR METRICS
    for (i = -padWidths(2); i < collocCoreSize(2) + padWidths(2); i++) {
        zColloc(i) = zLen*(zt(i) + zt(i + 1))/2.0;
    }

    xi_xStaggr = 1.0/xLen;
    xix2Staggr = pow(xi_xStaggr, 2.0);

    xi_xColloc = 1.0/xLen;
    xix2Colloc = pow(xi_xColloc, 2.0);

    et_yStaggr = 1.0/yLen;
    ety2Staggr = pow(et_yStaggr, 2.0);

    et_yColloc = 1.0/yLen;
    ety2Colloc = pow(et_yColloc, 2.0);

    zt_zStaggr = 1.0/zLen;
    ztz2Staggr = pow(zt_zStaggr, 2.0);

    zt_zColloc = 1.0/zLen;
    ztz2Colloc = pow(zt_zColloc, 2.0);
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to generate grid with tangent-hyperbolic stretching
 *
 *          The local collocated grids, \ref xColloc, \ref yColloc, and \ref zColloc are initialized from their corresponding
 *          transformed plane coordinates, \ref xi, \ref et, and \ref zt respectively using the tangent hyperbolic function.
 *          The corresponding grid derivative terms, \ref xi_xColloc, \ref xixxColloc, \ref ety2Colloc, etc are computed
 *          using analytical expressions for the tangent hyperbolic function.
 *
 *          Similarly, the staggered grids, \ref xStaggr, \ref yStaggr, and \ref zStaggr are initialized from the mid-point
 *          averaged values of the nodes in their corresponding transformed plane coordinates, \ref xi, \ref et, and \ref zt
 *          respectively using the tangent hyperbolic function.
 *          As before, the grid derivative terms for the staggered points are also computed using analytical expressions.
 *
 * \param   dim is an integer value that defines the direction along which tan-hyp grid is to be generated: 0 -> X, 1 -> Y, 2 -> Z
 ********************************************************************************************************************************************
 */
void grid::createTanHypGrid(int dim) {
    int i;

#ifndef TEST_RUN
    if (rankData.rank == 0) {
        switch (dim) {
            case 0: std::cout << "Generating tangent hyperbolic grid along X direction" << std::endl;
                    break;
            case 1: std::cout << "Generating tangent hyperbolic grid along Y direction" << std::endl;
                    break;
            case 2: std::cout << "Generating tangent hyperbolic grid along Z direction" << std::endl;
                    break;
        }
    }
#endif

    if (dim == 0) {
        // STAGGERED X-GRID POINTS FROM UNIFORM XI-GRID POINTS AND THEIR METRICS
        for (i = 0; i < staggrCoreSize(0); i++) {
            xStaggr(i) = xLen*(1.0 - tanh(thBeta[0]*(1.0 - 2.0*xi(i)))/tanh(thBeta[0]))/2.0;

            xi_xStaggr(i) = tanh(thBeta[0])/(thBeta[0]*xLen*(1.0 - pow((1.0 - 2.0*xStaggr(i)/xLen)*tanh(thBeta[0]), 2)));
            xixxStaggr(i) = -4.0*pow(tanh(thBeta[0]), 3)*(1.0 - 2.0*xStaggr(i)/xLen)/(thBeta[0]*xLen*xLen*pow(1.0 - pow(tanh(thBeta[0])*(1.0 - 2.0*xStaggr(i)/xLen), 2), 2));
            xix2Staggr(i) = pow(xi_xStaggr(i), 2.0);
        }

        // COLLOCATED X-GRID POINTS FROM UNIFORM XI-GRID POINTS AND THEIR METRICS
        for (i = 0; i < collocCoreSize(0); i++) {
            xColloc(i) = xLen*(1.0 - tanh(thBeta[0]*(1.0 - (xi(i) + xi(i + 1))))/tanh(thBeta[0]))/2.0;

            xi_xColloc(i) = tanh(thBeta[0])/(thBeta[0]*xLen*(1.0 - pow((1.0 - 2.0*xColloc(i)/xLen)*tanh(thBeta[0]), 2)));
            xixxColloc(i) = -4.0*pow(tanh(thBeta[0]), 3)*(1.0 - 2.0*xColloc(i)/xLen)/(thBeta[0]*xLen*xLen*pow(1.0 - pow(tanh(thBeta[0])*(1.0 - 2.0*xColloc(i)/xLen), 2), 2));
            xix2Colloc(i) = pow(xi_xColloc(i), 2.0);
        }
    }

#ifndef PLANAR
    if (dim == 1) {
        // STAGGERED Y-GRID POINTS FROM UNIFORM ETA-GRID POINTS AND THEIR METRICS
        for (i = 0; i < staggrCoreSize(1); i++) {
            yStaggr(i) = yLen*(1.0 - tanh(thBeta[1]*(1.0 - 2.0*et(i)))/tanh(thBeta[1]))/2.0;

            et_yStaggr(i) = tanh(thBeta[1])/(thBeta[1]*yLen*(1.0 - pow((1.0 - 2.0*yStaggr(i)/yLen)*tanh(thBeta[1]), 2)));
            etyyStaggr(i) = -4.0*pow(tanh(thBeta[1]), 3)*(1.0 - 2.0*yStaggr(i)/yLen)/(thBeta[1]*yLen*yLen*pow(1.0 - pow(tanh(thBeta[1])*(1.0 - 2.0*yStaggr(i)/yLen), 2), 2));
            ety2Staggr(i) = pow(et_yStaggr(i), 2.0);
        }

        // COLLOCATED Y-GRID POINTS FROM UNIFORM ETA-GRID POINTS AND THEIR METRICS
        for (i = 0; i < collocCoreSize(1); i++) {
            yColloc(i) = yLen*(1.0 - tanh(thBeta[1]*(1.0 - (et(i) + et(i + 1))))/tanh(thBeta[1]))/2.0;

            et_yColloc(i) = tanh(thBeta[1])/(thBeta[1]*yLen*(1.0 - pow((1.0 - 2.0*yColloc(i)/yLen)*tanh(thBeta[1]), 2)));
            etyyColloc(i) = -4.0*pow(tanh(thBeta[1]), 3)*(1.0 - 2.0*yColloc(i)/yLen)/(thBeta[1]*yLen*yLen*pow(1.0 - pow(tanh(thBeta[1])*(1.0 - 2.0*yColloc(i)/yLen), 2), 2));
            ety2Colloc(i) = pow(et_yColloc(i), 2.0);
        }
    }
#endif

    if (dim == 2) {
        // STAGGERED Z-GRID POINTS FROM UNIFORM ZETA-GRID POINTS AND THEIR METRICS
        for (i = 0; i < staggrCoreSize(2); i++) {
            zStaggr(i) = zLen*(1.0 - tanh(thBeta[2]*(1.0 - 2.0*zt(i)))/tanh(thBeta[2]))/2.0;

            zt_zStaggr(i) = tanh(thBeta[2])/(thBeta[2]*zLen*(1.0 - pow((1.0 - 2.0*zStaggr(i)/zLen)*tanh(thBeta[2]), 2)));
            ztzzStaggr(i) = -4.0*pow(tanh(thBeta[2]), 3)*(1.0 - 2.0*zStaggr(i)/zLen)/(thBeta[2]*zLen*zLen*pow(1.0 - pow(tanh(thBeta[2])*(1.0 - 2.0*zStaggr(i)/zLen), 2), 2));
            ztz2Staggr(i) = pow(zt_zStaggr(i), 2.0);
        }

        // COLLOCATED Z-GRID POINTS FROM UNIFORM ZETA-GRID POINTS AND THEIR METRICS
        for (i = 0; i < collocCoreSize(2); i++) {
            zColloc(i) = zLen*(1.0 - tanh(thBeta[2]*(1.0 - (zt(i) + zt(i + 1))))/tanh(thBeta[2]))/2.0;

            zt_zColloc(i) = tanh(thBeta[2])/(thBeta[2]*zLen*(1.0 - pow((1.0 - 2.0*zColloc(i)/zLen)*tanh(thBeta[2]), 2)));
            ztzzColloc(i) = -4.0*pow(tanh(thBeta[2]), 3)*(1.0 - 2.0*zColloc(i)/zLen)/(thBeta[2]*zLen*zLen*pow(1.0 - pow(tanh(thBeta[2])*(1.0 - 2.0*zColloc(i)/zLen), 2), 2));
            ztz2Colloc(i) = pow(zt_zColloc(i), 2.0);
        }
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to check the anisotropy of the grid
 *
 *          The function is called only if non-uniform grid is made.
 *          It scans through the entire grid cell-by-cell and checks the 2 aspect-ratios of the cell (one for 2D grids)
 *          The maximum value of aspect ratio is stored.
 *          Each MPI sub-domain checks within its limits and an MPI_Reduce call gets the global maximum.
 *
 ********************************************************************************************************************************************
 */
void grid::checkAnisotropy() {
    real xWidth, yWidth, zWidth;
    real localMax, globalMax;
    real xyRatio, yzRatio;
    real cellMaxAR;

    localMax = 0.0;
    for (int i = 0; i < collocCoreSize.ubound(0); i++) {
        for (int j = 0; j < collocCoreSize.ubound(0); j++) {
            for (int k = 0; k < collocCoreSize.ubound(0); k++) {
                xWidth = xColloc(i-1) - xColloc(i);
                yWidth = yColloc(j-1) - yColloc(j);
                zWidth = zColloc(k-1) - zColloc(k);
                xyRatio = std::max(xWidth/yWidth, yWidth/xWidth);
                yzRatio = std::max(yWidth/zWidth, zWidth/yWidth);
                cellMaxAR = std::max(xyRatio, yzRatio);
                if (cellMaxAR > localMax) localMax = cellMaxAR;
            }
        }
    }

    MPI_Allreduce(&localMax, &globalMax, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    if (globalMax > 5.0) {
        if (rankData.rank == 0) std::cout << "\nWARNING: Grid anisotropy exceeds limits. Finite-difference calculations will be inaccurate" << std::endl;
    } else {
        if (rankData.rank == 0) std::cout << "\nMaximum grid anisotropy is " << globalMax << std::endl;
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to gather global data about the grid into every rank
 *
 *          In certain cases, for example, implementation of non-homogeneous boundary conditions,
 *          the data regarding the global extents of the grid and the global coordinates are necessary.
 *          For this reason, the global grid data is made available to each rank through an MPI_Allgather call.
 *
 ********************************************************************************************************************************************
 */
void grid::gatherGlobal() {
    int i;
    int locSize, locDisp;
    int maxRank = std::max(rankData.npX, rankData.npY);

    int arrSize[maxRank];
    int arrDisp[maxRank];

    blitz::TinyVector<int, 3> collocGlobalSize;
    blitz::TinyVector<int, 3> staggrGlobalSize;
    blitz::TinyVector<int, 3> globalReIndexVal;

    staggrGlobalSize = globalSize + 2*padWidths;
    collocGlobalSize = globalSize + 2*padWidths - 1;
    globalReIndexVal = -padWidths;

    xCollocGlobal.resize(collocGlobalSize(0));     xCollocGlobal.reindexSelf(globalReIndexVal(0));        xCollocGlobal = 0.0;
    xStaggrGlobal.resize(staggrGlobalSize(0));     xStaggrGlobal.reindexSelf(globalReIndexVal(0));        xStaggrGlobal = 0.0;

#ifndef PLANAR
    yCollocGlobal.resize(collocGlobalSize(1));     yCollocGlobal.reindexSelf(globalReIndexVal(1));        yCollocGlobal = 0.0;
    yStaggrGlobal.resize(staggrGlobalSize(1));     yStaggrGlobal.reindexSelf(globalReIndexVal(1));        yStaggrGlobal = 0.0;
#endif

    zCollocGlobal.resize(collocGlobalSize(2));     zCollocGlobal.reindexSelf(globalReIndexVal(2));        zCollocGlobal = 0.0;
    zStaggrGlobal.resize(staggrGlobalSize(2));     zStaggrGlobal.reindexSelf(globalReIndexVal(2));        zStaggrGlobal = 0.0;

    // GATHERING THE STAGGERED GRID ALONG X-DIRECTION
    locSize = xStaggr.size() - 2*padWidths(0);
    locDisp = subarrayStarts(0);
    if (rankData.xRank == rankData.npX-1) {
        locSize += 2*padWidths(0);
    }
    MPI_Allgather(&locSize, 1, MPI_INT, arrSize, 1, MPI_INT, rankData.MPI_ROW_COMM);
    MPI_Allgather(&locDisp, 1, MPI_INT, arrDisp, 1, MPI_INT, rankData.MPI_ROW_COMM);
    MPI_Allgatherv(xStaggr.dataFirst(), locSize, MPI_FP_REAL, xStaggrGlobal.dataFirst(), arrSize, arrDisp, MPI_FP_REAL, rankData.MPI_ROW_COMM);

    // GATHERING THE COLLOCATED GRID ALONG X-DIRECTION
    locSize = xColloc.size() - 2*padWidths(0);
    locDisp = rankData.xRank*locSize;
    if (rankData.xRank == rankData.npX-1) {
        locSize += 2*padWidths(0);
    }
    MPI_Allgather(&locSize, 1, MPI_INT, arrSize, 1, MPI_INT, rankData.MPI_ROW_COMM);
    MPI_Allgather(&locDisp, 1, MPI_INT, arrDisp, 1, MPI_INT, rankData.MPI_ROW_COMM);
    MPI_Allgatherv(xColloc.dataFirst(), locSize, MPI_FP_REAL, xCollocGlobal.dataFirst(), arrSize, arrDisp, MPI_FP_REAL, rankData.MPI_ROW_COMM);

#ifndef PLANAR
    // GATHERING THE STAGGERED GRID ALONG Y-DIRECTION
    locSize = yStaggr.size() - 2*padWidths(1);
    locDisp = subarrayStarts(1);
    if (rankData.yRank == rankData.npY-1) {
        locSize += 2*padWidths(1);
    }
    MPI_Allgather(&locSize, 1, MPI_INT, arrSize, 1, MPI_INT, rankData.MPI_COL_COMM);
    MPI_Allgather(&locDisp, 1, MPI_INT, arrDisp, 1, MPI_INT, rankData.MPI_COL_COMM);
    MPI_Allgatherv(yStaggr.dataFirst(), locSize, MPI_FP_REAL, yStaggrGlobal.dataFirst(), arrSize, arrDisp, MPI_FP_REAL, rankData.MPI_COL_COMM);

    // GATHERING THE COLLOCATED GRID ALONG Y-DIRECTION
    locSize = yColloc.size() - 2*padWidths(1);
    locDisp = rankData.yRank*locSize;
    if (rankData.yRank == rankData.npY-1) {
        locSize += 2*padWidths(1);
    }
    MPI_Allgather(&locSize, 1, MPI_INT, arrSize, 1, MPI_INT, rankData.MPI_COL_COMM);
    MPI_Allgather(&locDisp, 1, MPI_INT, arrDisp, 1, MPI_INT, rankData.MPI_COL_COMM);
    MPI_Allgatherv(yColloc.dataFirst(), locSize, MPI_FP_REAL, yCollocGlobal.dataFirst(), arrSize, arrDisp, MPI_FP_REAL, rankData.MPI_COL_COMM);
#endif

    // GLOBAL AND LOCAL STAGGERED GRIDS ALONG Z-DIRECTION ARE SAME FOR ALL RANKS
    for (i = -padWidths(2); i < globalSize(2) + padWidths(2); i++) {
        zStaggrGlobal(i) = zStaggr(i);
    }

    // GLOBAL AND LOCAL COLLOCATED GRIDS ALONG Z-DIRECTION ARE SAME FOR ALL RANKS
    for (i = -padWidths(2); i < globalSize(2) + padWidths(2) - 1; i++) {
        zCollocGlobal(i) = zColloc(i);
    }
}
