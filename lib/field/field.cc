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
/*! \file field.cc
 *
 *  \brief Definitions for functions of class field
 *  \sa field.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "field.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the field class
 *
 *          The field class decides the limits necessary for a 3D array to store the data as per the specified grid staggering details.
 *          It initializes and stores necessary RectDomain objects for getting the core slice and various offset slices for performing
 *          finite difference operations.
 *          The upper and lower bounds necessary for the array are also calculated depending on the directions along which the mesh is
 *          staggered and those along which it is collocated.
 *          Finally, a blitz array to store the data of the field is resized according to the limits and initialized to 0.
 *
 * \param   gridData is a const reference to the global data contained in the grid class
 * \param   xStag is a const boolean value that is <B>true</B> when the grid is staggered along the x-direction and <B>false</B> when it is not
 * \param   yStag is a const boolean value that is <B>true</B> when the grid is staggered along the y-direction and <B>false</B> when it is not
 * \param   zStag is a const boolean value that is <B>true</B> when the grid is staggered along the z-direction and <B>false</B> when it is not
 ********************************************************************************************************************************************
 */
field::field(const grid &gridData, std::string fieldName, const bool xStag, const bool yStag, const bool zStag):
             gridData(gridData),
             xStag(xStag), yStag(yStag), zStag(zStag)
{
    this->fieldName = fieldName;

    fSize = gridData.collocFullSize;
    flBound = gridData.collocFullDomain.lbound();

    if (xStag) {
        fSize(0) = gridData.staggrFullSize(0);
        flBound(0) = gridData.staggrFullDomain.lbound()(0);
    }    

    if (yStag) {
        fSize(1) = gridData.staggrFullSize(1);
        flBound(1) = gridData.staggrFullDomain.lbound()(1);
     }

    if (zStag) {
        fSize(2) = gridData.staggrFullSize(2);
        flBound(2) = gridData.staggrFullDomain.lbound()(2);
    }

    F.resize(fSize);
    F.reindexSelf(flBound);

    mpiHandle = new mpidata(F, gridData.rankData);

    setCoreSlice();
    setBulkSlice();

    setWallSlices();

    setInterpolationSlices();

    mpiHandle->createSubarrays(fSize, cuBound + 1, gridData.padWidths, xStag, yStag);

    F = 0.0;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to shift a blitz RectDomain object a specified number of steps in specified dimension.
 *
 *          The RectDomain objects offer a view of the blitz arrays on which the shift function operates.
 *          These objects are shifted along the dimension specified in the argument, by <B>dim</B>, through a number of steps,
 *          to offer offset views.
 *
 * \param   dim is the integer input to specify the dimension (direction) of the shift. (x -> 0, y -> 1, z -> 2)
 * \param   core is the input RectDomain object which is to be shifted to get the new view
 * \param   steps is the integer value by which the input view must be offset along the dimension specified by <B>dim</B>
 *
 * \return  A RectDomain object that specifies the new offset view of the data
 ********************************************************************************************************************************************
 */
blitz::RectDomain<3> field::shift(int dim, blitz::RectDomain<3> core, int steps) {
    core.lbound()(dim) += steps;
    core.ubound()(dim) += steps;

    return core;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to create the core slice and its offset views
 *
 *          The core and full slices of the field differentiates the domain into the computational sub-domain and
 *          the ghost point/pad point regions which play only an auxilliary role in computing the derivatives.
 *          The core slice defined here also includes the walls of the domain.
 *          The core slice is also offset in all the directions for ease in computation of numerical derivatives.
 *
 ********************************************************************************************************************************************
 */
void field::setCoreSlice() {
    cuBound = gridData.collocCoreDomain.ubound();

    if (xStag) {
        cuBound(0) = gridData.staggrCoreDomain.ubound()(0);
    }

    if (yStag) {
        cuBound(1) = gridData.staggrCoreDomain.ubound()(1);
    }

    if (zStag) {
        cuBound(2) = gridData.staggrCoreDomain.ubound()(2);
    }

    // Following lines taken from Aether to correct periodic BCs for channel flow - test it thoroughly
    // They need to be commented when using Method 3 in setBulkSlice function below
    // Pushing the last point at the end of the domain inside by one unit of grid spacing for periodic domains
    //if (xStag and gridData.rankData.xRank == gridData.rankData.npX - 1 and gridData.inputParams.xPer) cuBound(0) -= 1;
    //if (yStag and gridData.rankData.yRank == gridData.rankData.npY - 1 and gridData.inputParams.yPer) cuBound(1) -= 1;
    //if (zStag and gridData.inputParams.zPer) cuBound(2) -= 1;

    fCore = blitz::RectDomain<3>(blitz::TinyVector<int, 3>(0, 0, 0), cuBound);

    // As of Dec 2019, the below slices are used only in the divergence calculation of vfield and gradient calculation of sfield and plainsf
    fCLft = shift(0, fCore, -1);
    fCRgt = shift(0, fCore,  1);

    fCFrt = shift(1, fCore, -1);
    fCBak = shift(1, fCore,  1);

    fCBot = shift(2, fCore, -1);
    fCTop = shift(2, fCore,  1);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to create the bulk slice and its offset views
 *
 *          The bulk and wall slices of the field differentiates the domain into the bulk of the fluid and the walls of the domain.
 *          The bulk slice is used in two places - i) to set wall slices - walls are considered to be the points just outside bulk,
 *                                                ii) the limits of all iterative solvers in src are set to bulk limits.
 *
 ********************************************************************************************************************************************
 */
void field::setBulkSlice() {
    blitz::TinyVector<int, 3> blBound;
    blitz::TinyVector<int, 3> buBound;

    blBound = gridData.collocCoreDomain.lbound();
    buBound = gridData.collocCoreDomain.ubound();

    if (xStag) {
        blBound(0) = gridData.staggrCoreDomain.lbound()(0);
        buBound(0) = gridData.staggrCoreDomain.ubound()(0);
    }

    if (yStag) {
        blBound(1) = gridData.staggrCoreDomain.lbound()(1);
        buBound(1) = gridData.staggrCoreDomain.ubound()(1);
    }

    if (zStag) {
        blBound(2) = gridData.staggrCoreDomain.lbound()(2);
        buBound(2) = gridData.staggrCoreDomain.ubound()(2);
    }

    // Bulk and core slices are differentiated only in the boundary sub-domains,
    // and that differentiation is imposed in the following lines
    // At all interior sub-domains after performing MPI domain decomposition,
    // the bulk and core slices are identical

    // Different ways of defining bulk are clubbed together
    // The correct method needs to be chosen through Swayamvar

    // Method 1: The method originally present in Saras
    /*
    if (xStag and gridData.rankData.xRank == 0) blBound(0) += 1;

    if (xStag and gridData.rankData.xRank == gridData.rankData.npX - 1) buBound(0) -= 1;

    if (yStag and gridData.rankData.yRank == 0) blBound(1) += 1;

    if (yStag and gridData.rankData.yRank == gridData.rankData.npY - 1) buBound(1) -= 1;

    if (zStag) {
        blBound(2) += 1;
        buBound(2) -= 1;
    }
    */

    // Method 2: The method implemented after long discussions - shift the entire bulk to one side
    /*
    if (xStag and gridData.rankData.xRank == 0 and not gridData.inputParams.xPer) blBound(0) += 1;

    if (xStag and gridData.rankData.xRank == gridData.rankData.npX - 1) {
        gridData.inputParams.xPer? buBound(0) -= 2: buBound(0) -= 1;
    }

    if (yStag and gridData.rankData.yRank == 0 and not gridData.inputParams.yPer) blBound(1) += 1;

    if (yStag and gridData.rankData.yRank == gridData.rankData.npY - 1) {
        gridData.inputParams.yPer? buBound(1) -= 2: buBound(1) -= 1;
    }

    if (zStag) {
        if (not gridData.inputParams.zPer) {
            blBound(2) += 1;
            buBound(2) -= 1;
        } else {
            buBound(2) -= 2;
        }
    }
    */

    // Method 3: Seemingly the oldest version that existed in Aether
    if (xStag and gridData.rankData.xRank == 0 and not gridData.inputParams.xPer) blBound(0) += 1;

    if (xStag and gridData.rankData.xRank == gridData.rankData.npX - 1 and not gridData.inputParams.xPer) buBound(0) -= 1;

    if (yStag and gridData.rankData.yRank == 0 and not gridData.inputParams.yPer) blBound(1) += 1;

    if (yStag and gridData.rankData.yRank == gridData.rankData.npY - 1 and not gridData.inputParams.yPer) buBound(1) -= 1;

    if (zStag) {
        if (not gridData.inputParams.zPer) {
            blBound(2) += 1;
            buBound(2) -= 1;
        }
    }

#ifdef PLANAR
    blBound(1) = 0;
    buBound(1) = 0;
#endif

    fBulk = blitz::RectDomain<3>(blBound, buBound);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to create the wall slices for the sub-domains
 *
 *          The wall slices of the sub-domain are for imposing the full-domain boundary conditions and hence are of importance
 *          only for the near-boundary sub-domains in parallel computations.
 *          Moreover, only the collocated grid has points on the boundary.
 *          The staggered grid points lie on either side of the domain boundaries.
 *          As a result the wall slices are defined only for those fields for which at least one of \ref field#xStag "xStag",
 *          \ref field#yStag "yStag" or \ref field#zStag "zStag" is false
 ********************************************************************************************************************************************
 */
void field::setWallSlices() {
    blitz::Array<blitz::TinyVector<int, 3>, 1> wlBound;
    blitz::Array<blitz::TinyVector<int, 3>, 1> wuBound;

    // 6 slices are stored in fWalls corresponding to the 6 faces of the 3D box
    fWalls.resize(6);

    wlBound.resize(6);
    wuBound.resize(6);

    // Wall slices are the locations where the BC (both Neumann and Dirichlet) is imposed.
    // In the places where these slices are being used, they should be on the LHS of equation.
    for (int i=0; i<6; i++) {
        wlBound(i) = F.lbound();
        wuBound(i) = F.ubound();
    }

    // The bulk slice corresponds to the part of the fluid within which all variables are computed at each time step.
    // Correspondingly, the boundary conditions are imposed on the layer just outside the bulk

    // UPPER BOUNDS OF LEFT WALL
    wlBound(0)(0) = wuBound(0)(0) = fBulk.lbound(0) - 1;

    // LOWER BOUNDS OF RIGHT WALL
    wuBound(1)(0) = wlBound(1)(0) = fBulk.ubound(0) + 1;

    // UPPER BOUNDS OF FRONT WALL
    wlBound(2)(1) = wuBound(2)(1) = fBulk.lbound(1) - 1;

    // LOWER BOUNDS OF BACK WALL
    wuBound(3)(1) = wlBound(3)(1) = fBulk.ubound(1) + 1;

    // UPPER BOUNDS OF BOTTOM WALL
    wlBound(4)(2) = wuBound(4)(2) = fBulk.lbound(2) - 1;

    // LOWER BOUNDS OF TOP WALL
    wuBound(5)(2) = wlBound(5)(2) = fBulk.ubound(2) + 1;

    for (int i=0; i<6; i++) {
        fWalls(i) = blitz::RectDomain<3>(wlBound(i), wuBound(i));
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to create the slices for interpolation while computing convective derivative
 *
 *          This function must be called before using the values contained in the arrays d2F_dx2, d2F_dy2 and d2F_dz2.
 ********************************************************************************************************************************************
 */
void field::setInterpolationSlices() {
    // INTERPOLATION SLICES FOR INTERPOLATING VALUES OF Vx FROM THE vfield
    // In all the below slices, we are considering interpolations between the following 8 variables
    //
    // Vx, Vy, Vz - these are face centered variables sitting on X, Y and Z planes respectively (like velocity)
    // Wx, Wy, Wz - these are edge centered variables sitting on X, Y and Z axes respectively (like vorticity)
    // Pc - this a cell centered variable (like temperature)
    // Qv - this a vertex centered variable (I don't know if anything sits here. But hey! completeness!)
    if (not xStag) {
        if (yStag) {
            if (zStag) {
                /* coll stag stag */
                /**X-Face centered configuration - Vx **/

                /* Interpolation of data - Vx ---> Vx
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => no change
                 *      - Z direction => no change
                 **/
                VxIntSlices.resize(1);
                VxIntSlices(0) = fCore;

                /* Interpolation of data - Vy ---> Vx
                 * Interpolation types:
                 *      - X direction => staggered to collocated
                 *      - Y direction => collocated to staggered
                 *      - Z direction => no change
                 **/
                VyIntSlices.resize(4);
                VyIntSlices(0) = fCore;
                VyIntSlices(1) = shift(1, VyIntSlices(0), -1);
                VyIntSlices(2) = shift(0, VyIntSlices(0), 1);
                VyIntSlices(3) = shift(1, VyIntSlices(2), -1);

                /* Interpolation of data - Vz ---> Vx
                 * Interpolation types:
                 *      - X direction => staggered to collocated
                 *      - Y direction => no change
                 *      - Z direction => collocated to staggered
                 **/
                VzIntSlices.resize(4);
                VzIntSlices(0) = fCore;
                VzIntSlices(1) = shift(2, VzIntSlices(0), -1);
                VzIntSlices(2) = shift(0, VzIntSlices(0), 1);
                VzIntSlices(3) = shift(2, VzIntSlices(2), -1);

                /* Interpolation of data - Pc ---> Vx
                 * Interpolation types:
                 *      - X direction => staggered to collocated
                 *      - Y direction => no change
                 *      - Z direction => no change
                 **/
                PcIntSlices.resize(2);
                PcIntSlices(0) = fCore;
                PcIntSlices(1) = shift(0, PcIntSlices(0), 1);

            } else {
                /* coll stag coll */
                /**Y-Edge centered configuration - Wy **/

                /* Interpolation of data - Vx ---> Wy
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => no change
                 *      - Z direction => staggered to collocated
                 **/
            }
        } else {
            if (zStag) {
                /* coll coll stag */
                /**Z-Edge centered configuration - Wz **/

                /* Interpolation of data - Vx ---> Wz
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => staggered to collocated
                 *      - Z direction => no change
                 **/
            } else {
                /* coll coll coll */
                /**Vertex centered configuration - Qv **/

                /* Interpolation of data - Vx ---> Qv
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => staggered to collocated
                 *      - Z direction => staggered to collocated
                 **/
            }
        }
    } else {
        if (yStag) {
            if (zStag) {
                /* stag stag stag */
                /**Cell centered configuration - Pc **/

                /* Interpolation of data - Vx ---> Pc
                 * Interpolation types:
                 *      - X direction => collocated to staggered
                 *      - Y direction => no change
                 *      - Z direction => no change
                 **/
                VxIntSlices.resize(2);
                VxIntSlices(0) = fCore;
                VxIntSlices(1) = shift(0, VxIntSlices(0), -1);

                /* Interpolation of data - Vy ---> Pc
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => collocated to staggered
                 *      - Z direction => no change
                 **/
                VyIntSlices.resize(2);
                VyIntSlices(0) = fCore;
                VyIntSlices(1) = shift(1, VyIntSlices(0), -1);

                /* Interpolation of data - Vz ---> Pc
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => no change
                 *      - Z direction => collocated to staggered
                 **/
                VzIntSlices.resize(2);
                VzIntSlices(0) = fCore;
                VzIntSlices(1) = shift(2, VzIntSlices(0), -1);

                /* Interpolation of data - Pc ---> Pc
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => no change
                 *      - Z direction => no change
                 **/
                PcIntSlices.resize(1);
                PcIntSlices(0) = fCore;

            } else {
                /* stag stag coll */
                /**Z-Face centered configuration - Vz **/

                /* Interpolation of data - Vx ---> Vz
                 * Interpolation types:
                 *      - X direction => collocated to staggered
                 *      - Y direction => no change
                 *      - Z direction => staggered to collocated
                 **/
                VxIntSlices.resize(4);
                VxIntSlices(0) = fCore;
                VxIntSlices(1) = shift(0, VxIntSlices(0), -1);
                VxIntSlices(2) = shift(2, VxIntSlices(0), 1);
                VxIntSlices(3) = shift(0, VxIntSlices(2), -1);

                /* Interpolation of data - Vy ---> Vz
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => collocated to staggered
                 *      - Z direction => staggered to collocated
                 **/
                VyIntSlices.resize(4);
                VyIntSlices(0) = fCore;
                VyIntSlices(1) = shift(1, VyIntSlices(0), -1);
                VyIntSlices(2) = shift(2, VyIntSlices(0), 1);
                VyIntSlices(3) = shift(1, VyIntSlices(2), -1);

                /* Interpolation of data - Vz ---> Vz
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => no change
                 *      - Z direction => no change
                 **/
                VzIntSlices.resize(1);
                VzIntSlices(0) = fCore;

                /* Interpolation of data - Pc ---> Vz
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => no change
                 *      - Z direction => staggered to collocated
                 **/
                PcIntSlices.resize(2);
                PcIntSlices(0) = fCore;
                PcIntSlices(1) = shift(2, PcIntSlices(0), 1);

            }
        } else {
            if (zStag) {
                /* stag coll stag */
                /**Y-Face centered configuration - Vy **/

                /* Interpolation of data - Vx ---> Vy
                 * Interpolation types:
                 *      - X direction => collocated to staggered
                 *      - Y direction => staggered to collocated
                 *      - Z direction => no change
                 **/
                VxIntSlices.resize(4);
                VxIntSlices(0) = fCore;
                VxIntSlices(1) = shift(0, VxIntSlices(0), -1);
                VxIntSlices(2) = shift(1, VxIntSlices(0), 1);
                VxIntSlices(3) = shift(0, VxIntSlices(2), -1);

                /* Interpolation of data - Vy ---> Vy
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => no change
                 *      - Z direction => no change
                 **/
                VyIntSlices.resize(1);
                VyIntSlices(0) = fCore;

                /* Interpolation of data - Vz ---> Vy
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => staggered to collocated
                 *      - Z direction => collocated to staggered
                 **/
                VzIntSlices.resize(4);
                VzIntSlices(0) = fCore;
                VzIntSlices(1) = shift(2, VzIntSlices(0), -1);
                VzIntSlices(2) = shift(1, VzIntSlices(0), 1);
                VzIntSlices(3) = shift(2, VzIntSlices(2), -1);

                /* Interpolation of data - Pc ---> Vy
                 * Interpolation types:
                 *      - X direction => no change
                 *      - Y direction => staggered to collocated
                 *      - Z direction => no change
                 **/
                PcIntSlices.resize(2);
                PcIntSlices(0) = fCore;
                PcIntSlices(1) = shift(1, PcIntSlices(0), 1);

            } else {
                /* stag coll coll */
                /**X-Edge centered configuration - Wx **/

                /* Interpolation of data - Vx ---> Wx
                 * Interpolation types:
                 *      - X direction => collocated to staggered
                 *      - Y direction => staggered to collocated
                 *      - Z direction => staggered to collocated
                 **/
            }
        }
    }

    // RESET INTERPOLATION SLICES FOR PLANAR GRID
#ifdef PLANAR
    if (fieldName == "Vy") {
        VxIntSlices.resize(1);
        VyIntSlices.resize(1);
        VzIntSlices.resize(1);

        VxIntSlices(0) = fCore;
        VyIntSlices(0) = fCore;
        VzIntSlices(0) = fCore;
    } else if (fieldName == "Vx") {
        VyIntSlices.resize(1);

        VyIntSlices(0) = fCore;
    } else if (fieldName == "Vz") {
        VyIntSlices.resize(1);

        VyIntSlices(0) = fCore;
    }
#endif
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to synchronise data across all processors when performing parallel computations
 *
 *          This function calls the \ref mpidata#syncData "syncData" function of mpidata class to perform perform data-transfer and thus update
 *          the sub-domain boundary pads.
 ********************************************************************************************************************************************
 */
void field::syncData() {
    mpiHandle->syncData();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to extract the maximum value from the field
 *
 *          The function uses the in-built blitz function to obtain the maximum value in an array.
 *          Note that this function *takes the maximum of the absolute value* of the field.
 *          While performing parallel computation, the function performs an <B>MPI_Allreduce()</B> to get
 *          the global maximum from the entire computational domain.
 *
 * \return  The real value of the maximum is returned (it is implicitly assumed that only real values are used)
 ********************************************************************************************************************************************
 */
real field::fieldMax() {
    real localMax, globalMax;

    localMax = blitz::max(blitz::abs(F));

    /***************************************************************************************************************
     * DID YOU KNOW?                                                                                               *
     * In the line above, most compilers will not complain even if you omitted the namespace specification blitz:: *
     * This behaviour wasted an hour of my development time (including the effort of making this nice box).        *
     * Check Ref. [4] in README for explanation.                                                                   *
     ***************************************************************************************************************/

    MPI_Allreduce(&localMax, &globalMax, 1, MPI_FP_REAL, MPI_MAX, MPI_COMM_WORLD);

    return globalMax;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given field
 *
 *          The unary operator += adds a given field to the field stored by the class and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another to be added to the member field
 *
 * \return  A pointer to itself is returned by the field class to which the operator belongs
 ********************************************************************************************************************************************
 */
field& field::operator += (field &a) {
    F += a.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given field
 *
 *          The unary operator -= subtracts a given field from the field stored by the class and returns
 *          a pointer to itself.
 *
 * \param   a is a reference to another field to be deducted from the member field
 *
 * \return  A pointer to itself is returned by the field class to which the operator belongs
 ********************************************************************************************************************************************
 */
field& field::operator -= (field &a) {
    F -= a.F;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to add a given scalar value
 *
 *          The unary operator += adds a given constant scalar value to the field stored by the class and returns
 *          a pointer to itself.
 *
 * \param   a is a real number to be added to the field
 *
 * \return  A pointer to itself is returned by the field class to which the operator belongs
 ********************************************************************************************************************************************
 */
field& field::operator += (real a) {
    F += a;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to subtract a given scalar value
 *
 *          The unary operator -= subtracts a given constant scalar value from the field stored by the class and returns
 *          a pointer to itself.
 *
 * \param   a is a real number to be subtracted from the field
 *
 * \return  A pointer to itself is returned by the field class to which the operator belongs
 ********************************************************************************************************************************************
 */
field& field::operator -= (real a) {
    F -= a;

    return *this;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a scalar value to the field
 *
 *          The operator = assigns a real value to the entire field.
 *
 * \param   a is a real number to be assigned to the field
 ********************************************************************************************************************************************
 */
void field::operator = (real a) {
    F = a;
}

/**
 ********************************************************************************************************************************************
 * \brief   Overloaded operator to assign a field to the field
 *
 *          The operator = copies the contents of the input field to itself.
 *
 * \param   a is the field to be assigned to the field
 ********************************************************************************************************************************************
 */
void field::operator = (field &a) {
    F = a.F;
}

field::~field() { }
