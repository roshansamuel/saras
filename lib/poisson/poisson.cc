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
/*! \file poisson.cc
 *
 *  \brief Definitions for functions of class poisson
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
 * \brief   Constructor of the base poisson class
 *
 *          The short base constructor of the poisson class merely assigns the const references to the grid and parser
 *          class instances being used in the solver.
 *          Moreover, it resizes and populates a local array of multi-grid sizes as used in the grid class.
 *          An array of strides to be used at different V-cycle levels is also generated and stored.
 *          Finally, the maximum allowable number of iterations for the Jacobi iterative solver being used at the
 *          coarsest mesh is set as \f$ N_{max} = N_x \times N_y \times N_z \f$, where \f$N_x\f$, \f$N_y\f$ and \f$N_z\f$
 *          are the number of grid points in the collocated grid at the local sub-domains along x, y and z directions
 *          respectively.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   solParam is a const reference to the user-set parameters contained in the parser class
 ********************************************************************************************************************************************
 */
poisson::poisson(const grid &mesh, const parser &solParam): mesh(mesh), inputParams(solParam) {
    int maxIndex = 15;

    mgSizeArray.resize(maxIndex);
    for (int i=0; i < maxIndex; i++) {
        mgSizeArray(i) = int(pow(2, i)) + 1;
    }

    mgSizeArray(0) = 1;

    strideValues.resize(inputParams.vcDepth + 1);
    for (int i=0; i<=inputParams.vcDepth; i++) {
        strideValues(i) = int(pow(2, i));
    }

    vLevel = 0;
    maxCount = mesh.collocCoreSize(0)*mesh.collocCoreSize(1)*mesh.collocCoreSize(2);

#ifdef TIME_RUN
    solveTimeComp = 0.0;
    solveTimeTran = 0.0;
    smothTimeComp = 0.0;
    smothTimeTran = 0.0;
#endif
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to initialize the arrays used in multi-grid
 *
 *          The memory required for various arrays in multi-grid solver are pre-allocated through this function.
 *          The function is called from within the constructor to perform this allocation once and for all.
 *          The arrays are initialized to 0.
 ********************************************************************************************************************************************
 */
void poisson::initializeArrays() {
    pressureData.resize(blitz::TinyVector<int, 3>(stagFull.ubound() - stagFull.lbound() + 1));
    pressureData.reindexSelf(stagFull.lbound());
    pressureData = 0.0;

    smoothedPres.resize(blitz::TinyVector<int, 3>(stagFull.ubound() - stagFull.lbound() + 1));
    smoothedPres.reindexSelf(stagFull.lbound());
    smoothedPres = 0.0;

    iteratorTemp.resize(blitz::TinyVector<int, 3>(stagFull.ubound() - stagFull.lbound() + 1));
    iteratorTemp.reindexSelf(stagFull.lbound());
    iteratorTemp = 0.0;

    inputRHSData.resize(blitz::TinyVector<int, 3>(stagFull.ubound() - stagFull.lbound() + 1));
    inputRHSData.reindexSelf(stagFull.lbound());
    inputRHSData = 0.0;

    residualData.resize(blitz::TinyVector<int, 3>(stagFull.ubound() - stagFull.lbound() + 1));
    residualData.reindexSelf(stagFull.lbound());
    residualData = 0.0;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to solve the poisson equation at the coarsest multi-grid level
 *
 *          This function operates exclusively at the lowest level of the multi-grid V-cycle.
 *          It uses the Jacobi iterative solver to solve the residual of the Poisson equation on the coarsest mesh.
 *          Note that the all calculations are performed assuming that the \ref vLevel variable is maximal when the function
 *          is being called.
 ********************************************************************************************************************************************
 */
void poisson::solve() { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to perform prolongation on the array being solved
 *
 *          Prolongation makes the grid finer by averaging values at two adjacent nodes onto an intermediate point between them
 *          As a result, the number of points in the domain increases from \f$ 2^N + 1 \f$ at the input level to \f$ 2^{N+1} + 1 \f$.
 *          The vLevel variable is accordingly reduced by 1 to reflect this ascent by one step up the V-Cycle.
 ********************************************************************************************************************************************
 */
void poisson::prolong() { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to perform smoothing operation on the input array
 *
 *          The smoothing operation is always performed on the data contained in the array \ref pressureData.
 *          The array \ref iteratorTemp is used to store the temporary data and it is continuously swapped with the
 *          \ref pressureData array at every iteration.
 *          This operation can be performed at any level of the V-cycle.
 *
 * \param   smoothCount is the integer value of the number of smoothing iterations to be performed
 ********************************************************************************************************************************************
 */
void poisson::smooth(const int smoothCount) { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to initialize the Range objects for accessing mesh derivatives in transformed plane
 *
 *          The Range objects defined here are used for reading the values of grid metrics at all the V-cycle levels.
 *          Since these values are required at all the grid levels, there are \ref parser#vcDepth "vcDepth" + 1 number of Range objects.
 ********************************************************************************************************************************************
 */
void poisson::initMeshRanges() { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to set the RectDomain variables for all future references throughout the poisson solver
 *
 *          The function sets the core and full domain staggered grid sizes for all the sub-domains.
 ********************************************************************************************************************************************
 */
void poisson::setStagBounds() { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate the local size indices of sub-domains after MPI domain decomposition
 *
 *          For the multi-grid solver, the number of grid nodes must be \f$ 2^N + 1 \f$ to perform V-Cycles
 *          The domain decomposition is also done in such a way that each sub-domain will also have \f$ 2^M + 1 \f$ points.
 *          As a result, the number of processors in each direction must be a power of 2.
 *          In this case, the the value \f$ M \f$ of \f$ 2^M + 1 \f$ can be computed as \f$ M = N - log_2(np) \f$
 *          where \f$ np \f$ is the number of processors along the direction under consideration.
 ********************************************************************************************************************************************
 */
void poisson::setLocalSizeIndex() { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to set the coefficients used for calculating laplacian and in smoothing
 *
 *          The function assigns values to the variables \ref hx, \ref hy etc.
 *          These coefficients are repeatedly used at many places in the Poisson solver.
 ********************************************************************************************************************************************
 */
void poisson::setCoefficients() { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to copy the staggered grid derivatives from the grid class to local arrays
 *
 *          Though the grid derivatives in the grid class can be read and accessed, they cannot be used directly
 *          along with the arrays defined in the poisson class as the local arrays have wide pads on both sides.
 *          Therefore, correspondingly wide arrays for grid derivatives are used, into which the staggered grid
 *          derivatives from the grid class are written and stored.
 *          This function serves this purpose of copying the grid derivatives.
 ********************************************************************************************************************************************
 */
void poisson::copyStaggrDerivs() { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to impose the boundary conditions of Poisson solver at different levels of the V-cycle
 *
 *          This function is called mainly during smoothing operations to impose the boundary conditions for the
 *          Poisson equation.
 *          The sub-domains close to the wall will have the Neumann boundary condition on pressure imposeed at the walls.
 *          Meanwhile at the interior boundaries at the inter-processor sub-domains, data is transferred from the neighbouring cells
 *          by calling the \ref updatePads function.
 ********************************************************************************************************************************************
 */
void poisson::imposeBC() { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to update the pad points of the local sub-domains at different levels of the V-cycle
 *
 *          This function is called mainly during smoothing operations by the \ref imposeBC function.
 *          It is also specifically called by the \ref solve function while solving the equation at the coarsest mesh level.
 *          At the interior boundaries at the inter-processor sub-domains, data is transferred from the neighbouring cells
 *          using a combination of MPI_Irecv and MPI_Send functions.
 ********************************************************************************************************************************************
 */
void poisson::updatePads() { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to create the MPI sub-array data types necessary to transfer data across sub-domains
 *
 *          The inter-domain boundaries of all the sub-domains at different V-cycle levels need data to be transfered at
 *          with different mesh strides.
 *          The number of sub-arrays along each edge/face of the sub-domains are equal to the number of V-cycle levels.
 *          Since this data transfer has to take place at all the mesh levels including the finest mesh, there will be
 *          vcDepth + 1 elements.
 ********************************************************************************************************************************************
 */
void poisson::createMGSubArrays() { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to perform one loop of V-cycle
 *
 *          The V-cycle of restrictions, prolongations and smoothings are performed within this function.
 *          First the input data contained in \ref pressureData is smoothed, after which the residual is computed and stored
 *          in the \ref residualData array.
 *          The restrictions, smoothing, and prolongations are performed on these two arrays subsequently.
 ********************************************************************************************************************************************
 */
void poisson::vCycle() { };

/**
 ********************************************************************************************************************************************
 * \brief   The core, publicly accessible function of poisson to compute the solution for the Poisson equation
 *
 *          The function calls the V-cycle as many times as set by the user.
 *          Before doing so, the input data is transferred into the data-structures used by the poisson class to
 *          perform restrictions and prolongations without copying.
 *          Finally, the computed solution is transferred back from the internal data-structures back into the
 *          scalar field supplied by the calling function.
 *
 * \param   inFn is a pointer to the plain scalar field (cell-centered) into which the computed soltuion must be transferred
 * \param   rhs is a const reference to the plain scalar field (cell-centered) which contains the RHS for the Poisson equation to solve
 ********************************************************************************************************************************************
 */
void poisson::mgSolve(plainsf &inFn, const plainsf &rhs) { };

/**
 ********************************************************************************************************************************************
 * \brief   Function to test whether strided data transfer is performing as expected
 *
 *          The function populates the arrays with predetermined values at all locations.
 *          It then calls updatePads at different vLevels and checks is the data is being transferred along x and y directions
 *          This done by printing the contents of the arrays for visual inspection for now.
 ********************************************************************************************************************************************
 */
real poisson::testTransfer() { return 0; };

/**
 ********************************************************************************************************************************************
 * \brief   Function to test whether prolongation operations are interpolating correctly
 *
 *          The function populates the arrays with predetermined values at all locations.
 *          It then calls prolong function at a lower vLevel and checks if the data is being interpolated correctly at higher vLevel
 *          This done by returning the average deviation from correct values as a real value
 ********************************************************************************************************************************************
 */
real poisson::testProlong() { return 0; };

/**
 ********************************************************************************************************************************************
 * \brief   Function to test whether periodic BC is being implemented properly
 *
 *          The function populates the arrays with predetermined values at all locations.
 *          It then calls imposeBC function at different vLevels and checks if the correct values of the functions are imposed at boundaries
 *          This done by printing the contents of the arrays for visual inspection for now.
 ********************************************************************************************************************************************
 */
real poisson::testPeriodic() { return 0; };

/**
 ********************************************************************************************************************************************
 * \brief   Function to test the solver used at the coarsest mesh level of V-Cycle
 *
 *          The function populates the arrays with predetermined values at all locations.
 *          It then calls solve function at lowest vLevel and compares the exact analytical values
 *          This done by printing the contents of the arrays for visual inspection for now.
 ********************************************************************************************************************************************
 */
real poisson::testSolve() { return 0; };

poisson::~poisson() {
#ifdef TIME_RUN
    if (mesh.rankData.rank == 0) {
        std::cout << std::left << std::setw(50) << "Time taken in computation within solve: "            << std::fixed << std::setprecision(6) << solveTimeComp << std::endl;
        std::cout << std::left << std::setw(50) << "Time taken in data-transfer within solve: "          << std::fixed << std::setprecision(6) << solveTimeTran << std::endl;
        std::cout << std::left << std::setw(50) << "Time taken in computation within smooth: "           << std::fixed << std::setprecision(6) << smothTimeComp << std::endl;
        std::cout << std::left << std::setw(50) << "Time taken in data-transfer within smooth: "         << std::fixed << std::setprecision(6) << smothTimeTran << std::endl;
    }
#endif
};
