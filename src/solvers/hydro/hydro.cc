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
/*! \file hydro.cc
 *
 *  \brief Definitions of common functions for both 2D and 3D runs of the solver class hydro - this class solves the basic Navier-Stokes equation.
 *  \sa hydro.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "hydro.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the base hydro class
 *
 *          The short base constructor of the hydro class merely assigns the const references to the grid and parser
 *          class instances being used in the solver.
 *          Also, the maximum allowable number of iterations for the Jacobi iterative solver being used to solve for the
 *          velocities implicitly is set as \f$ N_{max} = N_x \times N_y \times N_z \f$, where \f$N_x\f$, \f$N_y\f$ and \f$N_z\f$
 *          are the number of grid points in the collocated grid at the local sub-domains along x, y and z directions
 *          respectively.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   solParam is a const reference to the user-set parameters contained in the parser class
 ********************************************************************************************************************************************
 */
hydro::hydro(const grid &mesh, const parser &solParam, parallel &mpiParam):
            V(mesh, "V"),
            P(mesh, "P"),
            mesh(mesh),
            inputParams(solParam),
            inverseRe(1.0/inputParams.Re),
            mpiData(mpiParam),
            Pp(mesh, P),
            mgRHS(mesh, P),
            nseRHS(mesh, V),
            velocityLaplacian(mesh, V),
            pressureGradient(mesh, V),
            guessedVelocity(mesh, V)
{
    maxIterations = mesh.collocCoreSize(0)*mesh.collocCoreSize(1)*mesh.collocCoreSize(2);
}


/**
 ********************************************************************************************************************************************
 * \brief   The core publicly accessible function of the \ref hydro class to solve the Navier-Stokes equations
 *
 *          The NSE are integrated in time from within this function by calling \ref hydro#timeAdvance in a loop.
 *          The function keeps track of the non-dimensional time with \ref time and number of iterations with \ref iterCount.
 *          Both these values are continuously incremented from within the loop, and finally, when \ref time has reached the
 *          user-ser value in \ref parser#tMax "tMax", the time-integration loop is broken and the program exits.
 ********************************************************************************************************************************************
 */
void hydro::solvePDE() { };


/**
 ********************************************************************************************************************************************
 * \brief   The subroutine to solve the NS equations using the implicit Crank-Nicholson method
 *
 *          This function uses the values of velocity vector field and pressure scalar field, along with a specifed time-step
 *          to update the values of both fields by one time-step.
 *          Hence this function has to be repeatedly called in a loop from within the \ref solvePDE function to solve the equations.
 ********************************************************************************************************************************************
 */
void hydro::timeAdvance() { };


/**
 ********************************************************************************************************************************************
 * \brief   Function to solve the implicit equation for x-velocity
 *
 *          The implicit equation for \f$ u_x' \f$ of the implicit Crank-Nicholson method is solved using the Jacobi
 *          iterative method here.
 *
 *          The loop exits when the global maximum of the error in computed solution obtained using the \ref plainvf#vxMax "vxMax" function
 *          of scalar fields in sfield.h falls below the specified tolerance.
 *          If the solution doesn't converge even after an internally assigned maximum number for iterations, the solver
 *          aborts with an error message.
 *
 *          Note that this function uses the blitz index place holders firstIndex, secondIndex and thirdIndex.
 *          They are declared as i, j, and k respectively.
 *          Hence the variables i, j and k are not scalars in this function.
 ********************************************************************************************************************************************
 */
void hydro::solveVx() { };


/**
 ********************************************************************************************************************************************
 * \brief   Function to solve the implicit equation for y-velocity
 *
 *          The implicit equation for \f$ u_y' \f$ of the implicit Crank-Nicholson method is solved using the Jacobi
 *          iterative method here.
 *
 *          The loop exits when the global maximum of the error in computed solution obtained using the \ref plainvf#vyMax "vyMax" function
 *          of scalar fields in sfield.h falls below the specified tolerance.
 *          If the solution doesn't converge even after an internally assigned maximum number for iterations, the solver
 *          aborts with an error message.
 *
 *          Note that this function uses the blitz index place holders firstIndex, secondIndex and thirdIndex.
 *          They are declared as i, j, and k respectively.
 *          Hence the variables i, j and k are not scalars in this function.
 ********************************************************************************************************************************************
 */
void hydro::solveVy() { };


/**
 ********************************************************************************************************************************************
 * \brief   Function to solve the implicit equation for z-velocity
 *
 *          The implicit equation for \f$ u_z' \f$ of the implicit Crank-Nicholson method is solved using the Jacobi
 *          iterative method here.
 *
 *          The loop exits when the global maximum of the error in computed solution obtained using the \ref plainvf#vzMax "vzMax" function
 *          of scalar fields in sfield.h falls below the specified tolerance.
 *          If the solution doesn't converge even after an internally assigned maximum number for iterations, the solver
 *          aborts with an error message.
 *
 *          Note that this function uses the blitz index place holders firstIndex, secondIndex and thirdIndex.
 *          They are declared as i, j, and k respectively.
 *          Hence the variables i, j and k are not scalars in this function.
 ********************************************************************************************************************************************
 */
void hydro::solveVz() { };


/**
 ********************************************************************************************************************************************
 * \brief   Function to enable/disable periodic data transfer as per the problem
 *
 *          The function checks the xPer, yPer and zPer flags in the parser class
 *          and enables/disables MPI data transfer at boundaries accordingly
 *          By default, the MPI neighbours at boundaries are set for periodic data-transfer.
 *          This has to be disabled if the problem has non-periodic boundaries.
 *
 ********************************************************************************************************************************************
 */
void hydro::checkPeriodic() {
    // Disable periodic data transfer by setting neighbouring ranks of boundary sub-domains to NULL
    // Left and right walls
    if (not inputParams.xPer) {
        if (mpiData.rank == 0) {
            std::cout << "Using non-periodic boundary conditions along X Direction" << std::endl;
            std::cout << std::endl;
        }

        if (mpiData.xRank == 0)             mpiData.faceRanks(0) = MPI_PROC_NULL;
        if (mpiData.xRank == mpiData.npX-1) mpiData.faceRanks(1) = MPI_PROC_NULL;
    }

    // Front and rear walls
#ifdef PLANAR
    // Front and rear walls are by default non-periodic for 2D simulations
    if (mpiData.yRank == 0)             mpiData.faceRanks(2) = MPI_PROC_NULL;
    if (mpiData.yRank == mpiData.npY-1) mpiData.faceRanks(3) = MPI_PROC_NULL;

#else
    if (not inputParams.yPer) {
        if (mpiData.rank == 0) {
            std::cout << "Using non-periodic boundary conditions along Y Direction" << std::endl;
            std::cout << std::endl;
        }

        if (mpiData.yRank == 0)             mpiData.faceRanks(2) = MPI_PROC_NULL;
        if (mpiData.yRank == mpiData.npY-1) mpiData.faceRanks(3) = MPI_PROC_NULL;
    }
#endif

    // Inform user about BC along top and bottom walls
    if (not inputParams.zPer) {
        if (mpiData.rank == 0) {
            std::cout << "Using non-periodic boundary conditions along Z Direction" << std::endl;
            std::cout << std::endl;
        }
    }
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to set the coefficients used for solving the implicit equations of U, V and W
 *
 *          The function assigns values to the variables \ref hx, \ref hy, etc.
 *          These coefficients are repeatedly used at many places in the Poisson solver for implicit calculation of velocities.
 ********************************************************************************************************************************************
 */
void hydro::setCoefficients() {
    hx = mesh.dXi;
    hz = mesh.dZt;

    hz2hx2 = pow(mesh.dZt, 2.0)*pow(mesh.dXi, 2.0);

#ifdef PLANAR
    hx2 = pow(mesh.dXi, 2.0);
    hz2 = pow(mesh.dZt, 2.0);

#else
    hy = mesh.dEt;

    hx2hy2 = pow(mesh.dXi, 2.0)*pow(mesh.dEt, 2.0);
    hy2hz2 = pow(mesh.dEt, 2.0)*pow(mesh.dZt, 2.0);

    hx2hy2hz2 = pow(mesh.dXi, 2.0)*pow(mesh.dEt, 2.0)*pow(mesh.dZt, 2.0);
#endif
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to initialize the forcing terms for velocity
 *
 *          The forcing terms for the velocity field are initialized here.
 *          Out of the different forcings available in the force class,
 *          the appropriate forcing is chosen according to the parameters set by the user.
 ********************************************************************************************************************************************
 */
void hydro::initVForcing() {
    switch (inputParams.forceType) {
        case 0:
            if (mpiData.rank == 0) std::cout << "Running hydrodynamics simulation with zero velocity forcing" << std::endl << std::endl;
            vForcing = new zeroForcing(mesh, V);
            break;
        case 1:
            if (mpiData.rank == 0) std::cout << "Running hydrodynamics simulation with random velocity forcing" << std::endl << std::endl;
            vForcing = new randomForcing(mesh, V);
            break;
        case 2:
            if (mpiData.rank == 0) std::cout << "Running hydrodynamics simulation with rotation" << std::endl << std::endl;
            vForcing = new coriolisForce(mesh, V);
            break;
        case 5:
            if (mpiData.rank == 0) std::cout << "Running hydrodynamics simulation with constant pressure gradient along X-direction" << std::endl << std::endl;
            vForcing = new constantPGrad(mesh, V);
            break;
        default:
            if (mpiData.rank == 0) std::cout << "WARNING: Chosen velocity forcing is incompatible with hydrodynamics runs. Defaulting to zero forcing" << std::endl << std::endl;
            vForcing = new zeroForcing(mesh, V);
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to initialize the boundary conditions for velocity
 *
 *          The boundary conditions for all the 6 walls (4 in case of 2D simulations) are initialized here.
 *          Out of the different boundary conditions available in the boundary class,
 *          the appropriate BCs are chosen according to the type of problem being solved.
 ********************************************************************************************************************************************
 */
void hydro::initVBC() {
    if (inputParams.probType == 3) {
        // INFLOW AND OUTFLOW BCS
        uLft = new dirichletFC(mesh, V.Vx, 0, 1.0);
        uRgt = new neumannFC(mesh, V.Vx, 1, 0.0);
    } else {
        // NO-PENETRATION BCS
        uLft = new dirichletFC(mesh, V.Vx, 0, 0.0);
        uRgt = new dirichletFC(mesh, V.Vx, 1, 0.0);
    }

#ifndef PLANAR
    // NO-SLIP BCS
    uFrn = new dirichletCC(mesh, V.Vx, 2, 0.0);
    uBak = new dirichletCC(mesh, V.Vx, 3, 0.0);
#endif

    if (inputParams.zPer) {
        // PERIODIC BC
        uBot = new periodicCC(mesh, V.Vx, 4);
        uTop = new periodicCC(mesh, V.Vx, 5);
    } else {
        if (inputParams.probType == 1) {
            // NO-SLIP BCS FOR LDC
            uBot = new dirichletCC(mesh, V.Vx, 4, 0.0);
            uTop = new dirichletCC(mesh, V.Vx, 5, 1.0);
        } else {
            // NO-SLIP BCS
            uBot = new dirichletCC(mesh, V.Vx, 4, 0.0);
            uTop = new dirichletCC(mesh, V.Vx, 5, 0.0);
        }
    }

#ifndef PLANAR
    if (inputParams.probType == 3) {
        // INFLOW AND OUTFLOW BCS
        vLft = new dirichletCC(mesh, V.Vy, 0, 0.0);
        vRgt = new neumannCC(mesh, V.Vy, 1, 0.0);
    } else {
        // NO-SLIP BCS
        vLft = new dirichletCC(mesh, V.Vy, 0, 0.0);
        vRgt = new dirichletCC(mesh, V.Vy, 1, 0.0);
    }

    // NO-PENETRATION BCS
    vFrn = new dirichletFC(mesh, V.Vy, 2, 0.0);
    vBak = new dirichletFC(mesh, V.Vy, 3, 0.0);

    if (inputParams.zPer) {
        // PERIODIC BC
        vBot = new periodicCC(mesh, V.Vy, 4);
        vTop = new periodicCC(mesh, V.Vy, 5);
    } else {
        // NO-SLIP BCS
        vBot = new dirichletCC(mesh, V.Vy, 4, 0.0);
        vTop = new dirichletCC(mesh, V.Vy, 5, 0.0);
    }
#endif

    if (inputParams.probType == 3) {
        // INFLOW AND OUTFLOW BCS
        wLft = new dirichletCC(mesh, V.Vz, 0, 0.0);
        wRgt = new neumannCC(mesh, V.Vz, 1, 0.0);
    } else {
        // NO-SLIP BCS
        wLft = new dirichletCC(mesh, V.Vz, 0, 0.0);
        wRgt = new dirichletCC(mesh, V.Vz, 1, 0.0);
    }

#ifndef PLANAR
    // NO-SLIP BCS
    wFrn = new dirichletCC(mesh, V.Vz, 2, 0.0);
    wBak = new dirichletCC(mesh, V.Vz, 3, 0.0);
#endif

    if (inputParams.zPer) {
        // PERIODIC BC
        wBot = new periodicFC(mesh, V.Vz, 4);
        wTop = new periodicFC(mesh, V.Vz, 5);
    } else {
        // NO-SLIP BCS
        wBot = new dirichletFC(mesh, V.Vz, 4, 0.0);
        wTop = new dirichletFC(mesh, V.Vz, 5, 0.0);
    }
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to impose the boundary conditions for the X-component of velocity
 *
 *          The function first calls the \ref sfield#syncData "syncData" function of the Vx field to update the sub-domain pads.
 *          Then the boundary conditions are applied at the full domain boundaries by calling the imposeBC()
 *          of each boundary class object assigned to each wall.
 *          The order of imposing boundary conditions is - left, right, front, back, bottom and top boundaries.
 *          The corner values are not being imposed specifically and is thus dependent on the above order.
 *
 ********************************************************************************************************************************************
 */
void hydro::imposeUBCs() {
    V.Vx.syncData();

    if (not inputParams.xPer) {
        uLft->imposeBC();
        uRgt->imposeBC();
    }
#ifndef PLANAR
    if (not inputParams.yPer) {
        uFrn->imposeBC();
        uBak->imposeBC();
    }
#endif
    uTop->imposeBC();
    uBot->imposeBC();
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to impose the boundary conditions for the Y-component of velocity
 *
 *          The function first calls the \ref sfield#syncData "syncData" function of the Vy field to update the sub-domain pads.
 *          Then the boundary conditions are applied at the full domain boundaries by calling the imposeBC()
 *          of each boundary class object assigned to each wall.
 *          The order of imposing boundary conditions is - left, right, front, back, bottom and top boundaries.
 *          The corner values are not being imposed specifically and is thus dependent on the above order.
 *
 ********************************************************************************************************************************************
 */
void hydro::imposeVBCs() {
    V.Vy.syncData();

    if (not inputParams.xPer) {
        vLft->imposeBC();
        vRgt->imposeBC();
    }
#ifndef PLANAR
    if (not inputParams.yPer) {
        vFrn->imposeBC();
        vBak->imposeBC();
    }
#endif
    vTop->imposeBC();
    vBot->imposeBC();
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to impose the boundary conditions for the Z-component of velocity
 *
 *          The function first calls the \ref sfield#syncData "syncData" function of the Vz field to update the sub-domain pads.
 *          Then the boundary conditions are applied at the full domain boundaries by calling the imposeBC()
 *          of each boundary class object assigned to each wall.
 *          The order of imposing boundary conditions is - left, right, front, back, bottom and top boundaries.
 *          The corner values are not being imposed specifically and is thus dependent on the above order.
 *
 ********************************************************************************************************************************************
 */
void hydro::imposeWBCs() {
    V.Vz.syncData();

    if (not inputParams.xPer) {
        wLft->imposeBC();
        wRgt->imposeBC();
    }
#ifndef PLANAR
    if (not inputParams.yPer) {
        wFrn->imposeBC();
        wBak->imposeBC();
    }
#endif
    wTop->imposeBC();
    wBot->imposeBC();
};


/**
 ********************************************************************************************************************************************
 * \brief   Function to test whether periodic BC is being implemented properly
 *
 *          The function populates the arrays with predetermined values at all locations.
 *          It then calls imposeUBCs, imposeVBCs and imposeWBCs functions and checks if the correct values of the functions are imposed at boundaries
 ********************************************************************************************************************************************
 */
real hydro::testPeriodic() { return 0; };
