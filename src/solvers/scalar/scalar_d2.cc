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
/*! \file scalar_d2.cc
 *
 *  \brief Definitions of functions for 2D computations with the scalar class.
 *  \sa scalar.h
 *  \author Shashwat Bhattacharya, Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "scalar.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the scalar_d2 class derived from the base scalar class
 *
 *          The constructor passes its arguments to the base scalar class and then initializes all the scalar and
 *          vector fields necessary for solving the NS equations.
 *          The various coefficients for solving the equations are also set by a call to the \ref setCoefficients function.
 *          Based on the problem type specified by the user in the parameters file, and stored by the \ref parser class as
 *          \ref parser#probType "probType", the appropriate initial and boundary conditions are specified.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   solParam is a const reference to the user-set parameters contained in the parser class
 ********************************************************************************************************************************************
 */
scalar_d2::scalar_d2(const grid &mesh, const parser &solParam, parallel &mpiParam):
            scalar(mesh, solParam, mpiParam),
            mgSolver(mesh, inputParams)
{
    // SET VALUES OF COEFFICIENTS USED FOR COMPUTING LAPLACIAN
    setCoefficients();

    // INITIALIZE PRESSURE AND SCALAR
    if (inputParams.restartFlag) {
        // Fields to be read from HDF5 file are passed to reader class as a vector
        std::vector<field> readFields;

        // Populate the vector with required fields
        readFields.push_back(V.Vx);
        readFields.push_back(V.Vz);
        readFields.push_back(P.F);
        readFields.push_back(T.F);

        // Initialize reader object
        reader dataReader(mesh, readFields);

        time = dataReader.readData();

    } else {
        P = 1.0;
        T = 0.0;
        time = 0.0;
    }

    checkPeriodic();

    // Initialize velocity and temperature boundary conditions
    initVBC();
    initTBC();

    // Initialize velocity and temperature forcing fields
    initVForcing();
    initTForcing();

    imposeUBCs();
    imposeWBCs();

    imposeTBCs();
}


void scalar_d2::solvePDE() {

    real fwTime, prTime, rsTime;
    //set dt equal to input time step
    dt = inputParams.tStp;

    // Fields to be written into HDF5 file are passed to writer class as a vector
    std::vector<field> writeFields;

    // Populate the vector with required scalar fields
    writeFields.push_back(V.Vx);
    writeFields.push_back(V.Vz);
    writeFields.push_back(P.F);
    writeFields.push_back(T.F);

    // Initialize writer object
    writer dataWriter(mesh, writeFields);

    // Initialize probes
    if (inputParams.readProbes) {
        dataProbe = new probes(mesh, writeFields);
    }

    // Output file containing the time series of various variables
    tseries tsWriter(mesh, V, P, time, dt);

    // FILE WRITING TIME
    fwTime = time;

    // FIELD PROBING TIME
    prTime = time;

    // RESTART FILE WRITING TIME
    rsTime = time;

    timeStepCount = 0;

    // COMPUTE ENERGY AND DIVERGENCE FOR THE INITIAL CONDITION
    tsWriter.writeTSData(T, nu, kappa);

    // WRITE DATA AT t = 0
    if (not inputParams.restartFlag) {
        dataWriter.writeSolution(time);

        if (inputParams.readProbes) {
            dataProbe->probeData(time);
        }
    }

    // FOR RESTART RUNS, THE NEXT TIME FOR WRITING OUTPUT IS OBTAINED BY INCREMENTING WITH THE MOD OF RESTART TIME AND OUTPUT WRITE INTERVAL.
    // OTHERWISE THE WRITE INTERVAL IS ADDED DIRECTLY
    fwTime += inputParams.fwInt - std::fmod(time, inputParams.fwInt);
    prTime += inputParams.prInt - std::fmod(time, inputParams.prInt);
    rsTime += inputParams.rsInt - std::fmod(time, inputParams.rsInt);

    // TIME-INTEGRATION LOOP
    while (true) {
        // MAIN FUNCTION CALLED IN EACH LOOP TO UPDATE THE FIELDS AT EACH TIME-STEP
        computeTimeStep();
        if (inputParams.useCFL) {
            V.computeTStp(dt);
            if (dt > inputParams.tStp) {
                dt = inputParams.tStp;
            }
        }

        timeStepCount += 1;
        time += dt;

        if (timeStepCount % inputParams.ioCnt == 0) {
            tsWriter.writeTSData(T, nu, kappa);
        }

        if (inputParams.readProbes and std::abs(prTime - time) < 0.5*dt) {
            dataProbe->probeData(time);
            prTime += inputParams.prInt;
        }

        if (std::abs(fwTime - time) < 0.5*dt) {
            dataWriter.writeSolution(time);
            fwTime += inputParams.fwInt;
        }

        if (std::abs(rsTime - time) < 0.5*dt) {
            dataWriter.writeRestart(time);
            rsTime += inputParams.rsInt;
        }

        if (std::abs(inputParams.tMax - time) < 0.5*dt) {
            break;
        }
    }
}


void scalar_d2::computeTimeStep() {
    // BELOW FLAG MAY BE TURNED OFF FOR DEBUGGING/DIGNOSTIC RUNS ONLY
    // IT IS USED TO TURN OFF COMPUTATION OF NON-LINEAR TERMS
    // CURRENTLY IT IS AVAILABLE ONLY FOR THE 2D SCALAR SOLVER
    bool nlinSwitch = true;

    nseRHS = 0.0;
    tmpRHS = 0.0;

    // CALCULATE RHS OF NSE FROM THE NON LINEAR TERMS AND HALF THE VISCOUS TERMS
    V.computeDiff(nseRHS);
    nseRHS *= nu;

    // COMPUTE THE CONVECTIVE DERIVATIVE AND SUBTRACT IT FROM THE CALCULATED DIFFUSION TERMS OF RHS IN nseRHS
    V.computeNLin(V, nseRHS);

    //ADD VELOCITY FORCING TO THE RHS
    vForcing->addForcing(nseRHS);
    
    // RESET pressureGradient VFIELD AND CALCULATE THE PRESSURE GRADIENT
    pressureGradient = 0.0;
    P.gradient(pressureGradient, V);

    // ADD PRESSURE GRADIENT TO NON-LINEAR TERMS AND MULTIPLY WITH TIME-STEP
    nseRHS -= pressureGradient;
    nseRHS *= dt;

    // ADD THE CALCULATED VALUES TO THE VELOCITY AT START OF TIME-STEP
    nseRHS += V;

    // SYNCHRONISE THE RHS OF TIME INTEGRATION STEP THUS OBTAINED ACROSS ALL PROCESSORS
    nseRHS.syncData();

    // CALCULATE V IMPLICITLY USING THE JACOBI ITERATIVE SOLVER
    solveVx();
    solveVz();

    // CALCULATE THE RHS FOR THE POISSON SOLVER FROM THE GUESSED VALUES OF VELOCITY IN V
    V.divergence(mgRHS, P);
    mgRHS *= 1.0/dt;

    // USING THE CALCULATED mgRHS, EVALUATE Pp USING MULTI-GRID METHOD
    mgSolver.mgSolve(Pp, mgRHS);

    // SYNCHRONISE THE PRESSURE CORRECTION ACROSS PROCESSORS
    Pp.syncData();

    // ADD THE PRESSURE CORRECTION CALCULATED FROM THE POISSON SOLVER TO P
    P += Pp;

    // CALCULATE FINAL VALUE OF V BY SUBTRACTING THE GRADIENT OF PRESSURE CORRECTION
    Pp.gradient(pressureGradient, V);
    pressureGradient *= dt;
    V -= pressureGradient;

    // COMPUTE DIFFUSION AND NON-LINEAR TERMS FOR THE SCALAR EQUATION
    T.computeDiff(tmpRHS);
    tmpRHS *= kappa;

    if (nlinSwitch) {
        // COMPUTE THE CONVECTIVE DERIVATIVE AND SUBTRACT IT FROM THE CALCULATED DIFFUSION TERMS OF RHS IN tmpRHS
        T.computeNLin(V, tmpRHS);

    // EVEN WHEN NON-LINEAR TERM IS TURNED OFF, THE MEAN FLOW EFFECTS STILL REMAIN
    // HENCE THE CONTRIBUTION OF VELOCITY TO SCALAR EQUATION MUST BE ADDED
    // THIS CONTRIBUTION IS Uz FOR RBC AND SST, BUT Ux FOR VERTICAL CONVECTION
    } else {
        if (inputParams.probType == 5 || inputParams.probType == 6) {
            T.interTempF = 0.0;
            for (unsigned int i=0; i < T.F.VzIntSlices.size(); i++) {
                T.interTempF(T.F.fCore) += V.Vz.F(T.F.VzIntSlices(i));
            }

            tmpRHS.F += T.interTempF/T.F.VzIntSlices.size();

        } else if (inputParams.probType == 7) {
            T.interTempF = 0.0;
            for (unsigned int i=0; i < T.F.VxIntSlices.size(); i++) {
                T.interTempF(T.F.fCore) += V.Vx.F(T.F.VxIntSlices(i));
            }

            tmpRHS.F += T.interTempF/T.F.VxIntSlices.size();
        }
    }

    // ADD SCALAR FORCING TO THE TEMPERATURE EQUATION
    tForcing->addForcing(tmpRHS);

    // MULTIPLY WITH TIME-STEP AND ADD THE CALCULATED VALUE TO THE TEMPERATURE AT START OF TIME-STEP
    tmpRHS *= dt;
    tmpRHS += T;

    // SYNCHRONISE THE RHS OF TIME INTEGRATION STEP THUS OBTAINED ACROSS ALL PROCESSORS
    tmpRHS.syncData();

    // CALCULATE T IMPLICITLY USING THE JACOBI ITERATIVE SOLVER
    solveT();

    // IMPOSE BOUNDARY CONDITIONS ON V
    imposeUBCs();
    imposeWBCs();

    // IMPOSE BOUNDARY CONDITIONS ON T
    imposeTBCs();
}


void scalar_d2::solveVx() {
    int iterCount = 0;
    real maxError = 0.0;

    while (true) {
        int iY = 0;
#pragma omp parallel for num_threads(inputParams.nThreads) default(none) shared(iY)
        for (int iX = V.Vx.fBulk.lbound(0); iX <= V.Vx.fBulk.ubound(0); iX++) {
            for (int iZ = V.Vx.fBulk.lbound(2); iZ <= V.Vx.fBulk.ubound(2); iZ++) {
                guessedVelocity.Vx(iX, iY, iZ) = ((hz2 * mesh.xix2Colloc(iX) * (V.Vx.F(iX+1, iY, iZ) + V.Vx.F(iX-1, iY, iZ)) +
                                                     hx2 * mesh.ztz2Staggr(iZ) * (V.Vx.F(iX, iY, iZ+1) + V.Vx.F(iX, iY, iZ-1))) *
                        dt * nu / ( hz2hx2 * 2.0) + nseRHS.Vx(iX, iY, iZ))/
                    (1.0 + dt * nu * ((hz2 * mesh.xix2Colloc(iX) +
                                                     hx2 * mesh.ztz2Staggr(iZ)))/hz2hx2);
            }
        }

        V.Vx.F = guessedVelocity.Vx;

        imposeUBCs();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none) shared(iY)
        for (int iX = V.Vx.fBulk.lbound(0); iX <= V.Vx.fBulk.ubound(0); iX++) {
            for (int iZ = V.Vx.fBulk.lbound(2); iZ <= V.Vx.fBulk.ubound(2); iZ++) {
                velocityLaplacian.Vx(iX, iY, iZ) = V.Vx.F(iX, iY, iZ) - (
                                mesh.xix2Colloc(iX) * (V.Vx.F(iX+1, iY, iZ) - 2.0 * V.Vx.F(iX, iY, iZ) + V.Vx.F(iX-1, iY, iZ)) / (hx2) +
                                mesh.ztz2Staggr(iZ) * (V.Vx.F(iX, iY, iZ+1) - 2.0 * V.Vx.F(iX, iY, iZ) + V.Vx.F(iX, iY, iZ-1)) / (hz2)) *
                                                0.5 * dt * nu;
            }
        }

        velocityLaplacian.Vx(V.Vx.fBulk) = abs(velocityLaplacian.Vx(V.Vx.fBulk) - nseRHS.Vx(V.Vx.fBulk));

        maxError = velocityLaplacian.vxMax();

        if (maxError < inputParams.tolerance) {
            break;
        }

        iterCount += 1;

        if (iterCount > maxIterations) {
            if (mesh.rankData.rank == 0) {
                std::cout << "ERROR: Jacobi iterations for solution of Vx not converging. Aborting" << std::endl;
            }
            MPI_Finalize();
            exit(0);
        }
    }
}


void scalar_d2::solveVz() {
    int iterCount = 0;
    real maxError = 0.0;

    while (true) {
        int iY = 0;
#pragma omp parallel for num_threads(inputParams.nThreads) default(none) shared(iY)
        for (int iX = V.Vz.fBulk.lbound(0); iX <= V.Vz.fBulk.ubound(0); iX++) {
            for (int iZ = V.Vz.fBulk.lbound(2); iZ <= V.Vz.fBulk.ubound(2); iZ++) {
                guessedVelocity.Vz(iX, iY, iZ) = ((hz2 * mesh.xix2Staggr(iX) * (V.Vz.F(iX+1, iY, iZ) + V.Vz.F(iX-1, iY, iZ)) +
                                                     hx2 * mesh.ztz2Colloc(iZ) * (V.Vz.F(iX, iY, iZ+1) + V.Vz.F(iX, iY, iZ-1))) *
                        dt * nu / ( hz2hx2 * 2.0) + nseRHS.Vz(iX, iY, iZ))/
                    (1.0 + dt * nu * ((hz2 * mesh.xix2Staggr(iX) +
                                                     hx2 * mesh.ztz2Colloc(iZ)))/hz2hx2);
            }
        }

        V.Vz.F = guessedVelocity.Vz;

        imposeWBCs();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none) shared(iY)
        for (int iX = V.Vz.fBulk.lbound(0); iX <= V.Vz.fBulk.ubound(0); iX++) {
            for (int iZ = V.Vz.fBulk.lbound(2); iZ <= V.Vz.fBulk.ubound(2); iZ++) {
                velocityLaplacian.Vz(iX, iY, iZ) = V.Vz.F(iX, iY, iZ) - (
                                mesh.xix2Staggr(iX) * (V.Vz.F(iX+1, iY, iZ) - 2.0 * V.Vz.F(iX, iY, iZ) + V.Vz.F(iX-1, iY, iZ)) / (hx2) +
                                mesh.ztz2Colloc(iZ) * (V.Vz.F(iX, iY, iZ+1) - 2.0 * V.Vz.F(iX, iY, iZ) + V.Vz.F(iX, iY, iZ-1)) / (hz2)) *
                                                0.5 * dt * nu;
            }
        }

        velocityLaplacian.Vz(V.Vz.fBulk) = abs(velocityLaplacian.Vz(V.Vz.fBulk) - nseRHS.Vz(V.Vz.fBulk));

        maxError = velocityLaplacian.vzMax();

        if (maxError < inputParams.tolerance) {
            break;
        }

        iterCount += 1;

        if (iterCount > maxIterations) {
            if (mesh.rankData.rank == 0) {
                std::cout << "ERROR: Jacobi iterations for solution of Vz not converging. Aborting" << std::endl;
            }
            MPI_Finalize();
            exit(0);
        }
    }
}


void scalar_d2::solveT() {
    int iterCount = 0;
    real maxError = 0.0;

    while (true) {
        int iY = 0;
#pragma omp parallel for num_threads(inputParams.nThreads) default(none) shared(iY)
        for (int iX = T.F.fBulk.lbound(0); iX <= T.F.fBulk.ubound(0); iX++) {
            for (int iZ = T.F.fBulk.lbound(2); iZ <= T.F.fBulk.ubound(2); iZ++) {
                guessedScalar.F(iX, iY, iZ) = ((hz2 * mesh.xix2Staggr(iX) * (T.F.F(iX+1, iY, iZ) + T.F.F(iX-1, iY, iZ)) +
                                                hx2 * mesh.ztz2Colloc(iZ) * (T.F.F(iX, iY, iZ+1) + T.F.F(iX, iY, iZ-1))) *
                                          dt * kappa / ( hz2hx2 * 2.0) + tmpRHS.F(iX, iY, iZ))/
                                   (1.0 + dt * kappa * ((hz2 * mesh.xix2Staggr(iX) + hx2 * mesh.ztz2Colloc(iZ)))/hz2hx2);
            }
        }

        T = guessedScalar;

        imposeTBCs();

#pragma omp parallel for num_threads(inputParams.nThreads) default(none) shared(iY)
        for (int iX = T.F.fBulk.lbound(0); iX <= T.F.fBulk.ubound(0); iX++) {
            for (int iZ = T.F.fBulk.lbound(2); iZ <= T.F.fBulk.ubound(2); iZ++) {
                scalarLaplacian.F(iX, iY, iZ) = T.F.F(iX, iY, iZ) - (
                        mesh.xix2Staggr(iX) * (T.F.F(iX+1, iY, iZ) - 2.0 * T.F.F(iX, iY, iZ) + T.F.F(iX-1, iY, iZ)) / (hx2) +
                        mesh.ztz2Colloc(iZ) * (T.F.F(iX, iY, iZ+1) - 2.0 * T.F.F(iX, iY, iZ) + T.F.F(iX, iY, iZ-1)) / (hz2)) *
                                        0.5 * dt * kappa;
            }
        }

        scalarLaplacian.F(T.F.fBulk) = abs(scalarLaplacian.F(T.F.fBulk) - tmpRHS.F(T.F.fBulk));

        maxError = scalarLaplacian.fxMax();

        if (maxError < inputParams.tolerance) {
            break;
        }

        iterCount += 1;

        if (iterCount > maxIterations) {
            if (mesh.rankData.rank == 0) {
                std::cout << "ERROR: Jacobi iterations for solution of T not converging. Aborting" << std::endl;
            }
            MPI_Finalize();
            exit(0);
        }
    }
}


real scalar_d2::testPeriodic() {
    int iY = 0;
    real xCoord = 0.0;
    real zCoord = 0.0;

    nseRHS = 0.0;
    V = 0.0;

    for (int i=V.Vx.F.lbound(0); i <= V.Vx.F.ubound(0); i++) {
        for (int k=V.Vx.F.lbound(2); k <= V.Vx.F.ubound(2); k++) {
            V.Vx.F(i, 0, k) = sin(2.0*M_PI*mesh.xColloc(i)/mesh.xLen)*
                              cos(2.0*M_PI*mesh.zStaggr(k)/mesh.zLen);
            nseRHS.Vx(i, 0, k) = V.Vx.F(i, 0, k);
        }
    }
    for (int i=V.Vz.F.lbound(0); i <= V.Vz.F.ubound(0); i++) {
        for (int k=V.Vz.F.lbound(2); k <= V.Vz.F.ubound(2); k++) {
            V.Vz.F(i, 0, k) = -cos(2.0*M_PI*mesh.xStaggr(i)/mesh.xLen)*
                               sin(2.0*M_PI*mesh.zColloc(k)/mesh.zLen);
            nseRHS.Vz(i, 0, k) = V.Vz.F(i, 0, k);
        }
    }

    // EXPECTED VALUES IN THE PAD REGIONS OF Vx IF DATA TRANSFER HAPPENS WITH NO HITCH
    // X-VELOCITY IN LEFT AND RIGHT PADS
    for (int iX = 1; iX <= mesh.padWidths(0); iX++) {
        for (int iZ = V.Vx.fCore.lbound(2); iZ <= V.Vx.fCore.ubound(2); iZ += 1) {
            xCoord = mesh.xColloc(V.Vx.fCore.lbound(0)) - (mesh.xColloc(V.Vx.fCore.lbound(0) + iX) - mesh.xColloc(V.Vx.fCore.lbound(0)));
            nseRHS.Vx(V.Vx.fCore.lbound(0) - iX, iY, iZ) = sin(2.0*M_PI*xCoord/mesh.xLen)*
                                                           cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);

            xCoord = mesh.xColloc(V.Vx.fCore.ubound(0)) + (mesh.xColloc(V.Vx.fCore.ubound(0)) - mesh.xColloc(V.Vx.fCore.ubound(0) - iX));
            nseRHS.Vx(V.Vx.fCore.ubound(0) + iX, iY, iZ) = sin(2.0*M_PI*xCoord/mesh.xLen)*
                                                           cos(2.0*M_PI*mesh.zStaggr(iZ)/mesh.zLen);
        }
    }

    // X-VELOCITY IN BOTTOM AND TOP PADS
    for (int iZ = 1; iZ <= mesh.padWidths(2); iZ++) {
        for (int iX = V.Vx.fCore.lbound(0); iX <= V.Vx.fCore.ubound(0); iX += 1) {
            zCoord = mesh.zStaggr(V.Vx.fCore.lbound(2)) - (mesh.zStaggr(V.Vx.fCore.lbound(2) + iZ) - mesh.zStaggr(V.Vx.fCore.lbound(2)));
            nseRHS.Vx(iX, iY, V.Vx.fCore.lbound(2) - iZ) = sin(2.0*M_PI*mesh.xColloc(iX)/mesh.xLen)*
                                                           cos(2.0*M_PI*zCoord/mesh.zLen);

            zCoord = mesh.zStaggr(V.Vx.fCore.ubound(2)) + (mesh.zStaggr(V.Vx.fCore.ubound(2)) - mesh.zStaggr(V.Vx.fCore.ubound(2) - iZ));
            nseRHS.Vx(iX, iY, V.Vx.fCore.ubound(2) + iZ) = sin(2.0*M_PI*mesh.xColloc(iX)/mesh.xLen)*
                                                           cos(2.0*M_PI*zCoord/mesh.zLen);
        }
    }

    // Z-VELOCITY IN LEFT AND RIGHT PADS
    for (int iX = 1; iX <= mesh.padWidths(0); iX++) {
        for (int iZ = V.Vz.fCore.lbound(2); iZ <= V.Vz.fCore.ubound(2); iZ += 1) {
            xCoord = mesh.xStaggr(V.Vz.fCore.lbound(0)) - (mesh.xStaggr(V.Vz.fCore.lbound(0) + iX) - mesh.xStaggr(V.Vz.fCore.lbound(0)));
            nseRHS.Vz(V.Vz.fCore.lbound(0) - iX, iY, iZ) = -cos(2.0*M_PI*xCoord/mesh.xLen)*
                                                            sin(2.0*M_PI*mesh.zColloc(iZ)/mesh.zLen);

            xCoord = mesh.xStaggr(V.Vz.fCore.ubound(0)) + (mesh.xStaggr(V.Vz.fCore.ubound(0)) - mesh.xStaggr(V.Vz.fCore.ubound(0) - iX));
            nseRHS.Vz(V.Vz.fCore.ubound(0) + iX, iY, iZ) = -cos(2.0*M_PI*xCoord/mesh.xLen)*
                                                            sin(2.0*M_PI*mesh.zColloc(iZ)/mesh.zLen);
        }
    }

    // Z-VELOCITY IN BOTTOM AND TOP PADS
    for (int iZ = 1; iZ <= mesh.padWidths(2); iZ++) {
        for (int iX = V.Vz.fCore.lbound(0); iX <= V.Vz.fCore.ubound(0); iX += 1) {
            zCoord = mesh.zColloc(V.Vz.fCore.lbound(2)) - (mesh.zColloc(V.Vz.fCore.lbound(2) + iZ) - mesh.zColloc(V.Vz.fCore.lbound(2)));
            nseRHS.Vz(iX, iY, V.Vz.fCore.lbound(2) - iZ) = -cos(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                            sin(2.0*M_PI*zCoord/mesh.zLen);

            zCoord = mesh.zColloc(V.Vz.fCore.ubound(2)) + (mesh.zColloc(V.Vz.fCore.ubound(2)) - mesh.zColloc(V.Vz.fCore.ubound(2) - iZ));
            nseRHS.Vz(iX, iY, V.Vz.fCore.ubound(2) + iZ) = -cos(2.0*M_PI*mesh.xStaggr(iX)/mesh.xLen)*
                                                            sin(2.0*M_PI*zCoord/mesh.zLen);
        }
    }

    imposeUBCs();
    imposeWBCs();

    V -= nseRHS;

    return std::max(blitz::max(fabs(V.Vx.F)), blitz::max(fabs(V.Vz.F)));
}


scalar_d2::~scalar_d2() { }
