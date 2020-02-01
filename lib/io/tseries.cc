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
/*! \file tseries.cc
 *
 *  \brief Definitions for functions of class tseries
 *  \sa tseries.h
 *  \author Roshan Samuel, Ali Asad
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include <iostream>
#include "tseries.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the tseries class
 *
 *          The constructor initializes the variables and parameters for writing time-series data
 *          The header for the time-series data file is also written here.
 *          Depending on whether the run is for a scalar solver or hydro solver, the appropriate header is written.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   solverV is a reference to the velocity vector field whose data is used in calculating global quantities
 * \param   solverP is a const reference to the pressure scalar field whose dimensions are used to set array limits
 * \param   solverTime is a const reference to the real variable holding the value of current solver time
 * \param   timeStep is a const reference to the real variable holding the value of current time-step
 *
 ********************************************************************************************************************************************
 */
tseries::tseries(const grid &mesh, vfield &solverV, const sfield &solverP, const real &solverTime, const real &timeStep):
                 time(solverTime), tStp(timeStep), mesh(mesh), P(solverP), V(solverV), divV(mesh, P) {
    // Open TimeSeries file
    if (mesh.inputParams.restartFlag) {
        ofFile.open("output/TimeSeries.dat", std::fstream::out | std::fstream::app);
    } else {
        ofFile.open("output/TimeSeries.dat", std::fstream::out);
    }

    // UPPER AND LOWER LIMITS WHEN COMPUTING ENERGY IN STAGGERED GRID
    xLow = P.F.fCore.lbound(0);        xTop = P.F.fCore.ubound(0);
#ifndef PLANAR
    yLow = P.F.fCore.lbound(1);        yTop = P.F.fCore.ubound(1);
#endif
    zLow = P.F.fCore.lbound(2);        zTop = P.F.fCore.ubound(2);

    // STAGGERED GRIDS HAVE SHARED POINT ACROSS MPI-SUBDOMAINS - ACCORDINGLY DECREASE LIMITS
    if (mesh.rankData.xRank > 0) {
        xLow += 1;
    }

#ifndef PLANAR
    if (mesh.rankData.yRank > 0) {
        yLow += 1;
    }
#endif

    // TOTAL VOLUME FOR AVERAGING THE RESULT OF VOLUMETRIC INTEGRATION
    real localVol = 0.0;

    totalVol = 0.0;
#ifdef PLANAR
    for (int iX = xLow; iX <= xTop; iX++) {
        for (int iZ = zLow; iZ <= zTop; iZ++) {
            localVol += (mesh.dXi/mesh.xi_xColloc(iX))*(mesh.dZt/mesh.zt_zColloc(iZ));
        }
    }
#else
    for (int iX = xLow; iX <= xTop; iX++) {
        for (int iY = yLow; iY <= yTop; iY++) {
            for (int iZ = zLow; iZ <= zTop; iZ++) {
                localVol += (mesh.dXi/mesh.xi_xColloc(iX))*(mesh.dEt/mesh.et_yColloc(iY))*(mesh.dZt/mesh.zt_zColloc(iZ));
            }
        }
    }
#endif
    MPI_Allreduce(&localVol, &totalVol, 1, MPI_FP_REAL, MPI_SUM, MPI_COMM_WORLD);

    // WRITE THE HEADERS FOR BOTH STANDARD I/O AS WELL AS THE OUTPUT TIME-SERIES FILE
    if (mesh.rankData.rank == 0) {
        if (mesh.inputParams.probType <= 4) {
            std::cout << std::fixed << std::setw(6)  << "Time" << "\t" <<
                                    std::setw(16) << "Total energy" << "\t" << "Divergence" << std::endl;

            ofFile << "#VARIABLES = Time, Energy, Divergence, dt\n";
        } else {
            std::cout << std::fixed << std::setw(6)  << "Time" << "\t" <<
                                    std::setw(16) << "Re (Urms)" << "\t" << "Nusselt No" << "\t" << "Divergence" << std::endl;

            ofFile << "#VARIABLES = Time, Energy, NusseltNo, Divergence, dt\n";
        }
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Overloaded function to write the time-series data for hydro solver.
 *
 *          The function computes total energy and divergence for hydro solver runs.
 *          Only the root rank (rank 0) writes the output.
 *          One line is written to the standard I/O, while another is written to the time series dat file.
 *
 ********************************************************************************************************************************************
 */
void tseries::writeTSData() {
    V.divergence(divV, P);
    maxDivergence = divV.fxMax();

    if (maxDivergence > 1.0e5) {
        if (mesh.rankData.rank == 0) std::cout << std::endl << "ERROR: Divergence exceeds permissible limits. ABORTING" << std::endl << std::endl;
        MPI_Finalize();
        exit(0);
    }

    localEnergy = 0.0;
    totalEnergy = 0.0;
#ifdef PLANAR
    int iY = 0;
    for (int iX = xLow; iX <= xTop; iX++) {
        for (int iZ = zLow; iZ <= zTop; iZ++) {
            localEnergy += (pow((V.Vx.F(iX-1, iY, iZ) + V.Vx.F(iX, iY, iZ))/2.0, 2.0) +
                            pow((V.Vz.F(iX, iY, iZ-1) + V.Vz.F(iX, iY, iZ))/2.0, 2.0))*(mesh.dXi/mesh.xi_xColloc(iX))*(mesh.dZt/mesh.zt_zColloc(iZ));
        }
    }
#else
    for (int iX = xLow; iX <= xTop; iX++) {
        for (int iY = yLow; iY <= yTop; iY++) {
            for (int iZ = zLow; iZ <= zTop; iZ++) {
                localEnergy += (pow((V.Vx.F(iX-1, iY, iZ) + V.Vx.F(iX, iY, iZ))/2.0, 2.0) +
                                pow((V.Vy.F(iX, iY-1, iZ) + V.Vy.F(iX, iY, iZ))/2.0, 2.0) +
                                pow((V.Vz.F(iX, iY, iZ-1) + V.Vz.F(iX, iY, iZ))/2.0, 2.0))*(mesh.dXi/mesh.xi_xColloc(iX))*(mesh.dEt/mesh.et_yColloc(iY))*(mesh.dZt/mesh.zt_zColloc(iZ));
            }
        }
    }
#endif
    MPI_Allreduce(&localEnergy, &totalEnergy, 1, MPI_FP_REAL, MPI_SUM, MPI_COMM_WORLD);
    totalEnergy /= totalVol;

    if (mesh.rankData.rank == 0) {
        std::cout << std::fixed << std::setw(6)  << std::setprecision(4) << time << "\t" <<
                                    std::setw(16) << std::setprecision(8) << totalEnergy << "\t" << maxDivergence << std::endl;

        ofFile << std::fixed << std::setw(6)  << std::setprecision(4) << time << "\t" <<
                                std::setw(16) << std::setprecision(8) << totalEnergy << "\t" << maxDivergence << "\t" << tStp << std::endl;
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Overloaded function to write the time-series data for scalar solver.
 *
 *          The function computes total energy, Nusselt number and divergence for scalar solver runs.
 *          Only the root rank (rank 0) writes the output.
 *          One line is written to the standard I/O, while another is written to the time series dat file.
 *
 * \param   T is a const reference to the temperature scalar field whose data is used to compute Nusselt number
 * \param   nu is a const reference to the real variable containing the value of kinematic viscosity
 * \param   kappa is a const reference to the real variable containing the value of thermal diffusion
 *
 ********************************************************************************************************************************************
 */
void tseries::writeTSData(const sfield &T, const real nu, const real kappa) {
    // COMPUTE ENERGY AND DIVERGENCE FOR THE INITIAL CONDITION
    V.divergence(divV, P);
    maxDivergence = divV.fxMax();

    if (maxDivergence > 1.0e5) {
        if (mesh.rankData.rank == 0) std::cout << std::endl << "ERROR: Divergence exceeds permissible limits. ABORTING" << std::endl << std::endl;
        MPI_Finalize();
        exit(0);
    }

    localEnergy = 0.0;
    totalEnergy = 0.0;
    localUzT = 0.0;

#ifdef PLANAR
    int iY = 0;
    for (int iX = xLow; iX <= xTop; iX++) {
        for (int iZ = zLow; iZ <= zTop; iZ++) {
            localEnergy += (pow((V.Vx.F(iX-1, iY, iZ) + V.Vx.F(iX, iY, iZ))/2.0, 2.0) +
                            pow((V.Vz.F(iX, iY, iZ-1) + V.Vz.F(iX, iY, iZ))/2.0, 2.0))*(mesh.dXi/mesh.xi_xColloc(iX))*(mesh.dZt/mesh.zt_zColloc(iZ));

            localUzT += ((V.Vz.F(iX, iY, iZ) + V.Vz.F(iX, iY, iZ-1))/2.0)*T.F.F(iX, iY, iZ)*(mesh.dXi/mesh.xi_xColloc(iX))*(mesh.dZt/mesh.zt_zColloc(iZ));
        }
    }
#else
    for (int iX = xLow; iX <= xTop; iX++) {
        for (int iY = yLow; iY <= yTop; iY++) {
            for (int iZ = zLow; iZ <= zTop; iZ++) {
                localEnergy += (pow((V.Vx.F(iX-1, iY, iZ) + V.Vx.F(iX, iY, iZ))/2.0, 2.0) +
                                pow((V.Vy.F(iX, iY-1, iZ) + V.Vy.F(iX, iY, iZ))/2.0, 2.0) +
                                pow((V.Vz.F(iX, iY, iZ-1) + V.Vz.F(iX, iY, iZ))/2.0, 2.0))*(mesh.dXi/mesh.xi_xColloc(iX))*(mesh.dEt/mesh.et_yColloc(iY))*(mesh.dZt/mesh.zt_zColloc(iZ));

                localUzT += ((V.Vz.F(iX, iY, iZ) + V.Vz.F(iX, iY, iZ-1))/2.0)*T.F.F(iX, iY, iZ)*(mesh.dXi/mesh.xi_xColloc(iX))*(mesh.dEt/mesh.et_yColloc(iY))*(mesh.dZt/mesh.zt_zColloc(iZ));
            }
        }
    }
#endif

    MPI_Allreduce(&localEnergy, &totalEnergy, 1, MPI_FP_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&localUzT, &totalUzT, 1, MPI_FP_REAL, MPI_SUM, MPI_COMM_WORLD);
    totalEnergy /= totalVol;
    NusseltNo = 1.0 + (totalUzT/totalVol)/kappa;

    if (mesh.rankData.rank == 0) {
        std::cout << std::fixed << std::setw(6)  << std::setprecision(4) << time << "\t" <<
                                    std::setw(16) << std::setprecision(8) << sqrt(totalEnergy)/nu << "\t" << NusseltNo << "\t" << maxDivergence << std::endl;

        ofFile << std::fixed << std::setw(6)  << std::setprecision(4) << time << "\t" <<
                                std::setw(16) << std::setprecision(8) << sqrt(totalEnergy)/nu << "\t" << NusseltNo << "\t" << maxDivergence << "\t" << tStp << std::endl;
    }
}


tseries::~tseries() {
    ofFile.close();
}
