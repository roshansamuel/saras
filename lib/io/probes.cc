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
/*! \file probes.cc
 *
 *  \brief Definitions for functions of class probes
 *  \sa probes.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include <iostream>
#include "probes.h"
#include "mpi.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the probes class
 *
 *          The constructor initializes the variables and parameters for writing probed data
 *          The header for the ProbeData.dat file is also written here.
 *          Depending on whether the run is 2D or 3D, the appropriate header is written.
 *          Finally, the MPIStruct datatype for transferring probed data across processors is also initialized.
 *
 * \param   mesh is a const reference to the global data contained in the grid class.
 * \param   pFields is a vector containing a list of references to all the fields which need to be probed.
 *
 ********************************************************************************************************************************************
 */
probes::probes(const grid &mesh, std::vector<field> &pFields): mesh(mesh), pFields(pFields), numFields(pFields.size()) {
    placeProbes();

    // Open ProbesData file in append mode
    probeFile.open("output/ProbesData.dat", std::fstream::out | std::fstream::app);

    // If the solver is being started on a fresh run, add a header to the file
    if (not mesh.inputParams.restartFlag) {
        if (mesh.rankData.rank == 0) {
            std::ostringstream fileHeader;
            fileHeader.str(std::string());

#ifdef PLANAR
            fileHeader << "#VARIABLES = Time, X, Z, " << pFields[0].fieldName;
#else
            fileHeader << "#VARIABLES = Time, X, Y, Z, " << pFields[0].fieldName;
#endif
            for (unsigned int i = 1; i < numFields; i++) {
                fileHeader << ", " << pFields[i].fieldName;
            }
            probeFile << fileHeader.str() << std::endl;
        }
    }

    createMPIStruct();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to place all the probes in the MPI decomposed sub-domains
 *
 *          The indices of probes provided in probesList of parser are global indices.
 *          This function gives the local indices of the probes within the sub-domains after MPI decomposition.
 *          Note that the probes are at the staggered grid points - cell centers.
 ********************************************************************************************************************************************
 */
void probes::placeProbes() {
    blitz::TinyVector<int, 3> localIndices;

    for (unsigned int i = 0; i < mesh.inputParams.probesList.size(); i++) {
        localIndices = mesh.inputParams.probesList[i];
        localIndices(0) -= mesh.subarrayStarts(0);
        localIndices(1) -= mesh.subarrayStarts(1);

        if (mesh.inputParams.probesList[i](0) >= mesh.subarrayStarts(0) and mesh.inputParams.probesList[i](0) < mesh.subarrayEnds(0)) {
#ifndef PLANAR
            if (mesh.inputParams.probesList[i](1) >= mesh.subarrayStarts(1) and mesh.inputParams.probesList[i](1) < mesh.subarrayEnds(1)) {

                globalProbes.push_back(mesh.inputParams.probesList[i]);
                localProbes.push_back(localIndices);
            }
#else
            globalProbes.push_back(mesh.inputParams.probesList[i]);
            localProbes.push_back(localIndices);
#endif
        }
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to read the variables at the probe locations and return them
 *
 *          For each variable specified in the list of fields given to the constructor, obtain the value at the probe locations.
 ********************************************************************************************************************************************
 */
void probes::probeData(real time) {
    dataStruct *writeData;
    int arrayLen;

    if (mesh.rankData.rank == 0) {
        arrayLen = mesh.inputParams.probesList.size();
        writeData = new dataStruct[arrayLen];
    } else {
        arrayLen = localProbes.size();
        writeData = new dataStruct[arrayLen];
    }

    getData(writeData);

    gatherData(writeData);

    // Write data collected from all probes
    if (mesh.rankData.rank == 0) {
        for (int i = 0; i < arrayLen; i++) {
            probeFile << std::fixed << std::setw(6) << std::setprecision(4) << time << "\t";
#ifdef PLANAR
            probeFile << std::setw(4) << std::setprecision(0)
                      << writeData[i].x << "\t"
                      << writeData[i].z << "\t";
#else
            probeFile << std::setw(4) << std::setprecision(0)
                      << writeData[i].x << "\t"
                      << writeData[i].y << "\t"
                      << writeData[i].z << "\t";
#endif
            probeFile << std::setw(16) << std::setprecision(8)
                      << writeData[i].probeData[0];
            for (unsigned int j = 1; j < numFields; j++) {
                probeFile << "\t" << std::setw(16) << std::setprecision(8)
                          << writeData[i].probeData[j];
            }

            probeFile << std::endl;
        }
    }

    delete writeData;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to read the variables locally within each sub-domain
 *
 *          For each variable specified in the list of fields given to the constructor, obtain the value at the probe locations.
 *          Store the data in the reference pointer supplied to the function.
 ********************************************************************************************************************************************
 */
void probes::getData(dataStruct *outData) {
    for (unsigned int i = 0; i < localProbes.size(); i++) {
        outData[i].x = globalProbes[i](0);
#ifndef PLANAR
        outData[i].y = globalProbes[i](1);
#endif
        outData[i].z = globalProbes[i](2);

        for (unsigned int j = 0; j < numFields; j++) {
            outData[i].probeData[j] = pFields[j].F(localProbes[i]);
        }
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to gather the probed data to rank 0
 *
 *          For each variable specified in the list of fields given to the constructor, obtain the value at the probe locations.
 *          It returns a vector of strcts containing the variable values
 ********************************************************************************************************************************************
 */
void probes::gatherData(dataStruct *outData) {
    unsigned int localSize = localProbes.size();
    int *numProbes, *prbStarts;

    numProbes = new int[mesh.rankData.nProc];
    prbStarts = new int[mesh.rankData.nProc];

    dataStruct sendData[localSize];

    for (unsigned int i = 0; i < localSize; i++) {
        sendData[i].x = outData[i].x;
#ifndef PLANAR
        sendData[i].y = outData[i].y;
#endif
        sendData[i].z = outData[i].z;

        for (int j=0; j<10; j++) {
            sendData[i].probeData[j] = outData[i].probeData[j];
        }
    }

    MPI_Gather(&localSize, 1, MPI_INT, numProbes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&localSize, 1, MPI_INT, prbStarts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (mesh.rankData.rank == 0) {
        for (int i = mesh.rankData.nProc - 1; i > 0; --i) {
            prbStarts[i] = prbStarts[i-1];
        }
        prbStarts[0] = 0;
        for (int i = 1; i < mesh.rankData.nProc; ++i) {
            prbStarts[i] += prbStarts[i-1];
        }
    }

    // Gatherv instead of gather allows for non-uniformly distributed data from processes.
    MPI_Gatherv(&sendData, localSize, mpiStruct, outData, numProbes, prbStarts, mpiStruct, 0, MPI_COMM_WORLD);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to create the derived MPI datatype to be used for transferring probed data across processors in 3D runs
 *
 *          The probed data consists of coordinates (integers) and the field readings (real).
 *          An MPI_Struct datatype is defined accordingly to store and transmit this data across processors.
 ********************************************************************************************************************************************
 */
void probes::createMPIStruct() {
    dataStruct a;

    // Below method of defining MPI_struct is compatible with MPICH2 since it avoids use of deprecated functions and definitions
    int bLengths[4] = {1, 1, 1, 10};
    MPI_Datatype baseTypes[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_FP_REAL};
    MPI_Aint displacements[4];
    MPI_Aint intlb, intex;

    MPI_Type_get_extent(MPI_INT, &intlb, &intex);

    displacements[0] = (MPI_Aint) 0;

    displacements[1] = intex;
    displacements[2] = intex + intex;
    displacements[3] = intex + intex + intex + intex;   // Another intex has to be added due to padding of structs in C++. This could be compiler dependent

    MPI_Type_create_struct(4, bLengths, displacements, baseTypes, &mpiStruct);
    MPI_Type_commit(&mpiStruct);
}

probes::~probes() {
    probeFile.close();
}
