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
/*! \file parser.cc
 *
 *  \brief Definitions for functions of class parser
 *  \sa parser.h
 *  \author Roshan Samuel, Shashwat Bhattacharya
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include <iostream>
#include "parser.h"
#include "mpi.h"

parser::parser() {
    parseYAML();
    checkData();

    setGrids();
    setPeriodicity();

    if (readProbes) {
        parseProbes();

        testProbes();
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to open the yaml file and parse the parameters
 *
 *          The function opens the parameters.yaml file and parses the simulation parameters into its member variables that are publicly
 *          accessible.
 ********************************************************************************************************************************************
 */
void parser::parseYAML() {
    std::ifstream inFile;

    inFile.open("input/parameters.yaml", std::ifstream::in);

    YAML::Node yamlNode;
    YAML::Parser parser(inFile);

    parser.GetNextDocument(yamlNode);

    /********** Problem parameters **********/

    yamlNode["Program"]["Problem Type"] >> probType;
    yamlNode["Program"]["Initial Condition"] >> icType;
    yamlNode["Program"]["Domain Type"] >> domainType;
    yamlNode["Program"]["RBC Type"] >> rbcType;

    yamlNode["Program"]["Reynolds Number"] >> Re;
    yamlNode["Program"]["Rossby Number"] >> Ro;
    yamlNode["Program"]["Rayleigh Number"] >> Ra;
    yamlNode["Program"]["Prandtl Number"] >> Pr;
    yamlNode["Program"]["Taylor Number"] >> Ta;

    yamlNode["Program"]["X Length"] >> Lx;
    yamlNode["Program"]["Y Length"] >> Ly;
    yamlNode["Program"]["Z Length"] >> Lz;

    yamlNode["Program"]["Force"] >> forceType;

    yamlNode["Program"]["Heating Plate"] >> nonHgBC;
    yamlNode["Program"]["Plate Radius"] >> patchRadius;

    /********** Mesh parameters **********/

    yamlNode["Mesh"]["Mesh Type"] >> meshType;

    yamlNode["Mesh"]["X Beta"] >> betaX;
    yamlNode["Mesh"]["Y Beta"] >> betaY;
    yamlNode["Mesh"]["Z Beta"] >> betaZ;

    yamlNode["Mesh"]["X Index"] >> xInd;
    yamlNode["Mesh"]["Y Index"] >> yInd;
    yamlNode["Mesh"]["Z Index"] >> zInd;

    /********** Parallelization parameters **********/

    yamlNode["Parallel"]["Number of OMP threads"] >> nThreads;

    yamlNode["Parallel"]["X Number of Procs"] >> npX;
    yamlNode["Parallel"]["Y Number of Procs"] >> npY;

    /********** Solver parameters **********/

    yamlNode["Solver"]["Differentiation Scheme"] >> dScheme;
    yamlNode["Solver"]["Integration Scheme"] >> iScheme;
    yamlNode["Solver"]["Restart Run"] >> restartFlag;

    yamlNode["Solver"]["Use CFL Condition"] >> useCFL;
    yamlNode["Solver"]["Courant Number"] >> courantNumber;
    yamlNode["Solver"]["Time-Step"] >> tStp;
    yamlNode["Solver"]["Final Time"] >> tMax;

    yamlNode["Solver"]["I/O Count"] >> ioCnt;
    yamlNode["Solver"]["Solution Write Interval"] >> fwInt;
    yamlNode["Solver"]["Restart Write Interval"] >> rsInt;

    yamlNode["Solver"]["Record Probes"] >> readProbes;
    yamlNode["Solver"]["Probe Time Interval"] >> prInt;
    yamlNode["Solver"]["Probes"] >> probeCoords;

    /********** Multigrid parameters **********/

    yamlNode["Multigrid"]["Jacobi Tolerance"] >> tolerance;
    yamlNode["Multigrid"]["V-Cycle Depth"] >> vcDepth;
    yamlNode["Multigrid"]["V-Cycle Count"] >> vcCount;
    yamlNode["Multigrid"]["Pre-Smoothing Count"] >> preSmooth;
    yamlNode["Multigrid"]["Post-Smoothing Count"] >> postSmooth;
    yamlNode["Multigrid"]["Inter-Smoothing Count"] >> interSmooth;

    inFile.close();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to perform a check on the consistency of user-set parameters
 *
 *          In order to catch potential errors early on, a few basic checks are performed here to validate the paramters set
 *          by the user.
 *          Additional checks to be performed on the paramters can be added to this function if necessary.
 ********************************************************************************************************************************************
 */
void parser::checkData() {
    int gridSize, localSize, coarsestSize;

    // CHECK IF THE yInd VARIABLE IS SET CORRECTLY FOR A 2D/3D SIMULATION
#ifdef PLANAR
    if (yInd != 0) {
        std::cout << "WARNING: Y Index parameter of YAML file is non-zero although solver has been compiled with PLANAR flag. Setting Y Index to 0" << std::endl;
        yInd = 0;
    }
#else
    if (yInd == 0) {
        std::cout << "ERROR: Y Index parameter of YAML file is 0 for 3D simulation. ABORTING" << std::endl;
        MPI_Finalize();
        exit(0);
    }
#endif

    // CHECK IF LESS THAN 1 PROCESSOR IS ASKED FOR ALONG X-DIRECTION. IF SO, WARN AND SET IT TO DEFAULT VALUE OF 1
    if (npX < 1) {
        std::cout << "WARNING: Number of processors in X-direction is less than 1. Setting it to 1" << std::endl;
        npX = 1;
    }

    // CHECK IF LESS THAN 1 PROCESSOR IS ASKED FOR ALONG Y-DIRECTION. IF SO, WARN AND SET IT TO DEFAULT VALUE OF 1
    if (npY < 1) {
        std::cout << "WARNING: Number of processors in Y-direction is less than 1. Setting it to 1" << std::endl;
        npY = 1;
    }

    // CHECK IF DOMAIN TYPE STRING IS OF CORRECT LENGTH
    if (domainType.length() != 3) {
        std::cout << "ERROR: Domain type string is not correct. Aborting" << std::endl;
        MPI_Finalize();
        exit(0);
    }

    // CHECK IF MESH TYPE STRING IS OF CORRECT LENGTH
    if (meshType.length() != 3) {
        std::cout << "ERROR: Mesh type string is not correct. Aborting" << std::endl;
        MPI_Finalize();
        exit(0);
    }

    // CHECK IF THE TIME-STEP SET BY USER IS LESS THAN THE MAXIMUM TIME SPECIFIED FOR SIMULATION.
    if (tStp > tMax) {
        std::cout << "ERROR: Time step is larger than the maximum duration assigned for simulation. Aborting" << std::endl;
        MPI_Finalize();
        exit(0);
    }

    // CHECK IF MORE THAN 1 PROCESSOR IS ASKED FOR ALONG Y-DIRECTION FOR A 2D SIMULATION
    if (yInd == 0 and npY > 1) {
        std::cout << "ERROR: More than 1 processor is specified along Y-direction, but the yInd parameter is set to 0. Aborting" << std::endl;
        MPI_Finalize();
        exit(0);
    }

    // CHECK IF THE LENGTH OF ARRAY interSmooth IS LESS THAN vcDepth
    // THE SIZE OF interSmooth IS CONVERTED TO int TO AVOID -Wsign-compare WARNING
    // SIZE OF THIS ARRAY CAN NEVER BE TOO LARGE FOR THIS CONVERSION TO CAUSE ANY PROBLEMS ANYWAY
    if (int(interSmooth.size()) < vcDepth) {
        std::cout << "ERROR: The length of array of inter-smoothing counts is less than V-Cycle depths. Aborting" << std::endl;
        MPI_Finalize();
        exit(0);
    }

    // CHECK IF GRID SIZE SPECIFIED ALONG EACH DIRECTION IS SUFFICIENT ALONG WITH THE DOMAIN DIVIIONS TO REACH THE LOWEST LEVEL OF V-CYCLE DEPTH SPECIFIED
    // ALONG X-DIRECTION
    gridSize = int(pow(2, xInd));
    localSize = gridSize/npX;
    coarsestSize = int(pow(2, vcDepth+1));
    if (localSize < coarsestSize) {
        std::cout << "ERROR: The grid size and domain decomposition along X-direction results in sub-domains too coarse to reach the V-Cycle depth specified. Aborting" << std::endl;
        MPI_Finalize();
        exit(0);
    }

    // ALONG Y-DIRECTION
#ifndef PLANAR
    gridSize = int(pow(2, yInd));
    localSize = gridSize/npY;
    coarsestSize = int(pow(2, vcDepth+1));
    if (yInd > 0 and localSize < coarsestSize) {
        std::cout << "ERROR: The grid size and domain decomposition along Y-direction results in sub-domains too coarse to reach the V-Cycle depth specified. Aborting" << std::endl;
        MPI_Finalize();
        exit(0);
    }
#endif

    // ALONG Z-DIRECTION
    gridSize = int(pow(2, zInd));
    coarsestSize = int(pow(2, vcDepth+1));
    if (gridSize < coarsestSize) {
        std::cout << "ERROR: The grid size along Z-direction is too coarse to reach the V-Cycle depth specified. Aborting" << std::endl;
        MPI_Finalize();
        exit(0);
    }

#ifdef REAL_SINGLE
    if (tolerance < 5.0e-6) {
        std::cout << "ERROR: The specified tolerance for Jacobi iterations is too small for single precision calculations. Aborting" << std::endl;
        MPI_Finalize();
        exit(0);
    }
#endif
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to set the grid types along each direction based on meshType variable
 *
 *          The user specifies mesh type as a single string.
 *          This string has to be parsed to set the integer values xGrid, yGrid and zGrid.
 *          The values of these variables will determine the grid stretching along each direction appropriately.
 ********************************************************************************************************************************************
 */
void parser::setGrids() {
    // The integer values xGrid, yGrid and zGrid are set as below:
    // 0 - uniform grid
    // 1 - single-sided tangent-hyperbolic stretching
    // 2 - double-sided tangent-hyperbolic stretching
    xGrid = 0;
    yGrid = 0;
    zGrid = 0;

    char charMTypes[4];
    std::strcpy(charMTypes, meshType.c_str());

    switch (charMTypes[0]) {
        case 'S': xGrid = 1;
            break;
        case 'D': xGrid = 2;
            break;
    }

    switch (charMTypes[1]) {
        case 'S': yGrid = 1;
            break;
        case 'D': yGrid = 2;
            break;
    }

    switch (charMTypes[2]) {
        case 'S': zGrid = 1;
            break;
        case 'D': zGrid = 2;
            break;
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to set the periodicity of the domain based on domainType variable
 *
 *          The user specifies domain type as a single string.
 *          This string has to be parsed to set the boolean values xPer, yPer and zPer.
 *          The values of these variables will set the periodic boundary conditions appropriately.
 ********************************************************************************************************************************************
 */
void parser::setPeriodicity() {
    xPer = true;
    yPer = true;
    zPer = true;

    if (domainType[0] == 'N') xPer = false;
    if (domainType[1] == 'N') yPer = false;
    if (domainType[2] == 'N') zPer = false;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to parse the probeCoords string
 *
 *          The user specifies probeCoords in NumPy's linspace style with colons
 *          This function extracts the coordinates of the probes from the given string.
 ********************************************************************************************************************************************
 */
void parser::parseProbes() {
    std::string errorProbe;

    while (true) {
        std::vector<std::vector<int> > indexList;

        // Extract the leading set enclosed by square brackets
        std::string oneSet = probeCoords.substr(probeCoords.find('[') + 1, probeCoords.find(']') - 1);
        errorProbe = oneSet;

        oneSet.append(",");
        while (true) {
            std::vector<int> indexVector;
            std::istringstream iss;

            std::string indexData = oneSet.substr(oneSet.find_first_not_of(' '), oneSet.find(',') - oneSet.find_first_not_of(' '));
            indexData.erase(indexData.find_last_not_of(' ') + 1, indexData.length());

            // Erase the extracted set
            oneSet.erase(0, oneSet.find(',') + 1);

            if (indexData.find(":") == std::string::npos) {
                int indexVal;

                // Only a single integer
                iss.str(indexData);
                iss >> indexVal;
                indexVector.push_back(indexVal);
                indexList.push_back(indexVector);

            } else {
                unsigned int strIndex, endIndex, numIndex;

                // Further processing necessary to extract range
                std::string rangeStart = indexData.substr(0, indexData.find(':'));
                indexData.erase(0, indexData.find(':') + 1);
                iss.clear();
                iss.str(rangeStart);
                iss >> strIndex;

                std::string rangeEnd = indexData.substr(0, indexData.find(':'));
                indexData.erase(0, indexData.find(':') + 1);
                iss.clear();
                iss.str(rangeEnd);
                iss >> endIndex;

                iss.clear();
                iss.str(indexData);
                iss >> numIndex;

                if ((strIndex == endIndex) or (numIndex == 1)) {
                    indexVector.push_back(strIndex);

                } else {
                    for (unsigned int i = 0; i < numIndex; i++) {
                        real incIndex = ((real)endIndex - (real)strIndex)/((real)numIndex - 1);
                        int probeIndex = strIndex + (int)round((real)i*incIndex);
                        indexVector.push_back(probeIndex);
                    }
                }

                indexList.push_back(indexVector);
            }

            if (oneSet.length() < 2) {
                break;
            }
        }

        // Erase the extracted set
        probeCoords.erase(0, probeCoords.find(']') + 1);

        // Remove possible white-spaces between the sets
        probeCoords.erase(0, probeCoords.find_first_not_of(' '));

        // Add the coordinates from the extracted set to global coordinates list
#ifdef PLANAR
        if (indexList.size() == 2) {
            for (unsigned int iX = 0; iX < indexList[0].size(); iX++) {
                for (unsigned int iZ = 0; iZ < indexList[1].size(); iZ++) {
                    blitz::TinyVector<int, 3> probeLoc;
                    probeLoc = indexList[0][iX], 0, indexList[1][iZ];
                    probesList.push_back(probeLoc);
                }
            }
        } else if (indexList.size() == 3) {
            for (unsigned int iX = 0; iX < indexList[0].size(); iX++) {
                for (unsigned int iZ = 0; iZ < indexList[2].size(); iZ++) {
                    blitz::TinyVector<int, 3> probeLoc;
                    probeLoc = indexList[0][iX], 0, indexList[2][iZ];
                    probesList.push_back(probeLoc);
                }
            }
        } else {
            std::cout << "WARNING: Number of indices for the probe(s) " << errorProbe << " does not match dimensionality of problem." << std::endl;
        }
#else
        if (indexList.size() == 3) {
            for (unsigned int iX = 0; iX < indexList[0].size(); iX++) {
                for (unsigned int iY = 0; iY < indexList[1].size(); iY++) {
                    for (unsigned int iZ = 0; iZ < indexList[2].size(); iZ++) {
                        blitz::TinyVector<int, 3> probeLoc;
                        probeLoc = indexList[0][iX], indexList[1][iY], indexList[2][iZ];
                        probesList.push_back(probeLoc);
                    }
                }
            }
        } else {
            std::cout << "ERROR: Number of indices for the probe(s) [" << errorProbe << "] does not match dimensionality of problem. ABORTING" << std::endl;
            MPI_Finalize();
            exit(0);
        }
#endif

        // The number 3 is randomly chosen. Ideally if the string is smaller than that, it has no more info
        if (probeCoords.length() < 3) {
            break;
        }
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to test if the probes specified by user are valid
 *
 *          All the probe indices should lie within the domain limits.
 *          This function performs this check to avoid unpleasant surprises later on.
 ********************************************************************************************************************************************
 */
void parser::testProbes() {
    for (unsigned int i = 0; i < probesList.size(); i++) {
        if (probesList[i][0] < 0 or probesList[i][0] > int(pow(2, xInd)) + 1) {
            std::cout << "ERROR: The X index of the probe " << probesList[i] << " lies outside the bounds of the domain. ABORTING" << std::endl;
            MPI_Finalize();
            exit(0);
        }

#ifndef PLANAR
        if (probesList[i][1] < 0 or probesList[i][1] > int(pow(2, yInd)) + 1) {
            std::cout << "ERROR: The Y index of the probe " << probesList[i] << " lies outside the bounds of the domain. ABORTING" << std::endl;
            MPI_Finalize();
            exit(0);
        }
#endif

        if (probesList[i][2] < 0 or probesList[i][2] > int(pow(2, zInd)) + 1) {
            std::cout << "ERROR: The Z index of the probe " << probesList[i] << " lies outside the bounds of the domain. ABORTING" << std::endl;
            MPI_Finalize();
            exit(0);
        }
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to write all the parameter values to I/O
 *
 *          All the user set parameters have to be written to I/O so that the runlog or any out file
 *          contains all the relevant information about the case that was run.
 *          This public function has to be called from the solver by one rank only.
 ********************************************************************************************************************************************
 */
void parser::writeParams() {
    std::cout << std::endl << "Writing all parameters from the YAML input file for reference" << std::endl << std::endl;
    std::cout << "\t****************** START OF parameters.yaml ******************" << std::endl << std::endl;

    std::ifstream inFile;
    inFile.open("input/parameters.yaml", std::ifstream::in);
    std::string line;
    while (std::getline(inFile, line)) {
        std::cout << line << std::endl;
    }
    inFile.close();

    std::cout << std::endl << "\t******************* END OF parameters.yaml *******************" << std::endl;
    std::cout << std::endl;
}
