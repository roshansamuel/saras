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
/*! \file probes.h
 *
 *  \brief Class declaration of probes
 *
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#ifndef PROBES_H
#define PROBES_H

#include <vector>
#include <fstream>
#include <cstddef>
#include <blitz/array.h>
#include <mpi.h>

#include "field.h"
#include "grid.h"

/**
 ********************************************************************************************************************************************
 *  \struct dataStruct
 *  \brief The data obtained from the probes is stored in struct for quick transfer across processes.
 *
 *  The struct can store data from up to a maximum of 10 field variables.
 *  It also stores the global x, y and z indices of the probe location corresponding to each probe.
 ********************************************************************************************************************************************
 */
typedef struct dataStruct {
    /** Integer values of the global indices of the probe. */
    //@{
    int x, y, z;
    //@}

    /** Array of double/single precision numbers that contains the probed data from up to 10 field variables. */
    real probeData[10];

    /** Default constructor for the struct that initializes all values to 0. */
    dataStruct(): x(0), y(0), z(0) { for (int i=0; i<10; i++) probeData[i] = 0.0; }
} dataStruct;

class probes {
    public:
        probes(const grid &mesh, std::vector<field> &pFields);

        void probeData(real time);

        ~probes();

    private:
        const grid &mesh;

        std::vector<field> &pFields;

        const unsigned int numFields;

        std::vector<blitz::TinyVector<int, 3> > globalProbes, localProbes;

        std::ofstream probeFile;

        MPI_Datatype mpiStruct;

        void getData(dataStruct *outData);
        void gatherData(dataStruct *outData);

        void createMPIStruct();

        void placeProbes();
};
/**
 ********************************************************************************************************************************************
 *  \class probes probes.h "lib/io/probes.h"
 *  \brief Handles the writing of data from probes placed in the domain
 *
 *  The class places the probes in probesList provided by user through the parser class.
 *  It also provides an interface to the solver to read data from the probes.
 ********************************************************************************************************************************************
 */

#endif
