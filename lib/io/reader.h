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
/*! \file reader.h
 *
 *  \brief Class declaration of reader
 *
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#ifndef READER_H
#define READER_H

#include <blitz/array.h>
#include <sys/stat.h>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <dirent.h>

#include "field.h"
#include "grid.h"
#include "hdf5.h"

class reader {
    private:
        const grid &mesh;

        std::vector<field> &rFields;

#ifdef PLANAR
        blitz::Array<real, 2> fieldData;
#else
        blitz::Array<real, 3> fieldData;
#endif

        std::vector<hid_t> sourceDSpace, targetDSpace;

        std::vector< blitz::TinyVector<int, 3> > localSize;

        void initLimits();

        void copyData(field &outField);

    public:
        reader(const grid &mesh, std::vector<field> &rFields);

        real readData();

        ~reader();
};

/**
 ********************************************************************************************************************************************
 *  \class reader reader.h "lib/reader.h"
 *  \brief Class for all the global variables and functions related to reading input data for the solver.
 *
 *  The computational data for the solver can be read from HDF5 file.
 *  The class allows for both collocated and staggered grid data to be read from separate input files.
 ********************************************************************************************************************************************
 */

#endif
