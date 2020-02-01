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
/*! \file boundary.h
 *
 *  \brief Class declaration of boundary
 *
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <blitz/array.h>

#include "field.h"
#include "grid.h"

class boundary {
    public:
        boundary(const grid &mesh, field &inField, const int bcWall);

        virtual void imposeBC();

    protected:
        const grid &mesh;

        field &dField;

        bool rankFlag;
        const int wallNum;
        int shiftVal, shiftDim;
};

/**
 ********************************************************************************************************************************************
 *  \class boundary boundary.h "lib/boundary/boundary.h"
 *  \brief Contains all the global variables related to the imposing of boundary conditions, and functions to impose BCs
 *
 ********************************************************************************************************************************************
 */

class dirichletCC: public boundary {
    public:
        dirichletCC(const grid &mesh, field &inField, const int bcWall, const real bcValue);

        inline void imposeBC();
    private:
        const real fieldValue;
};

/**
 ********************************************************************************************************************************************
 *  \class dirichletCC boundary.h "lib/boundary/boundary.h"
 *  \brief The derived class from boundary to apply dirichlet boundary condition for a cell-centered variable.
 *
 ********************************************************************************************************************************************
 */

class dirichletFC: public boundary {
    public:
        dirichletFC(const grid &mesh, field &inField, const int bcWall, const real bcValue);

        inline void imposeBC();
    private:
        const real fieldValue;
};

/**
 ********************************************************************************************************************************************
 *  \class dirichletFC boundary.h "lib/boundary/boundary.h"
 *  \brief The derived class from boundary to apply dirichlet boundary condition for a face-centered variable.
 *
 ********************************************************************************************************************************************
 */

class periodicCC: public boundary {
    public:
        periodicCC(const grid &mesh, field &inField, const int bcWall);

        inline void imposeBC();
};

/**
 ********************************************************************************************************************************************
 *  \class periodicCC boundary.h "lib/boundary/boundary.h"
 *  \brief The derived class from boundary to apply periodic boundary condition for a cell-centered variable.
 *
 ********************************************************************************************************************************************
 */

class periodicFC: public boundary {
    public:
        periodicFC(const grid &mesh, field &inField, const int bcWall);

        inline void imposeBC();
};

/**
 ********************************************************************************************************************************************
 *  \class periodicFC boundary.h "lib/boundary/boundary.h"
 *  \brief The derived class from boundary to apply periodic boundary condition for a face-centered variable.
 *
 ********************************************************************************************************************************************
 */

class neumannCC: public boundary {
    public:
        neumannCC(const grid &mesh, field &inField, const int bcWall, const real bcValue);

        inline void imposeBC();
    private:
        const real fieldValue;
};

/**
 ********************************************************************************************************************************************
 *  \class neumannCC boundary.h "lib/boundary/boundary.h"
 *  \brief The derived class from boundary to apply neumann boundary condition for a cell-centered variable.
 *
 ********************************************************************************************************************************************
 */

class neumannFC: public boundary {
    public:
        neumannFC(const grid &mesh, field &inField, const int bcWall, const real bcValue);

        inline void imposeBC();
    private:
        const real fieldValue;
};

/**
 ********************************************************************************************************************************************
 *  \class neumannFC boundary.h "lib/boundary/boundary.h"
 *  \brief The derived class from boundary to apply neumann boundary condition for a face-centered variable.
 *
 ********************************************************************************************************************************************
 */

class hotPlateCC: public boundary {
    public:
        hotPlateCC(const grid &mesh, field &inField, const int bcWall, const real plateRad);

        void imposeBC();
    private:
        blitz::Array<bool, 3> wallMask;
        blitz::Array<real, 3> wallData;

        blitz::Array<real, 1> x, y, z;
        blitz::Array<real, 1> xGlo, yGlo, zGlo;

        const real patchRadius;

        void setXYZ();

        void createPatch(real patchRadius);
};

/**
 ********************************************************************************************************************************************
 *  \class hotPlateCC boundary.h "lib/boundary/boundary.h"
 *  \brief The derived class from boundary to apply mixed boundary condition involving a heated plate for a cell-centered variable.
 *
 ********************************************************************************************************************************************
 */

#endif
