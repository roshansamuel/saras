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
        /** A const reference to the global variables stored in the grid class to access mesh data. */
        const grid &mesh;

        /** Reference to the field onto which the boundary condition has to be applied. */
        field &dField;

        /** The flag is true for MPI ranks on which the boundary condition has to be applied. */
        bool rankFlag;

        /** The const integer denotes the wall at which the boundary condition is being applied. */
        const int wallNum;

        /** Denotes the dimension normal to the wall at which the boundary condition is applied. */
        int shiftDim;

        /** The number of points by which the view of the wall slice is shifted to when applying the boundary condition. */
        int shiftVal;
};

/**
 ********************************************************************************************************************************************
 *  \class boundary boundary.h "lib/boundary/boundary.h"
 *  \brief Contains all the global variables related to the imposing of boundary conditions, and functions to impose BCs
 *
 *         Depending on the location of the variable in the mesh (face-center, cell-center, etc.), the boundary condition (BC) may
 *         either be imposed directly, or through averaging.
 *         Moreover, the nature of the boundary condition (Neumann, Dirichlet, etc.) also affects the manner of imposing the BC.
 *         First of all, the wall on which the BC must be applied is passed to the constructor as an argument.
 *         This information is stored in the \ref wallNum variable, which takes the following values:
 *              - 0 = left wall along X direction
 *              - 1 = right wall along X direction
 *              - 2 = front wall along Y direction
 *              - 3 = rear wall along Y direction
 *              - 4 = bottom wall along Z direction
 *              - 5 = top wall along Z direction
 *
 *         Next the dimension normal to the wall has to be identifed.
 *         This is necessary for shifting the view of the wall slice when applying the BC within \ref imposeBC.
 *         This value is stored by the \ref shiftDim integer variable, which takes the following values:
 *              - 0 = X direction
 *              - 1 = Y direction
 *              - 2 = Z direction
 *
 *         The wall slice defined by the blitz RectDomain objects in the \ref field#fWalls has to be
 *         shifted to the left or right in order to access data at the points adjacent to the wall.
 *         This is done when averaging/interpolating the values at the wall.
 *         The integer value by which this shifting has to be done is specified by \ref shiftVal, which could be either +1 or -1.
 *         
 *         Finally not all subdomains need to have the BC applied to them.
 *         The boolean variable \ref rankFlag is set appropriately in the constructor so that \ref imposeBC updates only
 *         those MPI ranks whose boundaries need the BC applied to them.
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

class nullBC: public boundary {
    public:
        nullBC(const grid &mesh, field &inField, const int bcWall): boundary(mesh, inField, bcWall) { };

        inline void imposeBC() { };
};

/**
 ********************************************************************************************************************************************
 *  \class nullBC boundary.h "lib/boundary/boundary.h"
 *  \brief The derived class from boundary to impose null boundary condition that leaves the data unchanged.
 *
 ********************************************************************************************************************************************
 */

#endif
