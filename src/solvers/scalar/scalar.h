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
/*! \file scalar.h
 *
 *  \brief Class declaration of the scalar solver for both 2D and 3D cases.
 *
 *  \author Roshan Samuel, Shashwat Bhattacharya
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#ifndef SCALAR_H
#define SCALAR_H

#include <blitz/array.h>

#include "hydro.h"
#include <math.h>

class scalar: public hydro {
    public:
        /** The scalar field that stores the temperature field */
        sfield T;

        /** Instance of force class to handle temperature field forcing */
        force *tForcing;

        real nu, kappa; 

        scalar(const grid &mesh, const parser &solParam, parallel &mpiParam);

        virtual ~scalar() { };

    protected:
        plainsf tmpRHS;
        plainsf guessedScalar;
        plainsf scalarLaplacian;

        boundary *tLft, *tRgt, *tFrn, *tBak, *tTop, *tBot;

        void initTBC();
        void imposeTBCs();

        void initVForcing();
        void initTForcing();

        virtual void solveT();
};

/**
 ********************************************************************************************************************************************
 *  \class scalar scalar.h "lib/scalar.h"
 *  \brief The base class scalar to solve the incompressible Navier-Stokes equations with energy equation
 *
 *  The class initializes and stores the velocity vector field and the pressure scalar field along with a few auxilliary
 *  fields to solve the PDE.
 *  In addition to the fields in hydro class, the scalar equation is also solved here.
 *  It solves the NSE using the \ref solvePDE function from within which the implicit Crank-Nicholson method is used
 *  to solve the PDE.
 ********************************************************************************************************************************************
 */

class scalar_d2: public scalar {
    public:
        scalar_d2(const grid &mesh, const parser &solParam, parallel &mpiParam);

        void solvePDE();
        real testPeriodic();

        ~scalar_d2();

    private:
        multigrid_d2 mgSolver;

        void solveVx();
        void solveVz();

        void solveT();

        void timeAdvance();
};

/**
 ********************************************************************************************************************************************
 *  \class scalar_d2 scalar.h "lib/scalar.h"
 *  \brief The derived class from the scalar base class to solve the incompressible NSE in 2D with energy equation
 *
 *  Certain paramters to be used in the implicit calculation of velocity are defined separately from within the class.
 *  Since the class is instantiated when solving the NSE in 2D, the y-direction component of the grid is supressed.
 *  Consequently, the boundary conditions are imposed only on 4 sides of the domain.
 ********************************************************************************************************************************************
 */

class scalar_d3: public scalar {
    public:
        scalar_d3(const grid &mesh, const parser &solParam, parallel &mpiParam);

        void solvePDE();
        real testPeriodic();

        ~scalar_d3();

    private:
        multigrid_d3 mgSolver;

#ifdef TIME_RUN
        real visc_time, nlin_time, intr_time, impl_time, prhs_time, pois_time;
#endif

        void solveVx();
        void solveVy();
        void solveVz();

        void solveT();

        void timeAdvance();
};

/**
 ********************************************************************************************************************************************
 *  \class scalar_d3 scalar.h "lib/scalar.h"
 *  \brief The derived class from the scalar base class to solve the incompressible NSE in 3D with energy equation
 *
 *  Certain paramters to be used in the implicit calculation of velocity are defined separately from within the class.
 *  Moreover, it imposes boundary conditions on all the three faces of the computational domain.
 ********************************************************************************************************************************************
 */

#endif
