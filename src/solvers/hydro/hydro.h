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
/*! \file hydro.h
 *
 *  \brief Class declaration of the hydro solver for both 2D and 3D cases.
 *
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#ifndef HYDRO_H
#define HYDRO_H

#include <blitz/array.h>

#include "boundary.h"
#include "parallel.h"
#include "poisson.h"
#include "plainvf.h"
#include "tseries.h"
#include "writer.h"
#include "reader.h"
#include "probes.h"
#include "sfield.h"
#include "vfield.h"
#include "parser.h"
#include "force.h"
#include "grid.h"

class hydro {
    public:
        /** The vector field that stores the velocity field */
        vfield V;

        /** The scalar field that stores the pressure field */
        sfield P;

        /** Instance of force class to handle velocity field forcing */
        force *vForcing;

        hydro(const grid &mesh, const parser &solParam, parallel &mpiParam);

        virtual void solvePDE();
        virtual real testPeriodic();

        virtual ~hydro() { };

    protected:
        /** Integer value for the number of time-steps elapsed - it is incremented by 1 in each time-step. */
        int timeStepCount;

        /** Maximum number of iterations for the iterative solvers \ref hydro#solveVx, \ref hydro#solveVy and \ref hydro#solveVz */
        int maxIterations;

        real time, dt;

        real hx, hy, hz;
        real hx2, hz2, hz2hx2;
        real hx2hy2, hy2hz2, hx2hy2hz2;

        const grid &mesh;
        const parser &inputParams;

        const real inverseRe;

        /** Instance of the \ref probe class to collect data from probes in the domain. */
        probes *dataProbe;

        /** Instances of the \ref boundary class to impose boundary conditions on all the 6 walls for the 3 components of the velocity field. */
        //@{
        boundary *uLft, *uRgt, *uFrn, *uBak, *uTop, *uBot;
        boundary *vLft, *vRgt, *vFrn, *vBak, *vTop, *vBot;
        boundary *wLft, *wRgt, *wFrn, *wBak, *wTop, *wBot;
        //@}

        /** Instance of the \ref parallel class that holds the MPI-related data like rank, xRank, etc. */
        parallel &mpiData;

        /** Plain scalar field into which the pressure correction is calculated and written by the Poisson solver */
        plainsf Pp;
        /** Plain scalar field into which the RHS for pressure Poisson equation is written and passed to the Poisson solver */
        plainsf mgRHS;

        /** Plain vector field into which the RHS of the Navier-Stokes equation is written and stored */
        plainvf nseRHS;
        /** Plain vector field into which the RHS of the implicit equation for velocities are calculated during iterative solving. */
        plainvf velocityLaplacian;
        /** Plain vector field which stores the pressure gradient term. */
        plainvf pressureGradient;
        /** Plain vector field which serves as a temporary array during iterative solution procedure for velocity terms. */
        plainvf guessedVelocity;

        void checkPeriodic();
        void setCoefficients();

        void initVBC();
        void imposeUBCs();
        void imposeVBCs();
        void imposeWBCs();

        void initVForcing();

        virtual void solveVx();
        virtual void solveVy();
        virtual void solveVz();

        virtual void timeAdvance();
};

/**
 ********************************************************************************************************************************************
 *  \class hydro hydro.h "lib/hydro.h"
 *  \brief The base class hydro to solve the incompressible Navier-Stokes equations
 *
 *  The class initializes and stores the velocity vector field and the pressure scalar field along with a few auxilliary
 *  fields to solve the PDE.
 *  It solves the NSE using the \ref solvePDE function from within which the implicit Crank-Nicholson method is used
 *  to solve the PDE.
 ********************************************************************************************************************************************
 */

class hydro_d2: public hydro {
    public:
        hydro_d2(const grid &mesh, const parser &solParam, parallel &mpiParam);

        void solvePDE();
        real testPeriodic();

        ~hydro_d2();

    private:
        multigrid_d2 mgSolver;

        void solveVx();
        void solveVz();

        void timeAdvance();
};

/**
 ********************************************************************************************************************************************
 *  \class hydro_d2 hydro.h "lib/hydro.h"
 *  \brief The derived class from the hydro base class to solve the incompressible NSE in 2D
 *
 *  Certain paramters to be used in the implicit calculation of velocity are defined separately from within the class.
 *  Since the class is instantiated when solving the NSE in 2D, the y-direction component of the grid is supressed.
 *  Consequently, the boundary conditions are imposed only on 4 sides of the domain.
 ********************************************************************************************************************************************
 */

class hydro_d3: public hydro {
    public:
        hydro_d3(const grid &mesh, const parser &solParam, parallel &mpiParam);

        void solvePDE();
        real testPeriodic();

        ~hydro_d3();

    private:
        multigrid_d3 mgSolver;

#ifdef TIME_RUN
        real visc_time, nlin_time, intr_time, impl_time, prhs_time, pois_time;
#endif

        void solveVx();
        void solveVy();
        void solveVz();

        void timeAdvance();
};

/**
 ********************************************************************************************************************************************
 *  \class hydro_d3 hydro.h "lib/hydro.h"
 *  \brief The derived class from the hydro base class to solve the incompressible NSE in 3D
 *
 *  Certain paramters to be used in the implicit calculation of velocity are defined separately from within the class.
 *  Moreover, it imposes boundary conditions on all the three faces of the computational domain.
 ********************************************************************************************************************************************
 */

#endif
