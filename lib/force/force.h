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
/*! \file force.h
 *
 *  \brief Class declaration of force
 *
 *  \author Shashwat Bhattacharya, Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#ifndef FORCE_H
#define FORCE_H

#include <blitz/array.h>

#include "plainvf.h"
#include "sfield.h"
#include "vfield.h"

class force {
    public:
        force(const grid &mesh, vfield &U);

        virtual void addForcing(plainvf &Hv);
        virtual void addForcing(plainsf &Ht);

    protected:
        const grid &mesh;

        vfield &V;
};

/**
 ********************************************************************************************************************************************
 *  \class force force.h "lib/force/force.h"
 *  \brief Contains all the global variables related to the imposing of forcing, and associated functions
 *
 ********************************************************************************************************************************************
 */

class coriolisForce: public force {
    public:
        coriolisForce(const grid &mesh, vfield &U);

        inline void addForcing(plainvf &Hv);
        inline void addForcing(plainsf &Ht) { };
    private:
        real Fr;
};

/**
 ********************************************************************************************************************************************
 *  \class coriolisForce force.h "lib/force/force.h"
 *  \brief The derived class from force to add Coriolis forcing to the velocity field in rotating systems.
 *
 ********************************************************************************************************************************************
 */

class buoyantForce: public force {
    public:
        buoyantForce(const grid &mesh, vfield &U, const sfield &T);

        inline void addForcing(plainvf &Hv);
        inline void addForcing(plainsf &Ht) { };
    private:
        real Fb;

        const sfield &T;
};

/**
 ********************************************************************************************************************************************
 *  \class buoyantForce force.h "lib/force/force.h"
 *  \brief The derived class from force to add forcing due to buoyancy to the velocity field in convecting systems.
 *
 ********************************************************************************************************************************************
 */

class rotatingConv: public force {
    public:
        rotatingConv(const grid &mesh, vfield &U, const sfield &T);

        inline void addForcing(plainvf &Hv);
        inline void addForcing(plainsf &Ht) { };
    private:
        real Fb, Fr;

        const sfield &T;
};

/**
 ********************************************************************************************************************************************
 *  \class rotatingConv force.h "lib/force/force.h"
 *  \brief The derived class from force to add forcing due to both buoyancy and rotation to the velocity field in rotating convection simulations.
 *
 ********************************************************************************************************************************************
 */

class randomForcing: public force {
    public:
        randomForcing(const grid &mesh, vfield &U);

        inline void addForcing(plainvf &Hv);
        inline void addForcing(plainsf &Ht) { };
    private:
        blitz::Array<real, 3> Force_x, Force_y, Force_z;
};

/**
 ********************************************************************************************************************************************
 *  \class randomForcing force.h "lib/force/force.h"
 *  \brief The derived class from force to add random forcing to the velocity field.
 *
 ********************************************************************************************************************************************
 */

class constantPGrad: public force {
    public:
        constantPGrad(const grid &mesh, vfield &U): force(mesh, U) { };

        inline void addForcing(plainvf &Hv) {Hv.Vx += 1.0;};
        inline void addForcing(plainsf &Ht) { };
};

/**
 ********************************************************************************************************************************************
 *  \class constantPGrad force.h "lib/force/force.h"
 *  \brief The derived class from force to add forcing due to constant pressure gradient to the velocity field, specially in channel flow simulations.
 *
 ********************************************************************************************************************************************
 */

class zeroForcing: public force {
    public:
        zeroForcing(const grid &mesh, vfield &U): force(mesh, U) { };

        inline void addForcing(plainvf &Hv) { };
        inline void addForcing(plainsf &Ht) { };
};

/**
 ********************************************************************************************************************************************
 *  \class zeroForcing force.h "lib/force/force.h"
 *  \brief The derived class from force to add the default forcing of no forcing.
 *
 ********************************************************************************************************************************************
 */

#endif
