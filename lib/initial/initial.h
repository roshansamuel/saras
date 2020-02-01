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
/*! \file initial.h
 *
 *  \brief Class declaration of initial
 *
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#ifndef INITIAL_H
#define INITIAL_H

#include <blitz/array.h>

#include "vfield.h"
#include "grid.h"

class initial {
    public:
        initial(const grid &mesh);

        virtual void initializeField(vfield &uField);

    protected:
        const grid &mesh;
};

/**
 ********************************************************************************************************************************************
 *  \class initial initial.h "lib/initial/initial.h"
 *  \brief Contains all the global variables related to the imposing of initial conditions, and functions to impose them
 *
 ********************************************************************************************************************************************
 */

class taylorGreen: public initial {
    public:
        taylorGreen(const grid &mesh);

        void initializeField(vfield &uField);
};

/**
 ********************************************************************************************************************************************
 *  \class taylorGreen initial.h "lib/initial/initial.h"
 *  \brief The derived class from initial to impose initial condition of Taylor-Green vortices.
 *
 ********************************************************************************************************************************************
 */

class channelSine: public initial {
    public:
        channelSine(const grid &mesh);

        void initializeField(vfield &uField);
};

/**
 ********************************************************************************************************************************************
 *  \class channelSine initial.h "lib/initial/initial.h"
 *  \brief The derived class from initial to impose sinusoidal perturbation for channel flow.
 *
 ********************************************************************************************************************************************
 */

class channelRand: public initial {
    public:
        channelRand(const grid &mesh);

        void initializeField(vfield &uField);
};

/**
 ********************************************************************************************************************************************
 *  \class channelRand initial.h "lib/initial/initial.h"
 *  \brief The derived class from initial to impose random initial condition for channel flow.
 *
 ********************************************************************************************************************************************
 */

class zeroInitial: public initial {
    public:
        zeroInitial(const grid &mesh): initial(mesh) { };

        void initializeField(vfield &uField) {uField.Vx = 0.0; uField.Vy = 0.0; uField.Vz = 0.0;};
};

/**
 ********************************************************************************************************************************************
 *  \class zeroInitial initial.h "lib/initial/initial.h"
 *  \brief The derived class from initial to impose the default condition of 0 velocity.
 *
 ********************************************************************************************************************************************
 */

#endif
