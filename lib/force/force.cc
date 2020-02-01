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
/*! \file force.cc
 *
 *  \brief Definitions for functions of class force
 *  \sa force.h
 *  \author Shashwat Bhattacharya, Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "force.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the force class
 *
 *          The empty constructor merely initializes the local reference to the global mesh variable and vector field for velocity.
 *          The velocity vector field is used for its interpolation slices to be used in calculating forcing terms.
 *
 * \param   mesh is a const reference to the global data contained in the grid class
 * \param   U is a reference to the velocity vector field
 ********************************************************************************************************************************************
 */
force::force(const grid &mesh, vfield &U): mesh(mesh), V(U) { }


/**
 ********************************************************************************************************************************************
 * \brief   Prototype function to add the forcing field to a plain vector field
 *
 *          Based on the values of Fb, Fr, and other constants as applicable, the appropriate forcing field is calculated and
 *          added to the input plain field.
 *
 * \param   Hv is a reference to the plain vector field to which the forcing term is to be added (RHS of the NSE)
 ********************************************************************************************************************************************
 */
void force::addForcing(plainvf &Hv) { };


/**
 ********************************************************************************************************************************************
 * \brief   Prototype function to add the forcing field to a plain scalar field
 *
 *          Based on the values of Fb, Fr, and other constants as applicable, the appropriate forcing field is calculated and
 *          added to the input plain field.
 *
 * \param   Ht is a reference to the plain scalar field to which the forcing term is to be added (RHS of the temperature equation)
 ********************************************************************************************************************************************
 */
void force::addForcing(plainsf &Ht) { };
