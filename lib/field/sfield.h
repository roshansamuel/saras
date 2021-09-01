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
/*! \file sfield.h
 *
 *  \brief Class declaration of sfield - scalar field
 *
 *  \author Roshan Samuel, Ali Asad
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#ifndef SFIELD_H
#define SFIELD_H

#include "field.h"
#include "derivative.h"

// Forward declarations of relevant classes
class plainsf;
class plainvf;
class vfield;

class sfield {
    private:
        const grid &gridData;

        blitz::Array<real, 3> derivTempF;

    public:
        field F;

        /** derS is an instance of the derivative class used to compute derivatives */
        derivative derS;

        /** This string is used to identify the vector field, and is useful in file-writing */
        std::string fieldName;

        blitz::Array<real, 3> interTempF;

        sfield(const grid &gridData, std::string fieldName);

        void computeDiff(plainsf &H);
        void computeNLin(const vfield &V, plainsf &H);

        void gradient(plainvf &gradF, const vfield &V);

        void syncData();

        sfield& operator += (plainsf &a);
        sfield& operator -= (plainsf &a);

        sfield& operator += (sfield &a);
        sfield& operator -= (sfield &a);

        sfield& operator *= (real a);

        void operator = (plainsf &a);
        void operator = (sfield &a);
        void operator = (real a);

        ~sfield() { };
};

/**
 ********************************************************************************************************************************************
 *  \class sfield sfield.h "lib/sfield.h"
 *  \brief Scalar field class to store and operate on scalar fields
 *
 *  The class stores scalar fields in the form of an instance of the field class defined in <B>field.h</B>.
 *  While the class <B>field</B> merely stores data in the form of a blitz array and offers functions to compute derivatives
 *  over a uniform grid, the <B>sfield</B> class adds another layer of functionality along with the <B>grid</B> (<B>grid.h</B>)
 *  class to apply necessary grid transformation metrics and compute derivatives over a non-uniform grid.
 *  The scalar field is also equipped with a gradient operator, which returns a vector field (vfield).
 *  However, this operation is presently restricted to cell-centered scalar fields, i.e., those which are staggered in all
 *  the directions.
 *  Moreover, the \f$ (\mathbf{u}.\nabla)f \f$ operator is also provided as the function <B>computeNLin</B>
 ********************************************************************************************************************************************
 */

#endif
