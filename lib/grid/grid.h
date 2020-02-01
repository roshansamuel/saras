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
/*! \file grid.h
 *
 *  \brief Class declaration of grid
 *
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#ifndef GRID_H
#define GRID_H

#include <math.h>
#include <string>
#include <vector>
#include <blitz/array.h>

#include "parallel.h"

class grid {
    private:
        /** Grid stretching parameter for tangent-hyperbolic function along x, y and z directions */
        blitz::TinyVector<real, 3> thBeta;

        /** Array of the local values of \f$ \xi \f$ within the MPI-decomposed sub-domains in the transformed plane */
        blitz::Array<real, 1> xi;
        /** Array of the local values of \f$ \eta \f$ within the MPI-decomposed sub-domains in the transformed plane */
        blitz::Array<real, 1> et;
        /** Array of the local values of \f$ \zeta \f$ within the MPI-decomposed sub-domains in the transformed plane */
        blitz::Array<real, 1> zt;

        /** Array of the global values of \f$ \xi \f$ within the full domain in the transformed plane */
        blitz::Array<real, 1> xiGlo;
        /** Array of the global values of \f$ \eta \f$ within the full domain in the transformed plane */
        blitz::Array<real, 1> etGlo;
        /** Array of the global values of \f$ \zeta \f$ within the full domain in the transformed plane */
        blitz::Array<real, 1> ztGlo;

        void resizeGrid();
        void makeSizeArray();
        void setDomainSizes();
        void globalXiEtaZeta();

        void createUniformGrid();
        void createTanHypGrid(int dim);

        void checkAnisotropy();
        void gatherGlobal();

        void computeGlobalLimits();

    public:
        /** A const reference to the global variables stored in the parser class to access user set parameters */
        const parser &inputParams;

        /** A const reference to the global variables stored in the parallel class to access MPI related parameters */
        const parallel &rankData;

        /** The sizes of the core of MPI decomposed sub-domains without the pads (collocated points) - localNx, localNy, localNz */
        blitz::TinyVector<int, 3> collocCoreSize;

        /** The sizes of the MPI decomposed sub-domains including the pads on all sides (collocated points) - collocCoreSize + 2*padWidths */
        blitz::TinyVector<int, 3> collocFullSize;

        /** The sizes of the core of MPI decomposed sub-domains without the pads (staggered points) */
        blitz::TinyVector<int, 3> staggrCoreSize;

        /** The sizes of the MPI decomposed sub-domains including the pads on all sides (staggered points) - staggrCoreSize + 2*padWidths */
        blitz::TinyVector<int, 3> staggrFullSize;

        /** The sizes of the pad widths along the three directions - padX, padY, padZ */
        blitz::TinyVector<int, 3> padWidths;

        /** The size of the entire computational domain excluding the pads at the boundary of full domain - globalNx, globalNy, globalNz */
        blitz::TinyVector<int, 3> globalSize;

        /** The end indices of the MPI decomposed sub-domains within the global indexing of the full staggered grid - ztEn, etEn, xiEn */
        blitz::TinyVector<int, 3> subarrayEnds;

        /** The start indices of the MPI decomposed sub-domains within the global indexing of the full staggered grid - ztSt, etSt, xiSt */
        blitz::TinyVector<int, 3> subarrayStarts;

        /** Grid spacing in the transformed plane along the \f$ \xi \f$ direction */
        real dXi;

        /** Grid spacing in the transformed plane along the \f$ \eta \f$ direction */
        real dEt;

        /** Grid spacing in the transformed plane along the \f$ \zeta \f$ direction */
        real dZt;

        /** Length of the physical computational domain along the x direction */
        real xLen;

        /** Length of the physical computational domain along the y direction */
        real yLen;

        /** Length of the physical computational domain along the z direction */
        real zLen;

        /** Array of collocated grid sizes such that the corresponding staggered grid will be multi-grid compatible */
        blitz::Array<int, 1> sizeArray;

        /** Vector of indices pointing to the <B>sizeArray</B> that determines the global full domain size along the 3 directions */
        blitz::TinyVector<int, 3> sizeIndex;

        /** RectDomain object that defines the slice for the core of the local MPI decomposed sub-domain (collocated points) */
        blitz::RectDomain<3> collocCoreDomain;

        /** RectDomain object that defines the slice for the full extent of the local MPI decomposed sub-domain (collocated points) */
        blitz::RectDomain<3> collocFullDomain;

        /** RectDomain object that defines the slice for the core of the local MPI decomposed sub-domain (staggered points) */
        blitz::RectDomain<3> staggrCoreDomain;

        /** RectDomain object that defines the slice for the full extent of the local MPI decomposed sub-domain (staggered points) */
        blitz::RectDomain<3> staggrFullDomain;

        /*****************************************************************************************************************************************************/

        /** Collocated grid along the x-direction defined locally within MPI decomposed sub-domains */
        blitz::Array<real, 1> xColloc;

        /** Collocated grid along the y-direction defined locally within MPI decomposed sub-domains */
        blitz::Array<real, 1> yColloc;

        /** Collocated grid along the z-direction defined locally within MPI decomposed sub-domains */
        blitz::Array<real, 1> zColloc;

        /** Staggered grid along the x-direction defined locally within MPI decomposed sub-domains */
        blitz::Array<real, 1> xStaggr;

        /** Staggered grid along the y-direction defined locally within MPI decomposed sub-domains */
        blitz::Array<real, 1> yStaggr;

        /** Staggered grid along the z-direction defined locally within MPI decomposed sub-domains */
        blitz::Array<real, 1> zStaggr;

        /*****************************************************************************************************************************************************/

        /** Collocated grid along the x-direction for the full global domain */
        blitz::Array<real, 1> xCollocGlobal;

        /** Collocated grid along the y-direction for the full global domain */
        blitz::Array<real, 1> yCollocGlobal;

        /** Collocated grid along the z-direction for the full global domain */
        blitz::Array<real, 1> zCollocGlobal;

        /** Staggered grid along the x-direction for the full global domain */
        blitz::Array<real, 1> xStaggrGlobal;

        /** Staggered grid along the y-direction for the full global domain */
        blitz::Array<real, 1> yStaggrGlobal;

        /** Staggered grid along the z-direction for the full global domain */
        blitz::Array<real, 1> zStaggrGlobal;

        /*****************************************************************************************************************************************************/

        /** Array of the grid derivatives \f$ \frac{\partial\xi}{\partial x} \f$ at collocated grid points, defined locally within each sub-domain. */
        blitz::Array<real, 1> xi_xColloc;

        /** Array of the grid derivatives \f$ \frac{\partial\xi}{\partial x} \f$ at staggered grid points, defined locally within each sub-domain. */
        blitz::Array<real, 1> xi_xStaggr;

        /** Array of the grid derivatives \f$ \frac{\partial^2 \xi}{\partial x^2} \f$ at collocated grid points, defined locally within each sub-domain. */
        blitz::Array<real, 1> xixxColloc;

        /** Array of the grid derivatives \f$ \frac{\partial^2 \xi}{\partial x^2} \f$ at staggered grid points, defined locally within each sub-domain. */
        blitz::Array<real, 1> xixxStaggr;

        /** Array of the grid derivatives \f$ \left(\frac{\partial\xi}{\partial x}\right)^2 \f$ at collocated grid points, defined locally within each sub-domain. */
        blitz::Array<real, 1> xix2Colloc;

        /** Array of the grid derivatives \f$ \left(\frac{\partial\xi}{\partial x}\right)^2 \f$ at staggered grid points, defined locally within each sub-domain. */
        blitz::Array<real, 1> xix2Staggr;

        /*****************************************************************************************************************************************************/

        /** Array of the grid derivatives \f$ \frac{\partial\eta}{\partial y} \f$ at collocated grid points, defined locally within each sub-domain. */
        blitz::Array<real, 1> et_yColloc;

        /** Array of the grid derivatives \f$ \frac{\partial\eta}{\partial y} \f$ at staggered grid points, defined locally within each sub-domain. */
        blitz::Array<real, 1> et_yStaggr;

        /** Array of the grid derivatives \f$ \frac{\partial^2 \eta}{\partial y^2} \f$ at collocated grid points, defined locally within each sub-domain. */
        blitz::Array<real, 1> etyyColloc;

        /** Array of the grid derivatives \f$ \frac{\partial^2 \eta}{\partial y^2} \f$ at staggered grid points, defined locally within each sub-domain. */
        blitz::Array<real, 1> etyyStaggr;

        /** Array of the grid derivatives \f$ \left(\frac{\partial\eta}{\partial y}\right)^2 \f$ at collocated grid points, defined locally within each sub-domain. */
        blitz::Array<real, 1> ety2Colloc;

        /** Array of the grid derivatives \f$ \left(\frac{\partial\eta}{\partial y}\right)^2 \f$ at staggered grid points, defined locally within each sub-domain. */
        blitz::Array<real, 1> ety2Staggr;

        /*****************************************************************************************************************************************************/

        /** Array of the grid derivatives \f$ \frac{\partial\zeta}{\partial z} \f$ at collocated grid points, defined locally within each sub-domain. */
        blitz::Array<real, 1> zt_zColloc;

        /** Array of the grid derivatives \f$ \frac{\partial\zeta}{\partial z} \f$ at staggered grid points, defined locally within each sub-domain. */
        blitz::Array<real, 1> zt_zStaggr;

        /** Array of the grid derivatives \f$ \frac{\partial^2 \zeta}{\partial z^2} \f$ at collocated grid points, defined locally within each sub-domain. */
        blitz::Array<real, 1> ztzzColloc;

        /** Array of the grid derivatives \f$ \frac{\partial^2 \zeta}{\partial z^2} \f$ at staggered grid points, defined locally within each sub-domain. */
        blitz::Array<real, 1> ztzzStaggr;

        /** Array of the grid derivatives \f$ \left(\frac{\partial\zeta}{\partial z}\right)^2 \f$ at collocated grid points, defined locally within each sub-domain. */
        blitz::Array<real, 1> ztz2Colloc;

        /** Array of the grid derivatives \f$ \left(\frac{\partial\zeta}{\partial z}\right)^2 \f$ at staggered grid points, defined locally within each sub-domain. */
        blitz::Array<real, 1> ztz2Staggr;

        grid(const parser &solParam, parallel &parallelData);

        /**
        ********************************************************************************************************************************************
        * \brief   Function to check if a given set of global indices lie within a rank
        *
        *          Based on the data from rankData, the function checks if the given point lies within the subdomain of a processor
        *          The processor in whose sub-domain the point lies returns true.
        *
        * \param   gloIndex is a blitz TinyVector that contains the global indices to check if it lies within an MPI subdomain
        *
        * \return  A boolean value that evaluates to true if the point lies within the sub-domain
        ********************************************************************************************************************************************
        */
        inline bool pointInDomain(blitz::TinyVector<int, 3> gloIndex) const {
            if ((gloIndex(0) < subarrayEnds(0)) and (gloIndex(0) >= subarrayStarts(0)) and (gloIndex(1) < subarrayEnds(1)) and (gloIndex(1) >= subarrayStarts(1))) return true;

            return false;
        };

        /**
        ********************************************************************************************************************************************
        * \brief   Function to obtain local indices from global indices
        *
        *          Based on the data from rankData, the function computes the local indices
        *          in a manner that is consistent for both staggered and collocated indices.
        *
        * \param   gloIndex is a blitz TinyVector that contains the global indices for which local indices have to be found
        *
        * \return  A blitz TinyVector that contains the local indices computed from the given global indices
        ********************************************************************************************************************************************
        */
        inline blitz::TinyVector<int, 3> glo2loc(blitz::TinyVector<int, 3> gloIndex) const {
            blitz::TinyVector<int, 3> locIndex;

            if (pointInDomain(gloIndex)) {
                locIndex(0) = gloIndex(0) % collocCoreSize(0);
                locIndex(1) = gloIndex(1) % collocCoreSize(1);
                locIndex(2) = gloIndex(2);
            } else {
                locIndex = 0, 0, 0;
            }

            return locIndex;
        };

        /**
        ********************************************************************************************************************************************
        * \brief   Function to obtain global indices from local indices
        *
        *          Based on the data from rankData, the function computes the global indices
        *          in a manner that is consistent for both staggered and collocated indices.
        *
        * \param   locIndex is a blitz TinyVector that contains the local indices for which global indices have to be found
        *
        * \return  A blitz TinyVector that contains the global indices computed from the given local indices
        ********************************************************************************************************************************************
        */
        inline blitz::TinyVector<int, 3> loc2glo(blitz::TinyVector<int, 3> locIndex) const {
            blitz::TinyVector<int, 3> gloIndex;

            gloIndex(0) = rankData.xRank*collocCoreSize(0) + locIndex(0);
            gloIndex(1) = rankData.yRank*collocCoreSize(1) + locIndex(1);
            gloIndex(2) = locIndex(2);

            return gloIndex;
        };
};

/**
 ********************************************************************************************************************************************
 *  \class grid grid.h "lib/grid.h"
 *  \brief  Contains all the global variables related to the grid, its slices, limits, and grid derivatives used
 *          throughout the solver
 ********************************************************************************************************************************************
 */

#endif
