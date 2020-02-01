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
/*! \file mpidata.cc
 *
 *  \brief Definitions for functions of class mpidata
 *  \sa mpidata.h
 *  \author Roshan Samuel
 *  \date Nov 2019
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "mpidata.h"

/**
 ********************************************************************************************************************************************
 * \brief   Constructor of the mpidata class
 *
 *          The short constructor of mpidata class merely resizes the array of MPI_Status and MPI_Request datatypes.
 *          The former is used in non-blocking communication of MPI_Irecv, while the later is used in the MPI_Waitall
 *          function to complete the non-blocking communication call.
 *
 * \param   inputArray is the blitz array whose sub-arrays have to be created and synchronised across processors
 * \param   parallelData is a const reference to the global data contained in the parallel class
 ********************************************************************************************************************************************
 */
mpidata::mpidata(blitz::Array<real, 3> inputArray, const parallel &parallelData): dataField(inputArray), rankData(parallelData) {
    recvStatus.resize(4);
    recvRequest.resize(4);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to create the subarray MPI_Datatypes
 *
 *          Must be called only after the grid class has been initialized.
 *          The subarray data-types cannot be created within the constructor of the parallel class as it needs the grid parameters for
 *          setting the limits of the subarrays.
 *          For this, the grid class will have to be included in the parallel class.
 *
 *          However, the grid object cannot be passed to the parallel class as the grid class already includes the parallel object
 *          within itself, and a reverse include will raise cyclic dependency error.
 *          As a result, the mpidata class offers an additional layer over the parallel class for grid specific data transfer functions.
 *
 * \param   globSize stores the global size of a sub-domain - including core and pads
 * \param   coreSize stores the size of the core of the sub-domain and is similar to the collocCoreSize variable in the grid class
 * \param   padWidth contains the widths of pads along the 3 directions, namely padWidths TinyVector from the grid class
 * \param   xStag specifies whether the array to which the instance of \ref mpidata class is associated with has its data points staggered in x-direction or not
 * \param   yStag specifies whether the array to which the instance of \ref mpidata class is associated with has its data points staggered in y-direction or not
 ********************************************************************************************************************************************
 */
void mpidata::createSubarrays(const blitz::TinyVector<int, 3> globSize,
                              const blitz::TinyVector<int, 3> coreSize,
                              const blitz::TinyVector<int, 3> padWidth,
                              const bool xStag, const bool yStag) {
    /** The <B>loclSize</B> variable holds the local size of the sub-array slice to be sent/received within the sub-domain. */
    blitz::TinyVector<int, 3> loclSize;

    /** The <B>saStarts</B> variable holds the starting coordinates of the sub-array slice. */
    blitz::TinyVector<int, 3> saStarts;

    /** The <B>globCopy</B> variable holds a copy of the global size of the sub-domain. This keeps the original array safe*/
    blitz::TinyVector<int, 3> globCopy;

    globCopy = globSize;

    // CREATING SUBARRAYS FOR TRANSFER ACROSS THE 4 FACES OF EACH SUB-DOMAIN

    /************************************************************** NOTE ************************************************************\
     * MPI Subarrays assume that the starting index of arrays is 0, 0, 0                                                            *
     * But the arrays used here through Blitz start with the index (-padX, -padY, -padZ)                                            *
     * Hence the saStarts variable is shifted accordingly                                                                           *
    \********************************************************************************************************************************/

    //*****************************************************! ALONG XI-DIRECTION !***************************************************//
    // SEND SUB-ARRAY ON LEFT SIDE
    saStarts = padWidth;
    loclSize = coreSize;            loclSize(0) = padWidth(0);

    // STAGGERED GRID SHARE A POINT ACROSS SUB-DOMAIN BOUNDARIES AND HENCE SENDS A SLIGHTLY DIFFERENT DATA-SET
    // HOWEVER, IF THE DOMAIN IS PERIODIC (DEFAULT), THE DATA TO BE SENT ON THE LEFT IS FROM THE WALL POINT ITSELF
    // WHEN USING Method 3 OF setBulkSlice FUNCTION IN field.cc, THE BELOW MODIFICATION APPLIES TO ALL RANKS
    // ELSE, IT APPLIES ONLY TO xRank > 0
    //if (xStag and rankData.xRank > 0) {
    if (xStag) {
        saStarts(0) += padWidth(0);
    }

    MPI_Type_create_subarray(3, globCopy.data(), loclSize.data(), saStarts.data(), MPI_ORDER_C, MPI_FP_REAL, &sendSubarrayX0);
    MPI_Type_commit(&sendSubarrayX0);

    // RECEIVE SUB-ARRAY ON LEFT SIDE
    saStarts = padWidth;            saStarts(0) = 0;
    loclSize = coreSize;            loclSize(0) = padWidth(0);

    MPI_Type_create_subarray(3, globCopy.data(), loclSize.data(), saStarts.data(), MPI_ORDER_C, MPI_FP_REAL, &recvSubarrayX0);
    MPI_Type_commit(&recvSubarrayX0);


    // SEND SUB-ARRAY ON RIGHT SIDE
    saStarts = padWidth;            saStarts(0) = coreSize(0);
    loclSize = coreSize;            loclSize(0) = padWidth(0);

    // STAGGERED GRID SHARE A POINT ACROSS SUB-DOMAIN BOUNDARIES AND HENCE SENDS A SLIGHTLY DIFFERENT DATA-SET
    if (xStag) {
        saStarts(0) -= padWidth(0);
    }

    MPI_Type_create_subarray(3, globCopy.data(), loclSize.data(), saStarts.data(), MPI_ORDER_C, MPI_FP_REAL, &sendSubarrayX1);
    MPI_Type_commit(&sendSubarrayX1);

    // RECEIVE SUB-ARRAY ON RIGHT SIDE
    saStarts = padWidth;            saStarts(0) = coreSize(0) + padWidth(0);
    loclSize = coreSize;            loclSize(0) = padWidth(0);

    // FOR THE LAST RANK, IN PERIODIC BC (DEFAULT BECAUSE SEND NEIGHBOUR IS NULL RANK FOR NON-PERIODIC CASE) THE RECEIVED DATA IS WRITTEN INTO THE WALL POINT
    // BELOW LINE MUST BE COMMENTED WHEN USING Method 3 OF setBulkSlice FUNCTION IN field.cc
    //if (xStag and rankData.xRank == rankData.npX-1) {
    //    saStarts(0) -= padWidth(0);
    //}

    MPI_Type_create_subarray(3, globCopy.data(), loclSize.data(), saStarts.data(), MPI_ORDER_C, MPI_FP_REAL, &recvSubarrayX1);
    MPI_Type_commit(&recvSubarrayX1);


    //****************************************************! ALONG ETA-DIRECTION !***************************************************//
    // SEND SUB-ARRAY ON FRONT SIDE
    saStarts = padWidth;
    loclSize = coreSize;            loclSize(1) = padWidth(1);

    // STAGGERED GRID SHARE A POINT ACROSS SUB-DOMAIN BOUNDARIES AND HENCE SENDS A SLIGHTLY DIFFERENT DATA-SET
    // HOWEVER, FOR THE FIRST RANK, IN PERIODIC BC (DEFAULT) THE DATA FROM WALL POINT ITSELF IS SENT
    // WHEN USING Method 3 OF setBulkSlice FUNCTION IN field.cc, THE BELOW MODIFICATION APPLIES TO ALL RANKS
    // ELSE, IT APPLIES ONLY TO yRank > 0
    //if (yStag and rankData.yRank > 0) {
    if (yStag) {
        saStarts(1) += padWidth(1);
    }

    MPI_Type_create_subarray(3, globCopy.data(), loclSize.data(), saStarts.data(), MPI_ORDER_C, MPI_FP_REAL, &sendSubarrayY0);
    MPI_Type_commit(&sendSubarrayY0);

    // RECEIVE SUB-ARRAY ON FRONT SIDE
    saStarts = padWidth;            saStarts(1) = 0;
    loclSize = coreSize;            loclSize(1) = padWidth(1);

    MPI_Type_create_subarray(3, globCopy.data(), loclSize.data(), saStarts.data(), MPI_ORDER_C, MPI_FP_REAL, &recvSubarrayY0);
    MPI_Type_commit(&recvSubarrayY0);

    // SEND SUB-ARRAY ON REAR SIDE
    saStarts = padWidth;            saStarts(1) = coreSize(1);
    loclSize = coreSize;            loclSize(1) = padWidth(1);

    // STAGGERED GRID SHARE A POINT ACROSS SUB-DOMAIN BOUNDARIES AND HENCE SENDS A SLIGHTLY DIFFERENT DATA-SET
    if (yStag) {
        saStarts(1) -= padWidth(1);
    }

    MPI_Type_create_subarray(3, globCopy.data(), loclSize.data(), saStarts.data(), MPI_ORDER_C, MPI_FP_REAL, &sendSubarrayY1);
    MPI_Type_commit(&sendSubarrayY1);

    // RECEIVE SUB-ARRAY ON REAR SIDE
    saStarts = padWidth;            saStarts(1) = coreSize(1) + padWidth(1);
    loclSize = coreSize;            loclSize(1) = padWidth(1);

    // FOR THE LAST RANK, IN PERIODIC BC (DEFAULT BECAUSE SEND NEIGHBOUR IS NULL RANK FOR NON-PERIODIC CASE) THE RECEIVED DATA IS WRITTEN INTO THE WALL POINT
    // BELOW LINE MUST BE COMMENTED WHEN USING Method 3 OF setBulkSlice FUNCTION IN field.cc
    //if (yStag and rankData.yRank == rankData.npY-1) {
    //    saStarts(1) -= padWidth(1);
    //}

    MPI_Type_create_subarray(3, globCopy.data(), loclSize.data(), saStarts.data(), MPI_ORDER_C, MPI_FP_REAL, &recvSubarrayY1);
    MPI_Type_commit(&recvSubarrayY1);
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to send data across all sub-domain faces
 *
 *          This is the core function of the mpidata class.
 *          The end slices of each sub-domain recieves data from their corresponding neighbouring sub-domains,
 *          while the interior slices of each sub-domain sends data to their corresponding neighbouring sub-domains.
 *
 *          All the data slices are send as subarray MPI derived data-types created in the \ref createSubarrays function.
 *          As a result, \ref syncData must be called only after the subarrays have been created.
 *
 *          The data transfer is implemented here with a mixture of blocking and non-blocking communication calls.
 *          The recieves are non-blocking, while the sends are blocking. This combination prevents inter-processor deadlock.
 ********************************************************************************************************************************************
 */
void mpidata::syncData() {
    recvRequest = MPI_REQUEST_NULL;

    MPI_Irecv(dataField.dataFirst(), 1, recvSubarrayX0, rankData.nearRanks(0), 1, MPI_COMM_WORLD, &recvRequest(0));
    MPI_Irecv(dataField.dataFirst(), 1, recvSubarrayX1, rankData.nearRanks(1), 2, MPI_COMM_WORLD, &recvRequest(1));
    MPI_Irecv(dataField.dataFirst(), 1, recvSubarrayY0, rankData.nearRanks(2), 3, MPI_COMM_WORLD, &recvRequest(2));
    MPI_Irecv(dataField.dataFirst(), 1, recvSubarrayY1, rankData.nearRanks(3), 4, MPI_COMM_WORLD, &recvRequest(3));

    MPI_Send(dataField.dataFirst(), 1, sendSubarrayX0, rankData.nearRanks(0), 2, MPI_COMM_WORLD);
    MPI_Send(dataField.dataFirst(), 1, sendSubarrayX1, rankData.nearRanks(1), 1, MPI_COMM_WORLD);
    MPI_Send(dataField.dataFirst(), 1, sendSubarrayY0, rankData.nearRanks(2), 4, MPI_COMM_WORLD);
    MPI_Send(dataField.dataFirst(), 1, sendSubarrayY1, rankData.nearRanks(3), 3, MPI_COMM_WORLD);

    MPI_Waitall(4, recvRequest.dataFirst(), recvStatus.dataFirst());
}
