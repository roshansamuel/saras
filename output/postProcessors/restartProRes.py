#############################################################################################################################################
 # Saras
 # 
 # Copyright (C) 2019, Mahendra K. Verma
 #
 # All rights reserved.
 # 
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted provided that the following conditions are met:
 #     1. Redistributions of source code must retain the above copyright
 #        notice, this list of conditions and the following disclaimer.
 #     2. Redistributions in binary form must reproduce the above copyright
 #        notice, this list of conditions and the following disclaimer in the
 #        documentation and/or other materials provided with the distribution.
 #     3. Neither the name of the copyright holder nor the
 #        names of its contributors may be used to endorse or promote products
 #        derived from this software without specific prior written permission.
 # 
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 # ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 # WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 # DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 # ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 # (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 # LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 # ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 # (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 # SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 #
 ############################################################################################################################################
 ##
 ##! \file restartProRes.py
 #
 #   \brief Python script interpolate (prolong) or coarsen (reduce) a solution file.
 #
 #   \author Shashwat Bhattacharya and Ali Asad
 #   \date Nov 2019
 #   \copyright New BSD License
 #
 ############################################################################################################################################
 ##

import numpy as np
import h5py

def hdf5_reader(filename,dataset):
    file_V1_read = h5py.File(filename)
    dataset_V1_read = file_V1_read["/"+dataset]
    V1=dataset_V1_read[:,:,:]
    return V1


#Reduces the size of the array to a lower level, 2^(n-1)+1.
def restrict(function):
    global sInd, sLst

    restricted = np.zeros([rx + 1, ry + 1, rz + 1])

    for i in range(2, rx-1):
        for j in range(2, ry-1):
            for k in range(2, rz-1):
                restricted[i, j, k] = function[2*i - 1, 2*j - 1, 2*k - 1]

    return restricted


#Increases the size of the array to a higher level, 2^(n+1)+1.
def prolong(function):
    [lx, ly, lz] = np.shape(function)
    rx, ry, rz = 2*(lx-1), 2*(ly-1), 2*(lz-1)
    prolonged = np.zeros([rx + 1, ry + 1, rz + 1])
    for i in range(0, rx+1, 2):
        for j in range(0, ry+1, 2):
            for k in range(0, rz+1, 2):
                prolonged[i, j, k] = function[i/2, j/2, k/2]
    
    for i in range(1, rx, 2):
        for j in range(0, ry+1, 2):
            for k in range(0, rz+1, 2):
                prolonged[i, j, k] = (prolonged[i-1, j, k] + prolonged[i+1, j, k])/2.0

    for i in range(0, rx+1):
        for j in range(1, ry, 2):
            for k in range(0, rz+1):
                prolonged[i, j, k] = (prolonged[i, j-1, k] + prolonged[i, j+1, k])/2.0

    for i in range(0, rx+1):
        for j in range(0, ry+1):
            for k in range(1, rz, 2):
                prolonged[i, j, k] = (prolonged[i, j, k-1] + prolonged[i, j, k+1])/2.0
                
    return prolonged

def interU2FC(function):
	[lx, ly, lz] = np.shape(function)
	interpolated = np.zeros([lx-1,ly,lz])
	for i in range (lx-1):
		interpolated[i,:,:] = (function[i,:,:] + function[i+1,:,:])/2.0
	return interpolated
	
def interV2FC(function):
	[lx, ly, lz] = np.shape(function)
	interpolated = np.zeros([lx,ly-1,lz])
	for i in range (ly-1):
		interpolated[:,i,:] = (function[:,i,:] + function[:,i+1,:])/2.0
	return interpolated

def interW2FC(function):
	[lx, ly, lz] = np.shape(function)
	interpolated = np.zeros([lx,ly,lz-1])
	for i in range (lz-1):
		interpolated[:,:,i] = (function[:,:,i] + function[:,:,i+1])/2.0
	return interpolated
	
#=======================================================================
#========================Define solution time here======================
time = 500.0
fileNameIN = "Soln_0500.0000.h5" 
fileNameOUT = "SolnIP_0500.0000.h5"

T = hdf5_reader(fileNameIN, "T")
U = hdf5_reader(fileNameIN, "Vx")
V = hdf5_reader(fileNameIN, "Vy")
W = hdf5_reader(fileNameIN, "Vz")
P = hdf5_reader(fileNameIN, "P")
                        
#========================Prolong========================================
#'''
T_p = prolong(T)
U_p = prolong(U)
V_p = prolong(V)
W_p = prolong(W)
P_p = prolong(P)
print ("Done Prolonging from",T.shape[0],"X",T.shape[1],"X",T.shape[2])
print ("To",T_p.shape[0],"X",T_p.shape[1],"X",T_p.shape[2],"\n")
#'''

#========================Restrict========================================
'''
T_r = restrict(T)
U_r = restrict(U)
V_r = restrict(V)
W_r = restrict(W)
P_r = restrict(P)
print ("Done Restricting from",T.shape[0],"X",T.shape[1],"X",T.shape[2])
print ("To",T_r.shape[0],"X",T_r.shape[1],"X",T_r.shape[2])
'''

#========================Interpolate for restart file===================
#'''
U_i = interU2FC(U_p) 
V_i = interV2FC(V_p)
W_i = interW2FC(W_p)
print ("Done Interpolating \n")
#'''

#==================Writing Restart File for SARAS=======================
#'''
hw = h5py.File("restartFile.h5", "w")
dset1 = hw.create_dataset("Vx", data = U_i)
dset2 = hw.create_dataset("Vy", data = V_i)
dset3 = hw.create_dataset("Vz", data = W_i)
dset4 = hw.create_dataset("T", data = T_p)
dset5 = hw.create_dataset("P", data = P_p)
dset6 = hw.create_dataset("Time", data = time)
hw.close()
#'''
print ("Done Writing Restart file for SARAS\n")

#=================Restricted File (all data at cell-center)=============
'''
f1 = h5py.File(fileNameOUT, "w")
dset1 = f1.create_dataset("Vx", data = U_r)
dset2 = f1.create_dataset("Vy", data = V_r)
dset3 = f1.create_dataset("Vz", data = W_r)
dset4 = f1.create_dataset("T", data = T_r)
dset5 = f1.create_dataset("P", data = P_r)
f1.close()
'''

#=================Prolonged File (all data at cell-center)==============
'''
f1 = h5py.File(fileNameOUT, "w")
dset1 = f1.create_dataset("Vx", data = U_p)
dset2 = f1.create_dataset("Vy", data = V_p)
dset3 = f1.create_dataset("Vz", data = W_p)
dset4 = f1.create_dataset("T", data = T_p)
dset5 = f1.create_dataset("P", data = P_p)
f1.close()
'''
