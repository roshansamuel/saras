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
 ##! \file inter_Stag.py
 #
 #   \brief Python script to interpolate all the fields to staggered grid-points.
 #
 #   \author Shashwat Bhattacharya
 #   \date Nov 2019
 #   \copyright New BSD License
 #
 ############################################################################################################################################
 ##

#This code interpolates all the fields 

import h5py
import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt

def hdf5_reader(filename,dataset):
    file_V1_read = h5py.File(filename)
    dataset_V1_read = file_V1_read["/"+dataset]
    V1=dataset_V1_read[:,:,:]
    return V1

Ra = 6000.0
Pr=1.0

#kappa = np.sqrt(1.0/(Ra*Pr))
#nu = np.sqrt(Pr/Ra)

kappa = 1.0
nu = Pr


lx, ly, lz = 1.0, 1.0, 1.0

fileName = "Soln_0300.0000.h5"

U = hdf5_reader(fileName, "Vx")
V = hdf5_reader(fileName, "Vy")
W = hdf5_reader(fileName, "Vz")

T = hdf5_reader(fileName, "T")
P = hdf5_reader(fileName, "P")

U_st = np.zeros_like(P)
V_st = np.zeros_like(P)
W_st = np.zeros_like(P)

Nx_s = len(P[:,0,0])
Nx_c = len(U[:,0,0])
Ny_s = len(P[0,:,0])
Ny_c = len(V[0,:,0])
Nz_s = len(P[0,0,:])
Nz_c = len(W[0,0,:])

#print Nx_s, Nx_c

x = np.linspace(0,1,Nx_s)
y = np.linspace(0,1,Ny_s)
z = np.linspace(0,1,Nz_s)

U_st[1:Nx_c,:,:] = 0.5*(U[0:Nx_c-1,:,:]+U[1:Nx_c,:,:])
U_st[0,:,:] = 0.0
U_st[-1,:,:] = 0.0

V_st[:,1:Ny_c,:] = 0.5*(V[:,0:Ny_c-1,:]+V[:,1:Ny_c,:])
V_st[:,0,:] = 0.0
V_st[:,-1,:] = 0.0

W_st[:,:,1:Nz_c] = 0.5*(W[:,:,0:Nz_c-1]+W[:,:,1:Nz_c])
W_st[:,:,0] = 0.0
W_st[:,:,-1] = 0.0


f1 = h5py.File("SolnI_0300.0000.h5", "w")
dset1 = f1.create_dataset("Vx", data = U_st)
dset2 = f1.create_dataset("Vy", data = V_st)
dset3 = f1.create_dataset("Vz", data = W_st)
dset4 = f1.create_dataset("T", data = T)
dset5 = f1.create_dataset("P", data = P)
f1.close()
