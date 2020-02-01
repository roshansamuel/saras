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
 ##! \file cutter.py
 #
 #   \brief Python script to interchange indices and extract planes for 2D3C runs.
 #
 #   \author Shashwat Bhattacharya
 #   \date Nov 2019
 #   \copyright New BSD License
 #
 ############################################################################################################################################
 ##

import h5py
import numpy as np

def hdf5_reader(filename,dataset):
    file_V1_read = h5py.File(filename)
    dataset_V1_read = file_V1_read["/"+dataset]
    V1=dataset_V1_read[:,:,:]
    return V1

fileName = "Soln_0005.2000.h5"
d2cut = 0


U = hdf5_reader(fileName, "Vx")
V = hdf5_reader(fileName, "Vy")
W = hdf5_reader(fileName, "Vz")
T = hdf5_reader(fileName, "T")

#cutting U
[Nx, Ny, Nz] = U.shape
U_c = (U[:, 0:Ny-1, 0:Nz-1] + U[:, 0:Ny-1, 1:Nz] + U[:, 1:Ny, 0:Nz-1] + U[:, 1:Ny, 1:Nz])/4.0

#cutting V
[Nx, Ny, Nz] = V.shape
V_c = (V[0:Nx-1, :, 0:Nz-1] + V[0:Nx-1, :, 1:Nz] + V[1:Nx, :, 0:Nz-1] + V[1:Nx, :, 1:Nz])/4.0

#cutting W
[Nx, Ny, Nz] = W.shape
W_c = (W[0:Nx-1, 0:Ny-1, :] + W[0:Nx-1, 1:Ny, :] + W[1:Nx, 0:Ny-1, :] + W[1:Nx, 1:Ny, :])/4.0

#cutting T
[Nx, Ny, Nz] = T.shape
T_c = (T[0:Nx-1, 0:Ny-1, 0:Nz-1] + T[0:Nx-1, 0:Ny-1, 1:Nz] + 
        T[0:Nx-1, 1:Ny, 0:Nz-1] + T[0:Nx-1, 1:Ny, 1:Nz] +
        T[1:Nx, 0:Ny-1, 0:Nz-1] + T[1:Nx, 0:Ny-1, 1:Nz] + 
        T[1:Nx, 1:Ny, 0:Nz-1] + T[1:Nx, 1:Ny, 1:Nz])/8.0


if d2cut==1:
    [Nx, Ny, Nz] = T_c.shape
    U_t = np.zeros([Nx,1,Nz])
    V_t = np.zeros([Nx,1,Nz])
    W_t = np.zeros([Nx,1,Nz])
    T_t = np.zeros([Nx,1,Nz])
    
    U_t[:,0,:] = U_c[:,:,3]
    V_t[:,0,:] = W_c[:,:,3]
    W_t[:,0,:] = V_c[:,:,3]
    T_t[:,0,:] = T_c[:,:,3]

else:
    [Nx, Ny, Nz] = T_c.shape
    
    U_t = np.zeros_like(U_c)
    V_t = np.zeros_like(W_c)
    W_t = np.zeros_like(V_c)
    T_t = np.zeros_like(T_c)
    
    for i in range(Nx):
        print i+1, "out of", Nx
        for j in range(Ny):
            for k in range(Nz):
                U_t[i,j,k] = U_c[i,k,j]
                V_t[i,j,k] = W_c[i,k,j]
                W_t[i,j,k] = V_c[i,k,j]
                T_t[i,j,k] = T_c[i,k,j]
    

f1 = h5py.File("U.V1r.h5", "w")
dset1 = f1.create_dataset("U.V1r", data = U_t)
f1.close()

f2 = h5py.File("U.V2r.h5", "w")
dset2 = f2.create_dataset("U.V2r", data = V_t)
f2.close()

f3 = h5py.File("U.V3r.h5", "w")
dset3 = f3.create_dataset("U.V3r", data = W_t)
f3.close()

f4 = h5py.File("T.Fr.h5", "w")
dset4 = f4.create_dataset("T.Fr", data = T_t)
f4.close()
