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
 ##! \file time_series_plot_hydro.py
 #
 #   \brief Python script to plot time series output from the hydro solver.
 #
 #   \author Ali Asad
 #   \date Nov 2019
 #   \copyright New BSD License
 #
 ############################################################################################################################################
 ##

import numpy as np
import pylab as plt
from numpy import*
import matplotlib .pyplot as plt
import time

plt.ion()
plt.rcParams['xtick.major.size'] = 9

plt.rcParams['xtick.major.width'] = 1

plt.rcParams['xtick.minor.size'] = 5

plt.rcParams['xtick.minor.width'] = 1

plt.rcParams['ytick.major.size'] = 9

plt.rcParams['ytick.major.width'] = 1

plt.rcParams['ytick.minor.size'] = 5

plt.rcParams['ytick.minor.width'] = 1

font = {'family' : 'serif', 'weight' : 'normal', 'size' : 15}
plt.rc('font', **font)

#===============FIRST PLOT++++++++++++++++++++++++++++++++++++++++++++++

#Load Data==========================================================
data = np.loadtxt('../TimeSeries.dat',comments='#')
t = data[:,0]
Re = data[:,1]
Div = data[:,2]
dt = data[:,3]

#Define Figure and subplots=========================================
fig = plt.figure(1)
f_ax1 = fig.add_subplot(2,2,1)
f_ax2 = fig.add_subplot(2,2,2)
f_ax3 = fig.add_subplot(2,2,3)
f_ax4 = fig.add_subplot(2,2,4)
fig.suptitle(r"$Time\:Series\:Plots$" ,fontsize=30)


#Sub-Plot 1=============================================================
f_ax1.plot(t, dt, "b", lw=2.0)
    
f_ax1.set_xlabel(r"$t$", fontsize = 20)
f_ax1.set_ylabel(r"$dt$", fontsize = 20)
f_ax1.legend(loc = 0,fontsize=15)

#Sub-Plot 2=============================================================
f_ax2.plot(t, Re, "b", lw=2.0)

f_ax2.set_xlabel("$t$", fontsize = 20)
f_ax2.set_ylabel("$Re$", fontsize = 20)
f_ax2.legend(loc = 0,fontsize=15)

#Sub-Plot 3=============================================================
f_ax3.plot(t, Div, "b", lw=2.0)

f_ax3.set_xlabel("$t$", fontsize = 20)
f_ax3.set_ylabel(r" $\nabla \cdot u$", fontsize = 20)
f_ax3.legend(loc = 0,fontsize=15)


#fig.tight_layout()
figManager = plt.get_current_fig_manager()
figManager.resize(*figManager.window.maxsize())
plt.show()

#LOOP===================================================================
while (1) :
    #Load Data==========================================================
    data = np.loadtxt('../TimeSeries.dat',comments='#')
    t = data[:,0]
    Re = data[:,1]
    Div = data[:,2]
    dt = data[:,3]

    #Define Figure and subplots=========================================
    fig = plt.figure(1)
    f_ax1 = fig.add_subplot(2,2,1)
    f_ax2 = fig.add_subplot(2,2,2)
    f_ax3 = fig.add_subplot(2,2,3)
    f_ax4 = fig.add_subplot(2,2,4)
    
    #Sub-Plot 1=============================================================
    f_ax1.plot(t, dt, "b", lw=2.0)
    
    #Sub-Plot 2=============================================================
    f_ax2.plot(t, Re, "b", lw=2.0)

    #Sub-Plot 3=============================================================
    f_ax3.plot(t, Div, "b", lw=2.0)
    
    plt.pause(2)
    plt.show()

