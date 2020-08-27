# *** FEANBEAM *** A Finite Element Analysis Program for Beams
#
# .... Copyright (c) 1984-2020: Regents of the University of California
#                               All rights reserved
#
# -----------------------------------------------------------------------
#      Modification log                                Date (dd/mm/year)
#        Original version                                    27/08/2020
# -----------------------------------------------------------------------
#      Purpose: Generation of deformation plot for example cantilever.f90
#
#      Inputs:
#        step            - Step size of reading data (1 would read all data, 10 would read every 10th data set)
#
#      Outputs:
#        cantilever.png  - Image of the plot
# -----------------------------------------------------------------------
#      Instructions on how to use custom user mesh input command GETMesh
#        1. Generate data with cantilever.f90
#        2. Execute command python examples/cantilever.py
# -----------------------------------------------------------------------

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Parameters
skip = 1

data = []
wdir = os.getcwd()

test = 'cantilever'

datapath = wdir + '\\' + test

nsteps = len(os.listdir(datapath))
for i in range(0, nsteps, skip):
    data.append(np.loadtxt(datapath+'\\'+'step'+'{:03d}'.format(i)+'.dat', skiprows=1))
data = np.array(data)

fig = plt.figure(figsize=(6*1.3,6*0.27))
ax = fig.gca()
ax.set_xlim((-0.15, 1.05))
ax.set_ylim((-0.05, 0.22))

for i in range(0,nsteps):
    ax.plot(data[i,:,3], data[i,:,1])

plt.show()
#plt.savefig(test+'.png')  # to save image
