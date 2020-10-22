import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Parameters
skip = 10

data = []
wdir = os.getcwd()

test = 'arch-results'

datapath = wdir + '\\' + test

nsteps = len(os.listdir(datapath))
for i in range(0, nsteps, skip):
    data.append(np.loadtxt(datapath+'\\'+'step'+'{:03d}'.format(i)+'.dat', skiprows=1))
data = np.array(data)

fig = plt.figure(figsize=(6*1.3,6*1.3))
ax = fig.gca()
ax.set_xlim((-2.35, 2.05))
ax.set_ylim((-2.5, 2.22))

for i in range(len(data)):
    ax.plot(data[i,:,3], data[i,:,1])

plt.savefig(test+'.png')  # to save image
