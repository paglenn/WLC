import numpy as np
import os
from parameters import *

data_files = []
for fname in window_files:
    if not os.path.isfile(fname): exit('Data file(s) missing!')
    data_files.append( open(fname,'r') )


trajectory = [ list() for x in range(num_windows * (numsteps+1) ) ]
step = -1
for infile in data_files:

    for line in infile.readlines():
        line2 = line[:-1].split('\t')

        if 'step' in line2[0]:
            step += 1
            continue
        elif len(line2) == 1:
            continue
        elif 'finished' in line2[0]:
            break

        data = [float(x) for x in line2]
        trajectory[step].append(np.array(data))

T = trajectory[-1][-1]

# calculate net horizontal displacement and net orientation
def calculate_displacement():

    init_filament = trajectory[0]
    final_filament = trajectory[-1]
    _R = z = 0
    for i in range(len(final_filament)):
        scalar_projection = np.dot(final_filament[i],init_filament[i])/np.linalg.norm(init_filament[i])**2.
        z += scalar_projection * initial_filament[i]
        dR = final_filament[i] - scalar_projection*init_filament[i]
        if np.linalg.norm(dR) > 0:
            dR = dR / np.linalg.norm(dR)
        _R += dR * dx

    return z, _R

def calculate_displacements() :
    allR = []
    for filament in trajectory:
        R = np.zeros(filament[0].shape)

        for monomer in filament[:-1]:

            R += monomer

        allR.append(dx*R)

    return allR


def calculate_height_dist():
    # as tangential component of R to initial displacement
    t_init = trajectory[0][0]
    Z = []
    for R in calculate_displacements():

        tR = (np.dot(R,t_init)/np.linalg.norm(t_init)**2. ) * t_init
        Z.append(np.linalg.norm(tR) )

    return Z

# calculate net difference in tangent vectors
dT = trajectory[-1][-1] - trajectory[-1][0]

def free_energy():
    FE = - Ldx * (1./dx) * np.log(np.sinh(1./dx) )
    return FE

def mf_energy():
    return 1./2 * (1./L) * np.linalg.norm(dT) ** 2.


Z  = calculate_height_dist()

# convert from frequency distribution to probability dist

import matplotlib.pyplot as plt
binContents, bins  = np.histogram(Z,bins= 100,range=(0,L),density=True)
A = list() ; Z = list()
for i in range(len(bins[:-1])) :
    if binContents[i] != 0:
        Z.append(bins[i])
        A.append(-np.log(binContents[i]))

plt.plot(Z,A)
plt.title('Free energy')
plt.xlabel(r'$z$')
plt.ylabel(r'$\beta F(z)$')
plt.savefig('z_dist.png')
plt.show()


