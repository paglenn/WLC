import numpy as np
import os
from parameters import *
if not os.path.isfile(data_file): exit('Data file missing!')
infile = open(data_file,'r')

trajectory = [ list() for x in range(numsteps+1) ]
step = -1
for line in infile.readlines():
    line2 = line[:-1].split('\t')

    if 'step' in line2[0]:
        step += 1
        continue
    elif len(line2) == 1:
        continue

    data = [float(x) for x in line2]
    trajectory[step].append(np.array(data))

T = trajectory[-1][-1]


def calculate_displacements() :
    allR = []
    for filament in trajectory:
        R = np.zeros(filament[0].shape)

        for monomer in filament[:-1]:

            R += monomer

        allR.append(dx*R)

    return allR

def alt_calculate_height_dist():
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

#_R = calculate_displacements()
#_T = (np.dot(_T,_R)/np.linalg.norm(_R)**2.) * _R
#print("in-plane displacement: ", _R)
#F = free_energy()
#print("free energy: ", F)
Z  = alt_calculate_height_dist()

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


