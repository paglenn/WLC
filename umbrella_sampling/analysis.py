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

def calculate_height_dist() :
    init_filament = trajectory[0]
    Z = []
    for j in range(len(trajectory)):
        filament = trajectory[j]
        z = 0
        for i in range(len(trajectory[j])):
            scalar_projection = np.dot(filament[i],init_filament[i]) / np.linalg.norm(init_filament[i])**2.
            _t = scalar_projection * init_filament[i]
            z += np.linalg.norm(_t)
        Z.append(dx*z)
    return Z

# calculate net difference in tangent vectors
dT = trajectory[-1][-1] - trajectory[-1][0]

def free_energy():
    FE = - Ldx * (1./dx) * np.log(np.sinh(1./dx) )
    return FE

def mf_energy():
    return 1./2 * (1./L) * np.linalg.norm(dT) ** 2.

#_R = calculate_displacement()
#_T = (np.dot(_T,_R)/np.linalg.norm(_R)**2.) * _R
#print("in-plane displacement: ", _R)
#F = free_energy()
#print("free energy: ", F)
Z = calculate_height_dist()
import matplotlib.pyplot as plt
plt.hist(Z,bins= 100,range=(0,L),histtype='step',normed=True,log=True)
plt.show()


