import numpy as np
import os
infile = 'output.dat'
if not os.path.isfile(infile): exit('Data file missing!')
infile = open(infile,'r')

num_steps = int(1e5)
dx =1
L = 100
Ldx = int(L/dx)
num_actins = Ldx + 1

trajectory = [ list() for x in range(num_steps+1) ]
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

# calculate net horizontal displacement and net orientation
def calculate_displacement():

    init_filament = trajectory[0]
    final_filament = trajectory[-1]
    R = T = 0
    for i in range(len(final_filament)):
        scalar_projection = np.dot(final_filament[i],init_filament[i])/np.linalg.norm(init_filament[i])**2.
        dR = final_filament[i] - scalar_projection*init_filament[i]
        dR = dR / np.linalg.norm(dR)
        R += dR * dx

        T += final_filament[i]

    Tn = T / np.linalg.norm(T)
    return Tn, R

# calculate net difference in tangent vectors
dT = trajectory[-1][-1] - trajectory[-1][0]

def free_energy():
    FE = - Ldx * (1./dx) * np.log(np.sinh(1./dx) )
    return FE

def mf_energy():
    return 1./2 * (1./L) * np.linalg.norm(dT) ** 2.

Tn, R = calculate_displacement()
print("net orientation: ", Tn)
print("horizontal displacement: ", R)
F = free_energy()
print("free energy: ", F)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(trajectory[-1])
ax.plot(trajectory[0])

