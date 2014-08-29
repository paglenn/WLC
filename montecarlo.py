# actin filament monte carlo routine

import numpy as np
from actin import *

# choose initial configuration to be along z-axis
dx = 1
L =  100
nsteps = int(1e5)
N = int(L/dx) + 1
data_file = open('output.dat','w')

positions = np.zeros((N,3))
positions[:,2] = np.linspace(0,L,N)
A = [ actin_filament(positions[i],dx) for i in range(N) ]

# write initial config to file
write(A,0,data_file)

for j in range(1,nsteps+1):

    random_vec = np.random.randn(3)

    random_index= np.random.choice(range(N))
    #print(A[random_index].t)

    A[random_index].perturb(random_vec)
    #print(A[random_index].t)

    dE = delta_E(A,random_index)

    if dE > 0:
        if np.random.rand() > np.exp(-dE):
            A[random_index].revert()

    update(A,random_index)

    write(A,j,data_file)

data_file.close()





