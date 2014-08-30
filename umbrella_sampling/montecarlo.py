# single actin filament (WLC) monte carlo routine
import numpy as np
from actin import *
from parameters import *

# choose initial configuration to be along z-axis
'''
dx = 0.002
L =  0.2
numsteps = int(1e3)
N = int(L/dx) + 1
data_file = open('output.dat','w')

# umbrella sampling windows
num_windows = 10
windows = [(i*L/num_windows,(1+i)*L/num_windows) for i in range(num_windows) ]

'''
outfile = open(data_file, 'w')
N = num_actins
positions = np.zeros((N,3))
positions[:,2] = np.linspace(0,L,N)
A0 = [ actin_filament(positions[i],dx) for i in range(N) ]
# instead of rigid initial configuration, rotate all by their own angle
'''
for i in range(len(A)):

    phi = np.random.uniform(0,2*np.pi)

    A[i].rotate(phi)

    update(A,i)
'''
# write initial configuration to file
write_tangents(A,0,outfile)


for w in z_windows:

    A = adjust_z(A0,A0,w)

    for j in range(1,numsteps+1):

        random_vec = np.random.randn(3)

        random_index= np.random.choice(range(N))
        #print(A[random_index].t)

        A[random_index].perturb(random_vec)
        #print(A[random_index].t)

        if umbrella_bias :
            dE = deltaE_z_bias(A0,A,random_index,w)
        else:
            dE = deltaE(A,random_index)

        if dE > 0:
            if np.random.rand() > np.exp(-dE):
                A[random_index].revert()

        update(A,random_index)

        write_tangents(A,j,outfile)

outfile.close()





