# single actin filament (WLC) monte carlo routine
import numpy as np
from actin import *
from parameters import *

# choose initial configuration to be along z-axis
if umbrella_bias:
    out = dict()
    for i in range(num_windows):
        out[z_windows[i]] = open(window_files[i],'w')
else:
    outfile = open(data_file, 'w')

N = num_actins
positions = np.zeros((N,3))
positions[:,2] = np.linspace(0,L,N)
A0 = [ actin_filament(positions[i],dx) for i in range(N) ]

# write initial configuration to file

for dataFile in out.values():
    write_tangents(A0,0,dataFile)


for w in z_windows:

    A = list(A0)
    A = adjust_z(A0,A,w)

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

        write_tangents(A,j,out[w])

for dataFile in out.values():
    dataFile.write("program finished!!")
    dataFile.close()




