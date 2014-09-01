import numpy as np
import os
from parameters import *


data_array = None

def retrieve():
    data_files = []
    for fname in window_files:
        if not os.path.isfile(fname): exit('Data file(s) missing!')
        data_files.append( open(fname,'r') )
    data_array = [ list() for x in range(num_frames) ]
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
            data_array[step].append(np.array(data))

        infile.close()
    return data_array


#data_array =retrieve()
#T = data_array[-1][-1]
#Tp = T - np.dot(T,t0) * t0


def calculate_displacements(data_array) :
    allR = []
    for filament in data_array:
        R = np.zeros(filament[0].shape)

        for monomer in filament[:-1]:

            R += monomer

        allR.append(dx*R)

    return allR


def calculate_height_dist(data_array):
    # as tangential component of R to initial displacement
    #t0 = data_array[0][0]
    Z = []
    if os.path.isfile(z_file):
        zFile = open(z_file,'r')
        for line in zFile.readlines():
            Z.append(float(line[:-1]))
    else:
        zFile = open(z_file,'w')
        if data_array is None: data_array = retrieve()
        for R in calculate_displacements(data_array):
            Rt = np.dot(R,t0) * t0 # norm(t0) == 1
            z = np.linalg.norm(Rt)
            Z.append(z )
            zFile.write(str(z)+'\n')

    zFile.close()

    import matplotlib.pyplot as plt
    binContents, bins  = np.histogram(Z,bins= 100,range=(0,L),density=True)
    A = list() ; Z = list()
    for i in range(len(bins[:-1])) :
        if binContents[i] != 0:
            Z.append(bins[i])
            A.append(-np.log(binContents[i]))

    # shift windowed curves so they are continuous
    j = 0
    delta = dict()
    delta[j] = 0.
    for i in range(len(Z)):
        if Z[i] > z_window_edges[j]:
            j += 1
            A[i-1] = A[i-2]
            delta[j] = A[i] - A[i-1]
        A[i] -= delta[j]

    return Z,A

def calculate_horizon_dist(data_array):
    t_init = t0
    RP = list()
    if os.path.isfile(rp_file):
        rpFile = open(rp_file,'r')
        for line in rpFile.readlines():
            RP.append(float(line[:-1]))
    else:
        rpFile = open(rp_file,'w')
        if data_array is None: data_array = retrieve()
        for R in calculate_displacements(data_array):
            z = np.dot(R,t_init)  * t_init
            rp = np.linalg.norm(R - z)
            RP.append(rp)
            rpFile.write(str(rp) + '\n')

    rpFile.close()

    binContents, bins  = np.histogram(RP,bins= 100,range=(0,L),density=True)
    A = list() ; RP = list()
    for i in range(len(bins[:-1])) :
        if binContents[i] != 0:
            RP.append(bins[i])
            A.append(-np.log(binContents[i]))

    # shift windowed curves so they are continuous
    j = 0
    delta = dict()
    delta[j] = 0.
    for i in range(len(RP)):
        if RP[i] > rp_window_edges[j]:
            j += 1
            A[i-1] = A[i-2]
            delta[j] = A[i] - A[i-1]
        A[i] -= delta[j]

    return RP,A

F_rp = calculate_horizon_dist(data_array)
import matplotlib.pyplot as plt
plt.plot(F_rp[0],F_rp[1],'b-.')
plt.show()



'''
F_z = calculate_height_dist(data_array)
import matplotlib.pyplot as plt
plt.plot(F_z[0],F_z[1],'b-.')
#plt.savefig('z_dist.png')
#plt.vlines(window_edges,min(A),max(A))
plt.title('Free energy')
plt.xlabel(r'$z$')
plt.ylabel(r'$\beta F(z)$')
plt.show()
'''

