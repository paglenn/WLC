import numpy as np
import os
from parameters import *
import math

data_array = None

def retrieve():
    data_files = []
    for fname in window_files:
        if not os.path.isfile(fname): exit('Data file(s) missing!')
        data_files.append( open(fname,'r') )
    data_array = [ list() for x in range(total_frames) ]
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

def calculate_Tp(data_array):

    Tp = []
    for t in data_array:

        T = t[-1]
        tp_vec = T - np.dot(T,t0)
        tp= np.linalg.norm(tp_vec)
        Tp.append(tp)

    return Tp



def calculate_displacements(data_array) :
    allR = []
    for t in data_array:
        R = np.zeros(t[0].shape)
        for r in t[:-1]:
            R += r
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

    bins = [ list() for j in range(num_windows) ]
    binContents = [list() for j in range(num_windows) ]

    import matplotlib.pyplot as plt
    for j in range(num_windows):
        rng = (j*linesPerWindow,(j+1)*linesPerWindow)

        binContents[j],bins[j] = np.histogram(Z[rng[0]:rng[1]],bins=binsPerWindow,range=z_windows[j])
        #plt.semilogy(bins[j][:-1],binContents[j])
        print binContents[j]
        #plt.hist(Z,bins=binsPerWindow+1,range=z_windows[j],histtype='step')
    #plt.show()
    #quit()
    for j in range(num_windows-1):
        multiplier = 1.*binContents[j][-1]/binContents[j+1][0]
        for k in range(len(binContents[j+1])):
            binContents[j+1][k] *= multiplier
        #print binContents[j]
        #plt.plot(bins[j][:-1],binContents[j])
    #plt.show()
    #quit()
    for j in range(num_windows):
        for k in range(len(bins[j])-1):
            bins[j][k] = (bins[j][k+1] + bins[j][k])/2.
        bins[j] = np.delete(bins[j],bins[j][-1])
    for j in range(num_windows):
         print bins[j]
    maxFreq = max([max(x) for x in binContents])
    #print maxFreq
    A = list()
    Z = list()
    for j in range(num_windows):
        #print binContents[j]
        for k in range(binContents[j].shape[0]):
            if binContents[j][k] != 0:
                Z.append(bins[j][k])
                A.append(-math.log(binContents[j][k]/float(maxFreq)))
    print Z ,'\n', A
    plt.plot(Z,A)
    plt.savefig('z_dist.png')
    plt.show()
    quit()

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
            #A[i-1] = A[i-2]
            delta[j] = A[i] - A[i-1]
        A[i] -= delta[j]

    return RP,A

def plot_cosines() :
    cosFile = open(cos_file,'r')
    cos = []
    for line in cosFile.readlines():
        cos.append(float(line[:-1]))
    print np.mean(cos)
    binContents,bins = np.histogram(cos,range=(-1,1),bins=100,density=True)
    import matplotlib.pyplot as plt
    plt.semilogy(bins[:-1],binContents,'b.')

    x = np.linspace(0.95,1,100)
    P = 1./dx * np.exp(x/dx) /(2.*np.sinh(1./dx))
    plt.semilogy(x,P,'k')
    plt.savefig('cos_dist.png')
    #plt.show()
    return

#plot_cosines()

'''
F_rp = calculate_horizon_dist(data_array)
import matplotlib.pyplot as plt
plt.plot(F_rp[0],F_rp[1],'b-.')
plt.xlabel(r'$R_{\perp}$')
plt.ylabel(r'$\beta F(R_{\perp})$ ')
plt.title('Free energy' )
#plt.savefig('rp_dist.png')
plt.show()
'''


calculate_height_dist(data_array)
'''
F_z = calculate_height_dist(data_array)
import matplotlib.pyplot as plt
plt.plot(F_z[0],F_z[1],'b-.')
plt.savefig('z_dist.png')
#plt.vlines(window_edges,min(A),max(A))
plt.title('Free energy')
plt.xlabel(r'$z$')
plt.ylabel(r'$\beta F(z)$')
plt.show()
'''


