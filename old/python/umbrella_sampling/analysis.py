import numpy as np
import os
from parameters import *
import math
data_array = None

def calculate_height_dist(data_array):
    # as tangential component of R to initial displacement
    #t0 = data_array[0][0]
    Z = []
    if os.path.isfile(z_file):
        zFile = open(z_file,'r')
        for line in zFile.readlines():
            Z.append(float(line[:-1]))
    else:
        exit('%s missing!'%z_file)

    zFile.close()

    bins = [ list() for j in range(num_windows) ]
    binContents = [list() for j in range(num_windows) ]

    import matplotlib.pyplot as plt
    for j in range(num_windows):
        rng = (j*linesPerWindow,(j+1)*linesPerWindow)

        binContents[j],bins[j] = np.histogram(Z[rng[0]:rng[1]],bins=binsPerWindow,range=z_windows[j])
        #plt.semilogy(bins[j][:-1],binContents[j])
        print binContents[j]
        plt.hist(Z,bins=binsPerWindow+1,range=z_windows[j],histtype='step')
    plt.show()
    quit()
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
    #plt.savefig('z_dist.png')
    plt.show()
    quit()

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


