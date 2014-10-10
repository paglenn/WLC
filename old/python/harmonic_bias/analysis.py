# calculate z free energy given umbrella windows overlapping by one bin (EXP)
import numpy as np
import os
from parameters import *
import math
def calculate_height_dist():
    # as tangential component of R to initial displacement
    #t0 = [0][0]
    Z = []
    ctr = 0
    if os.path.isfile(z_file):
        zFile = open(z_file,'r')
        for line in zFile.readlines():
            ctr += 1
            try:
                Z.append(float(line[:-1]))
            except ValueError:
                print ctr
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
        plt.hist(Z,bins=binsPerWindow+1,range=z_windows[j],histtype='step')
    plt.show()
    quit()
    logP = [ [] for x in range(num_windows) ]
    for j in range(num_windows):
        logP[j] = np.log(binContents[j])

    for j in range(num_windows-1):
        difference = logP[j][-1] - logP[j+1][0]
        for k in range(len(logP[j+1])):
            logP[j+1][k] += difference
    maxP = max(max(x) for x in logP)
    for j in range(num_windows):
        for k in range(len(bins[j])-1):
            bins[j][k] = (bins[j][k+1] + bins[j][k])/2.
        bins[j] = np.delete(bins[j],[len(bins[j])-1])


    A = list()
    Z = list()
    for j in range(num_windows):
        #print binContents[j]
        for k in range(binContents[j].shape[0]):
            if binContents[j][k] != 0:
                Z.append(bins[j][k])
                A.append(-logP[j][k]+maxP )
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


calculate_height_dist()
'''
F_z = calculate_height_dist()
import matplotlib.pyplot as plt
plt.plot(F_z[0],F_z[1],'b-.')
plt.savefig('z_dist.png')
#plt.vlines(window_edges,min(A),max(A))
plt.title('Free energy')
plt.xlabel(r'$z$')
plt.ylabel(r'$\beta F(z)$')
plt.show()
'''


