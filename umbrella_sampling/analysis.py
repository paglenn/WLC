# calculate z free energy given umbrella windows overlapping by one bin (EXP)
import numpy as np
import os
import math
import matplotlib.pyplot as plt
cos_file = 'cos.dat'
dx = 0.1/100
k = 1./dx

def plot_cosines() :
    num_bins = 100

    cosFile = open(cos_file,'r')
    cos = []
    for line in cosFile.readlines():
        cos.append(float(line[:-1]))
    print np.mean(cos)
    binContents,bins = np.histogram(cos,range=(0.9,1),bins=100,density=True)
    cos = []; freq = []
    for j in range(num_bins):
        if binContents[j] != 0:
            cos.append(bins[j])
            freq.append(binContents[j])


    plt.plot(cos,np.log(freq),'b.')

    x = np.linspace(0.9,1,100)
    #P = 1./dx * np.exp(x/dx) /(2.*np.sinh(1./dx))
    logP = np.log(k) - np.log(1+math.exp(-2*k)) + k*(x-1)
    plt.plot(x,logP,'k')
    #plt.savefig('cos_dist.png')
    plt.show()
    return

plot_cosines()


#calculate_height_dist()
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


