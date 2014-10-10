import uwham
import numpy as np
import sys

Z = []
if len(sys.argv) == 1: exit("usage: python analysis_uwham.py num_windows")
numWindows = K = int( sys.argv[1] )
N_k = np.zeros(K)
u = list()
infile = open('uwham.dat','r')
for unformattedLine in infile.readlines():
    line = unformattedLine.split(' ')
    line = line[:-1] # remove \n
    #print line
    #print len(line)
    #quit()
    #zval = float(line[0])
    window = float(line[0])
    N_k[window] += 1.
    u.append( [ float(line[k]) for k in range(1,K+1) ])


u_kln = np.zeros((K,K,max(N_k)))
counter = 0
for k in range(K):
    for j in range(int(N_k[k])):
        u_kln[k,:,j] = u[counter]
        counter += 1

#print u_kln ; quit()

import pymbar
mbar = pymbar.MBAR(u_kln, N_k )
A = mbar.getFreeEnergyDifferences()[0]
print A
quit()
#results = uwham.UWHAM(u_kln,N_k)
#A = [a for a in results.f_k ]

import matplotlib.pyplot as plt
plt.plot([(i+0.5)/float(K) for i in range(K)],A)
plt.xlabel('z/L')
plt.ylabel('A(z/L)')
#plt.savefig('mbar.png')
plt.show()



