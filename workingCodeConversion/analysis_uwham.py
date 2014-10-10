import uwham
import numpy as np
import sys

Z = []
if len(sys.argv) == 1: exit("usage: python analysis_uwham.py num_windows")
numWindows = K = int( sys.argv[1] )
N_k = np.zeros(K)
u = list()
uklnFile = open('ukln.dat','r')
k = [2.1e2 for w in range(numWindows) ]
k[-1] *=
u_kn = np.zeros(
for unformattedLine in uklnFile.readlines():
    line = unformattedLine.split(' ')
    line = line[:-1] # remove \n
    #zval = float(line[0])
    window = float(line[0])
    N_k[window] += 1.
    u.append( [ float(line[k]) for k in range(1,K+1) ])

Nk_max = max(N_k)
u_kln = np.zeros((K,K,Nk_max))
counter = 0
for k in range(K):
    for j in range(int(N_k[k])):
        u_kln[k,:,j] = u[counter]
        counter += 1

u_kn = np.zeros([K,N_max],np.float64)
uknFile = open('ukn.dat','r')
for line in uknFile.readlines():



results = uwham.uwham(u_kln,N_k)


#A  = mbar.getFreeEnergyDifferences()[0]


import matplotlib.pyplot as plt
plt.plot([(i+0.5)/float(K) for i in range(K)],A)
plt.xlabel('z/L')
plt.ylabel('A(z/L)')
plt.title('%s windows'%numWindows)
#plt.savefig('mbar.png')
plt.show()



