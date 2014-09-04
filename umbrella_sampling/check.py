import numpy as np
import os
from parameters import *

def retrieve():
    data_files = []
    for fname in window_files:
        if not os.path.isfile(fname):
            print(fname)
            exit('Data file(s) missing!')
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

data_array = retrieve()

'''
for t in data_array:
    for actin in t:
        if abs(np.linalg.norm(actin) - 1.0) > 0.01*dx :
            print(np.linalg.norm(actin),"uh-oh")
            quit()
'''


# compute P( cos \theta  )

def compute_cos(data_array):
    import matplotlib.pyplot as plt
    cos = []
    for t in data_array:

        for i in range(num_actins - 1) :

            cos.append( np.dot(t[i],t[i+1]) )

    x = np.linspace(0.8,1,100)
    P = 1./dx * np.exp(x/dx) / (2*np.sinh(1./dx) )
    plt.hist(cos,bins=100,range=(0.8,1),normed=True,histtype='step',log=True)
    plt.semilogy(x,P,'k',lw=1.2)
    plt.ylabel('semilog cosine dist')
    plt.xlabel(r'$\cos \theta $')
    plt.ylabel(r'$P(\cos \theta ) $')
    plt.savefig('cos.png')

def compute_corr:
    # compute correlator
    R = range(0,num_actins//10)
    G = [0 for r in R ]
    num_samples = list(G)
    for t in data_array:

        for r in R:

            for i in range(num_actins - r):

                num_samples[r] += 1
                G[r] += np.dot(t[i],t[i+r])

    for r in R:

        if num_samples[r] != 0:

            G[r] = G[r] / num_samples[r]

    # rescale R
    R = np.array( [dx * r for r in R ] )

    # model prediction: G(s) ~ exp(-s/l_p)
    expR = np.exp(-R)
    plt.semilogy(R,np.array(G),'b.')
    plt.semilogy(R, expR,'k',lw=1.5 )
    plt.xlabel(r's/$l_p$')
    plt.ylabel('G(s) ')
    plt.savefig('corr.png')
