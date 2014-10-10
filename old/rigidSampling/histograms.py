import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) == 1: exit("usage: python analysis_uwham.py num_windows")
numWindows = int( sys.argv[1] )
windowFiles = [open("window_%i"%w,'r') for w in range(numWindows) ]

for f in windowFiles:
    Z = []
    for line in f.readlines():
        try:
            Z.append(float(line.split('\t')[1][:-1]))
        except ValueError:
            continue
    #plt.hist(Z)
    binContents,bins = np.histogram(Z,bins=10,density=True)
    bins = [0.5*sum(bins[i:i+2]) for i in range(len(bins)-1) ]
    print binContents
    plt.plot(bins,binContents)
#plt.savefig('k_5.png')
plt.show()


