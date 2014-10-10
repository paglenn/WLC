from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) == 1: exit("usage: py analysis.py numWindows")
numWindows = int( sys.argv[1] )
windowFiles = [open("hist_%i"%w,'r') for w in range(numWindows) ]

allZ = []
allBC = []
fs = 1
for f in windowFiles:
    Z,F = [], []
    for line in f.readlines():
        l = line.split('\t')
        Z.append(float( l[0] ) )
        F.append(float(l[1][:-1]))
    allZ.append(Z) # list of numpy arrays
    allBC.append(F)
    plt.plot(Z,F)
plt.show()

# combine histograms?
allZ = [np.array(z)/max(z) for z in allZ ]
allBC = [np.array(bc)/max(bc) for bc in allBC ]
for j in range(1,numWindows):
    allBC[j] = allBC[j] *allBC[j-1][-1]/allBC[j][0]
    plt.plot(allZ[j],allBC[j])
plt.show()

quit()




