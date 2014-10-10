from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
numWindows = 10
numBins = 100.
overlap = 0.01
itr = 0.
N = 100
windowFiles = [open("window_%i"%w,'r') for w in range(numWindows) ]

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

A = []
Z = []
# remove edge artifacts
for i in range(numWindows-1):
    allZ[i] = allZ[i][1:-1]
    allBC[i] = allBC[i][1:-1]
#allBC[i+1] = allBC[i+1][1:]

allBC = [np.array(x) for x in allBC ]
for i in range(numWindows-1):
    if allBC[i+1][0] != 0:
        allBC[i+1] = allBC[i+1] * allBC[i][-1]/allBC[i+1][0]
    else:
        exit('zero content bin')

for i in range(numWindows):
    for j in range(len(allBC[i])):
        if allBC[i][j] > 0:
            A.append(-np.log(allBC[i][j]))
            Z.append(allZ[i][j])

plt.plot(Z,A)
plt.xlabel('z/L')
plt.ylabel('A(z/L)')
plt.title('L = 0.1')
#plt.savefig('/Users/paulglen/Dropbox/point_1.png')
plt.show()





