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
itr = 0.
for f in windowFiles:
    Z = []
    mybins= np.linspace(itr/numWindows,(itr+1)/numWindows+overlap,numBins/numWindows+1)
    if windowFiles[-1] == f: mybins= np.linspace(itr/numWindows,(itr+1)/numWindows,numBins/numWindows)
    for line in f.readlines():
        try:
            Z.append(float(line.split('\t')[-1][:-1]))
        except:
            pass
    itr += 1.
    binContents,bins = np.histogram(Z,bins=mybins,density=True)
    #bins = [0.5*sum(bins[j:j+2]) for j in range(len(bins)-1) ]
    allZ.append(mybins) # list of numpy arrays
    allBC.append(binContents)

A = []
Z = []

# remove edge artifacts
for i in range(numWindows-1):
    allZ[i] = allZ[i][1:-1]
    allBC[i] = allBC[i][1:-1]

for i in range(numWindows-1):
    if allBC[i+1][0] != 0:
        allBC[i+1] = allBC[i+1] * allBC[i][-1]/allBC[i+1][0]
    else:
        exit('zero content bin')
'''
for i in range(numWindows-1):
    if allBC[i+1][0] != 0:
        allBC[i] = allBC[i]* allBC[i+1][0]/float(allBC[i][-1])
    else:
        exit('zero content bin')
'''





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





