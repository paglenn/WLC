import matplotlib.pyplot as plt
import numpy as np
import sys
if len(sys.argv) == 1: exit("usage: python analysis_uwham.py num_windows")
numWindows = int( sys.argv[1] )
windowFiles = [open("window_%i"%w,'r') for w in range(numWindows) ]
frameSkip = 10
itr = 0

for f in windowFiles:
    Z = []
    t = []
    for line in f.readlines():
        row = line.split('\t')
        if itr % frameSkip  == 0 :
            t.append(float(row[0]))
            Z.append(float(row[1][:-1]))
        itr += 1
    #plt.hist(Z)
    plt.plot(t,Z,'.')
#plt.savefig('k_5.png')
plt.xlabel('timesteps')
plt.ylabel('z')
plt.title('%s windows'%numWindows)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.show()

