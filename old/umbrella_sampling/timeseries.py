import matplotlib.pyplot as plt
import numpy as np
numWindows = 90
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
plt.show()

