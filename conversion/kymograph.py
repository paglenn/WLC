import matplotlib.pyplot as plt
import numpy as np
numWindows = 20
windowFiles = [open("window_%i"%w,'r') for w in range(numWindows) ]

for f in windowFiles:
    Z = []
    t = []
    for line in f.readlines():
        row = line.split('\t')
        try:
            t.append(float(row[0])/10)
            Z.append(float(row[1][:-1]))
        except ValueError:
            continue
    #plt.hist(Z)
    plt.plot(t,Z,'.')
#plt.savefig('k_5.png')
plt.show()

