from parameters import *
import matplotlib.pyplot as plt
import numpy as np

windowFiles = [open(fname,'r') for fname in window_files ]

for f in windowFiles:
    Z = []
    for line in f.readlines():
        try:
            Z.append(float(line.split('\t')[-1][:-1]))
        except ValueError:
            continue
    #plt.hist(Z)
    binContents,bins = np.histogram(Z,bins=10,density=True)
    bins = [0.5*sum(bins[i:i+2]) for i in range(len(bins)-1) ]
    plt.plot(bins,binContents)
plt.show()


