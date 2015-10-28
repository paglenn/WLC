import numpy as np
import matplotlib.pyplot as plt

a = np.genfromtxt("hist_0")
zvals = a[:,0]
logP = a[:,1]
l = 0.01

plt.plot(1- zvals/l, logP,'o')
plt.show()
