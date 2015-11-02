import numpy as np
import matplotlib.pyplot as plt

a = np.genfromtxt('moments.dat')
n = a[:,0]
M = a[:,1]

plt.plot(n,M)
plt.show()
