import numpy as np
import matplotlib.pyplot as plt

a = np.genfromtxt('hist_0')
zPz = a[:,0]*a[:,1]
zbar = np.sum(zPz)
print zbar
