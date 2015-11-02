import numpy as np
import math
def find_most_probable() :

    a = np.genfromtxt('hist_0')
    Z = a[:,0]
    A = a[:,1]
    return Z[np.argsort(A)[0]]


