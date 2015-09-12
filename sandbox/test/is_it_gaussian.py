import numpy as np
import math
import matplotlib.pyplot as plt

def is_it_gaussian(the_file):
    a = np.genfromtxt(the_file)
    binContents, binEdges = np.histogram(a, density=True)

    X = list()
    logP = list()
    for i in range(binContents.size) :
        if binContents[i] > 0 :
            X.append(0.5* (binEdges[i]+binEdges[i+1]) )
            logP.append(math.log(binContents[i]))

    X = np.array(X)
    logP = np.array(logP)
    testfn = - (X*X) / 2. - 0.5 * math.log(2*np.pi)
    plt.plot(X,logP, label = "test")
    plt.plot(X,testfn, label = "reference")
    plt.legend()
    plt.show()


the_file = "gaussian_output"
is_it_gaussian(the_file )
