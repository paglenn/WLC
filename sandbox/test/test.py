import numpy as np
import math

X = np.linspace(0,10,100) ;
Y = exp(-X*X/2 ) ; # standard normal

nsamples = 2000 ;
xp = np.zeros(2000) ;
yp = np.zeros(2000) ;
for i in range(nsamples):
    u = 100
    v = 100
    while u > math.exp(-v*v/2 ) :
        u = np.random.rand()
        v = np.random.rand()
