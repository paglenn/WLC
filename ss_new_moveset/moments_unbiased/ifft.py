import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
import math
import cmath

fin = 'moments.dat'
M = np.genfromtxt(fin)
N = M[:,0]
M = M[:,1]

#k = np.linspace(-100,100,10000)
#Gk = np.zeros_like(k,dtype=np.complex_)
#for n in N:
#    Gk += (1j*k)**n / math.factorial(n) * M[n-1]
#Gk += np.ones_like(Gk)

def G(k) :
    Gk = 0
    for n in N[:50]:
        if n > 150:
            lognfact = n*math.log(n/math.e) + 0.5 * math.log(2*math.pi*n)
            nfact = np.exp(lognfact)
        else:
            nfact = math.factorial(n)
        gk = cmath.exp(n*1j*math.pi/2.) * (k)**n / nfact * M[n-1]
        print 'gk: ' , gk
        Gk += gk
    Gk += 1.
    return Gk

def I(k,x):
    return G(k) * np.exp(1j*k*x)

def Re_I(k,x) :
    I_full = I(k,x)
    return I_full.real

def Im_I(k,x) :
    I_full = I(k,x)
    return I_full.imag

X = np.linspace(0,1,1001)
Px = np.zeros_like(X)
for i in range( X.size ):
    Px[i] = (1./2*math.pi) * scipy.integrate.quad(Re_I,-100,100,
            args=(X[i],))[0]
    #print 1j*(1./2*math.pi) * scipy.integrate.quad(Im_I,-100,100,
    #        args=(X[i],))[0]

plt.plot(X,Px)
plt.show()
