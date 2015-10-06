import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate

def w(x):
    #return 1.0
    return math.exp( - 1./(2.- x ) )

def f(x,u):
    ans = 0
    for n in range(u.size) :
        ans += u[n]* x**n
    return ans

def product(x,u,v) :
    return w(x) * f(x,u) * f(x, v)

def dot(u,v):
    a = 0
    b = 1
    ans, err  = scipy.integrate.quad(product, a,b,args=(u,v) )
    #print "erro: ", err
    if err > 1e-6: print "Warning: high error integration"
    return ans

def norm(u):
    return math.sqrt(dot(u,u))

def projection(u,v) :
    _v_ = norm(v)
    if _v_ == 0. :
        print "zero vector returned!"
        quit()
    return dot(u,v) / dot(v,v )

#N = 30

def genpoly(N) :
    U = np.eye(N)
    V = np.zeros((N,N))

    for i in range(N):
        #print U[:,i]
        V[:,i] = U[:,i]
        for j in range(i) :
            V[:,i] -= projection(U[:,i], V[:,j]) * V[:,j]

    for i in range(N):
        V[:,i] /= norm(V[:,i])
    V = V.T
    # now rows correspond to the Ln
    return V

#V = V.T # now rows correspond to the Ln
#np.savetxt("polynomials.dat", V)

