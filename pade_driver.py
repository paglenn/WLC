# test suite for Pade_multi.py

from Pade_multi import *
import Pade_multi as Pade
import numpy as np

def fn(coeffs , x , M, N) :
    # Pade function will have
    # x: 1-by-n array of variables
    ncols = x.size
    num = coeffs[0]
    den = 1.0
    itr = 0

    for m in range(ncols ) :

        itr = 1 + m * M
        for i in range(M) :
            num += coeffs[itr] * x[m] ** (i+1)
            itr +=1

        itr = 1 + ncols * M + m * N
        # ens: constant term of denominator is 1
        for j in range(N):
            print itr
            print coeffs[itr]
            den += coeffs[itr] *x[m] ** (j+1)
            itr += 1

    print 'num: ', num
    print 'den: ', den
    if den == 0.0 :
        den += 1e-15
    return num / den
# test function - Runge

#print ans
#print fn(ans, np.array([1.0]), 2, 2)
#quit()
# let's see how well it does with a polynomial
X1 = np.linspace(-1.0,1.0, 20)
X2 = np.copy(X1)
X1v,X2v = np.meshgrid(X1,X2)
Yv = X1v * X1v + X2v * X2v
X1 = np.ravel(X1v)
X2 = np.ravel(X2v)
Y = np.ravel(Yv)

XV = np.vstack((X1,X2,Y))
XV = XV.T
ans = Pade.Pade(XV,2,2)


YA = []
for i in range(XV.shape[0])  :
    print XV[i,:-1]
    YA.append(fn(ans, XV[i,:-1],2,2))
    print YA[-1]
YA = np.array(YA)
#YA = YA.reshape(X1v.shape)

print ans

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(X1,X2,Y,c=u'r')
ax.scatter(X1,X2,YA)
plt.show()
#ax.plot_surface(X1v,X2v,YA,color='r',label='approx')
#ax.plot_surface(X1v,X2v,Yv,label='exact')
#plt.plot(X,Y,'k.',label='exact')
#plt.plot(X,YA,label='approx')
#plt.legend()
plt.savefig('test.png')
#plt.show()

# single variable case
if 0 :
    X = np.linspace(0,1,100)
    Y = 1./ (1 + 25*X*X)
    XV = np.vstack((X,Y))
    XV = XV.T

    ans = Pade(XV,2,2)
    # test function - essentially singular exponential
    #X = np.linspace(0,0.9,100)
    #Y = np.exp(1/(1-X))
    #XV = np.vstack((X,Y))
    #XV = XV.T
    #print XV
    #ans = Pade(XV,2,2)

    def fn(coeffs , x ,M,N) :
        # Pade function will have
        # x: 1-by-n array of variables
        num = coeffs[0]
        den = 1.0
        itr = 0
        ncols = x.size


        m = 0
        itr = 1 + m * M
        for i in range(M) :
            print 'num', i+ 1, coeffs[itr]
            num += coeffs[itr] * x[m] ** (i+1)
            itr +=1

        itr = (ncols  ) * M + 1 + m * N
        # ens: constant term of denominator is 1
        for j in range(N):
            print 'num', j+ 1, coeffs[itr]
            den += coeffs[itr] *x[m] ** (j+1)
            itr += 1

        if den == 0 :
            quit()

        return num / den
