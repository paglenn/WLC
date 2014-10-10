# wham_driver.py
# plot FE from WHAM output file
from subprocess import call
import numpy as np
nb = 50
T = 100.
call("./wham"+" 0.0 1.0 {0} 1e-10 {1:.2f} 0 metadata.dat wham_output.dat".format(nb,T),shell=True)

infile = open('wham_output.dat','r')
Z = list()
A = list()
for line in infile.readlines():
    cols = line.split('\t')[:-1]
    if any('#' in x for x in cols): continue
    Z.append(float(cols[0]))
    A.append(float(cols[1]))

dA = [A[i+1] - A[i-1] for i in range(len(A[1:nb/4])) ]
z_lo = dA.index(min(dA))
z_hi = A.index(min(A))
z_target = Z[z_hi] - 2*np.std(Z)

for i in range(len(Z)):
    if Z[i] < z_target:
        z_hi = i
z_hi = 5*len(Z)//8
z_lo = 3*len(Z)//8
#z_hi = 62; z_lo = 24
import matplotlib.pyplot as plt
plt.plot(Z,A,lw=2.0)
import numpy as np
linearFit = np.polyfit(Z[z_lo:z_hi],A[z_lo:z_hi],1)
#linearFit = np.polyfit(Z[nb/4:3*nb/4+1],A[nb/4:3*nb/4+1],1)
fb = - 0.1*linearFit[0] # buckling force
plt.xlabel(r'$z/L$')
plt.ylabel(r'$A(z/L)$')
plt.plot([Z[z_lo],Z[z_hi]],[A[z_lo],A[z_hi]],'k',lw=0.5, label=r"bounds used for $F_{b}$")
plt.plot([Z[z_lo],Z[z_hi]],[A[z_lo],A[z_hi]],'k.',ms=12)
plt.title(r'$\Delta = %.3f\ell_{p}$ : $f_{b}\approx %.3f$'%(0.001,fb))
#plt.ylim(-0.1,max(A)+0.1)
#plt.xlim(0,1.1)
#plt.axvline(1.0,min(A),max(A),color='k',ls='--',label=r'$z=L$')
plt.legend()
#plt.savefig('fe_z.png')
plt.show()
