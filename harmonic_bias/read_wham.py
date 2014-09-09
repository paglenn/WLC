# wham_driver.py
# plot FE from WHAM output file
from subprocess import call
#call("./wham"+" 0.0 1.0 100 0.0001 200 0 metadata.dat wham_output.dat",shell=True)

infile = open('wham_output.dat','r')
Z = list()
A = list()
for line in infile.readlines():
    cols = line.split('\t')[:-1]
    if any('#' in x for x in cols): continue
    Z.append(float(cols[0]))
    A.append(float(cols[1]))

import matplotlib.pyplot as plt
plt.plot(Z,A)
plt.xlabel(r'$z/L$')
plt.ylabel(r'$F(z/L)$')
plt.title('Two-bin overlap, k = 100')
plt.savefig('fe_z.png')
plt.show()
