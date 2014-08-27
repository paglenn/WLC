import numpy as np
infile = open('output.dat','r')

num_actins = 101
num_steps = int(1e5)
dx = 0.1

trajectory = [ list() for x in range(num_steps+1) ]
step = -1
for line in infile.readlines():
    line2 = line[:-1].split('\t')

    if 'step' in line2[0]:
        step += 1
        continue
    elif len(line2) == 1:
        continue

    data = [float(x) for x in line2]
    trajectory[step].append(np.array(data))

# compute angles
cos = []
for filament in trajectory:

    for i in range(num_actins - 1) :

        cos.append( np.dot(filament[i],filament[i+1]) )

import matplotlib.pyplot as plt
plt.hist(cos,bins=100,range=(-1,1),normed=True,histtype='step')
x = np.linspace(-1,1,1000)
P = 1./dx * np.exp(x/dx) / (2*np.sinh(1./dx) )
plt.plot(x,P,'k',lw=1.5)
plt.show()





