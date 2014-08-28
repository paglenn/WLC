import numpy as np
infile = open('output.dat','r')

num_steps = int(1e5)
dx = 0.5
L = 500
num_actins = int(L/dx)  + 1

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


def free_energy():
    Ld = num_actins - 1
    FE = - Ld * (1./dx) * np.log(np.sinh(1./dx) )
    return FE

