import numpy as np
import os
infile = 'output.dat'
if not os.path.isfile(infile): exit('Data file missing!')
infile = open(infile,'r')

num_steps = int(1e6)
dx =1
L = 1000
Ldx = int(L/dx)
num_actins = Ldx + 1

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

# calculate tip deviation
filament = trajectory[-1]
end_actin = filament[-1]
T = [ filament[-1] for filament in trajectory ]

def calculate_displacement():

	original_filament = trajectory[0]
	final_filament = trajectory[-1]
	R = 0
	for i in range(len(final_filament)):
		scalar_projection = np.dot(final_filament[i],init_filament[i])/np.linalg.norm(init_filament[i])**2.
		dR = final_filament[i] - scalar_projection*final_filament[i]
		dR = dt / np.linalg.norm(dt)
		R += dR * dx
	return R

def free_energy():
    FE = - Ldx * (1./dx) * np.log(np.sinh(1./dx) )
    return FE

R = calculate_displacement()
print R
F = free_energy()
print F
