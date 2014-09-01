# simulation.py
# single actin filament (WLC) monte carlo routine
# umbrella sampling simulation adapted from that for Ising model
import numpy as np
from actin import *
from parameters import *
import copy

# choose initial configuration to be along z-axis
'''
out = dict()
for i in range(num_windows):
    out[rp_windows[i]] = open(window_files[i],'w')
'''

# write initial configuration to file
#zFile = open(z_file,'w')
rpFile = open(rp_file,'w')

positions = np.zeros((N,3))
positions[:,2] = np.linspace(0,L,N)
init_config = [ actin_filament(positions[i]) for i in range(N) ]
A = list(init_config)

for w in rp_windows:

    ################################
    # Passes over each window
    for itr in range(num_passes):

        ##############################################
        # Adjust horizontal displacement to target value
        A = adjust_rp(A,w)

        for j in range(1,numsteps+1):

            A = mc_step(A,w,rpFile)

            #write_tangents(A,j,out[w])

# close data files
'''
for dataFile in out.values():
    dataFile.write("program finished!!")
    dataFile.close()
'''
rpFile.close()
