# simulation.py
# single actin filament (WLC) monte carlo routine
# umbrella sampling simulation adapted from that for Ising model
import numpy as np
from actin import *
from parameters import *

# choose initial configuration to be along z-axis

out = dict()
for i in range(num_windows):
    out[z_windows[i]] = open(window_files[i],'w')


# write initial configuration to file
# rpFile = open(rp_file,'w')
zFile = open(z_file,'w')


for w in z_windows:
    t = [np.copy(t0) for x in range(num_actins) ]
    ################################
    # Passes over each window
    for itr in range(num_passes):


        ##############################################
        # Adjust horizontal displacement to target value
        t = adjust_z(t,w)
        #t = adjust_rp(t,w)
        #write_tangents(t,0,out[w])

        zFile.write("{0}\n".format(calculate_z(t) ) )

        for j in range(1,numsteps+1):

            umbrella_mc_step(t,w)

            zFile.write("{0}\n".format(calculate_z(t) ) )

            write_tangents(t,j,out[w])

# close data files

for dataFile in out.values():
    dataFile.write("program finished!!")
    dataFile.close()

#rpFile.close()
zFile.close()
