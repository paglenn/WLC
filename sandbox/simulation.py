# simulation.py
# single actin filament (WLC) monte carlo routine
# umbrella sampling simulation adapted from that for Ising model
import numpy as np
from actin import *
from parameters import *

# choose initial configuration to be along z-axis

out = open(data_file,'w')


# write initial configuration to file
rpFile = open(rp_file,'w')
zFile = open(z_file,'w')


t = [np.copy(t0) for x in range(num_actins) ]
zFile.write("{0}\n".format(calculate_z(t)/L ) )
rpFile.write("{0}\n".format(calculate_rp(t)/L ) )
for j in range(1,numsteps+1):

    mc_step(t)

    zFile.write("{0}\n".format(calculate_z(t)/L ) )
    rpFile.write("{0}\n".format(calculate_rp(t)/L ) )

    write_tangents(t,j,out)

# close data files
rpFile.close()
zFile.close()
out.close()
