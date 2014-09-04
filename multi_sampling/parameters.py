# simulation parameters
# note: l_p is the persistence length (~37 um for actin rod )

numsteps = int(1e4)

dx = 0.0027 # (units of l_p)
num_actins = N = 100
L = N * dx

import numpy as np

t0 = np.array([0,0,1])

nbr = dict()
for j in range(num_actins):

    if j == 0: nbr[j] = [j+1]
    elif j+1 == num_actins : nbr[j] = [j-1]
    else: nbr[j] = [j+1,j-1]

cos_file = 'cos.dat'
# for umbrella sampling
umbrella_bias = True

num_windows = 15

window_files = ["window_{0}.dat".format(i) for i in range(num_windows) ]

num_passes = 100 # passes per window

total_frames = num_windows * num_passes * numsteps

# natural units for z and rp are those of L/l_p
z_windows = [(i/num_windows,(1.+i)/num_windows) for i in range(num_windows) ]

z_windows.reverse()

z_window_edges = [ Z[1] for Z in z_windows]

z_window_edges.reverse()

z_file = 'zvals.dat'

rp_file = 'rpvals.dat'

rptp_file = 'rptpvals.dat'

tp_file = 'tpvals.dat'
