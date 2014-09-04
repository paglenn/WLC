# simulation parameters
# note: l_p is the persistence length (~37 um for actin rod )

numsteps = int(5e5)

L = 0.25 #(units of l_p) Whitelam et al 2008
dx = 0.005 # (units of l_p) calculated from proteopedia data
num_actins = N = int(L/dx)

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

num_windows = 1

window_files = ["window_{0}.dat".format(i) for i in range(num_windows) ]

num_passes = 1 # passes per window

total_frames = num_windows * num_passes * numsteps

# natural units for z and rp are those of L/l_p
z_windows = [(i/num_windows,(1.+i)/num_windows) for i in range(num_windows) ]

z_windows.reverse()

z_window_edges = [ Z[1] for Z in z_windows]

z_window_edges.reverse()

z_file = 'zvals.dat'

rp_windows = [(i/num_windows,(1+i)/num_windows) for i in range(num_windows) ]

rp_window_edges = [ rp[1] for rp in rp_windows]

rp_file = 'rpvals.dat'

rptp_file = 'rptpvals.dat'

tp_file = 'tpvals.dat'
