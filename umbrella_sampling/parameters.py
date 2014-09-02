# simulation parameters

numsteps = int(1e3)

#L = 0.5
L = 1.0

#dx = 0.002
dx = 0.01

Ldx = int(L/dx)

num_actins = N = Ldx + 1

data_file = 'output.dat'

import numpy as np

t0 = np.array([0,0,1])

nbr = dict()
for j in range(num_actins):

    if j == 0: nbr[j] = [j+1]
    elif j+1 == num_actins : nbr[j] = [j-1]
    else: nbr[j] = [j+1,j-1]





# for umbrella sampling
umbrella_bias = True

num_windows = 10

window_files = ["window_{0}.dat".format(i) for i in range(num_windows) ]

num_passes = 50

total_frames = num_windows * num_passes * numsteps

z_windows = [(i*L/num_windows,(1+i)*L/num_windows) for i in range(num_windows) ]

z_windows.reverse()

z_window_edges = [ Z[1] for Z in z_windows]

z_window_edges.reverse()

z_file = 'heights.dat'

rp_windows = [(i*L/num_windows,(1+i)*L/num_windows) for i in range(num_windows) ]

rp_window_edges = [ rp[1] for rp in rp_windows]

rp_file = 'horizons.dat'
