# simulation parameters

numsteps = int(1e2)

L = 0.5

dx = 0.002

Ldx = int(L/dx)

num_actins = N = Ldx + 1

data_file = 'output.dat'

import numpy as np

t0 = np.array([0,0,1])


# for umbrella sampling
umbrella_bias = True

num_windows = 10

window_files = ["window_{0}.dat".format(i) for i in range(num_windows) ]

num_passes = 500

num_frames = num_windows * num_passes * (numsteps+1)

z_windows = [(i*L/num_windows,(1+i)*L/num_windows) for i in range(num_windows) ]

z_windows.reverse()

z_window_edges = [ Z[1] for Z in z_windows]

z_window_edges.reverse()

z_file = 'heights.dat'

rp_windows = [(i*L/num_windows,(1+i)*L/num_windows) for i in range(num_windows) ]

rp_window_edges = [ rp[1] for rp in rp_windows]

rp_file = 'horizons.dat'

