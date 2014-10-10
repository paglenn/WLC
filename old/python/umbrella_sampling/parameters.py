# simulation parameters
# note: l_p is the persistence length (~37 um for actin rod )

numsteps = int(1e4)

dx = 0.0027 # (units of l_p -- from "Actin Based Motility" )
dx = 0.2 # better for testing
num_actins = N = 100
L = N * dx

import numpy as np

t0 = np.array([0,0,1])

nbr = dict()
for j in range(num_actins):

    if j == 0: nbr[j] = [j+1]
    elif j+1 == num_actins : nbr[j] = [j-1]
    else: nbr[j] = [j+1,j-1]

cos_file = 'cos.dat' #saves time for checking fidelity of simulation

# for umbrella sampling
umbrella_bias = True
num_passes = 50 # passes per window
num_windows = 20
num_bins = 200 # each window will have num_bins/num_windows+1 bins
binWidth = 1./ num_bins
overlap = binWidth
z_windows = [[1.*i/num_windows,(i+1)/float(num_windows)+overlap] for i in range(num_windows)]
linesPerWindow = num_passes*(numsteps+1)

z_windows[-1][1] -= overlap
binsPerWindow= num_bins/float(num_windows)+1
# determine z windows

window_files = ["window_{0}.dat".format(i) for i in range(num_windows) ]



total_frames = num_windows * num_passes * numsteps

# natural units for z and rp are those of L/l_p
# z_windows.reverse() # easier to start from higher z

z_window_edges = [ Z[1] for Z in z_windows]

#z_window_edges.reverse() # for matching in analysis/free energy computation

z_file = 'zvals.dat'

rp_windows = [(i/num_windows,(1+i)/num_windows) for i in range(num_windows) ]

rp_window_edges = [ rp[1] for rp in rp_windows]

rp_file = 'rpvals.dat'

rptp_file = 'rptpvals.dat'

tp_file = 'tpvals.dat'
