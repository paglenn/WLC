# simulation parameters
# note: l_p is the persistence length (~37 um for actin rod )

numsteps = int(600)
#dx = 0.0027 # (units of l_p -- from "Actin Based Motility" )
L = 4./15  # from "limits of filopodium stability"
num_actins = N = 10 # N == 100 and L == 4./15 corresponds well with dx above
dx = L/float(N)
gaus_var = 5*dx # variance in random step

# filament basis vector
import numpy as np
t0 = np.array([0,0,1])

nbrs = dict()
for j in range(num_actins):

    if j == 0: nbrs[j] = [j+1]
    elif j+1 == num_actins : nbrs[j] = [j-1]
    else: nbrs[j] = [j+1,j-1]


# for umbrella sampling
num_passes = 500 # passes per window
num_windows = 10
num_bins = 100 # each window will have num_bins/num_windows+1 bins
binWidth = 1./ num_bins
binOverlap = 5*binWidth
z_windows = [[1.*i/num_windows,(i+1)/float(num_windows)] for i in range(num_windows)]
for w in range(num_windows-1): z_windows[w][1] += binOverlap
linesPerWindow = num_passes*(numsteps+1)
K = 1./dx # umbrella sampling harmonic force constant
binsPerWindow= num_bins/float(num_windows)+binOverlap/binWidth
# determine z windows
Zmin = [0.5*sum(w) for w in z_windows ]
total_frames = num_windows * num_passes * numsteps
z_window_edges = [ Z[1] for Z in z_windows]

# files
window_files = ["window_{0}.dat".format(i) for i in range(num_windows) ]
z_file = 'zvals.dat'
cos_file = 'cos.dat'
rp_file = 'rpvals.dat'
rptp_file = 'rptpvals.dat'
tp_file = 'tpvals.dat'
metadata_file = 'metadata.dat'

