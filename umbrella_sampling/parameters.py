numsteps = int(1e3)
L = 0.5
dx = 0.01
Ldx = int(L/dx)
num_actins = Ldx + 1
data_file = 'output.dat'
import numpy as np
t0 = np.array([ 0,0,1])

# for umbrella sampling
umbrella_bias = True
num_windows = 10
num_passes = 100
total_runs = num_windows * num_passes * (numsteps+1)
z_windows = [(i*L/num_windows,(1+i)*L/num_windows) for i in range(num_windows) ]
z_windows.reverse()
window_files = ["window_{0}.dat".format(i) for i in range(num_windows) ]
window_edges = [ z[1] for z in z_windows]
window_edges.reverse()
z_file = 'heights.dat'

