numsteps = int(1e5)
L = 0.2
dx = 0.002
Ldx = int(L/dx)
num_actins = Ldx + 1
data_file = 'output.dat'

# for umbrella sampling
umbrella_bias = False
num_windows = 10
z_windows = [(i*L/num_windows,(1+i)*L/num_windows) for i in range(num_windows) ]

