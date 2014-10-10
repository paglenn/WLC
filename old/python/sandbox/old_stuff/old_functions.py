
''' # More robust but not necessary for now
def calculate_height_dist() :
    init_filament = trajectory[0]
    Z = []
    for j in range(len(trajectory)):
        filament = trajectory[j]
        z = np.zeros(filament[0].shape)
        for i in range(len(filament)-1):
            scalar_projection = np.dot(filament[i],init_filament[i]) / np.linalg.norm(init_filament[i])**2.
            _t = scalar_projection * init_filament[i]
            z = z +  _t
        Z.append(dx*np.linalg.norm(z))
    return Z
'''

# calculate net horizontal displacement and net orientation
'''
def calculate_displacement():

    init_filament = trajectory[0]
    final_filament = trajectory[-1]
    _R = z = 0
    for i in range(len(final_filament)):
        scalar_projection = np.dot(final_filament[i],init_filament[i])/np.linalg.norm(init_filament[i])**2.
        z += scalar_projection * initial_filament[i]
        dR = final_filament[i] - scalar_projection*init_filament[i]
        if np.linalg.norm(dR) > 0:
            dR = dR / np.linalg.norm(dR)
        _R += dR * dx

    return z, _R
'''
