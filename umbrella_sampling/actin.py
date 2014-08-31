import numpy as np
from parameters import dx,num_actins,t0

class actin_filament:

    def __init__(self,x,dx):

        self.x = x
        self.t = np.array([0,0,1])
        self.t_old = np.copy(self.t)
        dx = dx

    def perturb(self,vec):
        self.t_old = np.copy(self.t )
        self.t = self.t +  vec
        self.t = self.t / np.linalg.norm(self.t)

    def revert(self):
        self.t = self.t_old

    def shift(self,x,t):
        self.x = x + t*dx

    def rotate(self,phi):

        R = np.array([[1,0,0],[0,np.cos(phi),-np.sin(phi)],[0,np.sin(phi),np.cos(phi)]])
        self.t = np.dot(R,self.t)


def update(filament, index):

    for i in range(index+1,len(filament)):

        filament[i].shift(filament[i-1].x,filament[i-1].t)

def deltaE(filament,index) :

    E_old = 0.
    E_new = 0.

    if index + 1 != len(filament):
        E_old -= 1./dx * np.dot(filament[index+1].t,filament[index].t_old)
        E_new -= 1./dx * np.dot(filament[index+1].t,filament[index].t)

    if index != 0:
        E_old -= 1./dx * np.dot(filament[index-1].t,filament[index].t_old)
        E_new -= 1./dx * np.dot(filament[index-1].t,filament[index].t)


    return E_new - E_old

def calculate_z(A0,filament):
    z = 0
    N = len(filament)
    for i in range(N):
        scalar_projection = np.dot(filament[i].t,A0[i].t) / np.linalg.norm(A0[i].t)**2.
        proj = scalar_projection * A0[i].t
        z += proj * dx

    return np.linalg.norm(z)

def calculate_z(A0,filament):
    R = np.zeros(filament[0].t.shape)
    for monomer in filament[:-1]:
        R += dx*monomer.t
    z_vec = np.dot(R,t0) * t0  # norm(t0) == 1 so no division needed
    z = np.linalg.norm(z_vec)
    return z

def deltaE_z_bias(A0,filament,index,window,zFile):
    dE = deltaE(filament,index)

    z = calculate_z(A0,filament)

    if z < min(window) or z > max(window):

        dE = np.inf

    else:
        zFile.write(str(z)+'\n')

    return dE

def adjust_z(A0,filament,window):
    z = calculate_z(A0,filament)
    window_mean = 0.5*sum(window)
    width = max(window) - min(window)
    target = window_mean + width/2 * np.random.uniform(-1,1)
    tol = width / 4.

    while abs(z - target ) > tol :

        random_index= np.random.choice(range(num_actins))
        random_vec = np.random.randn(3)
        old_proj = np.dot(filament[random_index].t,t0)
        filament[random_index].perturb(random_vec)
        new_proj = np.dot(filament[random_index].t,t0)
        change = new_proj - old_proj

        if z > target and change > 0:
            filament[random_index].revert()
            change = 0
        if z < target and change < 0:
            filament[random_index].revert()
            change = 0

        if change != 0 :
            update(filament,random_index)
            z = calculate_z(A0,filament)

    return filament

def write_tangents(filament,step,outfile):

    outfile.write('step ' + str(step) +'\n')
    for actin in filament:

        outfile.write('{0}\t{1}\t{2}\n'.format(actin.t[0],actin.t[1],actin.t[2]) )

    outfile.write('\n')

def write_coordinates(filament,step,outfile):

    outfile.write('step' + str(step) + '\n')

    for actin in filament:

        outfile.write('{0}\t{1}\t{2}\n'.format(actin.x[0],actin.x[1],actin.x[2]) )

    outfile.write('\n')












