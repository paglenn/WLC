import numpy as np
from parameters import dx

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


def update(A, index):

    for i in range(index+1,len(A)):

        A[i].shift(A[i-1].x,A[i-1].t)

def deltaE(A,index) :

    E_old = 0.
    E_new = 0.

    if index + 1 != len(A):
        E_old -= 1./dx * np.dot(A[index+1].t,A[index].t_old)
        E_new -= 1./dx * np.dot(A[index+1].t,A[index].t)

    if index != 0:
        E_old -= 1./dx * np.dot(A[index-1].t,A[index].t_old)
        E_new -= 1./dx * np.dot(A[index-1].t,A[index].t)


    return E_new - E_old

def calculate_z(A0,A):
    z = 0
    N = len(A)
    for i in range(N):
        scalar_projection = np.dot(A[i].t,A0[i].t) / np.linalg.norm(A0[i].t)**2.
        proj = scalar_projection * A0[i].t
        z += np.linalg.norm(proj) * dx
    return z

def deltaE_z_bias(A0,A,index,window):
    dE = deltaE(A,index)

    z = calculate_z(A0,A)

    if z < min(window) or z > max(window):

        dE = np.inf

    return dE

def write_tangents(A,step,outfile):

    outfile.write('step ' + str(step) +'\n')
    for monomer in A:

        outfile.write('{0}\t{1}\t{2}\n'.format(monomer.t[0],monomer.t[1],monomer.t[2]) )

    outfile.write('\n')

def write_coordinates(A,step,outfile):

    outfile.write('step' + str(step) + '\n')

    for monomer in A:

        outfile.write('{0}\t{1}\t{2}\n'.format(monomer.x[0],monomer.x[1],monomer.x[2]) )

    outfile.write('\n')












