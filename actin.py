import numpy as np

class actin_filament:

    def __init__(self,x,dx):

        self.x = x
        self.t = np.array([0,0,1])
        self.t_old = np.copy(self.t)
        self.dx = dx

    def perturb(self,vec):
        self.t_old = np.copy(self.t )
        self.t = self.t +  vec
        self.t = self.t / np.linalg.norm(self.t)

    def revert(self):
        self.t = self.t_old

    def shift(self,x,t):
        self.x = x + t*self.dx

def update(A, index):

    for i in range(index+1,len(A)):

        A[i].shift(A[i-1].x,A[i-1].t)

def delta_E(A,index) :

    E_old = 0.
    E_new = 0.
    dx = A[0].dx

    if index + 1 != len(A):
        E_old -= 1./dx * np.dot(A[index+1].t,A[index].t_old)
        E_new -= 1./dx * np.dot(A[index+1].t,A[index].t)

    if index != 0:
        E_old -= 1./dx * np.dot(A[index-1].t,A[index].t_old)
        E_new -= 1./dx * np.dot(A[index-1].t,A[index].t)


    return E_new - E_old

def write(A,step,outfile):


    outfile.write('step ' + str(step) +'\n')
    for monomer in A:

        outfile.write('{0}\t{1}\t{2}\n'.format(monomer.t[0],monomer.t[1],monomer.t[2]) )

    outfile.write('\n')












