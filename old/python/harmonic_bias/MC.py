import numpy as np
import parameters as par
import math
import random
from numba import autojit

@autojit
def perturb(t,j):

    t_old = np.copy(t[j] )
    random_vec = par.gaus_var * np.random.randn(3)
    t[j] = t[j] +  random_vec
    t[j] = t[j] / np.linalg.norm(t[j])

    return t, t_old

def deltaE(t,t_old,index) :

    E_old = E_new = 0.
    q = 1./par.dx

    for nn in par.nbrs[index]:
        E_old -= q * np.dot(t[nn],t_old)
        E_new -= q * np.dot(t[nn],t[index])

    return E_new - E_old

@autojit
def calculate_r_vec(t):
    r_vec = np.zeros(t[0].shape)
    for j in range(par.num_actins):
        r_vec = r_vec +  t[j] * par.dx
    return r_vec


def sum_cosines(t):

    s = 0.
    for j in range(par.num_actins-1):

        s += np.dot(t[j],t[j+1])

    return s

def calculate_z(t):
    r = calculate_r_vec(t)
    z_vec = np.dot(r,par.t0) * par.t0  # norm(t0) == 1 so no division needed
    z = np.linalg.norm(z_vec)
    return z/par.L

def deltaE_bias(t,t_old,index,zmin):

    dE = deltaE(t,t_old,index)

    z = calculate_z(t)
    dE += 0.5 * par.K * (z - zmin)**2.

    return dE

def calculate_rp(t):

    r_vec = calculate_r_vec(t)
    z_vec = np.dot(r_vec,par.t0) * par.t0
    rp_vec = r_vec - z_vec
    rp = np.linalg.norm(rp_vec)
    return rp/par.L

def calculate_rptp(t):

    r_vec = calculate_r_vec(t)
    z_vec = np.dot(r_vec,par.t0) * par.t0
    rp_vec = r_vec - z_vec
    tN = t[-1]
    tp = tN - np.dot(tN, par.t0) * par.t0
    rptp = np.dot(tp,rp_vec)
    return rptp

def adjust_z(t,w_index):
    window = par.z_windows[w_index]
    z = 0.
    window_mean = 0.5*sum(window)
    w = window[1] - window[0]
    target = window_mean + w/2. * np.random.uniform(-1,1)
    tol = w/4.

    while abs(z - target ) > tol :

        random_index= np.random.choice(range(1,par.num_actins))
        t,t_old = perturb(t,random_index)
        #update(t,random_index)
        z_new = calculate_z(t)
        dz = z_new - z

        if z > target and dz > 0:
            t[random_index] = t_old
            #update(t,random_index)
            dz = 0
        if z < target and dz < 0:
            t[random_index] = t_old
            #update(t,random_index)
            dz = 0

        if dz != 0 :
            z = z_new

    return t

def write_tangents(t,step,outfile):

    outfile.write('step ' + str(step) +'\n')
    for j in range(par.num_actins):

        outfile.write('{0}\t{1}\t{2}\n'.format(t[j][0],t[j][1],t[j][2]) )

    outfile.write('\n')

def mc_step(t):

    random_index= np.random.choice(range(1,par.num_actins))
    #print(t[random_index].t)

    t,t_old = perturb(t,random_index)
    #print(t[random_index].t)

    dE = delta_E(t,random_index)

    if dE > 0:
        if np.random.rand() > np.exp(-dE):
            t[random_index] = t_old
            return False

    #update(t,random_index)
    return True

def umbrella_mc_step(t,w_index):

    random_index = np.random.choice(range(1,par.num_actins))

    t,t_old = perturb(t,random_index)

    dE = deltaE_bias(t,t_old,random_index,par.Zmin[w_index])

    if dE > 0 and random.uniform(0,1) > math.exp(-dE):
        t[random_index] = t_old
        return False

    t[0] = np.copy(par.t0)

    return True

def plot_positions(t):

    X,Y,Z = 3*[[0]]
    for x in range(1,len(t)):
        X.append(X[x-1] + t[x-1][0])
        Y.append(Y[x-1] + t[x-1][1])
        Z.append(Z[x-1] + t[x-1][2])

    return X,Y,Z

