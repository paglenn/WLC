import numpy as np
import parameters as par
from numba import autojit

def perturb(t,j):

    t_old = np.copy(t[j] )
    random_vec = np.random.randn(3)
    t[j] = t[j] +  random_vec
    t[j] = t[j] / np.linalg.norm(t[j])

    return t, t_old


def deltaE(t,t_old,index) :

    E_old = 0.
    E_new = 0.

    if index + 1 != par.num_actins:
        E_old -= 1./par.dx * np.dot(t[index+1],t_old)
        E_new -= 1./par.dx * np.dot(t[index+1],t[index])

    if index != 0:
        E_old -= 1./par.dx * np.dot(t[index-1],t_old)
        E_new -= 1./par.dx * np.dot(t[index-1],t[index])


    return E_new - E_old

@autojit
def calculate_r_vec(t):
    r_vec = np.zeros(t[0].shape)
    for j in range(par.num_actins):
        r_vec = r_vec +  t[j] * par.dx
    return r_vec

@autojit
def sum_cosines(t):

    s = 0.
    for j in range(par.num_actins-1):

        s += np.dot(t[j],t[j+1])

    return s

@autojit
def calculate_z(t):
    r = calculate_r_vec(t)
    z_vec = np.dot(r,par.t0) * par.t0  # norm(t0) == 1 so no division needed
    z = np.linalg.norm(z_vec)
    return z/par.L

def deltaE_z_bias(t,t_old,index,window):
    dE = deltaE(t,t_old,index)

    z = calculate_z(t)

    if z < min(window) or z > max(window):

        dE = np.inf

    return dE

@autojit
def calculate_rp(t):

    r_vec = calculate_r_vec(t)
    z_vec = np.dot(r_vec,par.t0) * par.t0
    rp_vec = r_vec - z_vec
    rp = np.linalg.norm(rp_vec)
    return rp/par.L

@autojit
def calculate_rptp(t):

    r_vec = calculate_r_vec(t)
    z_vec = np.dot(r_vec,par.t0) * par.t0
    rp_vec = r_vec - z_vec
    tN = t[-1]
    tp = tN - np.dot(tN, par.t0) * par.t0
    rptp = np.dot(tp,rp_vec)
    return rptp


def deltaE_rp_bias(t,t_old,index,window):
    dE = deltaE(t,t_old,index)

    rp = calculate_rp(t)

    if rp < min(window) or rp > max(window):

        dE = np.inf

    return dE

def adjust_z(t,window):
    z = calculate_z(t)
    window_mean = 0.5*sum(window)
    w = max(window) - min(window)
    target = window_mean + w/4. * np.random.uniform(-1,1)
    tol = w/4.

    while abs(z - target ) > tol :

        random_index= np.random.choice(range(par.num_actins))
        t,t_old = perturb(t,random_index)
        #update(t,random_index)
        z_new = calculate_z(t)
        change = z_new - z

        if z > target and change > 0:
            t[random_index] = t_old
            #update(t,random_index)
            change = 0
        if z < target and change < 0:
            t[random_index] = t_old
            #update(t,random_index)
            change = 0

        if change != 0 :
            z = z_new

    return t

def adjust_rp(t,window):
    rp = calculate_rp(t)
    window_mean = 0.5*sum(window)
    w = max(window) - min(window)
    target = window_mean + w/2. * np.random.uniform(-1,1)
    tol = w/4.

    while abs(rp - target ) > tol :
        #print('adjusting...',rp,'target: ',target)

        random_index= np.random.choice(range(par.num_actins))
        t,t_old = perturb(t,random_index)
        #update(t,random_index)
        rp_new = calculate_rp(t)

        change = rp_new - rp


        if rp > target and change > 0:
            t[random_index] = t_old
            #update(t,random_index)
            change = 0
        if rp < target and change < 0:
            t[random_index] = t_old
            #update(t,random_index)
            change = 0

        if change != 0 :
            rp = rp_new

    return t

def write_tangents(t,step,outfile):

    outfile.write('step ' + str(step) +'\n')
    for j in range(par.num_actins):

        outfile.write('{0}\t{1}\t{2}\n'.format(t[j][0],t[j][1],t[j][2]) )

    outfile.write('\n')

def mc_step(t):

    random_index= np.random.choice(range(par.num_actins))
    #print(t[random_index].t)

    t,t_old = perturb(t,random_index)
    #print(t[random_index].t)

    dE = delta_E(t,random_index)

    if dE > 0:
        if np.random.rand() > np.exp(-dE):
            t[random_index] = t_old
            return

    #update(t,random_index)
    return

def umbrella_mc_step(t,window):

    random_index= np.random.choice(range(par.num_actins))
    #print(t[random_index].t)

    t,t_old = perturb(t,random_index)
    #print(t[random_index].t)

    dE = deltaE_z_bias(t,t_old,random_index,window)

    if dE > 0:
        if np.random.rand() > np.exp(-dE):
            t[random_index] = t_old
            return False


    #update(t,random_index)

    return True















