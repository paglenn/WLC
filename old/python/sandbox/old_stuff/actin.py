import numpy as np
import parameters as par

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

def calculate_r_vec(t):
    r_vec = np.zeros(t[0].shape)
    for j in range(par.num_actins):
        r_vec += t[j] * par.dx
    return r_vec

def calculate_z(t):
    r = calculate_r_vec(t)
    z_vec = np.dot(r,par.t0) * par.t0  # norm(t0) == 1 so no division needed
    z = np.linalg.norm(z_vec)
    return z

def deltaE_z_bias(t,t_old,index,window):
    dE = deltaE(t,t_old,index)

    z = calculate_z(t)

    if z < min(window) or z > max(window):

        dE = np.inf

    return dE

def calculate_rp(t):

    r_vec = calculate_r_vec(t)
    z_vec = np.dot(r_vec,par.t0) * par.t0
    rp_vec = r_vec - z_vec
    rp = np.linalg.norm(rp_vec)
    return rp

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
    target = window_mean + w/2. * np.random.uniform(-1,1)
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
            z = calculate_z(t)

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
            rp = calculate_rp(t)

    return t

def write_tangents(t,step,outfile):

    outfile.write('step ' + str(step) +'\n')
    for j in range(par.num_actins):

        outfile.write('{0}\t{1}\t{2}\n'.format(t[j][0],t[j][1],t[j][2]) )

    outfile.write('\n')

def mc_step(t):

        random_index= np.random.choice(range(par.num_actins))
        #print(t[random_index].t)

        # define random axis
        axis = np.random.randn(3)
        axis = axis / np.linalg.norm(axis)

        # other way: define by (theta,phi pair)
        phi = np.random.uniform(0,2*np.pi)
        r = np.random.rand()
        theta = np.acos(1 - 2*r)
        axis = [np.cos(phi)*np.sin(theta), np.
        pocket  = [ random_index ]
        cluster = [ random_index ]

        while pocket != [] :
                j = np.random.choice(pocket)
                for nn in par.nbr[j]:
                    nn_para = np.dot(t[nn],axis)
                    j_para = np.dot(t[j],axis)
                    if nn_para*j_para > 0 and nn not in cluster:
                        p_add = 1 - np.exp(-2./par.dx * nn_para*j_para)
                        if np.random.rand() < p_add:
                            pocket.append(nn)
                            cluster.append(nn)
                pocket.remove(j)
        for j in cluster:
            proj = axis * np.dot(t[j],axis)
            t[j] = t[j] - 2*(t[j] - proj)

        return t

'''
def mc_step(t):

        random_index= np.random.choice(range(par.num_actins))
        #print(t[random_index].t)

        random_vec = np.random.randn(3)
        pocket  = [ random_index ]
        cluster = [ random_index ]

        while pocket != [] :
                j = np.random.choice(pocket)
                for nn in par.nbr[j]:
                    cos_old = np.dot(t[j],t[nn])
                    t_new = t[j] + random_vec
                    t_new = t_new / np.linalg.norm(t_new)
                    cos_new = np.dot(t_new,t[nn])
                    if np.dot(t[j],t[nn]) > 0 and nn not in cluster:
                        p_add = 1 - np.exp(-1./par.dx * (cos_old - cos_new))
                        if np.random.rand() < p_add:
                            pocket.append(nn)
                            cluster.append(nn)
                pocket.remove(j)

        for j in cluster:
            #print(len(cluster))
            t[j] = (t[j] + random_vec) / np.linalg.norm(t[j] + random_vec)

        return t
'''
