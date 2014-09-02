import numpy as np
import parameters as par

class actin_filament:

    def __init__(self,x):

        self.x = x
        self.t = np.copy(par.t0)
        self.t_old = np.copy(self.t)

    def perturb(self):

        self.t_old = np.copy(self.t )
        random_vec = np.random.randn(3)
        self.t = self.t +  random_vec
        self.t = self.t / np.linalg.norm(self.t)

    def revert(self):
        self.t = np.copy(self.t_old)

    def shift(self,x,t):
        self.x = x + t*par.dx

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
        E_old -= 1./par.dx * np.dot(filament[index+1].t,filament[index].t_old)
        E_new -= 1./par.dx * np.dot(filament[index+1].t,filament[index].t)

    if index != 0:
        E_old -= 1./par.dx * np.dot(filament[index-1].t,filament[index].t_old)
        E_new -= 1./par.dx * np.dot(filament[index-1].t,filament[index].t)


    return E_new - E_old

def calculate_r_vec(filament):
    r_vec = np.zeros(par.t0.shape)
    for actin in filament[:-1]:
        r_vec += par.dx * actin.t
    return r_vec

def calculate_z(filament):
    r = calculate_r_vec(filament)
    z_vec = np.dot(r,par.t0) * par.t0  # norm(t0) == 1 so no division needed
    z = np.linalg.norm(z_vec)
    return z

def deltaE_z_bias(filament,index,window,zFile):
    dE = deltaE(filament,index)

    z = calculate_z(filament)

    if z < min(window) or z > max(window):

        dE = np.inf

    else:
        zFile.write(str(z)+'\n')

    return dE

def calculate_rp(filament):

    r_vec = calculate_r_vec(filament)
    z_vec = np.dot(r_vec,par.t0) * par.t0
    rp_vec = r_vec - z_vec
    rp = np.linalg.norm(rp_vec)
    return rp

def deltaE_rp_bias(filament,index,window,rpFile):
    dE = deltaE(filament,index)

    rp = calculate_rp(filament)
    #print(rp)

    if rp < min(window) or rp > max(window):

        dE = np.inf

    else:
        rpFile.write(str(rp)+'\n')

    return dE

def adjust_z(filament,window):
    z = calculate_z(filament)
    window_mean = 0.5*sum(window)
    w = max(window) - min(window)
    target = window_mean + w/2. * np.random.uniform(-1,1)
    tol = w/4.

    while abs(z - target ) > tol :

        random_index= np.random.choice(range(par.num_actins))
        filament[random_index].perturb()
        #update(filament,random_index)
        z_new = calculate_z(filament)
        change = z_new - z

        if z > target and change > 0:
            filament[random_index].revert()
            #update(filament,random_index)
            change = 0
        if z < target and change < 0:
            filament[random_index].revert()
            #update(filament,random_index)
            change = 0

        if change != 0 :
            z = calculate_z(filament)

    return filament

def adjust_rp(filament,window):
    rp = calculate_rp(filament)
    window_mean = 0.5*sum(window)
    w = max(window) - min(window)
    target = window_mean + w/2. * np.random.uniform(-1,1)
    tol = w/4.

    while abs(rp - target ) > tol :
        #print('adjusting...',rp,'target: ',target)

        random_index= np.random.choice(range(par.num_actins))
        filament[random_index].perturb()
        #update(filament,random_index)
        rp_new = calculate_rp(filament)

        change = rp_new - rp


        if rp > target and change > 0:
            filament[random_index].revert()
            #update(filament,random_index)
            change = 0
        if rp < target and change < 0:
            filament[random_index].revert()
            #update(filament,random_index)
            change = 0

        if change != 0 :
            rp = calculate_rp(filament)

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


def mc_step(filament):

        random_index= np.random.choice(range(par.num_actins))
        #print(filament[random_index].t)

        filament[random_index].perturb()
        #print(filament[random_index].t)

        dE = delta_E(filament,random_index)

        if dE > 0:
            if np.random.rand() > np.exp(-dE):
                filament[random_index].revert()
                return filament

        #update(filament,random_index)
        return filament

def umbrella_mc_step(filament,window,outfile):

        random_index= np.random.choice(range(par.num_actins))
        #print(filament[random_index].t)

        filament[random_index].perturb()
        #print(filament[random_index].t)

        dE = deltaE_rp_bias(filament,random_index,window,outfile)

        if dE > 0:
            if np.random.rand() > np.exp(-dE):
                filament[random_index].revert()
                return filament

        #update(filament,random_index)
        return filament

#    return filament











