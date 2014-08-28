import numpy as np
infile = open('output.dat','r')

num_actins = 201
num_steps = int(1e5)
dx = 0.1

trajectory = [ list() for x in range(num_steps+1) ]
step = -1
for line in infile.readlines():
    line2 = line[:-1].split('\t')

    if 'step' in line2[0]:
        step += 1
        continue
    elif len(line2) == 1:
        continue

    data = [float(x) for x in line2]
    trajectory[step].append(np.array(data))

# compute P( cos \theta  )
import matplotlib.pyplot as plt
def plot_distribution():
    cos = []
    for filament in trajectory:

        for i in range(num_actins - 1) :

            cos.append( np.dot(filament[i],filament[i+1]) )

    plt.hist(cos,bins=100,range=(-1,1),normed=True,histtype='step')
    x = np.linspace(-1,1,100)
    P = 1./dx * np.exp(x/dx) / (2*np.sinh(1./dx) )
    plt.plot(x,P,'k',lw=1.2)
    plt.xlabel(r'$\cos \theta $')
    plt.ylabel(r'$P(\cos \theta ) $')
    plt.savefig('cos_hist.png')
    plt.show()

# compute correlator
def plot_correlator(mode='normal'):
    R = range(0,20)
    G = [0 for r in R ]
    num_samples = list(G)
    for filament in trajectory[1:]:

        for r in R:

            for i in range(num_actins - r):

                num_samples[r] += 1
                G[r] += np.dot(filament[i],filament[i+r])


    for r in R:

        G[r] = G[r] / num_samples[r]

    # rescale R
    R = [dx * r for r in R ]
    R = np.array(R)
    G = np.array(G)


    # model prediction: G(s) ~ exp(-s/l_p)
    expR = np.exp(-R)
    if mode == 'normal':
        plt.plot(np.array(R),G,'b.')
        plt.plot(R, expR,'k',lw=1.5 )
        plt.savefig('correlator.png')
    elif mode == 'log':
        plt.semilogy(np.array(R),G,'b.')
        plt.semilogy(R, expR,'k',lw=1.5 )
        plt.savefig('correlator_log.png')

    plt.xlabel(r's/$l_p$')
    plt.ylabel('G(s)')
    plt.show()

#plot_distribution()
plot_correlator()
#plot_correlator('log')

