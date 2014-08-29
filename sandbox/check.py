import numpy as np
from parameters import *
infile = open('output.dat','r')


trajectory = [ list() for x in range(num_steps+1) ]
step = -1
for line in infile.readlines():
    line2 = line[:-1].split('\t')

    if 'step' in line2[0]:
        step += 1
        continue
    elif len(line2) == 1:
        continue

    if step > num_steps: continue

    data = [float(x) for x in line2]
    trajectory[step].append(np.array(data))

for filament in trajectory:
    for actin in filament:
        if np.linalg.norm(actin) -1.0 > 1e-6 :
            print(np.linalg.norm(actin),"uh-oh")
            quit()


# compute P( cos \theta  )
import matplotlib.pyplot as plt
cos = []
for filament in trajectory:

    for i in range(num_actins - 1) :

        cos.append( np.dot(filament[i],filament[i+1]) )

x = np.linspace(-1,1,100)
P = 1./dx * np.exp(x/dx) / (2*np.sinh(1./dx) )
plt.subplot(221)
plt.hist(cos,bins=100,range=(-1,1),normed=True,histtype='step')
plt.plot(x,P,'k',lw=1.2)
plt.xlabel(r'$\cos \theta $')
plt.ylabel(r'$P(\cos \theta ) $')
plt.subplot(222)
plt.hist(cos,bins=100,range=(-1,1),normed=True,histtype='step',log=True)
plt.semilogy(x,P,'k',lw=1.2)
plt.xlabel('log scale')

# compute correlator
R = range(0,20)
G = [0 for r in R ]
num_samples = list(G)
for filament in trajectory[1:]:

    for r in R:

        for i in range(num_actins - r):

            num_samples[r] += 1
            G[r] += np.dot(filament[i],filament[i+r])


for r in R:

    if num_samples[r] != 0:

        G[r] = G[r] / num_samples[r]

# rescale R
R = [dx * r for r in R ]
R = np.array(R)
G = np.array(G)


# model prediction: G(s) ~ exp(-s/l_p)
expR = np.exp(-R)
plt.subplot(223)
plt.plot(np.array(R),G,'b.')
plt.plot(R, expR,'k',lw=1.5 )
plt.xlabel(r's/$l_p$')
plt.ylabel('G(s)')
#plt.savefig('correlator.png')
plt.subplot(224)
plt.semilogy(np.array(R),G,'b.')
plt.semilogy(R, expR,'k',lw=1.5 )
plt.xlabel('(semilog scale)')
#plt.ylabel('G(s) (log scale)')

plt.savefig('results.png')
#plt.show()
#plot_correlator('log')


