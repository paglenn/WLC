# get rp, tp, rptp values for a given z
import numpy as np
import os
#from parameters import *
import matplotlib.pyplot as plt
rp_file = 'rpvals.dat'
rptp_file = 'rptpvals.dat'
tp_file = 'tpvals.dat'
z_file = 'zvals.dat'

def calculate_TP():
    TP = []
    if not os.path.isfile(tp_file):
        print("data file ",tp_file,"missing!")
        exit()
    tpFile = open(tp_file,'r')
    for line in tpFile.readlines():
        TP.append(float(line[:-1]))

    tpFile.close()

    return TP

def calculate_RP():
    RP = [ ]
    if not os.path.isfile(rp_file):
        print("data file ",rp_file,"missing!")
        exit()
    rpFile = open(rp_file,'r')
    for line in rpFile.readlines():
        RP.append(float(line[:-1]))
    rpFile.close()
    return RP


def calculate_Z():
    Z = []
    if not os.path.isfile(z_file):
        print("data file ",z_file,"missing!")
        exit()
    zFile = open(z_file,'r')
    for line in zFile.readlines():
        Z.append(float(line[:-1]))
    zFile.close()
    return Z

def calculate_RPTP(RP,TP):
    RPTP = []
    if not os.path.isfile(rptp_file):
        print("data file ",rptp_file,"missing!")
        exit()
    rptpFile = open(rptp_file,'r')
    for line in rptpFile.readlines():
        RPTP.append(float(line[:-1]))
    rptpFile.close()
    return RPTP

def condition():
    TP = calculate_TP()
    RP = calculate_RP()
    RPTP = calculate_RPTP(RP,TP)

TP = calculate_TP()
mu = np.mean(TP)
var = np.var(TP)
plt.subplot(221)
binContents, bins = np.histogram(TP,bins=20,density=True)
gaus = np.sqrt(1/(2*np.pi*var)) * np.exp(-(bins[:-1]-mu)**2. / (2.*var) )
plt.semilogy(bins[:-1],binContents,marker='*',ls='',label='data')
plt.semilogy(bins[:-1],gaus,label='gaussian fit')
plt.title(r'$T_{\perp}$')
plt.ylabel(r'P($T_{\perp}$)')


RP = calculate_RP()
mu = np.mean(RP)
var = np.var(RP)
plt.subplot(222)
binContents, bins = np.histogram(RP,bins=20,density=True)
gaus = np.sqrt(1/(2*np.pi*var)) * np.exp(-(bins[:-1]-mu)**2. / (2.*var) )
plt.semilogy(bins[:-1],binContents,ls='.',marker='*',label='data')
plt.semilogy(bins[:-1],gaus,label='gaussian fit')
plt.title(r'$R_{\perp}$')
plt.ylabel('P($R_{\perp}$)')

RPTP = calculate_RPTP(RP,TP)
plt.subplot(223)
mu = np.mean(RPTP)
var = np.var(RPTP)
binContents, bins = np.histogram(RPTP,bins=20,density=True)
gaus = np.sqrt(1./(2.*np.pi*var)) * np.exp(-(bins[:-1]-mu)**2. / (2.*var) )
plt.semilogy(bins[:-1],binContents,ls='.',label='data',marker='*')
plt.semilogy(bins[:-1],gaus,label='gaussian fit')
plt.title(r'$R_{\perp}\cdot T_{\perp}$')
plt.ylabel(r'P($R_{\perp}\cdot T_{\perp}$)')

'''
Z = calculate_Z()
plt.subplot(224)
mu = np.mean(Z)
var = np.var(Z)
binContents, bins = np.histogram(Z,bins=20,density=True)
#gaus = np.sqrt(1/(2*np.pi*var)) * np.exp(-(bins[:-1]-mu)**2. / (var) )
plt.semilogy(bins[:-1],binContents,ls='.',label='data',marker='*')
#plt.semilogy(bins[:-1],gaus,label='gaussian fit')
plt.title('z/L')
'''
plt.suptitle('N=100,1e6 steps')

plt.savefig('dists.png')
print 'done!'



'''
z_X_joint = np.histogramdd(np.vstack((Z,TP,RP,RPTP)).T,normed=True)[0]
X_joint = np.histogramdd(np.vstack((TP,RP,RPTP)).T,normed=True)[0]
print(z_X_joint.shape,X_joint.shape)
#print(z_X_joint,X_joint)
#print(z_X_joint)
#z_cond = z_X_joint/X_joint
'''

