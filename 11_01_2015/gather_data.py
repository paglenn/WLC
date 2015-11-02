#!/usr/bin/env python
import os
import sys
import numpy as np

cwd = [ f for f in os.listdir('.') if 'test_' in f ]
fout = open('averages.dat','w')
fout.write('#RP\tTP\tRPTP\tZ_avg\tZmax\tZmp\n')
for direc in cwd:
    os.chdir(direc)
    if os.path.isfile('moments.dat') and os.path.isfile('progress.out'):

        # get Z
        fin = open('moments.dat','r')
        firstline  = fin.readlines()[0]
        if len(firstline) == 0 : continue
        fin.close()
        zbar = firstline.split()[1]
        zbar = float(zbar)


        #get RP, TP, RPTP
        fin = open('progress.out','r')
        lines = fin.readlines()
        firstline = lines[0].split()
        fin.close()
        rp_index = firstline.index('RP') + 2
        tp_index = firstline.index('TP') + 2
        rptp_index = firstline.index('RPTP') + 2
        rp =  firstline[rp_index]
        tp = firstline[tp_index]
        rptp = firstline[rptp_index]
        rp = round(float(rp),7)
        tp = round(float(tp),7)
        rptp = round(float(rptp),7)

        # read in max val. of z achieved
	fin = open('log.dat','r' )
	lines = fin.readlines()
	fin.close()
	lastline = lines[-1].split()
	zmax = float( lastline[-1] )

        # read in most probable z
        z_mp = find_most_probable()

        fout.write('%f\t%f\t%f\t%f\t%f\t\n'%(rp,tp,rptp,zbar,zmax,z_mp))
    #else :
        #print direc
    os.chdir('..')
