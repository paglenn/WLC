# simulation.py
# Kratky-Porod model MCMC simulation driver
# umbrella sampling simulation adapted from that for Ising model
import numpy as np
from MC import *
from parameters import *
import os

# choose initial configuration to be along z-axis

# write initial configuration to file
#rpFile = open(rp_file,'w')
#tpFile = open(tp_file,'w')
#rptpFile = open(rptp_file,'w')
progressFile = open(progress_file,'w')
zFile = open(z_file,'w')
windowFiles = dict(enumerate(open(f,'w') for f in window_files))

progressFile.write("N = {}\n".format(num_actins))
progressFile.write("delta = {0:.3f}\n".format(dx) )
progressFile.write("{} windows:{}\n".format(num_windows,z_windows) )
progressFile.write("{} passes/window\n".format(num_passes))
progressFile.write("{} steps/pass\n".format(numsteps) )
progressFile.write("k = {}\n".format(K) )
progressFile.write("bins = {}\n".format(num_bins))
progressFile.write("bin overlap = {}\n".format(binOverlap))

windowSeq = range(num_windows)
passSeq = range(num_passes)
stepSeq = range(1,numsteps+1)
filaSeq = range(num_actins)
num_acc = 0.

for wi in windowSeq:
    np.random.seed(wi)
    w = z_windows[wi]
    wmin = Zmin[wi]
    t = [np.copy(t0) for x in filaSeq ]
    ################################
    # Passes over each window
    for wpass in passSeq:

        ##############################################
        # Adjust height to random target value within window
        adjust_z(t,wi)
        z = calculate_z(t)
        #zFile.write('%f\n'%z)
        #tp = np.linalg.norm(t[-1] - np.dot(t[-1],t0)*t0 )

        windowFiles[wi].write("{0}\t{1}\n".format(j,z )  )
        #rpFile.write(   "{0}\n".format(calculate_rp(t) )  )
        #tpFile.write(   "{0}\n".format(tp)                  )
        #rptpFile.write( "{0}\n".format(calculate_rptp(t) )  )

        for step in stepSeq:

            # mc update
            acc = umbrella_mc_step(t,wi)
            # acc = mc_step(t)

            ######################################
            # Write reports of data to file

            z = calculate_z(t)
            '''
            zFile.write('{} {}'.format(z,wi))
            for n in windowSeq:
                u_kln = 0.5*K*(z-Zmin[n])**2.
                zFile.write(' {}'.format(u_kln ) )
            zFile.write('\n')
            '''


            windowFiles[wi].write("{0}\t{1}\n".format(step,z)      )
            #rpFile.write(   "{0}\n".format(calculate_rp(t) )  )
            #tp = np.linalg.norm(t[-1] - np.dot(t[-1],t0)*t0     )
            #tpFile.write(   "{0}\n".format(tp)                  )
            #rptpFile.write( "{0}\n".format(calculate_rptp(t) )  )
            ######################################
            if acc:
                num_acc += 1

        # write simulation progress
        # question: why doesn't this print after every pass? instead it bunches
        progressVars = [wi+1,num_windows,wpass+1,num_passes,numsteps]
        progressFile.write('window\t{0}/{1}\tpass\t{2}/{3}\tstep\t{4}\n'.format(*progressVars))

metaFile = open(metadata_file,'w')
metaFile.write('#window_file\tz_min\tk\n')
for j in windowSeq:
    data = [os.path.abspath(window_files[j]),Zmin[j],K ]
    metaFile.write('{0}\t{1}\t{2}\n'.format(*data))
metaFile.close()

# close data files
#rpFile.close()
#rptpFile.close()
#tpFile.close()
zFile.close()
for f in windowFiles:
    windowFiles[f].close()

# finally, progress file
progressFile.write("program finished!!\n")
progressFile.write("{0:.2%}\taccepted moves\n".format(num_acc/total_frames))
progressFile.close()
