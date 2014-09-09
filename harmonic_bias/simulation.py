# simulation.py
# single actin filament (WLC) monte carlo routine
# umbrella sampling simulation adapted from that for Ising model
import numpy as np
from montecarlo import *
from parameters import *
import os

# choose initial configuration to be along z-axis

# write initial configuration to file
rpFile = open(rp_file,'w')
tpFile = open(tp_file,'w')
rptpFile = open(rptp_file,'w')
progressFile = open("progress.out",'w')
zFile = open(z_file,'w')
#cosFile = open(cos_file,'w')
windowFiles = dict(enumerate(open(f,'w') for f in window_files))

num_acc = 0
t = [np.copy(t0) for x in range(num_actins) ]
for wi in range(num_windows):
    w = z_windows[wi]
    wmin = window_min[wi]
    ################################
    # Passes over each window
    for itr in range(num_passes):
        t = [np.copy(t0) for x in range(num_actins) ]

        ##############################################
        # Adjust horizontal displacement to target value
        t = adjust_z(t,w)
        tp = np.linalg.norm(t[-1] - np.dot(t[-1],t0)*t0 )

        z = calculate_z(t)
        zFile.write('%f\n'%z)
        windowFiles[wi].write("{0}\t{1}\n".format(j,z )  )
        rpFile.write(   "{0}\n".format(calculate_rp(t) )  )
        tpFile.write(   "{0}\n".format(tp)                  )
        rptpFile.write( "{0}\n".format(calculate_rptp(t) )  )

        for j in range(1,numsteps+1):

            # acc = mc_step(t)

            acc = umbrella_mc_step(t,wi)

            progressVars = [wi+1,num_windows,itr,j]
            progressFile.write('window\t{0}/{1}\tpass\t{2}\trun\t{3}\n'.format(*progressVars))

            z = calculate_z(t)
            zFile.write('%f\n'%z)
            windowFiles[wi].write("{0}\t{1}\n".format(j,z)      )
            rpFile.write(   "{0}\n".format(calculate_rp(t) )  )
            tp = np.linalg.norm(t[-1] - np.dot(t[-1],t0)*t0     )
            tpFile.write(   "{0}\n".format(tp)                  )
            rptpFile.write( "{0}\n".format(calculate_rptp(t) )  )
            #########################
            if acc:
                num_acc += 1



# close data files
rpFile.close()
rptpFile.close()
tpFile.close()
zFile.close()
for f in windowFiles:
    windowFiles[f].close()

mdFile = open('metadata.dat','w')
mdFile.write('#window_file\tz_min\tk\n')
for j in range(num_windows):
    data = [os.path.abspath(window_files[j]),window_min[j],K ]
    mdFile.write('{0}\t{1}\t{2}\n'.format(*data))
mdFile.close()

# finally, progress file
progressFile.write("program finished!!\n")
progressFile.write("{0}\taccepted moves\n".format(num_acc))
progressFile.close()
