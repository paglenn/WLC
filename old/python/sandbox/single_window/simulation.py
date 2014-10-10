# simulation.py
# single actin filament (WLC) monte carlo routine
# umbrella sampling simulation adapted from that for Ising model
import numpy as np
from montecarlo import *
from parameters import *

# choose initial configuration to be along z-axis

# write initial configuration to file
# rpFile = open(rp_file,'w')
zFile = open(z_file,'w')
rpFile = open(rp_file,'w')
tpFile = open(tp_file,'w')
rptpFile = open(rptp_file,'w')
progressFile = open("progress.out",'w')
cosFile = open(cos_file,'w')

num_acc = 0
t = [np.copy(t0) for x in range(num_actins) ]
for w in z_windows:
    ################################
    # Passes over each window
    for itr in range(num_passes):
        t = [np.copy(t0) for x in range(num_actins) ]

        ##############################################
        # Adjust horizontal displacement to target value
        #t = adjust_z(t,w)
        tp = np.linalg.norm(t[-1] - np.dot(t[-1],t0)*t0 )
        #t = adjust_rp(t,w)

        zFile.write(    "{0}\n".format(calculate_z(t)  )  )
        rpFile.write(   "{0}\n".format(calculate_rp(t) )  )
        tpFile.write(   "{0}\n".format(tp)                  )
        rptpFile.write( "{0}\n".format(calculate_rptp(t) )  )
        # write initial configuration
        for k in range(num_actins- 1):
            cosFile.write("{0:.8f}\n".format( np.dot(t[k],t[k+1]) ) )


        for j in range(1,numsteps+1):

            acc = mc_step(t)
            # acc = umbrella_mc_step(t,w)

            progressVars = [z_windows.index(w)+1,num_windows,itr,j]
            progressFile.write('window\t{0}/{1}\tpass\t{2}\trun\t{3}\n'.format(*progressVars))
            zFile.write(    "{0}\n".format(calculate_z(t)  )  )
            rpFile.write(   "{0}\n".format(calculate_rp(t) )  )
            tp = np.linalg.norm(t[-1] - np.dot(t[-1],t0)*t0     )
            tpFile.write(   "{0}\n".format(tp)                  )
            rptpFile.write( "{0}\n".format(calculate_rptp(t) )  )
            # and config at each step
            for k in range(num_actins- 1):
                cosFile.write("{0:.8f}\n".format( np.dot(t[k],t[k+1]) ) )

            if acc:
                num_acc += 1



# close data files
progressFile.write("program finished!!\n")
progressFile.write("{0}\taccepted moves\n".format(num_acc))
progressFile.close()
rpFile.close()
zFile.close()
rptpFile.close()
tpFile.close()
cosFile.close()
