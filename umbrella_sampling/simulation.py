# simulation.py
# single actin filament (WLC) monte carlo routine
# umbrella sampling simulation adapted from that for Ising model
import numpy as np
from montecarlo import *
from parameters import *

# choose initial configuration to be along z-axis

#out = dict()
#for i in range(num_windows):
#    out[z_windows[i]] = open(window_files[i],'w')


# write initial configuration to file
# rpFile = open(rp_file,'w')
zFile = open(z_file,'w')
rpFile = open(rp_file,'w')
tpFile = open(tp_file,'w')
rptpFile = open(rptp_file,'w')
progressFile = open("progress.out",'w')
cosFile = open(cos_file,'w')

for w in z_windows:
    t = [np.copy(t0) for x in range(num_actins) ]
    ################################
    # Passes over each window
    for itr in range(num_passes):


        ##############################################
        # Adjust horizontal displacement to target value
        t = adjust_z(t,w)
        tp = np.linalg.norm(t[-1] - np.dot(t[-1],t0)*t0 )
        #t = adjust_rp(t,w)
        #write_tangents(t,0,out[w])

        zFile.write(    "{0}\n".format(calculate_z(t)/L  )  )
        rpFile.write(   "{0}\n".format(calculate_rp(t)/L )  )
        tpFile.write(   "{0}\n".format(tp)                  )
        rptpFile.write( "{0}\n".format(calculate_rptp(t) )  )
        cosFile.write("{0}\n".format(sum_cosines(t)      )  )

        for j in range(1,numsteps+1):

            umbrella_mc_step(t,w)

            progressFile.write('window\t{0}\tpass\t{1}\trun\t{2}\n'.format(w,itr,j))
            zFile.write(    "{0}\n".format(calculate_z(t)/L  )  )
            rpFile.write(   "{0}\n".format(calculate_rp(t)/L )  )
            tp = np.linalg.norm(t[-1] - np.dot(t[-1],t0)*t0     )
            tpFile.write(   "{0}\n".format(tp)                  )
            rptpFile.write( "{0}\n".format(calculate_rptp(t) )  )
            cosFile.write("{0}\n".format(sum_cosines(t)      )  )

            #write_tangents(t,j,out[w])

# close data files

'''
for dataFile in out.values():
    dataFile.write("program finished!!")
    dataFile.close()
'''

progressFile.write("program finished!!")
progressFile.close()
rpFile.close()
zFile.close()
rptpFile.close()
tpFile.close()
cosFile.close()
