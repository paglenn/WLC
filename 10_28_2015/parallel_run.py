#!/usr/bin/env python

# stagger for different seed

import time
import os
import shutil
from multiprocessing import Pool
from itertools import product

numCondBins = 5
ranges = [ [1,2,3], [4,5,6] , [7,8,9], ]
theBins = range(numCondBins)
rptpBins = range(2*numCondBins)

thistime = 1
triplets = list( product(ranges[thistime], theBins, rptpBins))

def f(tpl) :
    i = tpl[0]
    j = tpl[1]
    k = tpl[2]
    new_dir_name = "test_%d_%d_%d"%(i,j,k)
    if os.path.exists(new_dir_name) :
            shutil.rmtree(new_dir_name)

    # make directory for output to prevent overwriiting
    os.mkdir(new_dir_name)
    shutil.copy('wlc',new_dir_name )
    shutil.copy('parameters.txt',new_dir_name )
    os.chdir(new_dir_name)
    t1 = i / float(numCondBins)
    t2 = j / float(numCondBins)
    t3 = (k- numCondBins) / float(numCondBins) * t1 * t2
    os.system("./wlc %f %f %f"%(t1, t2, t3))
    os.chdir('..')
    # change back to cwd
    '''
    new_dir_name = "test_%d_%d_%d"%(i,j,k)
    if os.path.exists(new_dir_name) :
            shutil.rmtree(new_dir_name)
    os.mkdir(new_dir_name)
    shutil.copy('wlc',new_dir_name )
    shutil.copy('parameters.txt',new_dir_name )
    shutil.copy('submit_wlc.sh',new_dir_name )
    os.chdir(new_dir_name)
    #os.system("sed 's@WORKDIR = /newhome/paulglen/test@WORKDIR=%s@' -i submit_wlc.sh"%os.getcwd())
    #os.system("sed 's@wlc\>@wlc %f %f %f@g' -i submit_wlc.sh"%(t1,t2,t3))
    #os.system("qsub -N %s -pe orte 1 ./submit_wlc.sh "%new_dir_name)
    os.system("./wlc %f %f %f"%(t1, t2, t3))
    os.chdir('..')
    time.sleep(1)
    '''

if __name__ == '__main__':

    pool = Pool(processes = 3* 2 * numCondBins **2)
    pool.map(f, triplets)
    #parsed =

