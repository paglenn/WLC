#!/usr/bin/env python

# stagger for different seed

import time
import os
import shutil

numCondBins = 10
ranges = [ [0,1,2], [3,4,5] , [6,7,8], [9] ]
ranges2 = range(10) 
for i in ranges2 :
	for j in range(numCondBins):
		for k in range(numCondBins) :
			t1 = i / float(numCondBins)
			t2 = j / float(numCondBins)
			t3 = k / float(numCondBins) 

			new_dir_name = "test_%d_%d_%d"%(i,j,k)
			if os.path.exists(new_dir_name) :
				shutil.rmtree(new_dir_name)
			os.mkdir(new_dir_name)
			shutil.copy('wlc',new_dir_name )
			shutil.copy('parameters.txt',new_dir_name )
			shutil.copy('submit_wlc.sh',new_dir_name )
			os.chdir(new_dir_name)
			os.system("sed 's@WORKDIR = /newhome/paulglen/test@WORKDIR=%s@' -i submit_wlc.sh"%os.getcwd())
			os.system("sed 's@wlc\>@wlc %f %f %f@g' -i submit_wlc.sh"%(t1,t2,t3))
			os.system("qsub -N %s -pe orte 1 ./submit_wlc.sh "%new_dir_name)
			#os.system("./wlc") 
			os.chdir('..')
			time.sleep(1)

