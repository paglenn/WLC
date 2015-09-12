#!/usr/bin/env python

# stagger for different seed

import time
import os
import shutil

numCondBins = 1
for i in range(numCondBins) :
	for j in range(numCondBins):
		for k in range(2*numCondBins) :
			t1 = i / float(numCondBins)
			t2 = j / float(numCondBins)
			t3 = (k- numCondBins) / float(numCondBins) * t1 * t2

			new_dir_name = "test_%d_%d_%d"%(i,j,k)
			if os.path.exists(new_dir_name) :
				shutil.rmtree(new_dir_name)
			os.mkdir(new_dir_name)
			shutil.copy('wlc',new_dir_name )
			shutil.copy('parameters.txt',new_dir_name )
			shutil.copy('submit_wlc.sh',new_dir_name )
			os.chdir(new_dir_name)
			os.system("sed 's@WORKDIR = /newhome/paulglen/test@WORKDIR=%s@' -i submit_wlc.sh"%os.getcwd()	)
			os.system("sed 's@wlc\>@wlc %f %f %f@g' -i submit_wlc.sh"%(t1,t2,t3))
			#os.system("qsub -pe orte 1 ./submit_wlc.sh")
			os.system("./wlc") 
			os.chdir('..')
			time.sleep(1)

