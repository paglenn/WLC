#!/usr/bin/env python

# stagger for different seed

import time
import os
import shutil

for i in range(10) :
    new_dir_name = "test_%d"%i
    if os.path.exists(new_dir_name) :
        shutil.rmtree(new_dir_name)
    os.mkdir(new_dir_name)
    shutil.copy('wlc',new_dir_name )
    os.chdir(new_dir_name)
    #os.system("./wlc")
    os.chdir('..')
    time.sleep(1)

