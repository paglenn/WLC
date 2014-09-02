import numpy as np
import os
from parameters import *

def calculate_TP():
    if not os.path.isfile(tp_file):
        print("data file ",tp_file,"missing!")
        exit()
    tpFile = open(tp_file,'w')
    for line in tpFile.readlines():
        tp.append(float(line[:-1]))

    tpFile.close()

    return tp

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

TP = calculate_TP()
RP = calculate_RP()
Z = calculate_Z()
RPTP = calculate_RPTP(RP,TP)

