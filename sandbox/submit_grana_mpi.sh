#!/bin/sh

# run this script by typing the following command,
# replacing $n with the actual number of processes:
# qsub -pe orte $n ./submit_mpi.sh
# $NSLOTS is initialized by -pe orte flag

#$ -S /bin/bash
#$ -N psi001-7575
#$ -q all.q
# alternatively, send to ibnet

# run from the current working directory
#$ -cwd
#$ -j n

# export environment variables
#$ -V

#echo This job is being run on $(hostname --short)

#echo Running $NSLOTS processes
date
# change this to your actual executable and arguments
#mpirun -np $NSLOTS ./GranaStack-orig/grana_monte_carlo/dist/Debug/GNU-Linux-x86/grana_monte_carlo  ./GranaStack-orig/grana_monte_carlo/test.cfg
#./GranaStack-orig/grana_monte_carlo/dist/Debug/GNU-Linux-x86/grana_monte_carlo  ./GranaStack-orig/grana_monte_carlo/test_psi001_7575.cfg
./wlc 

date
echo job completed
#echo job completed on $NSLOTS
