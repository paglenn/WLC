#!/bin/sh
## run this script by typing the following command,
## replacing $n with the actual number of processes:
## qsub -pe orte $n ./submit_mpi.sh
## $NSLOTS is initialized by -pe orte flag
##echo "error-free here" >  &2 

##$ -S /bin/sh
#$ -N test
#$ -q ibnet
### alternatively, send to ibnet

### run from the current working directory
#$ -cwd
#$ -j n

# export environment variables
#$ -V

echo This job is being run on $(hostname --short)
echo Running $NSLOTS processes

###echo "Environment:" 
### echo $ENV
date
# change this to your actual executable and arguments
#mpirun -np $NSLOTS ./GranaStack-orig/grana_monte_carlo/dist/Debug/GNU-Linux-x86/grana_monte_carlo  ./GranaStack-orig/grana_monte_carlo/test.cfg
#./GranaStack-orig/grana_monte_carlo/dist/Debug/GNU-Linux-x86/grana_monte_carlo  ./GranaStack-orig/grana_monte_carlo/test_psi001_7575.cfg
###WORKDIR = /newhome/paulglen/test
###cd $WORKDIR
##mkdir

scl enable devtoolset-1.1 bash ## this should solve environment issues ?  
/bin/env > ./env.out 
OLD_DIR=$PWD
echo $JOB_NAME > ./jobid  
WORKDIR=/state/partition1/paulglen/$JOB_NAME
rm -rf $WORKDIR 2> /dev/null  
mkdir -p $WORKDIR
cp  $OLD_DIR/* $WORKDIR/
cd $WORKDIR 
./wlc  &> $WORKDIR/output.dat
cp -r $WORKDIR/* $OLD_DIR/
date
echo job completed
#echo job completed on $NSLOTS
