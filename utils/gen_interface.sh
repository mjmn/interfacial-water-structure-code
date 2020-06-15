#!/bin/sh

# Perform the necessary steps to setup and submit a bunch of analysis jobs
# I make use of the fact that the analysis is proudly parallel; the trajectory
# is time-sliced and each slice is submitted as a seperate job. If doing this on 
# a cluster with slurm, array jobs could be used, but here I use PBS

# INPUT; can be from a loop
startT=${1:-1}
Nsteps=${2:-100}
Nframes=${3:-100}
opts=${4}
maxFrames=${5:-100000}

# SOURCE directory
SOURCE="$( dirname "${BASH_SOURCE[0]}" )"
echo "Running scripts found in ${SOURCE}"

for ((startT=1; startT<=${maxFrames}; startT+=${Nsteps}))
do
 echo ${startT}
 # Clear any excisting directory
 rm -r t${startT}
 mkdir t${startT}
 # The directory should excist at this point, if not, we can't do analysis 
 if [ ! -e t${startT} ]
 then
  echo "ERROR: no run directory found for ${startT}"
  stop
 fi
 # First use a python script to setup a directory with the batch script in it
 ${SOURCE}/gen_interface.py ${startT} ${Nframes} ${opts}
 cd t${startT}
 # Submit the batch job
 qsub job.batch
# sleep 2 # Just to not flood anything
 # Back to main directory
 cd ../
done
