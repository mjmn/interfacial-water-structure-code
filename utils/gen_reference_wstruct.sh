#!/bin/sh

# Perform the necessary steps to setup and submit a bunch of analysis jobs
# I make use of the fact that the analysis is proudly parallel; the trajectory
# is time-sliced and each slice is submitted as a seperate job. If doing this on 
# a cluster with slurm, array jobs could be used, but here I use PBS
# opts are the last argument and need to be provided between single quotes
# i.e., '--cylinder --rcut 6.0'

# INPUT; can be from a loop
startTI=${1:-1}
Nsteps=${2:-100}
Nframes=${3:-100}
prefix=${4:-"./traj/"}
opts="$5"

echo ${opts}

for ((startT=${startTI}; startT<100000; startT+=${Nsteps}))
do
 echo ${startT}
 if [ ! -e t${startT} ]
 then
  echo "ERROR: no run directory found for ${startT}"
  stop
 fi
 # First use a python script to setup a directory with the batch script in it
 /home/mniesen/Scripts/gen_reference_wstruct.py ${startT} ${Nframes} ${prefix} ${opts}
 # Make sure that the python script made the directory, then go there
 cd t${startT}
 # Submit the batch job
 qsub job_ref.batch
# sleep 2 # Just to not flood anything
 # Back to main directory
 cd ../
done
