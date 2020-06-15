#!/bin/bash

## USER INPUTS
stepT=${1:-1000} # Note that this is larger than usual, we skip frames to do time-averaging in Matlab
startT=${2:-1}
endT=${3:-100000}

# Tell the user what we're doing
echo "Looking for interface trajectory and bin_assignment (for the projection on polar angles) of each surface site using a sampling rate of ${stepT}" 

## Check for excisting files
if [ -e clr_trj.xyz ]
then
	rm clr_trj.xyz
fi
touch clr_trj.xyz
# bin assignment
if [ -e bin_assignment.txt ]
then
	rm bin_assignment.txt
fi
touch bin_assignment.txt

# Collect data
for (( i=${startT}; i<=${endT}; i+=${stepT} ))
do
	if [ -e t${i}/traj_iface.xyz ]
	then
	  # Get only the last frame
	  lframeL=$( echo "$( grep -B1 Surface t${i}/traj_iface.xyz | tail -n2 | head -n1 ) + 2" | bc ) # Number of lines for the last frame
	  tail -n${lframeL} t${i}/traj_iface.xyz >> clr_trj.xyz
        else
          echo "WARNING: Could not find file at t${i}/traj_iface.xyz; skipping!"
	fi
	if [ -e t${i}/bin_assignment.txt ]
	then
	  # Change format to be a column with the time step as the first column and the bin assignment in the second column
	  cat t${i}/bin_assignment.txt | awk -v var=${i} '{for (j=1;j<=NF;j++) print var" "$j}' >> bin_assignment.txt
	  #cat t${i}/bin_assignment.txt >> bin_assignment.txt
        else
          echo "WARNING: Could not find file at t${i}/bin_assignment.xyz; skipping!"
	fi
done
