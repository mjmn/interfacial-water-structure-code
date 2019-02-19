#!/bin/sh

# Parameters
ifexe=../IFP
qprexe=../QPR

# Determine the water surface
${ifexe} --T 10 --traj 1 --lx 64.578 --ly 64.578 --lz 64.578 > iface_output 2> iface_error

# Print timining
cat iface_output | grep seconds

# Compare the output to the reference output
# First use cmp, only try diff if there is a difference
if cmp -s coord_iface.txt Reference_output/coord_iface.txt; then
  echo "coord_iface.txt matches the reference!"
else
  echo "WARNING: Found differences in coord_iface.txt!"
  diff coord_iface.txt Reference_output/coord_iface.txt
fi

if cmp -s traj_iface.xyz Reference_output/traj_iface.xyz; then
  echo "traj_iface.xyz matches the reference!"
else
  echo "WARNING: Found differences in traj_iface.xyz!"
  diff traj_iface.xyz Reference_output/traj_iface.xyz 
fi
