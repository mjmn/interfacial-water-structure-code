#!/bin/sh

# Aggregate the data from these trajectories
## PARAMETERS
startT=${1:-1}
endT=${2:-100000}
stepT=${3:-100}

# The qlv score for each residue as a function of time
# Check if the file excist, if so overwrite
if [ -e qlv_traj.txt ]
then
 echo "Overwriting excisting qlv_traj.txt file!"
 rm qlv_traj.txt
fi
if [ -e qlv_traj_N.txt ]
then
 echo "Overwriting excisting qlv_traj_N.txt file!"
 rm qlv_traj_N.txt
fi
if [ -e wnear.txt ]
then
 echo "Overwriting excisting wnear.txt file!"
 rm wnear.txt
fi
if [ -e projected_height.txt ]
then
 echo "Overwriting excisting projected_height.txt and projected_height_N.txt file!"
 rm projected.txt
 rm projected_N.txt
fi
if [ -e projected.txt ]
then
 echo "Overwriting excisting projected.txt and projected_N.txt file!"
 rm projected.txt
 rm projected_N.txt
fi
# Create blank file
touch qlv_traj.txt
touch wnear.txt
touch projected_height.txt
touch projected_height_N.txt
touch projected.txt
touch projected_N.txt
touch qlv_traj_N.txt
# Now loop over the trajectory and add the frames to the output file
for (( t=${startT} ; t<=${endT} ; t+=${stepT} ))
do
 cat t${t}/qlv_avg_res_6.txt | awk '{print $2}' | paste -s -d" " >> qlv_traj.txt
 cat t${t}/qlv_avg_res_6.txt | awk '{print $4}' | paste -s -d" " >> qlv_traj_N.txt
 if [ -e t${t}/local_wstructure.txt ]
 then
  cat t${t}/local_wstructure.txt >> wnear.txt
 fi
 if [ -e t${t}/spherical_heightMap.txt ]
 then
  cat t${t}/spherical_heightMap.txt >> projected_height.txt
  cat t${t}/spherical_heightMap_N.txt >> projected_height_N.txt
 fi
 if [ -e t${t}/spherical_projection.txt ]
 then
  cat t${t}/spherical_projection.txt >> projected.txt
  cat t${t}/spherical_projection_N.txt >> projected_N.txt
 fi
done
echo "Wrote qlv trajectory to file qlv_traj.txt, and number of samples to qlv_traj_N.txt"
