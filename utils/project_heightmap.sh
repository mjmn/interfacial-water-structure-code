#!/bin/bash

# Wrapper script for project_heightmap.py

## User inputs
nbins=${1:-100}
t0=${2:-1}
tstep=${3:-100}
tsize=${4:-100}
tmax=${5:-100000}

# SOURCE directory
SOURCE="$( dirname "${BASH_SOURCE[0]}" )"
echo "Running scripts found in ${SOURCE}"
## Loop
for ((i=${t0}; i<=${tmax}; i+=${tstep}))
do
  cd t${i}
  ${SOURCE}/project_heightmap.py ${tsize} ${nbins} #0 0 90
  cd ../
done
