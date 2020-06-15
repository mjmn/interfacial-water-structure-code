#!/bin/bash

# Wrapper script for project_wdump.py

## User inputs
nbins=${1:-24}
nimage=${2:-4}
t0=${3:-1}
tstep=${4:-100}
tsize=${5:-100}
tmax=${6:-100000}

# SOURCE directory
SOURCE="$( dirname "${BASH_SOURCE[0]}" )"
echo "Running scripts found in ${SOURCE}"
## Loop
for ((i=${t0}; i<=${tmax}; i+=${tstep}))
do
  cd t${i}
  ${SOURCE}/project_wdump_linear.py ${tsize} ${nbins} ${nimage} #0 0 90
  cd ../
done
