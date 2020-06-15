#!/bin/bash

# Convert a histogram file created by output_water_structure.cc into a pjoint_sm2.txt input file.
# I use Bayesian Inference to convert histogram counts to probabilities; this results in the relative probability
# of each bin being given by (N[bin]+1)/(sum(N[:])+2)

# First determine sum(N[:])
Nsum=$( cat ${1} | awk '{ sum += $1 } END { print sum }' )
# Convert it to the normalization factor, the 225000 is to match the integrated relative probabilities of pbulk
#normfac=$( echo "225000/(${Nsum}+2)" | bc -l )
#normfac=$( echo "(900000-2)/(${Nsum})" | bc -l )
#normfac=$( echo "225000/(${Nsum})" | bc -l )
normfac=$( echo "400/(${Nsum})" | bc -l )
echo "${normfac}"

# Next add the 1 to each N[bin] and write pjoint_sm2.txt
#cat ${1} | awk -v var="${normfac}" '{ print "0.0 0.0 0.0 "($1+1)*var" 0.0" }' > pjoint_sm2.txt
#cat ${1} | awk -v var="${normfac}" '{ print "0.0 0.0 0.0 "$1*var+1.0" "$1 }' > pjoint_sm2.txt
cat ${1} | awk -v var="${normfac}" '{ print "0.0 0.0 0.0 "$1*var" "$1 }' > pjoint_sm2.txt
