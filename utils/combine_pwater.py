#!/home/mniesen/Software/Anaconda/anaconda3/bin/python

# Script that can combine water structure probability distribution files in Sucheols initial format into a single file in the new format.
# Specifically, take pofa.txt, pbulk_w2_tip3p.txt, and pjoint_sm2.txt --> pjoint.txt
# If the user wants to make a new reference by, i.e., not dividing by pbulk, the script can be modified to do that as well.

import numpy as np
import sys
import math

### PARAMETERS ###
astep = 0.1 # Step size in the a coordinate, typically 0.1
cstep = 0.02 # Step size for the cos(theta1) and cos(theta2) coordinates, typically 0.02

### Read the inputs as numpy arrays ###
# The p(a) 1-D distribution function
try:
 pofa=np.loadtxt('pofa.txt')
except:
 pofa=np.ones((140,2),dtype=float) # If it doesn't excist then it may already be in pjoint, so use all ones

# The p(cos(theta1),cos(theta2)) distribution function for bulk water (it is independent of a)
try:
 pbulk=np.loadtxt('pbulk_w2_tip3p.txt')
except:
 pbulk=np.ones((10000,3),dtype=float) # If it doesn't excist then don't normalize by bulk, so use all ones

# Joint distribution of the reference surface p(a,cos(theta1),cos(theta2)) is a required input
try:
 pjoint=np.loadtxt('pjoint_sm2.txt')
except:
 print("ERROR: Could not load pjoint_sm2.txt")
 sys.exit(1)


### Initialization ###
# The original pofa file has a larger range for collective variable a than the pjoint file. Make sure that the ranges are consistent.
min_pja = min(pjoint[:,0]); max_pja = max(pjoint[:,0]); # Min and Max values that are in the joint distribution
min_index = int((min_pja+4.0)*10); max_index = int((max_pja+4.0)*10)+1; # Convert to indeces for pofa
pofa = pofa[min_index:max_index,:]; # Slice out the part that is in the joint distribution as well
# Output array, its a manipulated version of pjoint
pjoint_o = np.copy(pjoint)
# Counters for array positions
index_a = 0; index_b = 0; index_c = 0;
# Number of values in each dimension
numa = max_index-min_index;
numc = int((2.0/cstep));
# The expectation value of \lambda when sampling over the reference distribution
lambda_ev = 0.0; normfac_lambda_ev = 0.0;


### Do the required math to get a single joint probability distribution out ###
for a in range(numa):
 index_b = 0; # Reset bulk array index every time we reset angles
 for c1 in range(numc):
  for c2 in range(numc): # Note that they run up to 0.98
   # Other than the probability density, the values in the output array are the same as the input, so don't update
   if (pbulk[index_b,2]>0.0) and (pjoint[index_c,3]>0.0): # Handle zeros in the probability distribution, avoiding division by zero
    pjoint_o[index_c,3] = pofa[index_a,1]*pjoint[index_c,3]/pbulk[index_b,2]; # Net probability density
    lambda_ev += -pjoint[index_c,3]*math.log(pjoint_o[index_c,3]); # Contribution to the expectation value of lambda, when sampling the bulk
    normfac_lambda_ev += pjoint[index_c,3];
   elif (pjoint[index_c,3]>0.0):
    print("WARNING: pbulk = 0.0 for a collective coordinate where pjoint > 0.0, approximating ")
    print("Pjoint: ",pjoint[index_c,:]," Pbulk: ",pbulk[index_b,:])
   # sys.exit(1)
    pjoint_o[index_c,3] = pofa[index_a,1] # Assume no angle bias, until I get to fix pbulk (*** Division by zero should not occur)
    lambda_ev += -pjoint[index_c,3]*math.log(pofa[index_a,1]); normfac_lambda_ev += pjoint[index_c,3];
   else:
    pjoint_o[index_c,3] = 1.0
   index_c += 1; # Update for every coordinate change
   index_b += 1; # Update for every angle coordinate change (see reset)
 index_a += 1; # Only update when the a coordinate changes

# Include subtraction of the expectation value of lambda when sampling the reference distribution, such that the generated pjoint.txt file calculates
# \lambda_ref - <\lambda_ref>_0. This is done by multiplying every element in the pjoint distribution with \exp{<\lambda_ref>_0}
lambda_ev /= normfac_lambda_ev
print("Calculated <\lambda_ref>_0 as: ",lambda_ev," correcting the pjoint appropriately.")
pjoint_o *= math.exp(lambda_ev)

### Write output to file ###
np.savetxt('pjoint.txt',pjoint_o,fmt='%12.6f')
