#!/usr/bin/python

# Normalize a given probability density such that \delta\lambda_ref = \lambda_ref - < \lambda_ref >_{0}
# This requires two probability density files of the same size, P_{ref}, the distribution by which we sample at the reference surface, and P_i, the modified probability density (i.e., divided by bulk probability density

import sys
import numpy as np

# Distribution by which we sample
pjointF = sys.argv[1]
pjointC = int(sys.argv[2])

# Distribution used to calculate the log-likelihood, this is the distribution we are going to normalize
piF = sys.argv[3]
piC = int(sys.argv[4])

# Read the data
psample = np.loadtxt(pjointF)
pref = np.loadtxt(piF)

# Calculate < \lambda_ref >_0
if (len(sys.argv)>5):
    minR = int(sys.argv[5])
    maxR = int(sys.argv[6])
    print("Considering only bins ",minR," to ",maxR," for normalization.")
    lambda_ev = -np.sum(psample[minR:maxR,pjointC]*np.log(pref[minR:maxR,piC]))/np.sum(psample[minR:maxR,pjointC])
else:
    lambda_ev = -np.sum(psample[:,pjointC]*np.log(pref[:,piC]))/np.sum(psample[:,pjointC])

print(" Calculated the expected value of lambda as: ",lambda_ev)
exp_lambda_ev = np.exp(lambda_ev)
print("Normalization factor :",exp_lambda_ev)

# Apply the normalization factor and save to file
pref[:,piC] = pref[:,piC]*exp_lambda_ev
# Save to file
np.savetxt('pjoint_normalized.txt',pref,fmt='%12.6f')
print("Normalized probability density written to file: pjoint_normalized.txt")
