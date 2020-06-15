#!/home/mniesen/Software/Anaconda/anaconda3/bin/python

import numpy as np
import sys

# Input
qlvfn = sys.argv[1]

# Read input info
qlv=np.loadtxt(qlvfn)
a=np.shape(qlv)

# Reshape
clr = np.reshape(qlv, (a[0]*a[1],))

# normalize
#maxv = max(clr); minv = min(clr);
#clr = (clr-minv)/(maxv-minv)

# Write
np.savetxt('colors.txt',clr,fmt='%12.6f')
