#!/usr/bin/python

import os
import sys

# EXECUTABLES
ows="/home/mniesen/Software/qref/OWS"
# Important environmental variables for running batch jobs on this cluster
exportlines = "export LD_LIBRARY_PATH=/home/mniesen/Software/gcc730/libexec/gcc/x86_64-pc-linux-gnu/7.3.0:$LD_LIBRARY_PATH\nexport LD_LIBRARY_PATH=/home/mniesen/Software/gcc730/lib64:$LD_LIBRARY_PATH"

# Python script to generate a input script for the analysis of water structure for a time-slice from an MD trajectory
startT=int(sys.argv[1])
duration=int(sys.argv[2])
prefix=sys.argv[3]
optL = []
for i in sys.argv[4:]:
  optL.append(i)
opts = ' '.join(optL) # Additional options
#opts = "--hist --rcut 6.0 --cylinder" # Hard-coded for now --cylinder --cylinder 

# Try to extract other required inputs from the pdb file
tfilen = prefix+"reference.pdb"
tlines = open(tfilen,'r').readlines()
# Intialize
nprot = 0

# Loop to extract info we need
for line in tlines:
  if line[0:5]=='CRYST':
    # Box-size info should be found at this defined location (GROMACS generated PDB)
    boxX = line[6:15].strip()
    boxY = line[15:24].strip()
    boxZ = line[24:33].strip()
  if line[0:4] == 'ATOM': # Atom info, could be protein
    resType=line[17:20]
    if ((resType!='SOL') and (resType!='TIP')): # This is not a water
      nprot += 1
    if ((resType=='SOL') or (resType=='TIP')): # Skip the remaining lines
      break

# We should have counted all the relevant information, generate the command-line argument that can be used to analyze this timeslice
cmd1 = "%s %s --prefix %s --T %i --lx %s --ly %s --lz %s --startT %i --nprot %i\n" % (ows, opts, prefix, duration, boxX, boxY, boxZ, startT, nprot)
# Will write a PBS batch file, write header info etc
header = "#!/bin/sh\n#PBS -l nodes=1:ppn=1\n#PBS -l walltime=999:00:00\n#PBS -N wn_t%s\n#PBS -q default\n\n%s\n\ncd $PBS_O_WORKDIR\n\n" % (startT,exportlines)

# Now setup a work-directory and go there
wdirn = 't%i' % startT;# os.mkdir(wdirn);
os.chdir(wdirn);

# Write the batch file
of = open('job_ref.batch','w')
of.write(header)
of.write(cmd1)
of.close()

# Submit the job to the cluster
# Done by the .sh script wrapper
