#!/usr/bin/python

import os
import sys

# Find the path to this script
path=os.path.dirname(sys.argv[0])
ifp=path+"/../IFP"

# EXECUTABLES
#ifp="/home/mniesen/Software/qref/IFP"
#qpr="/home/mniesen/Software/qref/QPR"
# Important environmental variables for running batch jobs on this cluster
exportlines = "export LD_LIBRARY_PATH=/home/mniesen/Software/gcc730/libexec/gcc/x86_64-pc-linux-gnu/7.3.0:$LD_LIBRARY_PATH\nexport LD_LIBRARY_PATH=/home/mniesen/Software/gcc730/lib64:$LD_LIBRARY_PATH"

# Python script to generate a input script for the analysis of water structure for a time-slice from an MD trajectory
startT=int(sys.argv[1])
duration=int(sys.argv[2])
optL = []
for i in sys.argv[3:]:
    optL.append(i)
opts = ' '.join(optL) # User-defined optional input parameters    
#opts="--cylinder"

# Try to extract other required inputs from the pdb file
tfilen = "./traj/reference.pdb"
tlines = open(tfilen,'r').readlines()
# Intialize
nwater = 0; nprot = 0; nres = 0;
prevRes = 0;

# Loop to extract info we need
for line in tlines:
  if line[0:5]=='CRYST':
    # Box-size info should be found at this defined location (GROMACS generated PDB)
    boxX = line[6:15].strip()
    boxY = line[15:24].strip()
    boxZ = line[24:33].strip()
  if line[0:4]=='ATOM':
    resType=line[17:20] # Residue type
    resID = line[22:26] # Residue ID
    if ((resType=='SOL') or (resType=='TIP')): # This is a water
      nwater += 1;
      if nprot == 0: # This is the first water, we can now set protein info
        nprot = int(line[4:11])-1; #nres = int(line[22:26])-1;
    if (nprot==0) and (resID!=prevRes):
        nres += 1; prevRes = resID; # Count unique protein residues

# We should have counted all the relevant information, generate the command-line argument that can be used to analyze this timeslice
cmd1 = "%s --T %i --traj 1 --lx %s --ly %s --lz %s --startT %i --nwater %i --nprot %i --nres %i %s\n" % (ifp, duration, boxX, boxY, boxZ, startT, nwater/3, nprot, nres, opts)
# Will write a PBS batch file, write header info etc
header = "#!/bin/sh\n#PBS -l nodes=1:ppn=1\n#PBS -l walltime=999:00:00\n#PBS -N wi_t%i\n#PBS -q default\n\n%s\necho $PBS_O_HOST\necho $PBS_O_NODEFILE\n\ncd $PBS_O_WORKDIR\n\n" % (startT, exportlines)

# Now setup a work-directory and go there
wdirn = 't%i' % startT; #os.mkdir(wdirn);
os.chdir(wdirn);

# Write the batch file
of = open('job.batch','w')
of.write(header)
of.write(cmd1)
of.close()

# Submit the job to the cluster
# Done by the .sh script wrapper
