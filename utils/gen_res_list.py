#!/usr/bin/python

# Python script that generates res_list.txt given a pdb file

import sys

pdbfn = sys.argv[1] # Pdb filename
maxres = int(sys.argv[2]) # Number of residues in protein

# Read input
lines = open(pdbfn,'r').readlines()

# Initialize
resnum = 0; atomnum = 1;
prevRes = ''
ol = []

# Get information from file
for line in lines[1:]:
 if line[0:4]=='ATOM':
  atomID = line[4:11] # Atom ID on this line
  resID = line[22:26] # Residue ID on this line
  if resID!=prevRes: # Update residue number if the residue ID changed
   resnum += 1
   prevRes = resID
  if "--splitChains" in sys.argv: # We want to assign different a different index to each residue, correcting for the numbering resetting in the PDB
      if resnum > maxres:
        break
      ol.append(str(atomnum)+" "+str(resnum))
  else: # Use residue numbers exactly as in the PDB file
      if int(resID) > maxres:
        break
      ol.append(str(atomnum)+" "+str(int(resID)))
  atomnum += 1

# Write output
fo = open('res_list.txt','w')
fo.write('\n'.join(ol))
fo.close()

