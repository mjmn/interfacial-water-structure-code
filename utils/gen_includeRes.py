#!/usr/bin/python

# Python script that generates res_list.txt given a pdb file

import sys

pdbfn = sys.argv[1] # Pdb filename
selection = []

for i in sys.argv[2:]: # Rest of the arguments are selection criteria
  selection.append((int(i.split('-')[0]),i.split('-')[1])) # Append a tupple with resid and chain

# Read input
lines = open(pdbfn,'r').readlines()

# Initialize
ol = []

# Get information from file
for line in lines[1:]:
 if line[0:4]=='ATOM':
  if line[17:20]=='SOL':
    break # Done with all protein lines
  resID = int(line[22:26]) # Residue ID on this line
  chainID = line[21]
  if (resID,chainID) in selection:
    ol.append('1')
  else:
    ol.append('0')

# Write output
fo = open('resInclude.txt','w')
fo.write('\n'.join(ol))
fo.close()

