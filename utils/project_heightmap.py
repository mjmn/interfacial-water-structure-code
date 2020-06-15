#!/usr/bin/python

# Convert output generated using the --wdump flag, for the analyze_water_structure.cc script into a 2D projection based on spherical coordinates. This 
# requires the following output files: traj_iface.xyz and qval_per_surface_site.txt.
#
# SYNTAX: ./project_wdump.py [number of frames in this time-window] [number of bins in the cos(theta) coordinate]
#
# The cos(phi)*sign[sin(phi)] coordinate will have twice as many bins as the cos(theta) coordinate. This script is typically called using the shell-wrapper
# project_wdump.sh


## Required modules
import sys
import numpy as np

## Read user input
if len(sys.argv)<3:
    print('ERROR: Incorrect number of input arguments.')
    print('Correct syntax: ./project_wdump.py [timesteps] [bins] [OPTIONAL: Rotation around each axis in degrees, i.e., 0 0 180 - flip z-axis]')
    sys.exit(1)
nframes = int(sys.argv[1])
nbins = int(sys.argv[2])
if len(sys.argv)>3:
    sx = np.pi*float(sys.argv[3])/180.0; sy = np.pi*float(sys.argv[4])/180.0; sz = np.pi*float(sys.argv[5])/180.0;
    print('Rotating coordinates by: ',sx,sy,sz)
    print('NOTE: Rotation matrices do not typically commute! Here we perform the rotation around the Z-axis first, then Y, then X: ovec = Rx*Ry*Rz*ivec')
    Rx = np.array([[1,0,0],[0,np.cos(sx),-np.sin(sx)],[0,np.sin(sx),np.cos(sx)]])
    Ry = np.array([[np.cos(sy),0,np.sin(sy)],[0,1,0],[-np.sin(sy),0,np.cos(sy)]])
    Rz = np.array([[np.cos(sz),-np.sin(sz),0],[np.sin(sz),np.cos(sz),0],[0,0,1]])
    Rmat = np.matmul(Rx,np.matmul(Ry,Rz)) # Complete rotation matrix
else:
    Rmat = np.eye(3) # Else, just use the identity matrix as a rotation matrix

## Initialize
# Files to read
trajLines = open('traj_iface.xyz','r').readlines() # Read interface trajectory
ctrajLine = 0 # Keep track of the lines we have already read
# Output data
qhist = np.zeros((nbins*nbins*2,1),dtype=float) # Contain q-values
nhist = np.zeros((nbins*nbins*2,1),dtype=int) # Contain counts in each bin
# For binning spherical coordinates
bwidth = 2.0/nbins # Width of a bin
# Needed to detect when the wdump file changes to the next frame

## Loop over time
for t in range(nframes):

    ## PART ONE: Convert the grid-coordinates into spherical coordinates, and assign each interfacial surface site to a bin in cos(theta),cos(phi) space
    natoms = int(trajLines[ctrajLine]) # Number of grid-sites at this time-frame
    # Initialize a numpy array to contain the cartesian coordinates at this time-frame
    ccoords = np.zeros((natoms,3),dtype=float)
    # Initialize a numpy array to contain the spherical coordinates, theta and phi, at this time-frame
    scoords = np.zeros((natoms,3),dtype=float)
    # Initialize a numpy array to contain the bin index to which each gridpoint contributes at this time-frame
    bcoords = np.zeros((natoms,1),dtype=int)
    for atom in range(natoms): # Fill the numpy array
        parts = trajLines[ctrajLine+2+atom].split()
        ccoords[atom,0] = float(parts[1]); ccoords[atom,1] = float(parts[2]); ccoords[atom,2] = float(parts[3]);
    # Now calculate COM and convert to spherical coordinates
    com = np.mean(ccoords, axis=0)
    ccoords -= com # Set the origin to [0,0,0]
    ccoords = np.matmul(Rmat,ccoords.T).T # Rotate using the user-defined rotation matrix
    # Conversion to spherical coordinates, we only care about the angles and not the value of 'r', so we don't calculate it here
    xy = ccoords[:,0]**2 + ccoords[:,1]**2
    scoords[:,0] = np.arctan2(np.sqrt(xy), ccoords[:,2]) # Theta
    scoords[:,1] = np.arctan2(ccoords[:,1], ccoords[:,0]) # Phi
    scoords[:,2] = np.sqrt(xy + ccoords[:,2]**2) # |r|
    # Further conversion, directly to bin coordinates to which this grid-point would contribute
    #bcoords[:,0] = 2*nbins*np.floor((np.cos(scoords[:,0])+1)/bwidth)  +  np.floor((2+np.sign(np.sin(scoords[:,1]))*(np.cos(scoords[:,1])+1))/bwidth)
    bcoords[:,0] = 2*nbins*np.floor((np.cos(scoords[:,0])+1)/bwidth)  + np.floor((2*(scoords[:,1]+np.pi)/np.pi)/bwidth) # Alternate, angle for phi and cos(theta), yields same volume grid-points 
    # Update the lines we have read
    ctrajLine += 2+natoms

    # Fill out height-map
    for atom in range(natoms):
        qhist[bcoords[atom]] += scoords[atom,2]
        nhist[bcoords[atom]] += 1

# To avoid weird behavior due to division by zero, replace all remaining zeros in nhist with 1, since qhist is always 0 where nhist is 0 this yields a 0 in qhist/nhist
nhist[nhist==0] = 1

# Divide the qhist by the counts, then write output
qhist = np.divide(qhist,nhist)
np.savetxt('spherical_heightMap.txt',qhist.T,fmt='%12.6f')
np.savetxt('spherical_heightMap_N.txt',nhist.T,fmt='%12.6f')
# Save the final bin assignment, useful for projecting onto a trajectory
np.savetxt('bin_assignment.txt',bcoords.T,fmt='%i')
