#!/usr/bin/python

# Convert output generated using the --wdump flag, for the analyze_water_structure.cc script into a 2D projection based on spherical coordinates. This 
# requires the following output files: traj_iface.xyz and qval_per_surface_site.txt.
#
# SYNTAX: ./project_wdump_linear.py [number of frames in this time-window] [number of images] [number of pixels per image dimension (image is NxN)]
#
# This script is typically called using the shell-wrapper project_wdump_linear.sh


## Required modules
import sys
import numpy as np

## Functions

def make_Rmat(sx, sy, sz):
    '''
    Create a rotation matrix that transforms coordinates by the given rotations.
    '''
    Rx = np.array([[1,0,0],[0,np.cos(sx),-np.sin(sx)],[0,np.sin(sx),np.cos(sx)]])
    Ry = np.array([[np.cos(sy),0,np.sin(sy)],[0,1,0],[-np.sin(sy),0,np.cos(sy)]])
    Rz = np.array([[np.cos(sz),-np.sin(sz),0],[np.sin(sz),np.cos(sz),0],[0,0,1]])
    Rmat = np.matmul(Rx,np.matmul(Ry,Rz)) # Complete rotation matrix
    return Rmat
    
def read_user_input():
    '''
    Read command-line input, for now using the sys module
    '''
    if len(sys.argv)<3:
        print('FATAL ERROR: Incorrect number of input parameters, required arguments are in [], optional arguments in {}.')
        print('Correct syntax: ./project_wdump_linear.py [timesteps] [# of images in azimuthal angle] [# of pixels per image side] {Rotation around X} {Rotation around Y} {Rotation around Z}')
        sys.exit(1)
    nframes = int(sys.argv[1])
    nimage = int(sys.argv[2])
    npixels = int(sys.argv[3])
    imgangle = np.pi/2.0 # Hardcoded for now, image size along 1 dimension in rads
    if (imgangle>np.pi*np.cos(imgangle/2.0)): # Check that this image angle is possible, if it's too large it will wrap around the sphere
        print('FATAL ERROR: Image angle is too large, the image would overlap with itself.')
        print('Viable image angles, theta_i, should follow the following criteria: theta_i <= pi*cos(theta_i/2)')
        sys.exit(1)
    # Further input arguments are optional rotation angles to be applied, this can be used to take images from random orientations
    if len(sys.argv)>4:
        sx = np.pi*float(sys.argv[4])/180.0; sy = np.pi*float(sys.argv[5])/180.0; sz = np.pi*float(sys.argv[6])/180.0;
        print('Rotating coordinates by: ',sx,sy,sz)
        print('NOTE: Rotation matrices do not typically commute! Here we perform the rotation around the Z-axis first, then Y, then X: ovec = Rx*Ry*Rz*ivec')
        Rmat = make_Rmat(sx, sy, sz)
    else:
        Rmat = np.eye(3) # Else, just use the identity matrix as a rotation matrix
    return nframes, nimage, npixels, imangle, Rmat

def load_data(traj_iface='traj_iface.xyz',traj_wstruct='qval_per_surface_site.txt'):
    '''
    Load in trajectory data
    '''
    trajLines = open(traj_iface,'r').readlines() # Read interface trajectory
    wdumpLines = open(traj_wstract,'r').readlines() # Read the water structure assigned to each interface grid-site at each time-frame
    return trajLines, wdumpLines

def read_frame(traj, cl, Rmat=np.eye(3)):
    ''' 
    Read a time-frame from a trajectory of the liquid-vapor interface, store in cartesian coordinates.
    '''
    natoms = int(traj[cl]) # Number of grid-sites at this time-frame (XYZ format)
    # Initialize
    ccoords = np.zeros((natoms,3),dtype=float)
    # Read the data
    for atom in range(natoms):
        parts = traj[cl+2+atom].split()
        ccoords[atom,0] = float(parts[1]); ccoords[atom,1] = float(parts[2]); ccoords[atom,2] = float(parts[3]);
    # Center the coordinates, and rotate by the optional rotation matrix
    com = np.mean(ccoords, axis=0)
    ccoords -= com
    ccoords = np.matmul(Rmat,ccoords.T).T
    # Update the number of lines we have read
    cl += natoms+2
    return ccoords, cl

def convert_to_polar(ccoords, Rmat=np.eye(3)):
    '''
    Given a numpy array of size (N,3) Cartesian coordinates and an optional rotation matrix, Rmat, return a set of Spherical coordinates (N,2) with the theta and phi angles
    '''
    # Initialize
    scoords = np.zeros((ccoords.shape[0],2), dtype=float)
    # First rotate
    rcoords = np.matmul(Rmat, ccoords.T).T
    # Now go to angle coordinates
    xy = rcoords[:,0]**2 + rcoords[:,1]**2
    scoords[:,0] = np.arctan2(np.sqrt(xy), rcoords[:,2]) # Theta, between 0 and pi
    scoords[:,1] = np.arctan2(rcoords[:,1], rcoords[:,0]) # Phi, between -pi and +pi
    return scoords

def create_image_centers(nimage):
    '''
    NOTE: The current method does not yield evenly distributed image centers! TODO: Replace with something that gets closer to an even distribution (i.e., spiral method).
    Create a list of the polar coordinates that define the center of an image, depends only on the number of images we want on the sphere.
    '''
    theta_step = np.pi/(nimage-1); phi_step = np.pi/nimage; # The difference with theta is due to non-periodicity in that coordinate (i.e., for nimage=2 we want centers at theta=0 and pi)
    icenter = []
    for i in range(nimage):
        for j in range(nimage*2):
            icenter.append(np.array([i*theta_step,j*phi_step - np.pi])) # Theta ranges in [0,pi], Phi ranges in [-pi,pi]
    return icenter

def find_pixel(coord, icenter, iangle, npixels):
    '''
    Determine the pixel to which this grid-point should be assigned.
    '''
    # NO LONGER NEEDED SINCE I START BY SETTING CENTERS TO [pi/2, 0]
    #dtheta, dphi = periodic_angles(coord, icenter) # Calculate the signed distance between this grid-point and the image center, considering periodicity
    # Calculate the bin in the theta coordinate
    theta_pixel = np.floor((coord[0]-(np.pi/2.0 - iangle/2.0))*npixels/iangle) # Get the pixel value
    # Calculate the bin in the phi-coordinate
    stheta = np.sin(coord[0]) # In order to maintain the same image-width, we need to bin over wider ranges in phi depending on the theta angle
    phi_range = iangle/stheta
    phi_pixel = np.floor((coord[1]+iangle/2.0)*npixels/phi_range) # Get the pixel value
    # Check if we're within the image
    if (theta_pixel<0) or (theta_pixel>=npixels) or (phi_pixel<0) or (theta_pixel>=npixels):
        return -1
    else:
        return theta_pixel*npixels + phi_pixel

def pixel_assignment(scoords, iangle, npixels):
    '''
    Assign grid-points to pixels.
    '''
    bcoords = np.zeros((scoords.shape[0],1), dtype=int)
    for gp in range(scoords.shape[0]):
        gp_coord = scoords[gp,:]
        bcoords[gp] = find_pixel(gp_coord, iangle, npixels)
    return bcoords

def set_polar_to_center(ic):
    '''
    Take a numpy array with [Theta, Phi] and output a rotation matrix such that this point will instead be at [pi/2,0].
    '''
    return make_Rmat(0, np.pi/2.0 - ic[0], -ic[1]) # Create a rotation matrix that sets a set of polar coordinates align with the X-axis 

def update_histograms(wtraj, cl, bc, q, n):
    '''
    Update histograms
    '''
    tfo = int(wtraj[cl].split()[0]) # Time-frame information in the wdump file is in the first column
    while(1):
        try:
            parts = wtraj[cl].split()
        except:
            print('Finished with the wdump file at time',t)
            break
        # Check if the time-frame changed
        tf = int(parts[0])
        if (tf!=tfo):
            tfo = tf
            break
        ### TODO: Insert histogram fill operations
        cl += 1
        return qhist, nhist, cl


        
## Main
if __name__ == "__main__":
    ### INPUT & INITIALIZE
    # Read user provided inputs
    nframes, nimage, npixels, imgangle, Rmat = read_user_input()
    # Read trajectory data
    trajLines, wdumpLines = load_data(traj_iface='traj_iface.xyz',traj_wstruct='qval_per_surface_site.txt')
    # Initialize containers for output + variables derrived from input
    qhist = np.zeros((npixels*npixels,2*nimage*nimage),dtype=float) # Contain q-values
    nhist = np.zeros((npixels*npixels,2*nimage*nimage),dtype=int) # Contain counts in each bin
    bwidth = imgangle/npixels # Size of a pixels side (used for binning)
    imgcenters = create_image_centers(nimage) # Create a list of image centers (in polar-coordinates)
    tfo = int(wdumpLines[0].split()[0]) # Number of waters in the first frame
    flag = 1 # Track whether we reached the end of a frame
    ctrajLine = 0; cwdumpLine = 0; # Lines already processed

    ### MAIN CALCULATION
    for t in range(nframes): # Loop over trajectory data
        ccoords, ctrajLine = read_frame(traj=trajLines, cl=ctrajLine, Rmat=Rmat) # Extract grid-point location in Cartesian coordinates (ccoords)
        bcoords = np.zeros((ccoords.shape[0],len(imgcenters)), dtype=int) # Initialize array for bin assignment
        for i in range(len(imgcenters)):
            ic = imgcenters[i]
            Rmat = set_polar_to_center(ic) # Create a rotation matrix that will set the image center at [theta=pi/2, phi=0], this makes further steps much easier
            scoords = convert_to_polar(ccoords=ccoords, Rmat=Rmat) # Convert from Cartesian to Polar coordinates, set image center to [pi/2, 0]
        #scoords, ctrajLine = read_frame_polar(traj=trajLines, cl=ctrajLine, Rmat=Rmat) # Extract grid-point location in angular coordinates (scoords)
            bcoords[:,i] = pixel_assigment(scoords=scoords, iangle=imgangle, npixels=npixels) # Assign grid-points to pixels, for this image
        qhist, nhist, wdumpLine = update_histograms(wtraj=wdumpLines, cl=cwdumpLine, bc=bcoords, q=qhist, n=nhist) # Update the output histograms with water data
        
    
    ## PART TWO: Assign the contribution of each water molecule to \delta\lambda_{ref} to the corresponding bin in cos(theta),cos(phi) space
    while (1): # I use break statements to exit this loop
        # Try read the line, checking that we didn't reach the end of the file
        try: 
            parts = wdumpLines[cwdumpLine].split()
        except:
            print('Finished with the wdump file at time',t)
            break
        # Check if we are still in the same time-frame
        tf = int(parts[0]) # Time-frame at which this \delta\lambda_{ref} contribution is calculated
        if ( tf != tfo ): # We're at the next frame
            tfo = tf
            break # Quit the while loop
        # Extract other information from this line
        gp = int(parts[1]) # Grid-point to which this contribution was assigned
        qval = float(parts[2]) # The \delta\lambda_{ref} contribution
        # Add to the corresponding output matrices
        qhist[bcoords[gp],icoords[gp]] += qval
        nhist[bcoords[gp],icoords[gp]] += 1
        # Move on to the next line
        cwdumpLine += 1

# To avoid weird behavior due to division by zero, replace all remaining zeros in nhist with 1, since qhist is always 0 where nhist is 0 this yields a 0 in qhist/nhist
nhist[nhist==0] = 1

# Divide the qhist by the counts, then write output
qhist = np.divide(qhist,nhist)
np.savetxt('image_projection.txt',qhist.T,fmt='%12.6f')
np.savetxt('image_projection_N.txt',nhist.T,fmt='%12.6f')
# Save the final bin assignment, useful for projecting onto a trajectory
np.savetxt('bin_assignment_linear.txt',bcoords.T,fmt='%i')
np.savetxt('image_assignment_linear.txt',icoords.T,fmt='%i')
