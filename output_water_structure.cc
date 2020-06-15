/* New script that reads a water structure trajectory, along with the corresponding protein trajectory, and outputs only the structure of water around a user provided list of residues. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <iomanip>

using namespace std;

/**********************/
/*  Define parameters */
/**********************/

/**********************/
/* Define variables   */
/**********************/

// Time-related
int Tstep = 0; // Duration of the trajectory slice to be analyzed (in timesteps)
int firstTraj = 0; // First trajectory frame to be read

// Simulated system related
int N_prot = 0; // Number of atoms in the protein
double L_box[3] = {0}; // Size of the periodic box (make sure it doesn't change over time, so run in a constant Volume ensemble)

// Parameters for analysis
double R_cut = 6.0; // Water needs to be within this cut-off of an atom to be included in output
char * prefix; // Name of the folder that contains the trajectory files
// Related to histogram spacing for output - now a user input, to deal with cases where less sample data is available
int Nabin = 90; // Number of bins in the coordinate that defines the distance from a water molecule to the W-C interface
int Nc1bin = 100; // Number of bins in the coordinate that defines the angle between the surface normal vector and the water OH1 bond
int Nc2bin = 100; // Number of bins in the coordinate that defines the angle between the surface normal vector and the water OH2 bond
// Histogram data range, currently not an input
double aRange = 9.0; // Distance from the W-C interface ranges between -1.0 and 8.0
double c1Range = 2.0; // Cosine values range between -1 and +1
double c2Range = 2.0;

// For output
bool verbose=false; // Write a lot more info, for debugging
bool raw=false; // Write an output file containing all selected waters
bool hist=false; // Write histogram data on the local water structure
bool cylinder=false; // Keep only the waters that are nearer to a tagged protein atom than any other protein atom (cylinder selection mode)

/**********************/
/* Declare functions  */
/**********************/
double calc_rsq(double * r1, double * r2); // Calculate the squared distance between two points

/**********************/
/* User input         */
/**********************/

void cl_parse(int argc, char * argv[]) {
  bool match;
  int i = 1;
  while ( i < argc ) {
  //for (int i=1; i<argc; i+=2) {
    match = false; // Track if we find a match for the provided input argument
    // Trajectory file related
    if(strcmp(argv[i],"--T")==0) {Tstep = atof(argv[i+1]);match=true;i+=2; continue;} // Number of frames to analyze 
    if(strcmp(argv[i],"--startT")==0) {firstTraj = atof(argv[i+1]);match=true;i+=2; continue;} // First frame to consider in analysis, useful for time-slicing
    // Analysis file related
    if(strcmp(argv[i],"--rcut")==0) {R_cut = atof(argv[i+1]);match=true;i+=2; continue;} // Cut-off distance for including waters in the output
    if(strcmp(argv[i],"--prefix")==0) {prefix = argv[i+1];match=true;i+=2; continue;} // File with atom index of protein atoms for which we want to output nearby waters
    // System related
    if(strcmp(argv[i],"--nprot")==0) {N_prot = atof(argv[i+1]);match=true;i+=2; continue;} // Number of atoms in the protein - avg_slice.py
    if(strcmp(argv[i],"--lx")==0) {L_box[0] = atof(argv[i+1]);match=true;i+=2; continue;} // Simulation box dimensions
    if(strcmp(argv[i],"--ly")==0) {L_box[1] = atof(argv[i+1]);match=true;i+=2; continue;}
    if(strcmp(argv[i],"--lz")==0) {L_box[2] = atof(argv[i+1]);match=true;i+=2; continue;}
    // Related to histogramming the water data
    if(strcmp(argv[i],"--nabin")==0) {Nabin = atof(argv[i+1]);match=true;i+=2; continue;} // Number of bins for the water-structure histogram in each CV dimension
    if(strcmp(argv[i],"--nc1bin")==0) {Nc1bin = atof(argv[i+1]);match=true;i+=2; continue;}
    if(strcmp(argv[i],"--nc2bin")==0) {Nc2bin = atof(argv[i+1]);match=true;i+=2; continue;}
    // Related to user desired output, these are flags that dont require further arguments
    if(strcmp(argv[i],"--raw")==0) {raw = true;match=true;i+=1; continue;}
    if(strcmp(argv[i],"--hist")==0) {hist = true;match=true;i+=1; continue;}
    if(strcmp(argv[i],"--verbose")==0) {verbose = true;match=true;i+=1; continue;}
    if(strcmp(argv[i],"--cylinder")==0) {cylinder = true;match=true;i+=1; continue;}
    // Check if the user provided a bad input argument
    if(match == false) {cout<<"Warning: Bad input parameter: "<<argv[i]<<endl; i+=1; continue;}
  }
  
  // Write out what the value of the parameters is
  cout << "Tstep:\t" << Tstep << "\n" << "startT:\t" << firstTraj << "\n"
       << "nprot:\t" << N_prot << "\tLx:\t" << L_box[0] << "\tLy:\t" << L_box[1] << "\tLz:\t" << L_box[2]
       << "\tNumber of bins in a:\t" << Nabin << "\tin c1:\t" << Nc1bin << "\tin c2:\t" << Nc2bin
       << "\nR_cut\t" << R_cut << "\nTrajectory file location:\t" << prefix << endl;
}

/**********************/
/* MAIN               */
/**********************/
int main(int argc, char * argv[]) {
  // Read user inputs
  cl_parse(argc, argv);
  cout.flush(); // Output information on read parameters
  
  cout << "Read User input" << endl;

/**********************/
/* INITIALIZE               */
/**********************/

  // Initialize; calculate derived variable and declare arrays used in the script
  double rc2 = R_cut*R_cut; // cut-off squared, used for calculating distances
  double rw[3] = {0}; // current water position
  double rp[N_prot][3] = {0}; // Protein positions
  memset( rp, 0, N_prot*3*sizeof(double) ); // Added to make the code compatible with C99 compilers, manually initialize the array because size is unknown upon compilation 
  char filename[1000]; // Name of the trajectory file
  int Nwater = 0; // Number of waters in the interface
  double r2 = 0.0; // Squared distance between a water and a protein atom
  double rmin2 = 0.0; // Minimum distance between a water and any protein atom

  // Variables that do not contain data to be processes, just help in loops and file IO
  string str1; // Will containt individual pdb file lines
  int counter=1; // Used to count in while loops etc
  int flagInclude=0; // Used to track if a water should be included in the output data
  int resInclude[N_prot] = {0}; // Read from input, contain a 1 for protein atoms that are being used to generate a new reference water structure file (waters will only be included if the nearest protein atom is one that is tagged by a 1 in this file)
  memset( resInclude, 0, N_prot*sizeof(int) );
  int tstr, n0, wID; // Info from coord_iface.txt file, to be copied to output (if water is within cutoff)
  double surfx, surfy, surfz, a0, c1, c2; // Info from coord_iface.txt file, to be copied to output (if water is within cutoff)

  // For histogramming local water structure
  int abin=1, c1bin=1, c2bin=1; // Discretized coordinates
  int whist[Nabin][Nc1bin][Nc2bin] = {0}; // Histogram with water structure
  memset( whist, 0, Nabin*Nc1bin*Nc2bin*sizeof(int) ); // Added to make the code compatible with C99 compilers, manually initialize the array because size is unknown upon compilation 
  

  double abinW = aRange/Nabin; // Binwdith = histogram_range/Nbins
  double c1binW = c1Range/Nc1bin;
  double c2binW = c2Range/Nc2bin;

  // Open input and output streams
  ifstream traj; // Protein trajectory
  ifstream wstruct("coord_iface.txt"); // Water structure
  // Open file for writing output
  FILE * wnear;
  wnear = fopen("local_wstructure.txt","w"); // Structure of nearby waters

/**********************/
/* READ INPUT FILE               */
/**********************/

  // Read user-provided file that contains which protein atoms we want to use to create the new reference water structure file (file with 1 or 0 for each protein atom, only water for which the nearest protein atom is tagged by a 1 will be retained)
  ifstream resIncludeF;
  sprintf(filename, "%sresInclude.txt", prefix);
  resIncludeF.open(filename);
  if (!resIncludeF) {cout << "ERROR: Failed to read information on which residues to include in the reference surface in: "<<filename<<endl; exit(1);} // Error if we can't read the file
  for ( int n = 0; n < N_prot; n++ ) {
    resIncludeF >> resInclude[n];
    if (verbose) {cout << "Protein atom " << n << " has include tag " << resInclude[n] << endl;}
  } 

  cout << "Finished initialization" << endl;

  // Iterate over the indicated time-slice
  for (int t = firstTraj; t<firstTraj+Tstep; t++)
    {
    cout << "Starting iteration "<<t<<endl;
/**********************/
/* READ PROTEIN INFO               */
/**********************/
    // Read the protein positions
    // File opening
    sprintf(filename, "%s%d.pdb", prefix, t); // Trajectory filename, if it is a pdb

    if (verbose) {cout << "Attempting to open trajectory at:\t" << filename << endl;}

    traj.open(filename);

    // Check for succesful file read
    if (!traj) { 
      cerr << "Failed to read trajectory file: "<<filename<<" trying .xyz "<<endl;
      sprintf(filename, "%s%d.xyz", prefix, t); // Trajectory filename, if it is an xyz
      traj.open(filename);
      if (!traj) {cerr << "ERROR: Failed to read trajectory file: "<<filename<<endl; exit(1);} // Fatal error, both fileformats didn't work
      counter = 0; // Reset counter
      // Skip comment lines
      getline(traj,str1); getline(traj,str1);
      //traj >> str1; traj >> str1;
      while (!traj.eof()) {
        if (counter+1 > N_prot) { cerr << "WARNING: Too many coordinates in file, "<<filename<<" skipping remaining lines"<<endl; break;}
        if (verbose) {cout << "Attempting to read line for atom"<<counter<<" from the provided xyz" << endl;}
        traj >> str1 >> rp[counter][0] >> rp[counter][1] >> rp[counter][2]; // XYZ is formatted by columns, not position
	if (verbose) {cout << rp[counter][0] << endl;}
        counter++; // Update the counter
      }
    } // Error out if we can't read the file as pdb, try xyz
    else { // The trajectory is pdb format
      // File reading
      counter = 0; // Reset counter
      while (!traj.eof()) {
        if (verbose) {cout << "Attempting to read a line from the provided pdb" << endl;}
        getline(traj,str1); // Get the entire content of this line
        if ( str1.substr(0,4) == "ATOM" ) { // This line contains coordinates we want
          // Check if we don't already have all coords we asked for
          if ( counter+1 > N_prot ) { cerr << "WARNING: Too many coordinates in file, "<<filename<<" skipping remaining lines"<<endl; break;}
          // Store coordinates
          if (verbose) {cout << "Trying to extract info from pdb ATOM line" << endl;}
          rp[counter][0] = std::stod (str1.substr(30,8));
          rp[counter][1] = std::stod (str1.substr(38,8));
          rp[counter][2] = std::stod (str1.substr(46,8));
          counter++; // Update the counter
        } // End of ATOM line if
      } // End of while

    } // If I pass this if else, I managed to read a protein trajectory
    if (verbose) {cout << "Succesfully opened trajectory file:\t" << filename << endl;}
    // Done with this file
    traj.close();

    if (verbose) {cout << "Finished reading protein coordinates" << endl;}

/**********************/
/* IDENTIFY NEARBY WATER               */
/**********************/
    // Distance calculation and file writing
    wstruct >> Nwater; // Number of interfacial waters in this frame

    cout << "Number of waters at time " << t << " are " << Nwater << endl;
    counter = 0;
    for ( int n = 0; n < Nwater; n++ ) { // For all the waters that could be included
      // Added important check; we didn't run out of water lines prematurely!
      if (wstruct.eof()) { cerr << "ERROR: unexpectedly ran out of water coordinates!" << endl; exit(1);}
      // I only need the water coordinates
      wstruct >> tstr >> n0 >> wID >> rw[0] >> rw[1] >> rw[2] >> surfx >> surfy >> surfz >> a0 >> c1 >> c2;
      // Check if the water is near enough to the surface to matter
      if ( (a0 > 8.0) or (a0 < -1.0) ) { continue; } // This water is too far away from, or too far inside(?), the surface to count
      // Initialize things for this water
      if (verbose) {cout << "Water "<< n << " is being considered" << endl;}
      if (verbose) {cout << "Water located at "<<rw[0] <<" "<< rw[1] <<" "<< rw[2] << endl;}
      flagInclude = 0; // The water will not be included in the output, unless this flag is set to 1 in the following loop
      rmin2 = rc2; // Initialized the cut-off distance, so that we do not include waters that are too far
      // Determine if any protein atom is within the cut-off of this water, and if so, which atom is nearest (new method) 
      for ( int i = 0; i < N_prot; i++ ) {
	if (verbose) { cout << "Protein at: "<< rp[i][0] << " " << rp[i][1] << " " << rp[i][2] << endl;}
        r2 = calc_rsq(rp[i], rw); // Distance between the water and this protein atom
        if (verbose) {cout << "This water is\t" << r2 << "\t away from protein atom " << i <<" Cut-off is "<<rmin2<< endl;}
        if (r2 < rmin2) { // New nearest water found
          if (cylinder) { // We want to only keep this water if it's nearest protein atom is a tagged residue
            rmin2 = r2;
	    if (verbose) {cout << "Using cylinder site assignment"<<endl;}
            //if (verbose) {cout <<i<<endl;}
	    if (verbose) {cout << "This water is tagged for inclusion based on "<<resInclude[i]<<endl;}
            flagInclude = resInclude[i]; // Whether the water is to be included depends on if this is one of the tracked protein residues
          } else { // We want to keep the water as long as it is near (within cut-off distance) of ANY tagged residue
            flagInclude += resInclude[i];
          }
        }
      }  // End of protein loop
      if ( flagInclude > 0 ) { // We want to keep info on this water, it's nearest protein residue is a tracked residue within the cut-off distance
        if (verbose) {cout << "Writing water info to output file" << endl;}
        if (raw) {fprintf(wnear, "%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", tstr, n0, wID, rw[0], rw[1], rw[2], surfx, surfy, surfz, a0, c1, c2);} // Write to file 
        if (hist) { // We want to output histograms of the local water structure
          abin = min(int((a0+1.0)/abinW),Nabin-1);
          c1bin = min(int((c1+1.0)/c1binW),Nc1bin-1);
          c2bin = min(int((c2+1.0)/c2binW),Nc2bin-1);  // Discretize internal coordinates for histogramming; Nbins/histogram_range = 1/binWidth
   	  whist[abin][c1bin][c2bin]++; // Add a count
        }
        counter++;
      }
    } // End of water loop

  if (verbose) {cout << "Finished iteration: " << t << endl;}
  cout << "Retained " << counter << " water molecules" << endl;
//  cout << "Current iteration is really "<<t<<endl;
//  cout << "Tstep is still known and is "<<Tstep<<endl;
//  cout << "firstTraj is still known and is "<<firstTraj<<endl;
//  cout << "Their sum is "<<firstTraj+Tstep<<endl;

  } // End of time loop
  cout << "Finished the loop" << endl;

/**********************/
/* WRITE               */
/**********************/

  // Write the histogram file, if required
  // Open file for writing local structure histogram
  if (hist) { // hist
  FILE * whistF;
  whistF = fopen("local_wstructure_histogram.txt","w");
    for ( int i = 0; i < Nabin; i++ ) {
      for ( int j = 0; j < Nc1bin; j++ ) {
        for ( int k = 0; k < Nc2bin; k++ ) {
          a0 = i*abinW - 1.0; c1 = j*c1binW - 1.0; c2 = k*c2binW -1.0; // Convert bin to internal coodinates; Nbins/histogram_range = 1/binWidth
          fprintf(whistF, "%f\t%f\t%f\t%d\n", a0, c1, c2, whist[i][j][k]); // Write coordinates and count to the histogram file
	  fflush(whistF);
        }
      }
    }
  fclose(whistF);
  }

  // Clean up before end
  wstruct.close();
  fclose(wnear);
  // MN:[XXX] Missing delete statements

  // Exit without error
  return 0;

} // Main end

/*************/
/* FUNCTIONS */
/*************/

/* calc_rsq(r1,r2); calculates the squared distance between two points in the simulation box, 
 * taking periodic boundaries into consideration. */
double calc_rsq ( double * r1, double * r2) {
  double r_sq = 0; // The output variable
  double dr = 0; // Used to store difference between vector components
  for (int i = 0; i < 3; i++ ) { // Loop over dimensions
    dr = fabs(r1[i] - r2[i]); // Difference between two points, absolute value
    if ( dr <= 0.5*L_box[i] ) { // Distances in a periodic box cannot be greater than 0.5*boxLength
      r_sq += dr*dr;
    }
    else { // Exception, need to correct dr
      dr = L_box[i]/2.0 - fmod(dr,L_box[i]/2.0); // Determine modulus, use that instead of dr
      r_sq += dr*dr;
    }
  }
  return r_sq;
} // End of calc_rsq function

