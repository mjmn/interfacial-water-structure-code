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
long int Tstep = 0; // Duration of the trajectory slice to be analyzed (in timesteps)
int firstTraj = 0; // First trajectory frame to be read

// Simulated system related
int N_prot = 0; // Number of atoms in the protein
double L_box[3] = {0}; // Size of the periodic box (make sure it doesn't change over time, so run in a constant Volume ensemble)

// Parameters for analysis
double R_cut = 10000.0; // Water needs to be within this cut-off of an atom to be included in output
char * prefix; // Name of the folder that contains the trajectory files

// For output
bool verbose=false; // Write a lot more info, for debugging
bool raw=false; // Write an output file containing all selected waters
bool hist=false; // Write histogram data on the local water structure

/**********************/
/* Declare functions  */
/**********************/
double calc_rsq(double * r1, double * r2); // Calculate the squared distance between two points

/**********************/
/* User input         */
/**********************/

void cl_parse(int argc, char * argv[]) {
  bool match;
  int i = 0;
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
    // Related to user desired output, these are flags that dont require further arguments
    if(strcmp(argv[i],"--raw")==0) {raw = true;match=true;i+=1; continue;}
    if(strcmp(argv[i],"--hist")==0) {hist = true;match=true;i+=1; continue;}
    if(strcmp(argv[i],"--verbose")==0) {verbose = true;match=true;i+=1; continue;}
    // Check if the user provided a bad input argument
    if(match == false) {cout<<"Warning: Bad input parameter: "<<argv[i]<<endl; i+=1; continue;}
  }
  
  // Write out what the value of the parameters is
  cout << "Tstep:\t" << Tstep << "\n" << "startT:\t" << firstTraj << "\n"
       << "nprot:\t" << N_prot << "\tLx:\t" << L_box[0] << "\tLy:\t" << L_box[1] << "\tLz:\t" << L_box[2]
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
  char filename[1000]; // Name of the trajectory file
  int Nwater = 0; // Number of waters in the interface
  double r2 = 0.0; // Squared distance between a water and a protein atom
  double rmin2 = 0.0; // Minimum distance between a water and any protein atom

  // Variables that do not contain data to be processes, just help in loops and file IO
  string str1; // Will containt individual pdb file lines
  int counter=1; // Used to count in while loops etc
  int flagInclude=0; // Used to track if a water should be included in the output data
  int resInclude[N_prot] = {0}; // Read from input, contain a 1 for protein atoms that are being used to generate a new reference water structure file (waters will only be included if the nearest protein atom is one that is tagged by a 1 in this file)
  int tstr, n0, wID; // Info from coord_iface.txt file, to be copied to output (if water is within cutoff)
  double surfx, surfy, surfz, a0, c1, c2; // Info from coord_iface.txt file, to be copied to output (if water is within cutoff)

  // For histogramming local water structure
  int abin=1, c1bin=1, c2bin=1; // Discretized coordinates
  int whist[91][101][101] = {0}; // Histogram with water structure

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
  if (!resIncludeF) {cerr << "ERROR: Failed to read information on which residues to include in the reference surface in: "<<filename<<endl; exit(1);} // Error if we can't read the file
  for ( int n = 0; n < N_prot; n++ ) {
    resIncludeF >> resInclude[n];
    if (verbose) {cout << "Protein atom " << n << " has include tag " << resInclude[n] << endl;}
  } 

  cout << "Finished initialization" << endl;

  // Iterate over the indicated time-slice
  for ( int t = firstTraj; t < firstTraj + Tstep; t++ ) {
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
      counter = 1; // Reset counter
      // Skip comment lines
      traj >> str1; traj >> str1;
      while (!traj.eof()) {
        if (verbose) {cout << "Attempting to read a line from the provided xyz" << endl;}
        if (counter > N_prot) { cerr << "WARNING: Too many coordinates in file, "<<filename<<" skipping remaining lines"<<endl; break;}
        traj >> str1 >> rp[counter][0] >> rp[counter][1] >> rp[counter][2]; // XYZ is formatted by columns, not position
	counter++; // Update the counter
      }
    } // Error out if we can't read the file as pdb, try xyz
    else { // The trajectory is pdb format
      // File reading
      counter = 1; // Reset counter
      while (!traj.eof()) {
        if (verbose) {cout << "Attempting to read a line from the provided pdb" << endl;}
        getline(traj,str1); // Get the entire content of this line
        if ( str1.substr(0,4) == "ATOM" ) { // This line contains coordinates we want
          // Check if we don't already have all coords we asked for
          if ( counter > N_prot ) { cerr << "WARNING: Too many coordinates in file, "<<filename<<" skipping remaining lines"<<endl; break;}
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
      flagInclude = 0; // The water will not be included in the output, unless this flag is set to 1 in the following loop
      rmin2 = rc2; // Initialized the cut-off distance, so that we do not include waters that are too far
      // Determine if any protein atom is within the cut-off of this water, and if so, which atom is nearest (new method) 
      for ( int i = 0; i < N_prot; i++ ) {
        r2 = calc_rsq(rp[i], rw); // Distance between the water and this protein atom
        if (verbose) {cout << "This water is\t" << r2 << "\t away" << endl;}
        if (r2 < rmin2) { // New nearest water found
          rmin2 = r2;
          flagInclude = resInclude[i]; // Whether the water is to be included depends on if this is one of the tracked protein residues
        }
      }  // End of protein loop
      if ( flagInclude == 1 ) { // We want to keep info on this water, it's nearest protein residue is a tracked residue within the cut-off distance
        if (verbose) {cout << "Writing water info to output file" << endl;}
        if (raw) {fprintf(wnear, "%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", tstr, n0, wID, rw[0], rw[1], rw[2], surfx, surfy, surfz, a0, c1, c2);} // Write to file 
        if (hist) { // We want to output histograms of the local water structure
          abin = int((a0+1.0)*10.0);
          c1bin = int((c1+1.0)*50.0);
          c2bin = int((c2+1.0)*50.0);  // Discretize internal coordinates for histogramming
          whist[abin][c1bin][c2bin]++; // Add a count
        }
        counter++;
      }
    } // End of water loop

  if (verbose) {cout << "Finished iteration:\t" << t << endl;}
  cout << "Retained " << counter << " water molecules" << endl;

  } // Time loop end

/**********************/
/* WRITE               */
/**********************/

  // Write the histogram file, if required
  // Open file for writing local structure histogram
  FILE * whistF;
  whistF = fopen("local_wstructure_histogram.txt","w");
  if (hist) {
    for ( int i = 0; i < 90; i++ ) {
      for ( int j = 0; j < 100; j++ ) {
        for ( int k = 0; k < 100; k++ ) {
          a0 = i/10.0 - 1.0; c1 = j/50.0 - 1.0; c2 = k/50.0 -1.0; // Convert bin to internal coodinates
          fprintf(whistF, "%f\t%f\t%f\t%d\n", a0, c1, c2, whist[i][j][k]); // Write coordinates and count to the histogram file
        }
      }
    }
  }

  // Clean up before end
  wstruct.close();
  fclose(wnear);
  fclose(whistF);

  // Exit without error
  exit(0);

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
      dr = fmod(dr,L_box[i]/2.0); // Determine modulus, use that instead of dr
      r_sq += dr*dr;
    }
  }
  return r_sq;
} // End of calc_rsq function

