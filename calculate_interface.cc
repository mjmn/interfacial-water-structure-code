/************ Locating the Intrinsic Water Interface for Protein (fixed)  *************/
/************ Taking the coordinates of water through the transformation ***********/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <iomanip>
//#include <mkl.h>
//#include <omp.h>

using namespace std;

//define parameters

#define PI 3.14159265358979
#define NUM_CORE 24 // number of processors to be used for parallelization

long int idum; // seed for random number generator

// Box-size related variables, will be re-set after reading user input
double L_box[3] = {150, 150, 150}; // {Lx, Ly, Lz}
double hL_box[3] = {75, 75, 75}; // 0.5*L_box (used in the calc_rsq function and elsewhere)

// MNiesen; 
double dx[3] = {1.0, 1.0, 1.0};                 // width of the lattice for computing coarse-grained water density
double cg_l = 2.4;                  // gaussian width for coarse-grained water density
double rho_bulk = 0.0336;           // bulk density (#/Angstrom^3)
double R_cut = 12.0; 	            // Include waters that are within this distance (A) of the proteins COM + Radius of Gyration

long int N_prot = 1231;             // Number of atoms in protein
long int N_water = 8537;            // Total number of water molecules
int N_res = 76;

long int Tstep = 20;
int sample = 100;
int N_sample;
int Ntraj = 1;
int firstTraj = 1; // MNiesen; introduced to allow for more easy splicing of trajectory frames (i.e., time-slicing) 

const double a_inc = 0.1;
const double a_range = 11.0;
const int a_bin = int(a_range/a_inc);
const double cos_inc = 0.02;
int cos_bin = int(2/cos_inc);
int a_bin_rho = 140;

//declare functions
double pbc_box(double x, int k);
int on_grid(int i, int k); // Same as pbc_box, but for the CG density grid
int safe_grid_index(int i, int j, int k); // Handle the fact that the grid does not necessarily span the entire simulation box
double pbc_disp(double x, int k);
double get_abs(double x);
double calc_rsq(double * v1, double * v2); 

// declaration of boolean variables
bool verbose=false; // Controls additional output
  

void cl_parse(int argc, char * argv[])
{
  bool match;
  int i = 1;
  while ( i < argc ) 
  //for(int i=1;i<argc;i+=2)
    {
      match = false;
      //if(strcmp(argv[i],"--N")==0){N = atoi(argv[i+1]);match=true;}
      if(strcmp(argv[i],"--rcut")==0){R_cut = atof(argv[i+1]);match=true;i+=2; continue;}
      if(strcmp(argv[i],"--T")==0){Tstep = atoi(argv[i+1]);match=true;i+=2; continue;}
      if(strcmp(argv[i],"--traj")==0){Ntraj = atoi(argv[i+1]);match=true;i+=2; continue;}
      if(strcmp(argv[i],"--dx")==0){dx[2] = atof(argv[i+1]);match=true;i+=2; continue;} // Desired lattice size, will be slightly adjusted later
      if(strcmp(argv[i],"--lx")==0){L_box[0] = atof(argv[i+1]);match=true;i+=2; continue;} // Box X-dimension - avg_slice.py
      if(strcmp(argv[i],"--ly")==0){L_box[1] = atof(argv[i+1]);match=true;i+=2; continue;} // Box Y-dimension - avg_slice.py
      if(strcmp(argv[i],"--lz")==0){L_box[2] = atof(argv[i+1]);match=true;i+=2; continue;} // Box Z-dimension - avg_slice.py
      if(strcmp(argv[i],"--cgl")==0){cg_l = atof(argv[i+1]);match=true;i+=2; continue;}
      if(strcmp(argv[i],"--bulk")==0){rho_bulk = atof(argv[i+1]);match=true;i+=2; continue;}
     // MNiesen; introduced several new input variables to make the program more flexible in dealing with varying systems
      if(strcmp(argv[i],"--startT")==0){firstTraj = atof(argv[i+1]);match=true;i+=2; continue;} // First trajectory frame to read
      if(strcmp(argv[i],"--nprot")==0){N_prot = atof(argv[i+1]);match=true;i+=2; continue;} // Number of protein ATOMS - avg_slice.py
      if(strcmp(argv[i],"--nres")==0){N_res = atof(argv[i+1]);match=true;i+=2; continue;} // Number of protein RESIDUES - avg_slice.py
      if(strcmp(argv[i],"--nwater")==0){N_water = atof(argv[i+1]);match=true;i+=2; continue;} // Number of water MOLECULES - avg_slice.py
      // MNiesen; flags that do not require other input arguments
      if(strcmp(argv[i],"--verbose")==0) {verbose = true;match=true;i+=1; continue;}
      // Catch for bad input
      if(match == false){cout<<"Warning: Bad input parameter: "<<argv[i]<<endl;i+=1; continue;}
    }

  // *** MN - Set L_box derivative variables
  hL_box[0] = 0.5*L_box[0]; // Half box length (used in PBC related functions, which are called often)
  hL_box[1] = 0.5*L_box[1];
  hL_box[2] = 0.5*L_box[2];

  // MNiesen - new step; adjust dx slightly, such that the simulation box can be divided into gridpoints of equal size
  // i.e., avoid weird boundary effects due to smaller voxels at the boundaries.
  for (int i=0; i < 3; i++) 
  {
    // Slightly change the lattice sizes in each dimension, to make the space discreetly divisible into lattice sites
    dx[i] = L_box[i]/round(L_box[i]/dx[2]); 
  }

  cout//<<"N:\t"<<N<<"\n"
      <<"Tstep:\t"<<Tstep<<"\n"
      <<"N_traj:\t"<<Ntraj<<"\n"
      <<"dx:\t"<<dx[0]<<"\n"
      <<"dy:\t"<<dx[1]<<"\n"
      <<"dz:\t"<<dx[2]<<"\n"
      <<"Lx:\t"<<L_box[0]<<"\n"
      <<"Ly:\t"<<L_box[1]<<"\n"
      <<"Lz:\t"<<L_box[2]<<"\n"
      <<"cg_l:\t"<<cg_l<<"\n"
      <<"rho_b:\t"<<rho_bulk<<"\n";
 if (verbose) { cout << "Verbose flag set, dumping info to standard output." << endl; }

}

int main(int argc, char * argv[]) {  
  // for timing
  clock_t startTime = clock();
  // parse for any new input parameters
  cl_parse(argc, argv);
  cout.flush();

  // declaration of parameters
  double r_p[3], r_c[3], grad[3], surf_norm[3], r2s[3];
  double solv_min[3], solv_max[3], center_of_mass[3];
  double R_g, R_hyd_sq, rho_p, disp_vec, r_sq, lx, ly, lz, d_o, d_n, d_OH, vec_OH, Norm;  
  int N_if_tot, N_iface, N_hyd_pre, N_nn, mol_id, count_hyd, count_iface, count_index, min_ind;
  int N_grid[3];               // Number of grid point (in each direction) for computing CG densities
  int nx, ny, nz, cnx, cny, cnz; // MN - New variables, discrete grid coordinate
  int bgrid[3]; // MN - New variable, buffer size around the grid that is needed to deal with PBC correctly

/* MNiesen; too soon, and Z-dimension ignored, moved down and corrected.
  N_grid[0] = int(L_box[0]/dx);
  N_grid[1] = int(L_box[1]/dx);

  N_if_tot = N_grid[0]*N_grid[1];
*/

  double R = 3.0*cg_l;           // Cutoff distance for Truncation
  double R2 = R*R;
  int maxind = ceil(R/dx[2])+1; // MN - New variable, used to optimize loop over the grid; the maximum number of neighboring grid points to consider
  double coeff = 0.5/(cg_l*cg_l); // Coefficient of Gaussian Tail
  double N_g = 0.00517322;        // Normalization Constant for CG Density function
  double rho_R = exp(-4.5);      // Truncation Constant for CG Density function
  double rho_c = rho_bulk*0.5;   // CG Density value for locating the interface

  // declaration of data arrays
  double ** H2O_pos;
  double ** Protein_pos;
  int * hyd_list;
  double *** cg_den;
  double ** iface;
  double ** coords;
  int * water_assignment;

  H2O_pos = new double * [3*N_water];	
  
  for (int n = 0; n < N_water; n++)
    {
      for (int j = 0; j < 3; j++)
	{
	  H2O_pos[3*n + j] = new double [3];
	}
    }

  Protein_pos = new double * [N_prot];	
  
  for (int n = 0; n < N_prot; n++)
    {
      Protein_pos[n] = new double [3];
    }

  hyd_list = new int [N_water];
  
  for (int n = 0; n < N_water; n++)
    {
      hyd_list[n] = 0;
    }

/* MNiesen; too soon and Z-dimension ignored, moved down and corrected
  iface = new double * [N_if_tot];

  for (int n = 0; n < N_if_tot; n++)
    {
      iface[n] = new double [6];   // interfacial position & the gradient in CG density (= -surf normal)
    }
*/
  
  // files to be written..
  ofstream logfile("log_iface.txt");

  logfile//<<"N:\t"<<N<<"\n"
	 <<"Tstep:\t"<<Tstep<<"\n"
	 <<"N_traj:\t"<<Ntraj<<"\n"
	 <<"dx:\t"<<dx[0]<<"\n"
	 <<"dy:\t"<<dx[1]<<"\n"
	 <<"dz:\t"<<dx[2]<<"\n"
	 <<"Lx:\t"<<L_box[0]<<"\n"
	 <<"Ly:\t"<<L_box[1]<<"\n"
	 <<"Lz:\t"<<L_box[2]<<"\n"
	 <<"cg_l:\t"<<cg_l<<"\n"
    	 <<"R:\t"<<R<<"\n"
	 <<"coeff:\t"<<coeff<<"\n"
	 <<"rho_b:\t"<<rho_bulk<<"\n"<<endl;
  
  logfile.flush();

  char file_prefix[50];
  int N_A, atom_type;
  int fformat=0; // MNiesen; new variable, tracks the format of the input trajectory file
  double t_num, x_p, y_p, z_p, N_p, tmp1, tmp2;
  double pos[3][3];
  string str1, str2, str3, str4, str5;  
  
  // open the trajectory file with a new stream	       
  ifstream traj;
  
  ofstream data_coord("coord_iface.txt");
  ofstream data_iface("traj_iface.xyz");

  
  /**** locate the interface & calculate the tranformed coordinate ****/
  for (int t = firstTraj; t < firstTraj+Tstep; t++)
    {
      // for debugging
      // cout << "reading frame: " <<  t << endL;
      // open and read the trajectory file (*.pdb)
      sprintf(file_prefix, "../traj/%d.pdb", t); // MNiesen; changed where to look for traj directory, this is temporary as it should not be hardcoded ***
      traj.open(file_prefix);
      
      if(!traj){
	cerr << "could not find file: " << file_prefix << " attempting xyz format instead." << endl;
        sprintf(file_prefix, "../traj/%d.xyz", t); // MNiesen; instead of immediately failing, see if the user is trying to read XYZ format instead
        traj.open(file_prefix);
        if(!traj){ // Also not XYZ, try GRO
	  cerr << "Also failed to find: " << file_prefix << ", trying gro format instead." << endl;
	  sprintf(file_prefix, "../traj/%d.gro",t);
	  traj.open(file_prefix);
	  if(!traj) { // Also not GRO, fail
	    cerr << "Also failed to find: "<<file_prefix << ", quitting.." << endl;
            exit(1);
	  } else {
	  fformat=2; // Update the file format to correctly read the GRO trajectory
	  }
        } else {
        fformat=1; // Update the file format to correctly read the XYZ trajectory
	}
      }
      else {
       fformat=0; // File is pdb format
      }

      if (verbose) { cout<< "File format code: " << fformat << endl;}

      // MNiesen; the following is meant to skip comment lines in the PDB. I don't like it, the 5 is hardcoded and assumes a specific number of commentlines
      // Ideally, commentlines need to be automatically recognized and skipped. Either use a available module to read PDB files, or pre-process so the number
      // of commentlines is known. ***
      if (fformat==0) {
	int istreampos; // For pdb file, we need to find where the atom coordinates start
        for (int i = 0 ; i < 5 ; i++) {  // PDB
	  getline(traj, str1);
	  if (str1.compare(0,4,"ATOM")==1) {
		  break;
		  traj.seekg(istreampos); // Reset where we are in the file to the first line before atoms
	  }
	  istreampos = traj.tellg(); // Keep track of where the atoms start in the pdb
        }
      } else { // XYZ or GRO
        for (int i = 0 ; i < 2 ; i++) {
	  getline(traj, str1);
        }
      }
      // End of comment skipping ***

      for (int i = 0; i < 3; i++)
	{	      
	  solv_min[i] = L_box[i];	  
	  solv_max[i] = 0;
	  center_of_mass[i] = 0;
	}
      
      // first loading protein atoms...
      for (int n = 0; n < N_prot; n++)
	{
          if (fformat==0) {
          // Read a formatted pdb file
            if (!getline(traj, str2)) {
	  	// Couldn't read line
		cerr << "Error reading protein coordinates" << endl;
		break;
 	    }
	    // The remainder of the pdb line is formatted, just get the coordinates and atom type
	    r_p[0] = std::stod (str2.substr(30, 8));
	    r_p[1] = std::stod (str2.substr(38, 8));
	    r_p[2] = std::stod (str2.substr(46, 8));
          } else if (fformat==1) { // XYZ file
            traj >> str1 >> r_p[0] >> r_p[1] >> r_p[2]; // XYZ is formatted by columns, not position
          } else { // GRO file, consider units are in nm, so convert to angstrom
            // Formatted consistently, so use getline
	    if (!getline(traj,str2)) {
	      cerr << "Error reading protein coordinates" << endl;
	      break;
	    }
	    // Get coordinates from the expected positions in the GRO file
	    r_p[0] = 10.0*std::stod (str2.substr(20, 8));
	    r_p[1] = 10.0*std::stod (str2.substr(28, 8));
	    r_p[2] = 10.0*std::stod (str2.substr(36, 8));
	  }

	  if (verbose) { cout << "Adding protein atom with coordinates: " << r_p[0] << " " << r_p[1] << " " << r_p[2] << endl; }

	  // check for the min or max position of atoms in the simulation box
	  for (int i = 0; i < 3; i++)
	    {
	      Protein_pos[n][i] = r_p[i]; 

	      center_of_mass[i] += r_p[i];
		
	      if (r_p[i] < solv_min[i])
		solv_min[i] = r_p[i];
	      if (r_p[i] > solv_max[i])
		solv_max[i] = r_p[i];
	    }
	}

      cout << "loaded in protein coordinates" << endl;
      if (verbose) { cout << "Solvent region currently set to " << solv_min[0] << " " << solv_max[0] << " " << solv_min[1] << " " << solv_max[1] << " " << solv_min[2] << " " << solv_max[2] << endl; }

      // Determine the center of mass      
      for (int i = 0; i < 3; i++)
	{
	  center_of_mass[i] /= (double)N_prot;
	}
	
      
      // calculate the radius of gyration
      R_g = 0;
      #pragma omp parallel for reduction(+:R_g)
      for (int n = 0; n < N_prot; n++)
	{
	  R_g += calc_rsq(Protein_pos[n], center_of_mass);	  
	}

      R_g /= N_prot;
      R_g = sqrt(R_g);       

      R_hyd_sq = (R_g + R_cut)*(R_g + R_cut);

      if (verbose) { cout << "R_hyd_sq set to: " << R_hyd_sq << endl; }
      
      count_hyd = 0;
      // *** count_index = atom_type;
      // then loading water molecules...
      for (int n = 0; n < N_water; n++)
	{
	  // *** count_index++;
	  for (int j = 0; j < 3; j++)
	    {
          // Check file type and read accordingly
          if (fformat==0) { // PDB file
            if (!getline(traj, str2)) {
		// Couldn't read line
		cerr << "Error reading water coordinates" << endl;
		break;
 	    }
	    // The remainder of the pdb line is formatted, just get the coordinates and atom type
	    // atom_type = std::stoi (str2.substr(22, 4)); // This is going to be unnecessary ***
	    pos[j][0] = std::stod (str2.substr(30, 8));
	    pos[j][1] = std::stod (str2.substr(38, 8));
	    pos[j][2] = std::stod (str2.substr(46, 8));
          } else if (fformat==1) { // XYZ file
            traj >> str1 >> pos[j][0] >> pos[j][1] >> pos[j][2]; // Simple column formatting
          } else { // GRO file, convert nm to Angstrom
            // Formatted consistently, so use getline
	    if (!getline(traj,str2)) {
	      cerr << "Error reading water coordinates" << endl;
	      break;
	    }
	    // Get coordinates from the expected positions in the GRO file
	    pos[j][0] = 10.0*std::stod (str2.substr(20, 8));
	    pos[j][1] = 10.0*std::stod (str2.substr(28, 8));
	    pos[j][2] = 10.0*std::stod (str2.substr(36, 8));
	  }
	}
	  
	  
	  for (int j = 0; j < 3; j++)
	    {
	      for (int k = 0; k < 3; k++)
		{
		  H2O_pos[3*n + j][k] = pos[j][k];
		}
	      
	    }

	  if (calc_rsq(H2O_pos[3*n], center_of_mass) < R_hyd_sq)	    
	    {
	      hyd_list[count_hyd] = n;
	      count_hyd++;
	    }	  
	}
      
      cerr << "loaded in water coordinates" << endl;
      N_hyd_pre = count_hyd;    // set the number of waters potentially included in the hydration shell
      cerr << "loaded in water coordinates: Using " << N_hyd_pre << " waters" << endl;
      
      // set the region of solvation and the number of grid points for that region
      for (int i = 0; i < 3; i++)
	{
	  // MN - Updated these to define a grid that is guaranteed to fit within the hydration water sphere, this avoids boundary effects
          solv_min[i] = int((solv_min[i]-R_cut/sqrt(2))/dx[i])*dx[i];
	  solv_max[i] = (int((solv_max[i]+R_cut/sqrt(2))/dx[i])+ 1)*dx[i];
          // Old
	  //solv_min[i] = int(solv_min[i]/dx[i])*dx[i] - R_cut;
	  //solv_max[i] = (int(solv_max[i]/dx[i])+ 1)*dx[i] + R_cut;
	  
	  // MN - Added a check that the solvent region does not go beyond the box boundary
	  if (solv_max[i] - solv_min[i] > L_box[i] ) { solv_min[i] = 0; solv_max[i] = L_box[i]; cerr << "Warning: Cut-off radius " << R_cut << " results in a water-grid that is larger than the simulation box, setting water-grid boundary to equal the box boundary." << endl; }

	  N_grid[i] = ceil((solv_max[i] - solv_min[i])/dx[i]);
          bgrid[i] = floor((L_box[i] + solv_min[i] - solv_max[i])/dx[i]); // Buffer, difference in size b/w grid and periodic box
          //bgrid[i] = 0; // ONLY UNCOMMENT TO TEST WITHOUT PROPER PBC BUFFER REGION (should give same result, but much slower)
	}
      
      cerr << solv_min[0] << "  " << solv_max[0] << endl;
      cerr << "system size is: " << N_grid[0]*dx[0] << "  " << N_grid[1]*dx[1] << "  " << N_grid[2]*dx[2] << endl;
      
      // MNiesen; iface array definitions moved here, since I now know how big the grid is
      N_if_tot = N_grid[0]*N_grid[1]*N_grid[2]; // Maximum number of grid points
      iface = new double * [N_if_tot];

      for (int n = 0; n < N_if_tot; n++)
      {
        iface[n] = new double [6]; // interfacial position & the gradient in CG density (=  -surf normal)
      }

      // declare the array for CG density
      cg_den = new double ** [N_grid[0]];
      for (int i = 0; i < N_grid[0]; i++)
	{
	  cg_den[i] = new double * [N_grid[1]];
	  for (int j = 0; j < N_grid[1]; j++)
	    {
	      cg_den[i][j] = new double [N_grid[2]];
	      for (int k = 0; k < N_grid[2]; k++)
		{
		  cg_den[i][j][k] = 0;
		}
	    }
	}

      // Report on timing
      cout << "Before: " << double( clock() - startTime ) / CLOCKS_PER_SEC << " seconds." << endl; 

      // *** 2018/11/16 - MN Reworking the cg_den loop; switching indeces to loop over all waters and consider
      // only a grid of defined size around each water. This should drastically reduce the size of the loops
      // and avoid several IF statements
      for (int n = 0; n < N_water; n++)
       { // First determine where in the grid this water is
	nx = floor( (H2O_pos[3*n][0]-solv_min[0]) / dx[0] );
	ny = floor( (H2O_pos[3*n][1]-solv_min[1]) / dx[1] );
	nz = floor( (H2O_pos[3*n][2]-solv_min[2]) / dx[2] );
        // Now loop over only those grid points that could be near enough to this water
        for (int idx = -maxind; idx <= maxind; idx++ )
   	 {
     	  cnx = safe_grid_index(nx+idx,N_grid[0],bgrid[0]); // Return a safe index, if possible
	  if (cnx < 0 ) { // There was no safe index possible
	    continue; // Skip considering all grid-points with this x-value
          }
	  for (int idy = -maxind; idy <= maxind; idy++ )
	   {
            cny = safe_grid_index(ny+idy,N_grid[1],bgrid[1]); // Return a safe index, if possible
	    if (cny < 0 ) { // There was no safe index possible
	      continue; // Skip considering grid-points with this y-value
	    }
	    for (int idz = -maxind; idz <= maxind; idz++ )
	     {
              cnz = safe_grid_index(nz+idz,N_grid[2],bgrid[2]); // Return a safe index, if possible
	      if (cnz < 0 ) { // There was no safe index possible
		continue; // Skip considering grid-points with this z-value
	      }
	      // Possible grid point to evaluate, calculate the distance between this point and the water
	      r_p[0] = pbc_box(solv_min[0] + cnx*dx[0], 0);
	      r_p[1] = pbc_box(solv_min[1] + cny*dx[1], 1);
	      r_p[2] = pbc_box(solv_min[2] + cnz*dx[2], 2);
	      // rsq calculation
	      r_sq = calc_rsq(r_p, H2O_pos[3*n]);
	      if (r_sq < R2) {
	       cg_den[cnx][cny][cnz] += N_g*(exp(-r_sq*coeff) - rho_R);
	      }
	     }
	   }
	 }
       }
     // Report on timing
     cout << "CG density: " << double( clock() - startTime ) / CLOCKS_PER_SEC << " seconds." << endl; 
        

      count_iface = 0;
      // loop for locating interface within the region of solvation
      for (int i = 1; i < N_grid[0] - 1; i++)
	{
	  for (int j = 1; j < N_grid[1] - 1; j++)
	    {
	      for (int k = 1; k < N_grid[2] - 1; k++)
		{
		  r_p[0] = pbc_box(solv_min[0] + i*dx[0], 0);
		  r_p[1] = pbc_box(solv_min[1] + j*dx[1], 1);
		  r_p[2] = pbc_box(solv_min[2] + k*dx[2], 2);

		  rho_p = cg_den[i][j][k];
		  if (rho_p < rho_c)
		    {
		      N_nn = 0;
		      // check whether cg_den(nn(i,j,k)) > rho_c
		      if (cg_den[i+1][j][k] > rho_c)
			N_nn++;
		      if (cg_den[i-1][j][k] > rho_c)
			N_nn++;
		      if (cg_den[i][j+1][k] > rho_c)
			N_nn++;
		      if (cg_den[i][j-1][k] > rho_c)
			N_nn++;
		      if (cg_den[i][j][k+1] > rho_c)
			N_nn++;
		      if (cg_den[i][j][k-1] > rho_c)
			N_nn++;
		      
		      if (N_nn > 0)
			{
			  // calculation of gradient
			  grad[0] = (cg_den[i+1][j][k] - cg_den[i-1][j][k])*0.5/dx[0];
			  grad[1] = (cg_den[i][j+1][k] - cg_den[i][j-1][k])*0.5/dx[1];
			  grad[2] = (cg_den[i][j][k+1] - cg_den[i][j][k-1])*0.5/dx[2];

			  r_sq = 0; 
			  for (int n = 0; n < 3; n++)
			    {
			      r_sq += grad[n]*grad[n];
			    }
			  
			  disp_vec = (rho_c - rho_p)/r_sq;

			  r_sq = sqrt(r_sq); // magnitude of gradient
			  for (int n = 0; n < 6; n++)
			    {
			      if (n < 3)				
				iface[count_iface][n] = pbc_box(r_p[n] + disp_vec*grad[n], n);
			      else				
				iface[count_iface][n] = -grad[n-3]/r_sq;  // surface normal vector
			    }
			  
			  count_iface++;
			}
		    }
		}
	    }
	}

      N_iface = count_iface;


      /******* Calculating the new coordinates *******/     

      // declaring the coordinate array for the hydration-shell waters
      coords = new double * [N_hyd_pre];
      // declaring the array that tracks for each water to which surface grid-site it was assigned
      water_assignment = new int [N_hyd_pre];
      for (int n = 0; n < N_hyd_pre; n++)
	{
	  coords[n] = new double [6]; 
	}

      // locate the interfacial position nearest from each molecule
      #pragma omp parallel for private(r_p, mol_id, surf_norm, r2s, d_o, d_n, min_ind, d_OH, vec_OH)
      for (int n = 0; n < N_hyd_pre; n++)
	{
	  mol_id = hyd_list[n];
	  d_o = 3.0*L_box[2]*L_box[2];
	  for (int i = 0; i < N_iface; i++)
	    {
	      for (int j = 0; j < 3; j++)
		{
		  r_p[j] = iface[i][j];
		}
	      
	      d_n = calc_rsq(r_p, H2O_pos[3*mol_id]);

	      if (d_o > d_n)
		{
		  d_o = d_n;
		  min_ind = i;
		}		
	    }
	  
	  // MN - Store which surface site the water was assigned to
	  water_assignment[n] = min_ind;

	  for (int k = 0; k < 3; k++)
	    {
	      surf_norm[k] = iface[min_ind][3+k];
	      coords[n][k] = surf_norm[k];
	      r2s[k] = iface[min_ind][k] - H2O_pos[3*mol_id][k];  // vector from water to the interfacial point
	    }

	  // get the proximity of a particle to the interface
	  coords[n][3] = 0;
	  for (int k = 0; k < 3; k++)
	    {
	      coords[n][3] += pbc_disp(r2s[k], k)*surf_norm[k];
	    }

	  // get the angular coordinates of a particle w.r.t. surface normal
	  for (int j = 1; j < 3; j++)
	    {
	      coords[n][3+j] = 0;
	      d_OH = 0;
	      for (int k = 0; k < 3; k++)
		{
		  vec_OH = pbc_disp(H2O_pos[3*mol_id + j][k] - H2O_pos[3*mol_id][k], k);
		  coords[n][3+j] += surf_norm[k]*vec_OH;
		  d_OH += vec_OH*vec_OH;      // calculation for the actual bond length
		}
	      d_OH = sqrt(d_OH);
	      coords[n][3+j] /= d_OH;
	    }
	      
	}
	  
      /******* save the interfacial properties into files *******/
      
      // write the surface normal and the new coordinate
      data_coord << N_hyd_pre << endl;
      for (int n = 0; n < N_hyd_pre; n++)
	{
	  mol_id = hyd_list[n];
	  // MN - Output modified to contain the surface grid site to which the water was assigned
	  // data_coord << t << " " << n << " " << mol_id + N_res <<" " << H2O_pos[3*mol_id][0] << " " << H2O_pos[3*mol_id][1] << " " << H2O_pos[3*mol_id][2] << " " << coords[n][0] << " " << coords[n][1] << " " << coords[n][2] << " " << coords[n][3] << " " << coords[n][4] << " " << coords[n][5] << endl;
	  data_coord << t << " " << n << " " << water_assignment[n] <<" " << H2O_pos[3*mol_id][0] << " " << H2O_pos[3*mol_id][1] << " " << H2O_pos[3*mol_id][2] << " " << coords[n][0] << " " << coords[n][1] << " " << coords[n][2] << " " << coords[n][3] << " " << coords[n][4] << " " << coords[n][5] << endl;
	}
      data_coord << '\n';      

      // write the coordinates of instantaneous interface
      // MNiesen; changed to make the trajectories much smaller. VMD can read xyz files with changing number of atoms per frame using readvarxyz
      data_iface << N_iface << endl;
      //data_iface << N_if_tot << endl;
      data_iface << "Surface. Timestep: " << t << endl;
	
      for (int i = 0; i < N_iface; i++)
	{
	  data_iface << 0 << " " << iface[i][0] << " " << iface[i][1] << " " << iface[i][2] << endl;	    
	}
      // MNiesen; removed the placement of dummy atoms, just use readvarxyz in vmd
      // putting dummy surface positions (for vmd-readable traj)
      /*for (int i = N_iface; i < N_if_tot; i++)
	{
	  data_iface << "0 0 0 0" << endl;
	}
      */
      // record on the log file
      if (t % 2 == 0)
	{
	  logfile << "Calculation is currently  " << fixed << setprecision(3) << ((float)t-(float)firstTraj)/(float)Tstep*100 << " % done.\n";
	  logfile.flush();
	}

      // de-allocation of memory
      for (int n = 0; n < N_hyd_pre; n++)
	{
	  delete[] coords[n];
	}
      delete[] coords;      

      for (int i = 0; i < N_grid[0]; i++)
	{
	  for (int j = 0; j < N_grid[1]; j++)
	    {
	      delete[] cg_den[i][j];
	    }
	  delete[] cg_den[i];
	}
      delete[] cg_den;

      
      traj.close();
      
    }
      

  data_coord.close();
  data_iface.close();
      
  logfile << "Complete! \n";
  logfile.close();
  
  for (int i = 0; i < N_if_tot; i++)
    {
      delete[] iface[i];
    }

  delete[] iface;

  for (int n = 0; n < N_water; n++)
    {
      for (int j = 0; j < 3; j++)
	{
	  delete[] H2O_pos[3*n+j];
	}
    }

  for (int n = 0; n < N_prot; n++)
    {
      delete[] Protein_pos[n];
    }
  delete[] Protein_pos;
  delete[] H2O_pos;
  delete[] hyd_list;
  
  // Report on timing
  cout << "Completed in: " << double( clock() - startTime ) / CLOCKS_PER_SEC << " seconds." << endl; 
  
  
}

// Function to handle periodic boundary effects; atoms are placed in a box bounded by
// 0 and L_box[k]. Here, k is the dimension of the coordinate (0,1, or 2)
double pbc_box(double x, int k) {

// New implementation, can handle multiple box sizes displaced
  x = fmod( ( L_box[k] + fmod( x , L_box[k] ) ) , L_box[k] );

// Original implementation
/*  if (x >= L_box[k])
    x -= L_box[k];
  else if (x < 0)
    x += L_box[k];
*/    
  return x;
}

// New function, given integer i and maximum value k. Return an integer that is on the line [0,k] (using periodic boundary conditions)
int on_grid(int i, int k) {
  i = ( k + ( i % k ) ) % k;
return i;
}

// New function, given a grid-index value return a "safe" grid index value.
// In effect, this handles grid-indices that are negative or that go beyond the maximum
// grid index. Not strictly by PBC, but considering the fact that the grid does not span 
// the entire simulation box; i.e., it is periodic with a buffer region (k)
// i is the un-safe grid index, j is the grid dimension to consider, and k the size of the buffer
// |---k---|-------j-------|.... with i a point on that line
int safe_grid_index(int i, int j, int k) {
  if ( ( i < 0 ) || ( i >= j ) ) { // If unsafe, we need to process
    if ( i < -k ) { // It is negative enough to wrap around, considering the buffer region
      return on_grid(i+k,j); // Return safe index
    } else if ( i >= (j+k) ) { // It is positive enough to wrap around
      return on_grid(i-k,j); // Return safe index
    } else {
      return -1; // We are not on the grid!
    }
  } else { // it was already a safe index
    return i;
  }
}

// Need comment
double pbc_disp(double x, int k) {

  if (x > 0.5*L_box[k])
    x -= L_box[k];
  else if (x < -0.5*L_box[k])
    x += L_box[k];

  return x;
}

// Get the absolute value of double x
double get_abs(double x) {
  
  if (x < 0)
    x *= -1.0;
  
  return x;
  
}

// Calculate the squared distance between two vectors, this function does
// not allow for distances greater than half a box length in any direction
// if you just want a squared distance, without PBC consideration, a new 
// function is needed!
double calc_rsq (double * v1, double * v2) {

  // New version, considers PBC in all dimensions and not just XY
  // also uses standard min template, which is preferable
  double r_sq = 0.0; // Initialize
  double dv; // Used to store absolute distance between two points
  for ( int i = 0; i < 3; i++) // Loop over dimensions
  { // Add the minimum distance, either between points or through box boundary
    dv = get_abs(v1[i]-v2[i]);
    if ( dv < hL_box[i] ) // MN - use hL_box here, which is 0.5*L_box (used often so globally defined)
     {
      r_sq += dv*dv;
     }
    else
     {
      dv = hL_box[i] - fmod(dv, hL_box[i]);
      r_sq += dv*dv;
      //r_sq += ( L_box[i] - dv )*( L_box[i] - dv );
     }
  }
  return r_sq;

  // Old version
  /*double r_sq = (H2O_pos1[2] - r_p[2])*(H2O_pos1[2] - r_p[2]);

  for (int i = 0; i < 2; i++)
    {
      if (get_abs(r_p[i] - H2O_pos1[i]) < L_box[i]*0.5)
	{
	  r_sq += (H2O_pos1[i] - r_p[i])*(H2O_pos1[i] - r_p[i]);
	}
      else
	{
	  r_sq += pow(L_box[i] - get_abs(r_p[i] - H2O_pos1[i]), 2);
	}
    }
  return r_sq;
  */
}
