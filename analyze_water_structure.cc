/************ Computation of Hydrophilicity Order Parameter, q_lv for Protein Surface *************/
/************ based on the residue id and the distance between water and each residue *************/

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

long int idum; // seed for random number generator

double L_box[3] = {90, 90, 90}; // {Lx, Ly, Lz}

long int Tstep = 1;
int sample = 1;
int N_sample;
int firstTraj = 1; // MNiesen; new input argument to allow time-slicing

long int N_prot = 9060;             // Number of atoms in protein
int N_res = 578;                    // Number of residues

double R_cut = 6.0;            // cutoff distance for partitioning water in q_lv

// MN - Restructured the histogramming to allow user to change the number of bins
// As part of a larger structural change, we now only load one file wrt water-structure probability
// it combines all other factors (i.e., normalization by isotropic angle probabilities) in that one input file via pre-processing. This increases generality and flexibility.
double aRange = 9.0; // Distnace from the W-C interface ranges between -1.0 and 8.0
double c1Range = 2.0; // Cosine values range between -1.0 and 1.0
double c2Range = 2.0;
// Number of bins in each dimension, default values
int Nabin = 90; // Number of bins in the coordinate that describes the distance to the W-C interface
int Nc1bin = 100; // Number of bins in the coordinate that describes the O-H1 bond angle
int Nc2bin = 100; // Number of bins in the coordinate that describes the O-H2 bond angle
// Derived variables are handled during initialization in main

/*double a_range = 7.0;
double a_inc = 0.1;

int a_bin = int(a_range/a_inc);
double cos_inc = 0.02;
int cos_bin = 100;
int a_bin_rho = 140;
*/

//declare functions
double pbc(double x, int k);
double get_abs(double x);
int int_round(double x);
double calc_rsq(double * r1, double * r2);
/* MN - Removed; guessing water structure probabilities by interpolation. This should also be handled by pre-processing the water structure probability file, or by simply not having too many bins.
double guess_plv(double *** p_cos, int a_index, int c1_index, int c2_index);
double guess_pbulk(double ** p_cos, int c1_index, int c2_index);
*/

// Defining some user options
bool verbose = false; // Print a bunch of additional information on runtime (mostly for debugging)
bool cylinder = false; // A different mode for assigning waters to residues, each water is only assigned to the nearest residue in the cylinder mode
bool water_dump = false; // For each frame, write out the q-value contribution of each water molecule + the grid surface index to which that water was assigned

void cl_parse(int argc, char * argv[])
{
  int i=1;
  bool match;
 // for(int i=1;i<argc;i+=2)
    while ( i < argc ) // New way for doing user inputs, better for handling flags and wrong inputs
    {
      match = false;
      //if(strcmp(argv[i],"--N")==0){N = atoi(argv[i+1]);match=true;}
      if(strcmp(argv[i],"--T")==0){Tstep = atoi(argv[i+1]);match=true;i+=2;continue;}
      if(strcmp(argv[i],"--ts")==0){sample = atoi(argv[i+1]);match=true;i+=2;continue;}
      if(strcmp(argv[i],"--rcut")==0){R_cut = atof(argv[i+1]);match=true;i+=2;continue;}
      if(strcmp(argv[i],"--lx")==0){L_box[0] = atof(argv[i+1]);match=true;i+=2;continue;}
      if(strcmp(argv[i],"--ly")==0){L_box[1] = atof(argv[i+1]);match=true;i+=2;continue;}
      if(strcmp(argv[i],"--lz")==0){L_box[2] = atof(argv[i+1]);match=true;i+=2;continue;}
      // MNiesen; new input arguments, to allow time-splicing and more flexibility in the system definitions
      if(strcmp(argv[i],"--startT")==0){firstTraj = atof(argv[i+1]);match=true;i+=2;continue;} // Start time for analysis
      if(strcmp(argv[i],"--nprot")==0){N_prot = atof(argv[i+1]);match=true;i+=2;continue;} // Number of protein ATOMS - avg_slice.py
      if(strcmp(argv[i],"--nres")==0){N_res = atof(argv[i+1]);match=true;i+=2;continue;} // Number of protein RESIDUES - avg_slice.py
      // Related to histogramming the water data
      if(strcmp(argv[i],"--nabin")==0) {Nabin = atof(argv[i+1]);match=true;i+=2; continue;} // Number of bins for the water-structure histogram in each CV dimension
      if(strcmp(argv[i],"--nc1bin")==0) {Nc1bin = atof(argv[i+1]);match=true;i+=2; continue;}
      if(strcmp(argv[i],"--nc2bin")==0) {Nc2bin = atof(argv[i+1]);match=true;i+=2; continue;}
      // Boolean flags
      if(strcmp(argv[i],"--verbose")==0){verbose=true;match=true;i+=1;continue;} // Write more output to screen, mostly for debugging
      if(strcmp(argv[i],"--cylinder")==0){cylinder=true;match=true;i+=1;continue;} // Assign each water only to the nearest protein residue
      if(strcmp(argv[i],"--wdump")==0){water_dump=true;match=true;i+=1;continue;} // Write the q-value contribution of each water for each frame
      if(match == false){cout<<"Warning: Bad input parameter: "<<argv[i]<<endl;i+=1;continue;}
    }

  cout<<"Tstep:\t"<<Tstep<<"\n"
      <<"sample_time:\t"<<sample<<"\n"
      <<"R_cut:\t"<<R_cut<<"\n"
      << "\tNumber of bins in a:\t" << Nabin << "\tin c1:\t" << Nc1bin << "\tin c2:\t" << Nc2bin
      <<"\nLx:\t"<<L_box[0]<<"\n"
      <<"Ly:\t"<<L_box[1]<<"\n"
      <<"Lz:\t"<<L_box[2]<<"\n";
  
}

int main(int argc, char * argv[]) { 
  // For timing
  clock_t startTime = clock();
  // parse for any new input parameters
  cl_parse(argc, argv);
  cout.flush();

  // initialize - MN [XXX] Clean this up, removing variables that are no longer used
  double r_p[3], pos_tmp[3], surf_norm[3], int_coord[3];
  double tmp1, tmp2, r_sq, rsq_min, a0, c1, c2, N0, rho_val, a_tmp, q_val, plv_tmp, pbulk_tmp;
  int N_hyd, count_surf, min_ind, atom_id, mol_id, a_index, c1_index, c2_index, a_ind_rho;

  double Rcut_sq = R_cut*R_cut;
  
  char file_prefix[50];

  // Hisogram initialization, replacing previous method with something more straightforward
  double whist[Nabin][Nc1bin][Nc2bin] = {0}; // to contain the histogram w water structure, double incase it was normalized beforehand etc.
  memset( whist, 0, Nabin*Nc1bin*Nc2bin*sizeof(int) ); // Initialize memory for the variable size hisogram, compatible w C99
  // Derived properties
  double abinW = aRange/Nabin; // Bindiwth = histogram_range/Nbins
  double c1binW = c1Range/Nc1bin;
  double c2binW = c2Range/Nc2bin;

  //N_sample = (Tstep - 1)/sample;      

  // declaration of arrays
  double ** Protein_pos;
  double ** qlv_res;
  int * res_prot;

  Protein_pos = new double * [N_prot];
  qlv_res = new double * [N_res];
  res_prot = new int [N_prot];
 
  // load the residue information
  ifstream data_res("res_list.txt");  // containing the columns of atom id and corresponding residue id
  
  if(!data_res){
    cerr << "could not find file: res_list.txt" << endl;
    exit(1);
  }

  for (int n = 0; n < N_prot; n++)
    {
      Protein_pos[n] = new double [3];

      data_res >> atom_id >> mol_id;
      
      res_prot[n] = mol_id;
    }

  if (verbose) {cout << "Read in residue information"<<endl;}
  data_res.close();
  
  for (int n = 0; n < N_res; n++)
    {
      qlv_res[n] = new double [3];

      for (int k = 0; k < 3; k++)
	{
	  qlv_res[n][k] = 0;
	}
    }

  // MN - New read of water-structure histogram

  // load the joint prob distibutions of liq-vap interface...
  ifstream p_lv ("../ref_structure/pjoint.txt"); // MNiesen; changed the source directory location, should not be hardcoded ***
  if(!p_lv){
    cerr << "could not find file: ../ref_structure/pjoint.txt" << endl;
    exit(1);
   }
  
  for (int a = 0; a < Nabin; a++) {      
      for (int y1 = 0; y1 < Nc1bin; y1++) {
	  for (int y2 = 0; y2 < Nc2bin; y2++) {
	      p_lv >> a0 >> c1 >> c2 >> whist[a][y1][y2] >> N0;
	      // Since I went with flexible histograms, I should add a check here that my file matches user specifications [XXX]
	    }
	}	  
    }
  
  p_lv.close();
    
  // MN - Commented out old initialization of water structure histograms, due for deletion [XXX]
  /*
  // load the joint prob distibutions of bulk...
  ifstream p_bulk ("../ref_structure/pbulk_w2_tip3p.txt"); // MNiesen; changed the source directory location, should not be hardcoded ***
  if(!p_bulk){
    cerr << "could not find file: ../tip3p/pbulk_w2_tip3p.txt" << endl;
    exit(1);
  }
  

  double ** pcos_bulk;

  pcos_bulk = new double * [cos_bin];
  
  for (int y1 = 0; y1 < cos_bin; y1++)
    {
      pcos_bulk[y1] = new double [cos_bin];
	  
      for (int y2 = 0; y2 < cos_bin; y2++)
	{
	  p_bulk >> c1 >> c2 >> pcos_bulk[y1][y2];
	}
    }
  p_bulk.close();

// load the joint prob distibutions of liq-vap interface...
  ifstream p_lv ("../ref_structure/pjoint_sm2.txt"); // MNiesen; changed the source directory location, should not be hardcoded ***
  if(!p_lv){
    cerr << "could not find file: ../tip3p/pjoint_sm2.txt" << endl;
    exit(1);
  }

  
  double *** pcos_lv;
  double * pofa;
  
  pcos_lv = new double ** [a_bin];
  
  for (int a = 0; a < a_bin; a++)
    {      
      pcos_lv[a] = new double * [cos_bin];

      for (int y1 = 0; y1 < cos_bin; y1++)
	{
	  pcos_lv[a][y1] = new double [cos_bin];
	  
	  for (int y2 = 0; y2 < cos_bin; y2++)
	    {
	      p_lv >> a0 >> c1 >> c2 >> pcos_lv[a][y1][y2] >> N0;
	    }
	}	  
    }
  
  p_lv.close();

  ifstream pa_lv ("../ref_structure/pofa.txt"); // MNiesen; changed the source directory location, should not be hardcoded ***
  if(!p_lv){
    cerr << "could not find file: ../tip3p/pofa.txt" << endl;
    exit(1);
  }  

  pofa = new double [a_bin_rho];
  
  for (int a = 0; a < a_bin_rho; a++)
    {
     // MNiesen; changed this line due to changes in the pofa file format (it now has only 2 columns)
     // pa_lv >> a0 >> N0 >> pofa[a];
     pa_lv >> a0 >> pofa[a];
    }

  pa_lv.close();
  */

  if (verbose) {cout << "Read in reference water structure information" << endl;}
  
  // files to be written..
  sprintf(file_prefix, "qlv_avg_res_%1.f.txt", R_cut);
  ofstream data_qlv(file_prefix);
  ofstream logfile("log_qlv.txt");
  ofstream wdump_file("qval_per_surface_site.txt");

  logfile<<"Tstep:\t"<<Tstep<<"\n"
	 <<"sample_time:\t"<<sample<<"\n"
	 <<"R_cut:\t"<<R_cut<<"\n"
	 <<"Lx:\t"<<L_box[0]<<"\n"
	 <<"Ly:\t"<<L_box[1]<<"\n"
	 <<"Lz:\t"<<L_box[2]<<"\n"<<endl;
  
  logfile.flush();

  int atom_type;
  int fformat=0; // MNiesen; new variable to allow different filetypes to be recognized
  double N_A, t_num;
  double pos[3][3];
  string str1, str2, str3, str4, str5;
  

  // open the file that contains waters' interfacial coordinates
  ifstream traj;
  ifstream data_coord("coord_iface.txt");        
  
  /************************ Start the Analysis *************************/
  for (int t = firstTraj; t < firstTraj+Tstep; t++)
    {
      // open and read the trajectory file (*.pdb)
      sprintf(file_prefix, "../traj/%d.pdb", t); // MNiesen; changed the source directory location, should not be hardcoded ***
      traj.open(file_prefix);

      // MNiesen; Check if the file was read succesfully, if not try XYZ format
      if (!traj) {
        sprintf(file_prefix, "../traj/%d.xyz", t);
        traj.open(file_prefix);
        if (!traj) { // It isn't an XYZ file, try GRO
          sprintf(file_prefix, "../traj/%d.gro", t);
	  traj.open(file_prefix);
	  if (!traj) { // Not a GRO file either, give up
            cerr << "Failed to read trajectory file: " << file_prefix << endl;
	    exit(1);
	  } else {
	  fformat=2; // We're reading a GRO file
	  }
        } else {
        fformat=1; // We're reading a XYZ file
	}
      }
      else {
        fformat=0; // We're reading a PDB file
      }

      if (verbose) {cout << "File format code: " << fformat << endl;}
     
      // [XXX] MNiesen; hardcoded the number of commentlines here, appropriate for GROMACS (trjconv) generated pdb. Should be automated ***
      if (fformat==0) {
        for(int i = 0 ; i < 5 ; i++){
          getline(traj, str1);
        }
      }
      else { // Comment skipping for XYZ or GRO
        for(int i = 0; i < 2; i++) {
          getline(traj, str1);
        }
      }
      // *** End comment skipping
      
      // loading protein atoms...
      for (int n = 0; n < N_prot; n++)
	{
          // First check for the type of file we are reading
          if (fformat==0) { // PDB
            if (!getline(traj, str2)) {
		// Couldn't read line
		cerr << "Error reading protein coordinates from the PDB file" << endl;
		break;
 	    }
	    // The remainder of the pdb line is formatted, just get the coordinates and atom type
	    r_p[0] = std::stod (str2.substr(30, 8));
	    r_p[1] = std::stod (str2.substr(38, 8));
	    r_p[2] = std::stod (str2.substr(46, 8));
	  } else if (fformat==1) { // XYZ format
            traj >> str1 >> r_p[0] >> r_p[1] >> r_p[2]; // Column format
          } else { // GRO format
            // Use standard GRO formatting, and convert nm to Angstrom
	    if (!getline(traj,str2)) {
	      cerr << "Error reading protein coordinates from the GRO file" << endl;
	      break;
	    }
	    r_p[0] = 10.0*std::stod (str2.substr(20,8));
	    r_p[1] = 10.0*std::stod (str2.substr(28,8));
	    r_p[2] = 10.0*std::stod (str2.substr(36,8));
	  }

          // Store protein coordinates
	  for (int i = 0; i < 3; i++)
	    {
	      Protein_pos[n][i] = r_p[i];
	    }
	}

      traj.close();
      if (verbose) {cout << "Read protein coordinates" << endl;}
 
      // read the intrinsic coordinate file
      data_coord >> N_hyd;      
      for (int n = 0; n < N_hyd; n++)
	{
	  data_coord >> tmp1 >> tmp2 >> mol_id;

	  for (int j = 0; j < 3; j++)
	    {
	      data_coord >> pos_tmp[j];
	    }
	  for (int j = 0; j < 3; j++)
	    {
	      data_coord >> surf_norm[j];
	    }
	  for (int j = 0; j < 3; j++)
	    {
	      data_coord >> int_coord[j];
	    }
	  
	  // MNiesen: The following block of code has been significantly updated, to allow for flexible numbers of histogram bins and to only take in a single water-structure histogram instead of three separate files
	  // Update the histograms for interfacial structure, q(x,y)
	  if (int_coord[0] < 3.5 && int_coord[0] >= -1.0)
	    {
              if (verbose) {cout<<"Adding info for water "<<n<<" with a "<<int_coord[0]<<" to calculate local water structure"<<endl;}
	      rsq_min = Rcut_sq;
	      min_ind = N_prot + 1;
              //#pragma omp parallel for private(r_sq)
	      for (int i = 0; i < N_prot; i++)
		{
		  r_sq = calc_rsq(Protein_pos[i], pos_tmp);

		  if (r_sq <= rsq_min)
		    {			     
		      a_index = min(int((int_coord[0] + 1.0)/abinW),Nabin-1);
		
		      //a_ind_rho = int((int_coord[0] + 4.0)*10.0);
			  
		      c1_index = min(int((int_coord[1] + 1.0)/c1binW),Nc1bin-1);
		      c2_index = min(int((int_coord[2] + 1.0)/c2binW),Nc2bin-1);

		      /*
		      if (c1_index==cos_bin)
			c1_index--;
		      if (c2_index==cos_bin)
			c2_index--;
		      */

		      plv_tmp = whist[a_index][c1_index][c2_index];
		      //pbulk_tmp = pcos_bulk[c1_index][c2_index];
		      if (verbose) {cout<<"Succesful look up of water structure probability in the reference "<<plv_tmp<<endl;}
		      if (verbose) {cout<<"Bin index "<<a_index<<" "<<c1_index<<" "<<c2_index<<endl;}
		      if (verbose) {cout<<"Raw value "<<int_coord[0]<<" "<<int_coord[1]<<" "<<int_coord[2]<<endl;}
	              /* MN: Guessing the population of a bin by interpolation should not be necessary	
		      if (plv_tmp == 0)
			{
			  plv_tmp = guess_plv(whist, a_index, c1_index, c2_index);
			}*/
		      /*if (pbulk_tmp == 0)
			{
			  pbulk_tmp = guess_pbulk(pcos_bulk, c1_index, c2_index);
			}*/		

		      if (plv_tmp > 0)
			{
                          if (verbose) {cout<<"Trying to include a contribution for for residue "<<res_prot[i]-1<<endl;}	
			  // include the probability of corresponding config.
			  if (cylinder) { // In this mode, only include the water to the nearest residue
			    q_val = -log(plv_tmp);
                            rsq_min = r_sq; // This protein atom is the nearest, store it
                            min_ind = i; // Store the atom index
                          } else { // In this mode, add the water to all residues within the cut-off distance
			    q_val = -log(plv_tmp);
			    qlv_res[res_prot[i] - 1][0] += q_val;
			    qlv_res[res_prot[i] - 1][1] += q_val*q_val;
			    qlv_res[res_prot[i] - 1][2]++;
			  } // End output mode check
		          if (water_dump) {wdump_file<<t<<" "<<mol_id<<" "<<q_val<<endl;} // MN - Useful for coloring the surface and debugging
			} // End check for finite probabilities			
		    } // End distance check      		  
		} // End protein residue loop
		if (cylinder) { // If we store a water to only the nearest protein residue, do that here
		  if (min_ind <= N_prot) {
		    qlv_res[res_prot[min_ind] -1][0] += q_val;
		    qlv_res[res_prot[min_ind] -1][1] += q_val*q_val;
		    qlv_res[res_prot[min_ind] -1][2]++;
		  }
	        }
	    } // End water distance from surface check	    
	} // End water molecule loop
 
      // record on the log file
      logfile << "Calculation is currently  " << fixed << setprecision(3) << ((double)t - (double)firstTraj)/(double)Tstep*100 << " % done.\n";
      logfile.flush();
    	   
    }

  // average the quantity over the given length of time
  for (int i = 0; i < N_res; i++)
    {	      
      if (qlv_res[i][2] > 0)
	{
	  qlv_res[i][0] /= qlv_res[i][2];
	  qlv_res[i][1] /= qlv_res[i][2];
	  qlv_res[i][1] -= qlv_res[i][0]*qlv_res[i][0];
	}	     
	     
      // save the q-parameter
      data_qlv << i+1 << " " << qlv_res[i][0] << " " << sqrt(qlv_res[i][1]) << " " << qlv_res[i][2] << endl;

      delete[] qlv_res[i];
    }
  
  data_coord.close();
  data_qlv.close();
  wdump_file.close();

  logfile << "Complete! \n";
  logfile.close();

  for (int n = 0; n < N_prot; n++)
    {
      delete[] Protein_pos[n];
    }

  delete[] Protein_pos;
  delete[] res_prot;
  delete[] qlv_res;

  /* MN; With new initialization, delete doesn't work for memset array
  for (int a = 0; a < Nabin; a++)
    {
      for (int y1 = 0; y1 < Nc1bin; y1++)
	{
	  delete[] whist[a][y1];
	}
      delete[] whist[a];
    }
  */ 
  //delete[] whist;
  //delete[] pofa;
  
  /*
  for (int y1 = 0; y1 < cos_bin; y1++)
    {
      delete[] pcos_bulk[y1];
    }
  delete[] pcos_bulk;*/    

  // Report timing
  cout << "Completed in: " << double( clock() - startTime ) / CLOCKS_PER_SEC << " seconds." << endl; 
  
}

double pbc(double x, int k) {

  if (x > 0.5*L_box[k])
    x -= L_box[k];
  else if (x < -0.5*L_box[k])
    x += L_box[k];

  return x;
}

int int_round(double x) {

  int n = int(x);

  if (x - n >= 0.5)
    n++;

  return n;

}

double get_abs(double x) {
  
  if (x < 0)
    x *= -1.0;
  
  return x;
  
}

/*
double guess_plv(double *** p_cos, int a_index, int c1_index, int c2_index) {

  double pcos_val = 0;
  double count_ind = 0;

  if (c1_index - 1 >= 0)
    {
      if (p_cos[a_index][c1_index - 1][c2_index] > 0)
	{
	  pcos_val += p_cos[a_index][c1_index - 1][c2_index];
	  count_ind++;
	}      
    }
  if (c1_index + 1 < cos_bin)
    {
      if (p_cos[a_index][c1_index + 1][c2_index] > 0)
	{
	  pcos_val += p_cos[a_index][c1_index + 1][c2_index];
	  count_ind++;
	}      
    }
  if (c2_index - 1 >= 0)
    {
      if (p_cos[a_index][c1_index][c2_index - 1] > 0)
	{
	  pcos_val += p_cos[a_index][c1_index][c2_index - 1];
	  count_ind++;
	}      
    }
  if (c2_index + 1 < cos_bin)
    {
      if (p_cos[a_index][c1_index][c2_index + 1] > 0)
	{
	  pcos_val += p_cos[a_index][c1_index][c2_index + 1];
	  count_ind++;
	}      
    }  

  if (count_ind > 0)
    {
      pcos_val /= count_ind;
    }
  else
    {
      if (c1_index - 1 >= 0 && c2_index - 1 >= 0)
	{
	  if (p_cos[a_index][c1_index - 1][c2_index - 1] > 0)
	    {
	      pcos_val += p_cos[a_index][c1_index - 1][c2_index - 1];
	      count_ind++;
	    }      
	}
      if (c1_index + 1 < cos_bin && c2_index - 1 >= 0)
	{
	  if (p_cos[a_index][c1_index + 1][c2_index - 1] > 0)
	    {
	      pcos_val += p_cos[a_index][c1_index + 1][c2_index - 1];
	      count_ind++;
	    }      
	}
      if (c2_index - 1 >= 0 && c2_index + 1 < cos_bin)
	{
	  if (p_cos[a_index][c1_index - 1][c2_index + 1] > 0)
	    {
	      pcos_val += p_cos[a_index][c1_index - 1][c2_index + 1];
	      count_ind++;
	    }      
	}
      if (c1_index + 1 < cos_bin && c2_index + 1 < cos_bin)
	{
	  if (p_cos[a_index][c1_index + 1][c2_index + 1] > 0)
	    {
	      pcos_val += p_cos[a_index][c1_index + 1][c2_index + 1];
	      count_ind++;
	    }      
	}

      if (count_ind > 0)
	pcos_val /= count_ind;
    }
  return pcos_val;
}
  
double guess_pbulk(double ** p_cos, int c1_index, int c2_index) {

  double pcos_val = 0;
  double count_ind = 0;

  if (c1_index - 1 >= 0)
    {
      if (p_cos[c1_index - 1][c2_index] > 0)
	{
	  pcos_val += p_cos[c1_index - 1][c2_index];
	  count_ind++;
	}      
    }
  if (c1_index + 1 < cos_bin)
    {
      if (p_cos[c1_index + 1][c2_index] > 0)
	{
	  pcos_val += p_cos[c1_index + 1][c2_index];
	  count_ind++;
	}      
    }
  if (c2_index - 1 >= 0)
    {
      if (p_cos[c1_index][c2_index - 1] > 0)
	{
	  pcos_val += p_cos[c1_index][c2_index - 1];
	  count_ind++;
	}      
    }
  if (c2_index + 1 < cos_bin)
    {
      if (p_cos[c1_index][c2_index + 1] > 0)
	{
	  pcos_val += p_cos[c1_index][c2_index + 1];
	  count_ind++;
	}      
    }      

  if (count_ind > 0)
    {
      pcos_val /= count_ind;
    }
  else
    {
      if (c1_index - 1 >= 0 && c2_index - 1 >= 0)
	{
	  if (p_cos[c1_index - 1][c2_index - 1] > 0)
	    {
	      pcos_val += p_cos[c1_index - 1][c2_index - 1];
	      count_ind++;
	    }      
	}
      if (c1_index + 1 < cos_bin && c2_index - 1 >= 0)
	{
	  if (p_cos[c1_index + 1][c2_index - 1] > 0)
	    {
	      pcos_val += p_cos[c1_index + 1][c2_index - 1];
	      count_ind++;
	    }      
	}
      if (c2_index - 1 >= 0 && c2_index + 1 < cos_bin)
	{
	  if (p_cos[c1_index - 1][c2_index + 1] > 0)
	    {
	      pcos_val += p_cos[c1_index - 1][c2_index + 1];
	      count_ind++;
	    }      
	}
      if (c1_index + 1 < cos_bin && c2_index + 1 < cos_bin)
	{
	  if (p_cos[c1_index + 1][c2_index + 1] > 0)
	    {
	      pcos_val += p_cos[c1_index + 1][c2_index + 1];
	      count_ind++;
	    }      
	}

      if (count_ind > 0)
	pcos_val /= count_ind;
    }

  return pcos_val;
}
*/

/* calc_rsq(r1,r2); calculates the squared distance between two points in the simulation box,
 *  * taking periodic boundaries into consideration. */
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

/*
double calc_rsq (double * r_p, double * pos1) {

  double r_sq = 0;  
  
  for (int i = 0; i < 3; i++)
    {
      if (get_abs(r_p[i] - pos1[i]) < L_box[i]*0.5)
	{
	  r_sq += (pos1[i] - r_p[i])*(pos1[i] - r_p[i]);
	}
      else
	{
	  r_sq += pow(L_box[i] - get_abs(r_p[i] - pos1[i]), 2);
	}
    }
  
  return r_sq;
}
*/
