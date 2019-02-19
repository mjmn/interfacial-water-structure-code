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

double a_range = 7.0;
double a_inc = 0.1;

int a_bin = int(a_range/a_inc);
double cos_inc = 0.02;
int cos_bin = 100;
int a_bin_rho = 140;

//declare functions
double pbc(double x, int k);
double get_abs(double x);
int int_round(double x);
double calc_rsq(double * r_p, double * pos1);
double guess_plv(double *** p_cos, int a_index, int c1_index, int c2_index);
double guess_pbulk(double ** p_cos, int c1_index, int c2_index);

void cl_parse(int argc, char * argv[])
{
  bool match;
  for(int i=1;i<argc;i+=2)
    {
      match = false;
      //if(strcmp(argv[i],"--N")==0){N = atoi(argv[i+1]);match=true;}
      if(strcmp(argv[i],"--T")==0){Tstep = atoi(argv[i+1]);match=true;}
      if(strcmp(argv[i],"--ts")==0){sample = atoi(argv[i+1]);match=true;}
      if(strcmp(argv[i],"--rcut")==0){R_cut = atof(argv[i+1]);match=true;}
      if(strcmp(argv[i],"--lx")==0){L_box[0] = atof(argv[i+1]);match=true;}
      if(strcmp(argv[i],"--ly")==0){L_box[1] = atof(argv[i+1]);match=true;}
      if(strcmp(argv[i],"--lz")==0){L_box[2] = atof(argv[i+1]);match=true;}
      // MNiesen; new input arguments, to allow time-splicing and more flexibility in the system definitions
      if(strcmp(argv[i],"--startT")==0){firstTraj = atof(argv[i+1]);match=true;} // Start time for analysis
      if(strcmp(argv[i],"--nprot")==0){N_prot = atof(argv[i+1]);match=true;} // Number of protein ATOMS - avg_slice.py
      if(strcmp(argv[i],"--nres")==0){N_res = atof(argv[i+1]);match=true;} // Number of protein RESIDUES - avg_slice.py
      if(match == false){cout<<"Warning: Bad input parameter: "<<argv[i]<<endl;}
    }

  cout<<"Tstep:\t"<<Tstep<<"\n"
      <<"sample_time:\t"<<sample<<"\n"
      <<"R_cut:\t"<<R_cut<<"\n"
      <<"Lx:\t"<<L_box[0]<<"\n"
      <<"Ly:\t"<<L_box[1]<<"\n"
      <<"Lz:\t"<<L_box[2]<<"\n";
  
}

int main(int argc, char * argv[]) { 
  // For timing
  clock_t startTime = clock();
  // parse for any new input parameters
  cl_parse(argc, argv);
  cout.flush();

  // initialize
  double r_p[3], pos_tmp[3], surf_norm[3], int_coord[3];
  double tmp1, tmp2, r_sq, rsq_min, a0, c1, c2, N0, rho_val, a_tmp, q_val, plv_tmp, pbulk_tmp;
  int N_hyd, count_surf, min_ind, atom_id, mol_id, a_index, c1_index, c2_index, a_ind_rho;

  double Rcut_sq = R_cut*R_cut;
  
  char file_prefix[50];

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

  data_res.close();
  
  for (int n = 0; n < N_res; n++)
    {
      qlv_res[n] = new double [3];

      for (int k = 0; k < 3; k++)
	{
	  qlv_res[n][k] = 0;
	}
    }
  
  // load the joint prob distibutions of bulk...
  ifstream p_bulk ("../tip3p/pbulk_w2_tip3p.txt"); // MNiesen; changed the source directory location, should not be hardcoded ***
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
  ifstream p_lv ("../tip3p/pjoint_sm2.txt"); // MNiesen; changed the source directory location, should not be hardcoded ***
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

  ifstream pa_lv ("../tip3p/pofa.txt"); // MNiesen; changed the source directory location, should not be hardcoded ***
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
  
  
  // files to be written..
  sprintf(file_prefix, "qlv_avg_res_%1.f.txt", R_cut);
  ofstream data_qlv(file_prefix);
  ofstream logfile("log_qlv.txt");

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
        if (!traj) { // Can't find any file, exit
          cerr << "Failed to read trajectory file: " << file_prefix << endl;
          exit(1);
        }
        fformat=1; // We're reading a XYZ file
      }
      else {
        fformat=0; // We're reading a PDB file
      }
     
      // MNiesen; hardcoded the number of commentlines here, appropriate for GROMACS (trjconv) generated pdb. Should be automated ***
      if (fformat==0) {
        for(int i = 0 ; i < 5 ; i++){
          getline(traj, str1);
        }
      }
      else {
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
	  }
          else { // XYZ format
            traj >> str1 >> r_p[0] >> r_p[1] >> r_p[2]; // Column format
          }

          // Store protein coordinates
	  for (int i = 0; i < 3; i++)
	    {
	      Protein_pos[n][i] = r_p[i];
	    }
	}

      traj.close();
      
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
	      
	  // Update the histograms for interfacial structure, q(x,y)
	  if (int_coord[0] < 3.0 && int_coord[0] >= -1.0)
	    {
	      rsq_min = 50.0;
	      min_ind = 0;
              //#pragma omp parallel for private(r_sq)
	      for (int i = 0; i < N_prot; i++)
		{
		  r_sq = calc_rsq(Protein_pos[i], pos_tmp);

		  if (r_sq <= Rcut_sq)
		    {			     
		      a_index = int((int_coord[0] + 1.0)*10.0);
		
		      a_ind_rho = int((int_coord[0] + 4.0)*10.0);
			  
		      c1_index = int((int_coord[1] + 1.0)*50.0);
		      c2_index = int((int_coord[2] + 1.0)*50.0);

		      if (c1_index==cos_bin)
			c1_index--;
		      if (c2_index==cos_bin)
			c2_index--;
		
		      plv_tmp = pcos_lv[a_index][c1_index][c2_index];
		      pbulk_tmp = pcos_bulk[c1_index][c2_index];
		
		      if (plv_tmp == 0)
			{
			  plv_tmp = guess_plv(pcos_lv, a_index, c1_index, c2_index);
			}
		      if (pbulk_tmp == 0)
			{
			  pbulk_tmp = guess_pbulk(pcos_bulk, c1_index, c2_index);
			}		

		      if (plv_tmp > 0 && pbulk_tmp > 0)
			{	
			  // include the probability of corresponding config.
			  q_val = -log(pofa[a_ind_rho]*plv_tmp/pbulk_tmp);
			  qlv_res[res_prot[i] - 1][0] += q_val;
			  qlv_res[res_prot[i] - 1][1] += q_val*q_val;
			  qlv_res[res_prot[i] - 1][2]++;
			}
			
		    }      		  
		}
	    }	    

	}
 
      // record on the log file
      logfile << "Calculation is currently  " << fixed << setprecision(3) << (double)t/(double)Tstep*100 << " % done.\n";
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

  logfile << "Complete! \n";
  logfile.close();

  for (int n = 0; n < N_prot; n++)
    {
      delete[] Protein_pos[n];
    }

  delete[] Protein_pos;
  delete[] res_prot;
  delete[] qlv_res;

  for (int a = 0; a < a_bin; a++)
    {
      for (int y1 = 0; y1 < cos_bin; y1++)
	{
	  delete[] pcos_lv[a][y1];
	}
      delete[] pcos_lv[a];
    }
     
  delete[] pcos_lv;
  delete[] pofa;
  
  for (int y1 = 0; y1 < cos_bin; y1++)
    {
      delete[] pcos_bulk[y1];
    }
  delete[] pcos_bulk;    

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
