Code to analyze the water structure at a surface (i.e., protein) during an MD simulation. Currently the code can be used as a post-processing analysis tool, provided one or multiple trajectory frames in PDB or XYZ format. This file contains instructions for compiling and some of the basic usage of the code.

###  COMPILATION ###
At the moment, the code is compiled by a simple shell script that is included. Simply type ./compile.sh in the main directory and the code should compile. If you want to customize the compilation, the compile.sh script can be edited. Succesful complation should yield three executable:
IFP - Given a protein trajectory; calculate the location of the protein/water dividing surface during the trajectory, and a trajectory of the water molecules near the surface in terms of the collective variables described in S. Shin et al., JCTC 2018 (distance from the surface and angle between the OH bonds and the surface normal). These outputs are needed to use the other two executables.
OWS - Used to generate new reference water structure files. This allows to generate arbitrary new reference files, such as the water structure at a specific part of a protein surface.
QPR - Given the output from IFP, and a reference water structure; calculate the log-likelihood over trajectory time that a site on the protein surface has the same water structure as in the given reference.

Further description of the executables, and instructions on how to use them (with an example workflow) are provided below.

### TESTS ###
To ensure that the code has compiled correctly and to test that any modifications made to the code have not broken some functionality, a test folder is provided. To perform the tests, go to the ./Test/ directory and execute the included shell script: ./run_test.sh 

The resulting output will either state that the test outputs are identical to what is expected (the code works as intended), or list the differences. The test makes use of an included short trajectory, and uses a water reference structure for TIP3P water at an ideal planar hydrophobic surface.

### USING THE CODE ###
This section describes a typical workflow for the analysis of a MD trajectory with a protein solvated in water. Not all steps are required and optional steps are indicated, the workflow should be used as a general example and can be adjusted depending on what the user is trying to do. The workflow uses all three of the main executables, and some of the utility scripts that are included in the ./utils/ folder. A more complete description of what each script does is provided in the UTILITIES and EXECUTABLES sections.

Step 1) Calculate the location of the Willard-Chandler water interface throughout the MD trajectory and output the structure of waters near the interface in the reduced coordinates described in S. Shin et al., JCTC 2018:

./IFP --lx [box size X] --ly [box size Y] --lz [box size Z] --startT [first frame to analyze] --T [number of frames to analyze] --nwater [number of water molecules] --nprot [number of protein atoms] --nres [number of protein residues]

This requires the MD trajectory to be in the default location (../traj), with each frame saved as a .pdb or .xyz file with filename [frame].pdb or [frame].xyz. The --T and --startT options enable time-slicing of long MD trajectories, to paralellize the analysis. This step will generate the following output files:
- traj_iface.xyz : A trajectory of the location of the calculated W-C interface. Due to the changing size of the interface, this is an xyz trajectory with varying number of particles per frame. It can be displayed in VMD using "readvarxyz", further information on this VMD utility can be found on the VMD website.
- coord_iface.txt : Information on the water molecules that are near the W-C interface, including their structure in terms of the distance to the interface, and the angle between the OH bonds and a vector perpendicular to the interface.
- log_iface.txt : Information dump on the calculation, can be expanded by using the '--verbose' flag when calling ./IFP.

Step 2) [OPTIONAL] Generate a reference water-structure distribution. The purpose of our code is to compare the structure of water molecules at a protein surface to a well-defined reference. In some cases, it may be useful to make such a reference distribution from a protein in water trajectory (i.e., if we want to create a reference that describes a functional site on a protein such as a ligand binding site or PPI site). To do this we use the '/utils/gen_includeRes.py' script and the 'OWS' executable:

./utils/gen_includeRes.py [pdbfilename with the simulated protein] [user selection of residues and chains]

This step defines which part of the protein makes up the reference surface, in terms of specific amino acid residues. For example, if we want to generate a reference water-distribution for a PPI that consists of residues 10 11 and 12 in chain A of the protein contained in protein.pdb, we would execute the following command: "./utils/gen_includeRes.py protein.pdb 10-A 11-A 12-A"
The output of this utility script is a single file "resInclude.txt" that contains an entry for each atom in the protein; the entry is "1" for protein atoms that make up the reference surface, and "0" otherwise. This file should be places in the MD trajectory file for the next step:

./OWS --prefix [path to the MD trajectory] --lx [box size X] --ly [box size Y] --lz [box size Z] --startT [first frame to analyze] --T [number of frames to analyze] --nprot [number of protein atoms] --rcut [cut-off distance beyond which waters are not considered; default 6 Angstrom] --hist --cylinder

The reference water structure is calculated only using waters that are near a set of user-defined residues, this enables us to generate water structures for specific functional sites on the protein surface. This requires MD frames saved as a .pdb or .xyz file with filename [frame].pdb or [frame].xyz. The --T and --startT options enable time-slicing of long MD trajectories, to paralellize the calculation, outputs can be combined afterwards using the ./utils/aggregate_local_structure_histograms.py script. The '--hist' flag tells the script to output a histogram of the water structure at the reference site, and the 'cylinder' flag tells the script that each water molecule should only be assigned to the nearest protein residue (alternate options described in the EXECUTABLES section). This step will generate the following output files:
- local_wstructure_histogram.txt : contains the histogram of water structure at the reference surface. The user can control the number of bins in each of the three collective coordinates using the --nabin [number of bins in the distance to the interface] --nc1bin [number of bins in the OH1/surface normal angle] --nc2bin [number of bins in the OH2/surface normal angle].
The executable can generate other outputs that, along with additional input options, are described in the EXECUTABLES section.

Step 3) In the final step, we track the log-likelihood that the water-structure near each protein residue matches that of a user-provided reference water-structure distribution. The user can generate their own reference water-structure distribution using the protocol described in Step 2, or use a default reference water-structure distribution (i.e., the one provided in the ./tip3p folder describes water at a planar hydrophobic surface). The default location for the reference water structure files is "../ref_structure/", in the most recent version of the code this folder should contain a single file "pjoint_sm2.txt" that contains 5 columns; [distance to the interface] [OH1/surface normal angle] [OH2/surface normal angle] [probability (can be normalized histogram count] [raw histogram count]. With these files available, we can use the "QPR" executable to perform the analysis:

./QPR --lx [box size X] --ly [box size Y] --lz [box size Z] --starT [first frame to analyze] --T [number of frames to analyze] --nprot [number of protein atoms] --nres [number of protein residues]

NOTE: if using a custom water-structure reference file, add "--nabin --nc1bin and nc2bin" flags to match what was used in Step 2 to generate the reference file. Also make sure that your reference file has the correct number of columns (5), if you just want to use histogram count as probability, the 4th and 5th columns can be identical.

This will output the log-likelihood for each protein residue that the water structure near that residue, R, in timewindow, startT to startT+T, matches the provided reference water-structure. The user can define which protein atom belong to which residue (allowing additional flexibility) using the "res_list.txt" file, this file contains two columns; atom IDs + residue IDs. It can be automatically generated using the script "./utils/gen_res_list.py protein.pdb [number of residues in protein.pdb]", in which case it just uses standard amino acid definitions. Running the "QPR" executable will generate the following output files:
- qlv_avg_res_6.txt : Contains the output of the analysis, for each user-defined residue this has; [resID] [log-likelihood] [standard deviation of the log-likelihood] [number of waters analyzed to do the analysis].
- log_qlv.txt : Information dump on the analysis, can be expanded using the '--verbose' flag when calling ./QPR.

In case the user wants to track the log-likelihood over time, use time slicing with the '--startT' and '--T' flags. Resulting outputs can be combined into a matrix with time on one axis and residueID on the other axis using the script "./utils/collect_data.sh". This will generate separate files for the log-likelihood and the number of waters used to calculate those values.
