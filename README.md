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


