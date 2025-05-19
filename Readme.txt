This is a repository to reproduce the calculations and plots for the Paper: 
"How Thermostats Influence Dynamics Across Time Scales: A Systematic Study from Fast Motions to Slow Transitions"

The DOI and Archive Link will follow once they are done.

The Code is split in three parts Gromacs, AIMD, MSM
The AIMD part is currently missing since the contributer is on parental leave

In Gromacs, are the codes for calculating VACF and PACF from gromacs simulations. 
This Part is written in Julia
Please read the paper for simulastion details. The Calculate codes can be copied into the simulation folder and executed there.
Each Calculate Code saves the results in .JLD2 files. 

From these JLD2 files the plots can be created with the Plot_Figure_X files.
For the VACF the TRR is necessary, for the PACF a .xvg file obtained from the edr file containing the Pressure Tensor
For the Energy distribution a .xvg file with the total energy (also obtained from the edr) is needed.
The Plot_Figure_X files contain more code, that create additional plots used in the SI or where no longer used.

The MSM code is compatible with Gromacs simulations and the folder contains a more detailed readme on how to Calculate the MSMs for Pentane.

