# Making initial inputs

Initial inputs to solution phase lammps files involves putting some number of molecules you are interested in and then put them in a solvent.

- Initial goal is to be able to feed in a base xyz file in to this project, turn it in to a molecule and specify how many of them I would like in the box
- Then I would like to solvate it in tip3p water.

This will be acheived by using the typical equilibration work that is implemented then bringing the simulation down to near 0K.
The .cpt file will be useful as it is a restart file.
1. use a .ctp file as an input file
2. ramp the temperature down using NVT ensemble
3. output .trr file with positions and velocities


# Running Initial Equilibration Runs

- NVP ensemble
- NVT ensemble
- Check system equilibration
- Bring tmeperature down to near 0K

-  final output files need coords and velocities
    trr files have coordinates tpr files have both, but cannot be read by MDAnalysis 

# Switch to ReaxFF forcefield

I want this to be as easy as specifying a file to look at.

# Data Analysis

I will likely use ben_sci_tools stuff that I have already made for species out data. 

I would also like to zip data for easy transfer.
