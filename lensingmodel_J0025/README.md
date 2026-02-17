# Scripts to compute the maximum possible lensing magnification with only one image

Inside each fsubdirectory (containing the script for each filter):

- compute_mu_thetaE.py: the main function. You need to have an executable glafic file to run the script. The executable is not included here but can be downloaded at https://github.com/oguri/glafic2/tree/main?tab=readme-ov-file
- compute_positions.py: a script for computing the position of the point source relative to the lensing galaxy.
- muthetaE.txt: contains the derived mu - thetaE relation.
- plot_mu_thetaE.py: a script to plot the dependency of magnification (mu) on theta_E. 
- point.dat: contains the position of the lensed image (from compute_positions.py). This file will be read by glafic.

All other files are intermediate products (output files of glafic runs).
