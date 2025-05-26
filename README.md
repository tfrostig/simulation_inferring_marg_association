# Simulation 
replication code for the simulation in the manuscript: "Inferring on joint associations from marginal associations and a reference sample" 

# Installation 
run the script `installation.R`. Notice, that the `PSAT` package requires `tmg`, which requires a manual installation, install the `tmg` from 
 - https://cran.r-project.org/src/contrib/Archive/tmg/

# R Version 
Tested R version 4.0.4 

 # Run Code 
 - Open Project 
 - Update parameters (number of iterations and cores) and run `run_simulations.R`
 - Update the relevant `RESULTS_DIR` directory in `generate_plots.R`
   

# Notes 
The simulation is extremely slow when attempting to replicate exactly the results in the manuscript (due to the computational complexity and the number of scenarios) and might take several days, to replicate the results exactly increase the number of iteration to 1000. 
