# Simulation 
replication code for the simulation in the manuscript: "Inferring on joint associations from marginal associations and a reference sample" 

# Installation 
run the script `installation.R`. Notice, that the `PSAT` package requires `tmg`, which requires a manual installation, install the `tmg` from 
 - https://cran.r-project.org/src/contrib/Archive/tmg/

# R Version 
Tested R version 4.0.4 

 # Run Code 
 - Open Project 
 - Open the `run_simulations.R` script. Update the parameters:  There are three parameters: `n_iter` - number of iterations used in the simulation (default 50),
 `n_cores` - number of cores used in parallel (default 10) and `parallel` - Boolean flag whether to use parallelization at all. 
 - Run `run_simulations.R`
 - The results will be saved in the `RESULTS_DIR` directory, which is defined in `run_simulations.
 - Update the relevant `RESULTS_DIR` directory in `generate_plots.R`
   

# Notes 
The simulation is extremely slow when attempting to replicate exactly the results in the manuscript (due to the computational complexity and the number of scenarios) and might take several days, to replicate the results exactly increase the number of iteration to 1000. 
