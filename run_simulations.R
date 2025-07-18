source('simulation/forward_selection_sim.R')
source('simulation/genetic_no_selection_sim.R')
source('simulation/normal_dist_sim.R')
source('simulation/genetic_w_selection_sim.R')

## Parameters
n_iter = 50
n_cores = 10
parallel = TRUE
output_dir = paste0("output/results", Sys.Date())

# create directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Set seed for reproducibility
set.seed(123)

# Run simulations
runGeneticPSATSim(output_dir = output_dir, n_iter = n_iter, parallel = parallel, n_cores = n_cores)
runGeneticStandardSim(output_dir = output_dir, n_iter = n_iter, parallel = parallel, n_cores = n_cores)
runNormalSim(output_dir = output_dir, n_iter = n_iter, parallel = parallel, n_cores = n_cores)
runForwardSelectionSim(output_dir = output_dir, n_iter = n_iter, parallel = parallel, n_cores = n_cores)
