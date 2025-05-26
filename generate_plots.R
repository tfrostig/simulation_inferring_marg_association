library(glue)
source("plotting/create_plots.R")

# Global constants
RESULTS_DIR <- "output/results2025-05-25"
PLOTS_DIR <- file.path(RESULTS_DIR, "plots")

# Ensure plots directory exists
ensure_directory(PLOTS_DIR)

# Define file paths
result_files <- list(
  psat = file.path(RESULTS_DIR, "FULL_RESULTS_PSAT.csv"),
  normal = file.path(RESULTS_DIR, "FULL_RESULTS_NORMAL.csv"),
  standard = file.path(RESULTS_DIR, "FULL_RESULTS_STANDARD.csv"),
  fs = file.path(RESULTS_DIR, "FULL_RESULTS_FS.csv")
)

# Generate all plots
generate_plots <- function() {
  createNormalPlots(results_file = result_files$normal, figure_dir = PLOTS_DIR)
  createForwardSelectionPlots(results_file = result_files$fs, figure_dir = PLOTS_DIR)
  createPSATPlots(results_file = result_files$psat, figure_dir = PLOTS_DIR)
  createStandardSimPlots(results_file = result_files$standard, figure_dir = PLOTS_DIR)

  message("All plots successfully generated in: ", PLOTS_DIR)
}

# Run the main function
generate_plots()
