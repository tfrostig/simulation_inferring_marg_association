# install_dependencies.R
# Script to install all required dependencies for genetic simulation and plotting

# Function to check and install CRAN packages
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    message("Installing missing packages: ", paste(new_packages, collapse = ", "))
    install.packages(new_packages, repos = "https://cloud.r-project.org")
  }
}

# Function to check and install GitHub packages
install_github_if_missing <- function(packages, github_urls) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    message("Installing GitHub packages: ", paste(new_packages, collapse = ", "))
    if(!require(remotes)) install.packages("remotes", repos = "https://cloud.r-project.org")
    for(i in seq_along(new_packages)) {
      remotes::install_github(github_urls[i])
    }
  }
}

# Required CRAN packages
cran_packages <- c(
  "MASS",        # For multivariate normal distribution
  "dplyr",       # Data manipulation
  "magrittr",    # Pipe operators
  "doParallel",  # Parallel processing
  "doRNG",       # Reproducible parallel random numbers
  "CovTools",    # Covariance matrix tools
  "tidyverse",   # Data manipulation and visualization
  "cowplot",     # Plot composition
  "latex2exp",   # LaTeX expressions in plots
  "gridExtra",   # Grid-based plotting
  "unglue",      # String parsing
  "ggsci",       # Scientific journal color palettes
  "ggplot2",     # Plotting
  "parallel",    # Parallel computation
  "foreach"      # Foreach looping construct
)


# NOTE: tmg must be installed for PSATInference, the package is not available on CRAN
# must be downloaded and installed as tar
# https://cran.r-project.org/src/contrib/Archive/tmg/


# Required GitHub packages
github_packages <- c("ECCCM", "PSATinference")
github_urls <- c("tfrostig/ECCCM", "tfrostig/PSAT")

# Create directory for outputs
dir.create("output", showWarnings = FALSE)

# Install required CRAN packages


message("Checking and installing required CRAN packages...")
install_if_missing(cran_packages)

# Install required GitHub packages
message("Checking and installing required GitHub packages...")
install_github_if_missing(github_packages, github_urls)

message("All dependencies installed successfully!")
message("Run 'main.R' to execute the simulations and generate plots.")
