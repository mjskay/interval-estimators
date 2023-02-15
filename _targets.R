# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(dplyr)
library(rlang)
library(purrr)

# Set target options:
tar_option_set(
  packages = c("dplyr", "tidyr", "distributional", "rlang", "purrr", "sets"),
  format = "qs"
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multiprocess")

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# Replace the target list below with your own:
list(
  sims_targets,
  intervals_targets,
  eti_intervals_targets_combined,
  si_intervals_targets_combined,
  hdi_intervals_targets_combined
)
