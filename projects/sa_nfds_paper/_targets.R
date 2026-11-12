library(targets)

pkgload::load_all("nfdsmodelling", quiet=TRUE)
tar_option_set(
  packages = c("ggplot2", "dplyr", "gt"),
  imports  = "nfdsmodelling"
)

source("projects/sa_nfds_paper/config.R")
tar_source("projects/sa_nfds_paper/R")
tar_source("projects/sa_nfds_paper/targets")

c(
  targets_data,
  targets_fits,
  targets_fit_diagnostics,
  targets_figures,
  targets_sensitivity,
  targets_projected_ipd,
  targets_vaccine_design
)
