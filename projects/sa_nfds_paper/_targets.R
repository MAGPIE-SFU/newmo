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
  targets_vaccines,
  targets_data,
  targets_fits
)
