library(targets)

pkgload::load_all("nfdsmodelling", quiet=TRUE)
tar_option_set(
  package = c("ggplot2"),
  imports = "nfdsmodelling"
)

source("projects/sa_nfds_paper/config.R")
tar_source("projects/sa_nfds_paper/R")
tar_source("projects/sa_nfds_paper/targets")

c(
  targets_vaccines
)
