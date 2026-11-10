FROM rocker/r-ver:4.5.3

# System libraries needed by R package binaries (igraph, used by targets)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libglpk40 libxml2 \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled \
    targets pkgload here \
    ggplot2 dplyr tidyr stringr readr patchwork \
    deSolve

WORKDIR /repo
ENV TAR_PROJECT=sa_nfds_paper
CMD ["Rscript", "-e", "targets::tar_make()"]
