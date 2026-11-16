############################################################
## 00_setup.R
## Run once, interactively, from the project root. Creates the renv
## library and installs the packages. INLA is not on CRAN; it comes
## from the INLA repository and is large (a few hundred MB).
############################################################

if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::init(bare = TRUE)   # new, empty project library

cran_pkgs <- c("dplyr", "tidyr", "readr", "purrr", "haven", "sf", "terra",
               "survey", "spdep", "Matrix", "ggplot2", "patchwork",
               "emdi",    # emdi: frequentist Fay-Herriot cross-check (04)
               "sn",      # sn: INLA suggests it, but inla.posterior.sample() needs it (03)
               "leaflet", "htmlwidgets", "htmltools")   # interactive map (06)
renv::install(cran_pkgs)

## INLA (stable). On Apple silicon the binary is installed by default;
## if it fails, see https://www.r-inla.org/download-install
## Required dependencies only: INLA's suggested packages include
## Bioconductor packages (graph, Rgraphviz) that are not needed here.
options(repos = c(INLA = "https://inla.r-inla-download.org/R/stable",
                  CRAN = "https://cloud.r-project.org"))
install.packages("INLA", dependencies = c("Depends", "Imports", "LinkingTo"))
renv::snapshot(prompt = FALSE)

## Quick check
library(INLA)
message("Setup finished: INLA ", inla.version("version"))
