############################################################
## sens_real.R
## Sensitivity run on the real data: refits 03 to 05 with one or more
## switches away from their defaults (see functions.R), for example
##   Rscript R/sens_real.R MOZSAE_ZERO=include
## Run from the project root after run_all.R. Every output goes to
## data_work/sensitivity/<MOZSAE_VAR>_<MOZSAE_ZERO>/ (private, never
## published): output/ and docs/ are not touched, and 06 is not run.
############################################################

args <- commandArgs(trailingOnly = TRUE)
stopifnot("Give switches as NAME=value" = length(args) >= 1, all(grepl("^MOZSAE_(VAR|ZERO)=", args)))
stopifnot("Not in simulation mode" = !identical(Sys.getenv("MOZSAE_SIM"), "1"))
sw <- sub("^[^=]*=", "", args); names(sw) <- sub("=.*$", "", args)
do.call(Sys.setenv, as.list(sw))
source("R/functions.R")
stopifnot("Every switch is at its default: use run_all.R" = sens_mode)

## ---- Inputs from 01 and 02 ---------------------------------------------
inputs <- c("direct_district.rds", "direct_province.rds", "covariates_district.rds")
stopifnot("Run run_all.R first" = all(file.exists(file.path("data_work", inputs))))
file.copy(file.path("data_work", inputs), dir_work, overwrite = TRUE)

## ---- Models and outputs ----------------------------------------------
logs <- file.path(dir_out, "logs"); dir.create(logs, showWarnings = FALSE)
for (s in c("R/03_models_inla.R", "R/04_fh_emdi.R", "R/05_outputs.R")) {
  log <- file.path(logs, sub("\\.R$", ".log", basename(s)))
  message("Running ", s, " ...")
  if (system2("Rscript", s, stdout = log, stderr = log, env = args) != 0) stop(s, " failed; see ", log)
}
message("sens_real.R finished: ", dir_root)
