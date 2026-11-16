############################################################
## run_sim.R: test the pipeline on SIMULATED data (no DHS data read).
## Run from the project root after 02_covariates.R:
##   Rscript run_sim.R
## Everything is written to sim/ (excluded from version control).
## Sampling variances follow MOZSAE_VAR (default "smoothed"; see
## R/functions.R). To compare the variants over 20 seeds:
##   Rscript R/sim_variance_study.R
############################################################
scripts <- c("R/sim_01_synthetic.R", "R/03_models_inla.R", "R/04_fh_emdi.R",
             "R/05_outputs.R", "R/06_web.R", "R/sim_check.R")
dir.create("sim/output/logs", showWarnings = FALSE, recursive = TRUE)
Sys.setenv(MOZSAE_SIM = "1")
for (s in scripts) {
  log <- file.path("sim/output/logs", sub("\\.R$", ".log", basename(s)))
  message("Running ", s, " (simulation) ...")
  status <- system2("Rscript", s, stdout = log, stderr = log, env = "MOZSAE_SIM=1")
  if (status != 0) stop(s, " failed; see ", log)
}
message("Simulation finished. Recovery check: sim/output/sim_recovery_districts.csv")
