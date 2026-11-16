############################################################
## run_all.R: run the pipeline from a clean session, from the
## project root. 00_setup.R is run once, interactively, beforehand.
## Each script runs in its own R process; console output goes to
## output/logs/.
############################################################
scripts <- c("R/01_direct.R", "R/02_covariates.R", "R/03_models_inla.R",
             "R/04_fh_emdi.R", "R/05_outputs.R", "R/06_web.R")
dir.create("output/logs", showWarnings = FALSE, recursive = TRUE)
for (s in scripts) {
  log <- file.path("output/logs", sub("\\.R$", ".log", basename(s)))
  message("Running ", s, " ...")
  status <- system2("Rscript", s, stdout = log, stderr = log)
  if (status != 0) stop(s, " failed; see ", log)
}
message("All scripts finished.")
