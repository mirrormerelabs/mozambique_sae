############################################################
## sim_variance_study.R
## Compares variants of the sampling model in 03 and 04 on the same
## synthetic data, over several simulation seeds. SIMULATED DATA ONLY:
## no DHS data are read. Two studies (see functions.R and 03, section 1):
##   "variance" MOZSAE_VAR: "direct", "smoothed", "smoothed1"
##   "zero"     MOZSAE_ZERO with smoothed variances: "smoothed" (zero-share
##              districts excluded), "smoothed_zero" (included with the
##              empirical logit)
##   "zero_plasmode" the same, in the plasmode simulation scenario
## Every variant sets both switches, so defaults do not affect a study.
##
## For each seed and variant, runs sim_01, 03, 04 and sim_check with
## MOZSAE_SIM=1 in sim/runs/seed_<k>/<variant>/ (sim_01 is rerun in each
## folder with the same seed, so the variants see identical data). The
## first seed is the one used by run_sim.R.
##
## Run from the project root after 02_covariates.R:
##   Rscript R/sim_variance_study.R [seeds, default 20] [parallel jobs, default 4] [study, default variance]
## Writes to sim/output/ (<study> = variance or zero):
##   sim_<study>_comparison.csv          mean, minimum and maximum over seeds
##   sim_<study>_comparison_by_seed.csv  one row per seed, variant and model
## DIC, WAIC and the log score compare models within a variant; they are
## not comparable across variants, whose likelihoods use different
## sampling variances. RMSE, coverage and interval width are, since they
## are scored against the truth.
############################################################

Sys.setenv(MOZSAE_SIM = "1")
source("R/functions.R")
suppressPackageStartupMessages(library(parallel))

args <- commandArgs(trailingOnly = TRUE)
n_seeds <- if (length(args) >= 1) as.integer(args[1]) else 20L
n_jobs  <- if (length(args) >= 2) as.integer(args[2]) else 4L
study   <- if (length(args) >= 3) args[3] else "variance"
studies <- list(
  variance = list(direct = c(MOZSAE_VAR = "direct", MOZSAE_ZERO = "exclude"),
                  smoothed = c(MOZSAE_VAR = "smoothed", MOZSAE_ZERO = "exclude"),
                  smoothed1 = c(MOZSAE_VAR = "smoothed1", MOZSAE_ZERO = "exclude")),
  zero = list(smoothed = c(MOZSAE_VAR = "smoothed", MOZSAE_ZERO = "exclude"),
              smoothed_zero = c(MOZSAE_VAR = "smoothed", MOZSAE_ZERO = "include")),
  ## as "zero", in the plasmode scenario (sim_01_synthetic.R)
  zero_plasmode = list(plasmode_smoothed = c(MOZSAE_VAR = "smoothed", MOZSAE_ZERO = "exclude", MOZSAE_SIM_SCENARIO = "plasmode"),
                       plasmode_smoothed_zero = c(MOZSAE_VAR = "smoothed", MOZSAE_ZERO = "include", MOZSAE_SIM_SCENARIO = "plasmode")))
## the base scenario unless a variant says otherwise
studies <- map(studies, \(st) map(st, \(v) if ("MOZSAE_SIM_SCENARIO" %in% names(v)) v else c(v, MOZSAE_SIM_SCENARIO = "base")))
stopifnot("study must be variance, zero or zero_plasmode" = study %in% names(studies))
## Plasmode truth: snapshot of the public M3 estimates of the first real
## run, made once (later real runs do not change it)
if (study == "zero_plasmode" && !file.exists("sim/plasmode_truth.rds")) {
  pm <- read_csv("output/s03_district_estimates.csv", show_col_types = FALSE) %>%
    filter(model == "M3") %>% select(district_id, p_true = p_mean)
  stopifnot(nrow(pm) == n_districts)
  saveRDS(pm, "sim/plasmode_truth.rds")
}
variants <- names(studies[[study]])
scripts <- c("R/sim_01_synthetic.R", "R/03_models_inla.R", "R/04_fh_emdi.R", "R/sim_check.R")
jobs <- expand.grid(k = seq_len(n_seeds), variance = variants, stringsAsFactors = FALSE)
jobs$seed <- 20261116L + jobs$k - 1L
jobs$run <- sprintf("runs/seed_%02d/%s", jobs$k, jobs$variance)

## ---- 1. Runs -----------------------------------------------------------
run_job <- function(j) {
  logs <- file.path("sim", j$run, "output", "logs")
  dir.create(logs, showWarnings = FALSE, recursive = TRUE)
  sw <- studies[[study]][[j$variance]]
  env <- c("MOZSAE_SIM=1", paste0("MOZSAE_SIM_RUN=", j$run), paste0("MOZSAE_SIM_SEED=", j$seed),
           paste0(names(sw), "=", sw))
  for (s in scripts) {
    log <- file.path(logs, sub("\\.R$", ".log", basename(s)))
    if (system2("Rscript", s, stdout = log, stderr = log, env = env) != 0) return(paste(j$run, s, "failed"))
  }
  "ok"
}
t0 <- Sys.time()
status <- unlist(mclapply(split(jobs, seq_len(nrow(jobs))), run_job, mc.cores = n_jobs, mc.preschedule = FALSE))
message("Runs finished in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " minutes")
if (any(status != "ok")) stop(paste(status[status != "ok"], collapse = "\n"))

## ---- 2. Collect --------------------------------------------------------
read_run <- function(j) {
  out <- file.path("sim", j$run, "output")
  rec <- read_csv(file.path(out, "sim_recovery_districts.csv"), show_col_types = FALSE) %>%
    mutate(model = case_when(grepl("^Direct", estimator) ~ "Direct", grepl("emdi", estimator) ~ "FH emdi",
                             TRUE ~ sub(" .*", "", estimator)),
           group = case_when(grepl("in likelihood", estimator) ~ "in", grepl("out of sample", estimator) ~ "out",
                             grepl("zero-share", estimator) ~ "zero", TRUE ~ "all")) %>%
    select(model, group, districts, rmse, coverage_95, mean_interval_width) %>%
    pivot_wider(names_from = group, values_from = c(districts, rmse, coverage_95, mean_interval_width))
  cmp <- read_csv(file.path(out, "s03_model_comparison.csv"), show_col_types = FALSE) %>%
    select(model, dic, waic, log_score, pit_share_outside_05_95, cpo_failures)
  fx <- read_csv(file.path(out, "sim_recovery_fixed.csv"), show_col_types = FALSE) %>%
    group_by(model) %>%
    summarise(fixed_contain_truth = sum(interval_contains_truth),
              fixed_contain_realised = sum(interval_contains_realised), .groups = "drop")
  de <- read_csv(file.path(out, "s03_design_effects.csv"), show_col_types = FALSE)
  rec %>% left_join(cmp, by = "model") %>% left_join(fx, by = "model") %>%
    mutate(variance = j$variance, seed = j$seed,
           deff_pooled = as.numeric(de$value[grepl("^Pooled", de$item)]), .before = 1)
}
by_seed <- map_dfr(split(jobs, seq_len(nrow(jobs))), read_run) %>%
  mutate(model = factor(model, levels = c("Direct", "M1", "M2", "M3", "M4", "FH emdi")),
         variance = factor(variance, levels = variants)) %>%
  arrange(variance, model, seed)

summary_tab <- by_seed %>%
  select(-seed) %>%
  pivot_longer(-c(variance, model), names_to = "metric") %>%
  filter(!is.na(value)) %>%
  group_by(variance, model, metric) %>%
  summarise(seeds = n(), mean = mean(value), min = min(value), max = max(value), .groups = "drop") %>%
  arrange(metric, model, variance)

rnd <- \(x) mutate(x, across(where(is.double), \(v) round(v, 4)))
write_out(rnd(by_seed), paste0("sim_", study, "_comparison_by_seed"))
write_out(rnd(summary_tab), paste0("sim_", study, "_comparison"))
message("sim_variance_study.R finished (SIMULATED DATA).")
