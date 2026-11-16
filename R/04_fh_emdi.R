############################################################
## 04_fh_emdi.R
## Frequentist cross-check: standard Fay-Herriot model (REML) with the
## arcsine transformation for proportions and bias-corrected back-
## transformation, MSE by parametric bootstrap (emdi; Kreutzmann et
## al., 2019). Same covariates, districts and sampling variances as
## M2/M3 (see 03, section 1). Districts not in the model
## data are out-of-sample (synthetic prediction).
##
## Purpose: show that the Bayesian M2 reproduces the established
## frequentist estimator, and give MSE-based CVs in the form used in
## official SAE poverty maps (e.g. Corral et al., 2022).
############################################################

source("R/functions.R")
suppressPackageStartupMessages(library(emdi))
set.seed(20261116)

d <- readRDS(file.path(dir_work, "model_data.rds"))
est <- read_csv(file.path(dir_out, "s03_district_estimates.csv"), show_col_types = FALSE)

fh_data <- d %>%
  transmute(district_id, z_rwi, z_log_tt, z_log_dens,
            p_direct = ifelse(in_model, p_lik, NA_real_),   # empirical share for 0/1 districts (MOZSAE_ZERO)
            ## sampling variances as in 03 (direct or smoothed; MOZSAE_VAR);
            ## smoothed: n_eff = n_households / pooled design effect
            var_direct = ifelse(in_model, var_p, NA_real_),
            n_eff = ifelse(in_model, n_eff_fh, NA_real_)) %>%
  as.data.frame()

fh_fit <- fh(fixed = p_direct ~ z_rwi + z_log_tt + z_log_dens,
             vardir = "var_direct", combined_data = fh_data, domains = "district_id",
             method = "reml", transformation = "arcsin", backtransformation = "bc",
             eff_smpsize = "n_eff", MSE = TRUE, mse_type = "boot", B = 500)

fh_est <- estimators(fh_fit, indicator = "FH", MSE = TRUE, CV = TRUE)$ind %>%
  as_tibble() %>%
  transmute(district_id = as.integer(as.character(Domain)), p_fh = FH, mse_fh = FH_MSE, cv_fh = FH_CV)

compare <- est %>% filter(model %in% c("M2", "M3")) %>%
  select(model, district_id, p_mean, cv) %>%
  left_join(select(d, district_id, in_model), by = "district_id") %>%   # not in the published CSV
  inner_join(fh_est, by = "district_id") %>%
  group_by(model, in_model) %>%
  summarise(districts = n(), pearson_r = cor(p_mean, p_fh), spearman_rho = cor(p_mean, p_fh, method = "spearman"),
            mean_abs_diff = mean(abs(p_mean - p_fh)), median_cv_bayes = median(cv), median_cv_fh = median(cv_fh),
            .groups = "drop")

model_info <- tibble(item = c("REML variance of the random effect (arcsine scale)", "Adjusted R2 (Lahiri and Suntornchost)"),
                     value = c(fh_fit$model$variance, fh_fit$model$model_select$AdjR2))

rnd <- \(x) mutate(x, across(where(is.double), \(v) round(v, 4)))
write_out(rnd(fh_est), "s04_fh_emdi_estimates")
write_out(rnd(compare), "s04_fh_vs_inla")
write_out(rnd(model_info), "s04_fh_model_info")
message("04_fh_emdi.R finished.")
