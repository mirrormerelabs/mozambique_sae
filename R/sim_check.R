############################################################
## sim_check.R
## Recovery check for the simulation (run after 03 and 04 with
## MOZSAE_SIM=1). Compares each model's district estimates and fixed
## effects with the known truth from sim_01_synthetic.R.
##
## What a working pipeline should show:
##   - M2 and M3 have lower RMSE than the direct estimates, and M3
##     lower than M1 (covariates and spatial structure help);
##   - 95% intervals cover the true share in roughly 90% to 98% of
##     districts (exact coverage is not expected from one draw);
##   - fixed-effect intervals of M2 and M3 contain the true values,
##     and the BYM2 mixing parameter phi is not far from 0.6;
##   - the frequentist Fay-Herriot (04) gives estimates close to M2.
##
## Districts are scored in two fixed groups, by the 01 rule (2+ clusters
## and 0 < p < 1): "in likelihood" and "out of sample". With MOZSAE_VAR =
## "smoothed1", single-cluster districts enter the likelihood but are
## still scored with the second group, so every variant is scored on the
## same districts.
############################################################

source("R/functions.R")
stopifnot("Set MOZSAE_SIM=1 (use run_sim.R)" = sim_mode)

truth <- readRDS(file.path(dir_work, "truth_district.rds"))
direct <- readRDS(file.path(dir_work, "direct_district.rds"))
est <- read_csv(file.path(dir_out, "s03_district_estimates.csv"), show_col_types = FALSE) %>%
  left_join(truth, by = "district_id")
truth_par <- read_csv(file.path(dir_out, "sim_truth_parameters.csv"), show_col_types = FALSE)

## ---- District estimates ---------------------------------------------
score <- function(df, label) df %>% summarise(
  estimator = label, districts = n(),
  rmse = sqrt(mean((p_hat - p_true)^2)), mean_abs_error = mean(abs(p_hat - p_true)),
  ## NA where predictions are constant (M1 out of sample: intercept only)
  correlation = if (length(p_hat) > 2 && sd(p_hat) > 0) cor(p_hat, p_true) else NA_real_,
  coverage_95 = if (all(c("lo", "hi") %in% names(df))) mean(p_true >= lo & p_true <= hi) else NA_real_,
  mean_interval_width = if (all(c("lo", "hi") %in% names(df))) mean(hi - lo) else NA_real_)

in_model_ids <- direct$district_id[direct$in_model]
## Also scored on their own (they belong to "out of sample"): 2+ clusters
## and a direct share of exactly 0 or 1 (MOZSAE_ZERO)
zero_ids <- direct$district_id[coalesce(direct$n_clusters >= 2 & (direct$p_direct == 0 | direct$p_direct == 1), FALSE)]
dir_df <- direct %>% filter(in_model) %>% left_join(truth, by = "district_id") %>%
  transmute(p_hat = p_direct, p_true, lo = p_direct - 1.96 * se_direct, hi = p_direct + 1.96 * se_direct)

recovery <- bind_rows(
  score(dir_df, "Direct (districts in likelihood)"),
  map_dfr(unique(est$model), function(mm) bind_rows(
    score(est %>% filter(model == mm, district_id %in% in_model_ids) %>%
            transmute(p_hat = p_mean, p_true, lo = p_q025, hi = p_q975), paste(mm, "(districts in likelihood)")),
    score(est %>% filter(model == mm, !district_id %in% in_model_ids) %>%
            transmute(p_hat = p_mean, p_true, lo = p_q025, hi = p_q975), paste(mm, "(out of sample)")),
    score(est %>% filter(model == mm) %>%
            transmute(p_hat = p_mean, p_true, lo = p_q025, hi = p_q975), paste(mm, "(all districts)")),
    if (length(zero_ids) > 0)
      score(est %>% filter(model == mm, district_id %in% zero_ids) %>%
              transmute(p_hat = p_mean, p_true, lo = p_q025, hi = p_q975), paste(mm, "(zero-share, 2+ clusters)")))))

fh_file <- file.path(dir_out, "s04_fh_emdi_estimates.csv")
if (file.exists(fh_file)) {
  fh <- read_csv(fh_file, show_col_types = FALSE) %>% left_join(truth, by = "district_id")
  ## 95% interval from the bootstrap MSE (normal approximation)
  fh <- transmute(fh, district_id, p_hat = p_fh, p_true, lo = p_fh - 1.96 * sqrt(mse_fh), hi = p_fh + 1.96 * sqrt(mse_fh))
  recovery <- bind_rows(recovery,
    score(filter(fh, district_id %in% in_model_ids), "Frequentist FH, emdi (districts in likelihood)"),
    score(filter(fh, !district_id %in% in_model_ids), "Frequentist FH, emdi (out of sample)"),
    score(fh, "Frequentist FH, emdi (all districts)"),
    if (length(zero_ids) > 0) score(filter(fh, district_id %in% zero_ids), "Frequentist FH, emdi (zero-share, 2+ clusters)"))
}

## ---- Parameters -----------------------------------------------------
## true_value: the coefficients of the data-generating process.
## realised_value: least squares of the realised theta_true on the
## covariates. With one draw of u, part of u is collinear with the
## (strongly correlated) covariates and is absorbed by the coefficients,
## so the models can only recover the realised values; a miss against
## true_value alone is therefore not evidence of bias.
md <- readRDS(file.path(dir_work, "model_data.rds")) %>% left_join(truth, by = "district_id")
realised <- coef(lm(theta_true ~ z_rwi + z_log_tt + z_log_dens, data = md))
fixed <- read_csv(file.path(dir_out, "s03_fixed_effects.csv"), show_col_types = FALSE) %>%
  filter(model %in% c("M2", "M3")) %>%
  left_join(rename(truth_par, term = parameter, true_value = value), by = "term") %>%
  mutate(realised_value = unname(realised[term]),
         interval_contains_truth = q025 <= true_value & q975 >= true_value,
         interval_contains_realised = q025 <= realised_value & q975 >= realised_value)
## phi in the simulation is the variance share of a smoothed-noise field,
## not of the scaled ICAR in INLA's BYM2, so only rough agreement is
## expected. Precision: 1 / sigma_total^2.
hyper <- read_csv(file.path(dir_out, "s03_hyperparameters.csv"), show_col_types = FALSE) %>%
  filter(model == "M3") %>%
  mutate(true_value = case_when(grepl("^Phi", parameter) ~ truth_par$value[truth_par$parameter == "phi"],
                                grepl("^Precision", parameter) ~ 1 / truth_par$value[truth_par$parameter == "sigma_total"]^2),
         interval_contains_truth = q025 <= true_value & q975 >= true_value)

rnd <- \(x) mutate(x, across(where(is.double), \(v) round(v, 4)))
recovery <- mutate(as_tibble(recovery), variance = var_method, zero_share = zero_method, .before = 1)
write_out(rnd(recovery), "sim_recovery_districts")
write_out(rnd(fixed), "sim_recovery_fixed")
write_out(rnd(hyper), "sim_recovery_hyper")
print(rnd(recovery), n = Inf, width = Inf)
message("sim_check.R finished.")
