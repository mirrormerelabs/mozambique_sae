############################################################
## 03_models_inla.R
## Area-level (Fay-Herriot) models on the logit scale, fitted in R-INLA.
##
##   y_i = logit(p_direct_i),  V_i = var_i / (p_i (1 - p_i))^2  (delta method)
##   y_i | theta_i ~ N(theta_i, V_i)        (sampling model, V_i known)
## var_i is the direct sampling variance or, by default, a smoothed one
## from a design-effect generalised variance function (section 1;
## MOZSAE_VAR in functions.R).
##   theta_i = beta0 + x_i' beta + u_i      (linking model)
##
## Four linking models:
##   M1  intercept + IID district effect (no covariates)
##   M2  covariates + IID effect          (Bayesian Fay-Herriot)
##   M3  covariates + BYM2 effect         (spatial Fay-Herriot; main model)
##   M4  intercept + BYM2 effect          (spatial smoothing only)
## Districts not in the likelihood (fewer than 2 clusters, or a direct
## estimate of 0 or 1; with MOZSAE_VAR = "smoothed1", fewer than 1)
## enter with y = NA and are predicted.
##
## The Gaussian likelihood has its log precision fixed at 0 and is
## scaled by 1/V_i, so each observation has precision exactly 1/V_i.
## Priors (PC): BYM2 total SD, P(sigma > 1) = 0.01; mixing parameter
## phi, P(phi < 0.5) = 2/3 (Riebler et al., 2016); same SD prior for IID.
##
## Model checks: DIC, WAIC, log score from CPO (leave-one-district-out
## predictive density of the direct estimate, which includes sampling
## variance), PIT values. Outputs are district-level only.
############################################################

source("R/functions.R")
suppressPackageStartupMessages({ library(INLA); library(spdep); library(Matrix) })
## sn: only suggested by INLA, but inla.posterior.sample() (section 5) fails without it
stopifnot("Package sn is needed: see R/00_setup.R" = requireNamespace("sn", quietly = TRUE))
set.seed(20261116)

adm2 <- read_districts()
direct <- readRDS(file.path(dir_work, "direct_district.rds"))
covs <- readRDS(file.path(dir_work, "covariates_district.rds"))

d <- direct %>%
  left_join(select(covs, district_id, pop, rwi, log_tt, log_dens), by = "district_id") %>%
  arrange(district_id) %>%
  mutate(z_rwi = as.numeric(scale(rwi)), z_log_tt = as.numeric(scale(log_tt)),
         z_log_dens = as.numeric(scale(log_dens)),
         idx = district_id,
         direct_usable = in_model)   # the 01 rule: 2+ clusters, 0 < p < 1, se > 0
d$district_id <- as.integer(d$district_id)
stopifnot(nrow(d) == n_districts, all(d$district_id == seq_len(n_districts)))

## ---- 1. Sampling variances (MOZSAE_VAR; see functions.R) ---------------
## Direct variances from 2 or 3 clusters are very noisy, yet the model
## treats them as known: a district whose variance is underestimated is
## over-trusted, and intervals under-cover. "smoothed" replaces them
## with a design-effect generalised variance function (GVF; Corral et
## al., 2022):
##   V_i = deff_pooled * p~_i (1 - p~_i) / n_households_i
## p~_i is a stabilised share that does not use district i's own direct
## estimate, so the variance is not correlated with the sampling error:
## the leave-one-out synthetic prediction from a least-squares
## regression of logit(p_direct) on the covariates, fitted to the
## districts with a usable direct estimate (exact leave-one-out from the
## hat values). The province direct estimate was not used: it includes
## district i's own clusters and does not vary within a province.
## deff_pooled is the mean of deff_i = V_direct_i / (p~_i (1 - p~_i) / n_i)
## weighted by degrees of freedom (clusters - 1), i.e. a pooled variance,
## unbiased because each V_direct_i is. The median was not used: a
## variance from 2 or 3 clusters is right-skewed (median about 0.7 times
## the mean at 2 df), so the median would bring back the
## underestimation the smoothing is meant to remove. With the DHS, the
## design effect also absorbs the person weights and household size.
u <- which(d$direct_usable)
pre <- lm(qlogis(p_direct) ~ z_rwi + z_log_tt + z_log_dens, data = d[u, ])
eta_tilde <- unname(predict(pre, newdata = d))
eta_tilde[u] <- qlogis(d$p_direct[u]) - unname(residuals(pre) / (1 - hatvalues(pre)))   # leave one out
d <- d %>% mutate(p_tilde = plogis(eta_tilde),
                  deff = ifelse(direct_usable, se_direct^2 / (p_tilde * (1 - p_tilde) / n_households), NA_real_))
deff_pooled <- with(d[u, ], sum((n_clusters - 1) * deff) / sum(n_clusters - 1))

## "smoothed1": single-cluster districts (no direct variance) also enter
## the likelihood, with the smoothed variance
if (var_method == "smoothed1")
  d$in_model <- coalesce(d$n_clusters >= 1 & d$p_direct > 0 & d$p_direct < 1, FALSE)
smooth <- var_method != "direct"

## Districts with 2+ clusters and a direct share of exactly 0 or 1
## (MOZSAE_ZERO; see functions.R). Excluding them drops the least (or
## most) poor districts from the likelihood and predicts them from
## covariates and neighbours, which pulls them towards the middle. With
## "include", they enter with the empirical logit, y = logit((x + 0.5) /
## (n_eff + 1)), x = p_direct * n_eff, n_eff = n_households / deff_pooled,
## and the smoothed variance, which does not depend on their own estimate.
d$zero_share <- coalesce(d$n_clusters >= 2 & (d$p_direct == 0 | d$p_direct == 1), FALSE)
if (zero_method == "include") {
  stopifnot("MOZSAE_ZERO = include needs smoothed sampling variances" = smooth)
  d$in_model <- d$in_model | d$zero_share
}
d <- d %>% mutate(
  ## share used in the likelihood (empirical logit for 0 and 1 shares)
  p_lik = case_when(!in_model ~ NA_real_,
                    zero_share ~ (p_direct * n_households / deff_pooled + 0.5) / (n_households / deff_pooled + 1),
                    TRUE ~ p_direct),
  var_p = case_when(!in_model ~ NA_real_, smooth ~ deff_pooled * p_tilde * (1 - p_tilde) / n_households,
                    TRUE ~ se_direct^2),                                   # proportion scale (04)
  V = case_when(!in_model ~ NA_real_, smooth ~ var_p / (p_tilde * (1 - p_tilde))^2,
                TRUE ~ se_direct^2 / (p_direct * (1 - p_direct))^2),      # logit scale, delta method
  n_eff_fh = case_when(!in_model ~ NA_real_, smooth ~ n_households / deff_pooled, TRUE ~ n_eff),   # arcsine FH (04)
  y = ifelse(in_model, qlogis(p_lik), NA_real_),
  prec_scale = ifelse(in_model, 1 / V, 1))   # NA rows: scale unused
stopifnot(!anyNA(d$V[d$in_model]), all(d$V[d$in_model] > 0))

q <- quantile(d$deff[u], c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1))
v_direct_logit <- with(d[u, ], se_direct^2 / (p_direct * (1 - p_direct))^2)
v_smooth_logit <- with(d[u, ], deff_pooled / (n_households * p_tilde * (1 - p_tilde)))
design_effects <- tibble(
  item = c("Sampling variances used (MOZSAE_VAR)", "Districts with a direct variance",
           "Pooled design effect (mean weighted by clusters - 1)", "Unweighted mean of district design effects",
           "District design effect: minimum", "10th percentile", "25th percentile", "Median",
           "75th percentile", "90th percentile", "Maximum",
           "Median ratio, smoothed to direct variance (logit scale)",
           "R2 of the preliminary regression for the stabilised share",
           "Districts with 2+ clusters and a direct share of 0 or 1 (MOZSAE_ZERO)",
           "These districts in the likelihood"),
  value = c(var_method, length(u), round(c(deff_pooled, mean(d$deff[u]), q,
                                           median(v_smooth_logit / v_direct_logit), summary(pre)$r.squared), 3),
            paste0(sum(d$zero_share), " (", zero_method, ")"), sum(d$zero_share & d$in_model)))

## ---- 1b. Neighbourhood graph (queen contiguity) -------------------------
nb <- poly2nb(adm2, queen = TRUE)
n_islands <- sum(card(nb) == 0)
W <- nb2mat(nb, style = "B", zero.policy = TRUE)
g <- inla.read.graph(Matrix(W, sparse = TRUE))
## INLA scales the BYM2 structured effect by connected component and
## treats islands as IID (adjust.for.con.comp, default TRUE).

## ---- 2. Models -----------------------------------------------------
pc_sd  <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
pc_bym <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
               phi  = list(prior = "pc", param = c(0.5, 2/3)))
xs <- "z_rwi + z_log_tt + z_log_dens"
forms <- list(
  M1 = as.formula("y ~ 1 + f(idx, model = 'iid', hyper = pc_sd)"),
  M2 = as.formula(paste("y ~", xs, "+ f(idx, model = 'iid', hyper = pc_sd)")),
  M3 = as.formula(paste("y ~", xs, "+ f(idx, model = 'bym2', graph = g, scale.model = TRUE, constr = TRUE, hyper = pc_bym)")),
  M4 = as.formula("y ~ 1 + f(idx, model = 'bym2', graph = g, scale.model = TRUE, constr = TRUE, hyper = pc_bym)"))

fit_fh <- function(form) {
  inla(form, family = "gaussian", data = d, scale = d$prec_scale,
       control.family = list(hyper = list(prec = list(initial = 0, fixed = TRUE))),
       control.predictor = list(compute = TRUE),
       ## return.marginals.predictor: off by default in recent INLA; needed
       ## for the district posteriors in section 4
       control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE,
                              return.marginals.predictor = TRUE),
       control.inla = list(strategy = "laplace"))
}
fits <- map(forms, fit_fh)

## ---- 3. Model comparison -------------------------------------------
obs <- which(d$in_model)
comparison <- imap_dfr(fits, function(m, name) {
  cpo <- m$cpo$cpo[obs]; fail <- m$cpo$failure[obs]
  tibble(model = name, dic = m$dic$dic, waic = m$waic$waic,
         log_score = -mean(log(cpo)),            # lower is better
         cpo_failures = sum(fail > 0, na.rm = TRUE),
         pit_share_outside_05_95 = mean(m$cpo$pit[obs] < 0.05 | m$cpo$pit[obs] > 0.95))
})

fixed <- imap_dfr(fits, \(m, name) m$summary.fixed %>% as_tibble(rownames = "term") %>%
                   transmute(model = name, term, mean, sd, q025 = `0.025quant`, q975 = `0.975quant`))
hyper <- imap_dfr(fits, \(m, name) m$summary.hyperpar %>% as_tibble(rownames = "parameter") %>%
                   transmute(model = name, parameter, mean, sd, q025 = `0.025quant`, q975 = `0.975quant`))

## ---- 4. District estimates, all models ---------------------------
## Posterior of p_i = expit(theta_i) from the marginal of the linear
## predictor; exceedance P(p_i > 0.5): more than half the population in
## the national poorest 40%.
district_est <- function(m, name) {
  map_dfr(seq_len(nrow(d)), function(i) {
    marg <- m$marginals.linear.predictor[[i]]
    zm <- inla.zmarginal(inla.tmarginal(plogis, marg), silent = TRUE)
    tibble(district_id = d$district_id[i], model = name,
           p_mean = zm$mean, p_sd = zm$sd, p_q025 = zm$quant0.025, p_q975 = zm$quant0.975,
           pr_above_50 = 1 - inla.pmarginal(qlogis(0.5), marg))
  })
}
estimates <- imap_dfr(fits, district_est) %>%
  left_join(select(d, district_id, district, province, n_clusters, in_model, direct_usable), by = "district_id") %>%
  mutate(cv = p_sd / p_mean)

## Precision gain for districts in the likelihood: CV of the direct
## estimate against the posterior CV
gain <- estimates %>% filter(direct_usable) %>%
  left_join(select(d, district_id, p_direct, se_direct), by = "district_id") %>%
  mutate(cv_direct = se_direct / p_direct) %>%
  group_by(model) %>%
  summarise(districts = n(), median_cv_direct = median(cv_direct), median_cv_model = median(cv),
            median_ratio = median(cv / cv_direct), share_cv_below_20pct_direct = mean(cv_direct < 0.2),
            share_cv_below_20pct_model = mean(cv < 0.2), .groups = "drop")

## ---- 5. Posterior draws of p for M3 (province aggregation in 05) ----
n_draws <- 1000
draws <- inla.posterior.sample(n_draws, fits$M3)
pred_rows <- grep("^Predictor:", rownames(draws[[1]]$latent))
stopifnot(length(pred_rows) == nrow(d))
p_draws <- plogis(sapply(draws, \(s) s$latent[pred_rows, 1]))   # n_districts x n_draws
saveRDS(p_draws, file.path(dir_work, "M3_p_draws.rds"))
saveRDS(fits, file.path(dir_work, "inla_fits.rds"))
saveRDS(d, file.path(dir_work, "model_data.rds"))

## ---- 6. Write (district-level model output only) ---------------------
rnd <- \(x) mutate(x, across(where(is.double), \(v) round(v, 4)))
write_out(rnd(comparison), "s03_model_comparison")
write_out(rnd(fixed), "s03_fixed_effects")
write_out(rnd(hyper), "s03_hyperparameters")
write_out(rnd(gain), "s03_precision_gain")
write_out(design_effects, "s03_design_effects")
## Published: cluster counts grouped; in_model and direct_usable are left
## out, since with the cluster group they would show which named
## districts have a direct estimate of exactly 0 or 1 (both in
## model_data.rds in data_work/)
write_out(rnd(estimates %>% mutate(clusters = cluster_group(n_clusters)) %>%
                select(model, district_id, district, province, clusters,
                                   p_mean, p_sd, p_q025, p_q975, cv, pr_above_50)),
          "s03_district_estimates")
write_out(tibble(item = c("Districts", "Districts in likelihood", "Districts predicted out of sample",
                          "Islands in the contiguity graph"),
                 value = c(nrow(d), length(obs), nrow(d) - length(obs), n_islands)),
          "s03_model_data_summary")
message("03_models_inla.R finished.")
