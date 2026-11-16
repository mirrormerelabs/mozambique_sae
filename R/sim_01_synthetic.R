############################################################
## sim_01_synthetic.R
## Synthetic stand-in for 01_direct.R, for testing 03 to 06 before the
## DHS data may be used. NO DHS DATA ARE READ.
##
## Run with MOZSAE_SIM=1 (see run_sim.R): all outputs go to sim/, and
## figures and the web map are labelled as simulated.
##
## Data-generating process, on the real districts and real public
## covariates (02_covariates.R):
##   theta_i = beta0 + beta' z_i + u_i                 (true district logit)
##   u_i     = sigma * (sqrt(phi) * s_i + sqrt(1 - phi) * v_i)
##             s_i: spatially smoothed noise on the district graph,
##             v_i: IID noise (both standardised)
##   clusters per district ~ Poisson(617 x population share)
##   cluster share p_c = expit(m_i + e_c), e_c ~ N(0, tau^2), with m_i
##             set so that E[p_c] = expit(theta_i), the true share
##   households per cluster = 23; poor ~ Bernoulli(p_c)
## Direct estimates are cluster-mean estimators with between-cluster
## variance, mimicking a two-stage design, and the same in_model rule
## as 01. Province "direct" estimates are population-weighted over
## clusters. The truth is saved for the recovery check (sim_check.R).
############################################################

source("R/functions.R")
stopifnot("Set MOZSAE_SIM=1 (use run_sim.R)" = sim_mode)
suppressPackageStartupMessages(library(spdep))
set.seed(as.integer(Sys.getenv("MOZSAE_SIM_SEED", "20261116")))   # other seeds: sim_variance_study.R

## ---- Parameters (known truth) ---------------------------------------
## MOZSAE_SIM_SCENARIO:
##   "base"       (default) moderate shares; few districts with a direct
##                share of 0
##   "plasmode"   true district shares are the M3 posterior means of the
##                first real run (16 November 2026; public estimates,
##                snapshot in sim/plasmode_truth.rds, made by
##                sim_variance_study.R), with a cluster SD of 1.05,
##                matching its pooled design effect of about 6.5. Many
##                very low shares where clusters are, as in the south,
##                for testing MOZSAE_ZERO. The coefficients have no true
##                value (NA).
scenario <- Sys.getenv("MOZSAE_SIM_SCENARIO", "base")
stopifnot("MOZSAE_SIM_SCENARIO must be base or plasmode" = scenario %in% c("base", "plasmode"))
truth_par <- tibble(
  parameter = c("(Intercept)", "z_rwi", "z_log_tt", "z_log_dens", "sigma_total", "phi", "tau_cluster"),
  value     = if (scenario == "base") c(qlogis(0.40), -0.6, 0.4, -0.2, 0.5, 0.6, 0.8)
              else c(NA, NA, NA, NA, NA, NA, 1.05))
tp <- setNames(truth_par$value, truth_par$parameter)

## ---- Districts, provinces and covariates (public data) -------------
adm2 <- read_districts()
covs <- readRDS(file.path("data_work", "covariates_district.rds"))   # real covariates from 02
saveRDS(covs, file.path(dir_work, "covariates_district.rds"))          # copy into sim/data_work

adm1_tmp <- tempfile("adm1_"); dir.create(adm1_tmp)
unzip(file_adm1_zip, exdir = adm1_tmp)
adm1 <- st_read(list.files(adm1_tmp, pattern = "json$", full.names = TRUE), quiet = TRUE) %>%
  st_make_valid() %>% transmute(province = norm_province(NAME_1))
unlink(adm1_tmp, recursive = TRUE)
district_province <- st_join(st_point_on_surface(adm2), adm1, join = st_within) %>%
  st_drop_geometry() %>% select(district_id, province)

d <- adm2 %>% st_drop_geometry() %>%
  left_join(district_province, by = "district_id") %>%
  left_join(select(covs, district_id, pop, rwi, log_tt, log_dens), by = "district_id") %>%
  arrange(district_id) %>%
  mutate(z_rwi = as.numeric(scale(rwi)), z_log_tt = as.numeric(scale(log_tt)),
         z_log_dens = as.numeric(scale(log_dens)))   # same standardisation as 03
stopifnot(nrow(d) == n_districts, !anyNA(d))

## ---- True district effects --------------------------------------------
nb <- poly2nb(adm2, queen = TRUE)
W <- nb2mat(nb, style = "W", zero.policy = TRUE)
s <- rnorm(n_districts)
for (k in 1:4) s <- 0.5 * s + 0.5 * as.numeric(W %*% s)   # smooth over neighbours
s <- as.numeric(scale(s)); v <- as.numeric(scale(rnorm(n_districts)))
u <- tp[["sigma_total"]] * (sqrt(tp[["phi"]]) * s + sqrt(1 - tp[["phi"]]) * v)
d <- d %>% mutate(theta_true = tp[["(Intercept)"]] + tp[["z_rwi"]] * z_rwi + tp[["z_log_tt"]] * z_log_tt +
                    tp[["z_log_dens"]] * z_log_dens + u,
                  p_true = plogis(theta_true))
if (scenario == "plasmode") {
  f_plasmode <- file.path("sim", "plasmode_truth.rds")
  stopifnot("Missing sim/plasmode_truth.rds: run sim_variance_study.R" = file.exists(f_plasmode))
  d <- d %>% select(-theta_true, -p_true) %>%
    left_join(readRDS(f_plasmode), by = "district_id") %>%
    mutate(p_true = pmin(pmax(p_true, 0.001), 0.999), theta_true = qlogis(p_true))
  stopifnot(!anyNA(d$p_true))
}

## ---- Synthetic survey -------------------------------------------------
## The direct estimator targets the district share, the mean of p_c over
## all clusters. expit is not linear, so clusters centred on theta_true
## would give a share E[expit(theta_true + e_c)], pulled towards 0.5
## (by up to 0.03 here), and not p_true. Cluster logits are therefore
## centred on m_i, solved so that E[expit(m_i + e_c)] = p_true: the
## estimand is p_true, and the coefficients above are its true values.
share_given_centre <- function(m) integrate(\(e) plogis(m + e) * dnorm(e, 0, tp[["tau_cluster"]]), -Inf, Inf)$value
d$centre <- vapply(d$p_true, \(p) uniroot(\(m) share_given_centre(m) - p, c(-15, 15), tol = 1e-10)$root, 0)
n_hh <- 23
d$n_clusters <- rpois(n_districts, 617 * d$pop / sum(d$pop))
clusters <- d %>% select(district_id, province, pop, centre, n_clusters) %>%
  filter(n_clusters > 0) %>%
  uncount(n_clusters, .id = "cluster_in_district") %>%
  mutate(p_cluster = plogis(centre + rnorm(n(), 0, tp[["tau_cluster"]])),
         share_poor = rbinom(n(), n_hh, p_cluster) / n_hh)

direct_d <- clusters %>% group_by(district_id) %>%
  summarise(n_clusters = n(), p_direct = mean(share_poor),
            se_direct = ifelse(n() >= 2, sd(share_poor) / sqrt(n()), NA_real_), .groups = "drop")

direct <- d %>% select(district_id, district, province) %>%
  left_join(direct_d, by = "district_id") %>%
  mutate(n_clusters = coalesce(n_clusters, 0L), n_households = n_clusters * n_hh,
         in_model = coalesce(n_clusters >= 2 & p_direct > 0 & p_direct < 1 & se_direct > 0, FALSE),
         n_eff = ifelse(in_model, p_direct * (1 - p_direct) / se_direct^2, NA_real_))

## Province estimates: clusters weighted by district population / clusters
province_direct <- clusters %>%
  left_join(count(clusters, district_id, name = "n_cl"), by = "district_id") %>%
  mutate(w = pop / n_cl) %>%
  group_by(province) %>%
  summarise(p_direct = sum(w * share_poor) / sum(w),
            se_direct = sqrt(n() / (n() - 1) * sum((w * (share_poor - p_direct))^2)) / sum(w),
            .groups = "drop") %>%
  mutate(lo = p_direct - 1.96 * se_direct, hi = p_direct + 1.96 * se_direct)

## ---- Write (sim/ only) ---------------------------------------------------
saveRDS(direct, file.path(dir_work, "direct_district.rds"))
saveRDS(province_direct, file.path(dir_work, "direct_province.rds"))
saveRDS(select(d, district_id, theta_true, p_true), file.path(dir_work, "truth_district.rds"))
write_out(truth_par, "sim_truth_parameters")
write_out(tibble(item = c("Districts", "Synthetic clusters", "Districts with 0 clusters", "Districts in likelihood",
                          "Median clusters per sampled district", "Population-weighted mean true share"),
                 value = c(n_districts, nrow(clusters), sum(direct$n_clusters == 0), sum(direct$in_model),
                           median(direct$n_clusters[direct$n_clusters > 0]), round(weighted.mean(d$p_true, d$pop), 3))),
          "sim_summary")
message("sim_01_synthetic.R finished (SIMULATED DATA).")
