############################################################
## 05_outputs.R
## Figures for the repository README, and the province benchmark.
##
## Run after 03 and 04. (Figure S1, the covariate maps, is written by
## 02_covariates.R, since it needs public data only.)
##
##   Figure 1  Lead map: M3 posterior mean (diverging scale centred on
##             the national share, 40%) and the width of its 95%
##             interval, side by side; districts without survey clusters
##             outlined
##   Figure 2  Precision gain: width of the 95% interval of the direct
##             estimate and of M3 for the same districts, by cluster
##             group; districts without a usable direct estimate in a
##             separate panel (M3 only)
##   Figure 3  Districts ranked by M3 posterior mean, with 95% intervals,
##             marked where P(share > 50%) >= 0.9
##   Figure 4  Province benchmark: M3 aggregated with population weights
##             against the DHS direct province estimates
##   Figure 5  Model check: direct against M3 for districts in the
##             likelihood. Written to data_work/figures/ only: its y
##             values are the published M3 estimates, so its points
##             could be matched to named districts.
##   Figure S2 Survey clusters per district, and P(share > 50%)
##
## DATA PROTECTION: no cluster location is drawn, and direct district
## estimates are never mapped or labelled with a district name (some
## rest on two clusters). Figure 2 shows direct interval widths as unnamed,
## unlinked points only, by cluster group; Figure 5 is not published.
############################################################

source("R/functions.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("R/figures.R")

adm2 <- read_districts()
stopifnot("Run 01, 03 and 04 first" =
            all(file.exists(file.path(dir_work, c("model_data.rds", "M3_p_draws.rds", "direct_province.rds")))))

## ---- Model outputs -----------------------------------------------------
d <- readRDS(file.path(dir_work, "model_data.rds"))
p_draws <- readRDS(file.path(dir_work, "M3_p_draws.rds"))
prov_direct <- readRDS(file.path(dir_work, "direct_province.rds"))
est <- read_csv(file.path(dir_out, "s03_district_estimates.csv"), show_col_types = FALSE)
## Exact cluster counts from data_work/ (the published CSV has groups only)
m3 <- est %>% filter(model == "M3") %>% left_join(select(d, district_id, n_clusters), by = "district_id")
map_df <- adm2 %>% left_join(select(m3, district_id, p_mean, p_q025, p_q975, cv, pr_above_50, n_clusters),
                             by = "district_id")

## ---- 1. Province benchmark ---------------------------------------------
stopifnot(nrow(p_draws) == nrow(d))
benchmark <- map_dfr(unique(d$province), function(pv) {
  i <- which(d$province == pv)
  w <- d$pop[i] / sum(d$pop[i])
  agg <- colSums(p_draws[i, , drop = FALSE] * w)
  tibble(province = pv, districts = length(i), model_mean = mean(agg),
         model_q025 = quantile(agg, 0.025), model_q975 = quantile(agg, 0.975))
}) %>%
  left_join(prov_direct %>% rename(direct = p_direct, direct_lo = lo, direct_hi = hi), by = "province") %>%
  mutate(difference = model_mean - direct, z = difference / se_direct,
         model_mean_in_direct_ci = model_mean >= direct_lo & model_mean <= direct_hi)
write_out(benchmark %>% mutate(across(where(is.double), \(v) round(v, 4))), "s05_province_benchmark")

## ---- Figure 1: lead map ------------------------------------------------
cap1 <- paste0("Share of the population in the national poorest 40% (DHS wealth quintiles 1 and 2), ", nrow(d),
               " districts. Colours are centred on the national share of 40%:",
               " orange districts are less poor than the country as a whole, blue districts poorer.",
               "\nSpatial Fay-Herriot model (BYM2, R-INLA). Sources: ", src_dhs, "; ",
               sub("Sources: ", "", src_open))
## Uncertainty as the width of the 95% interval (percentage points): the
## CV is unstable for very low shares, and dominated the map in the
## south. The CV stays in the CSVs and the README summary.
map_df <- map_df %>% mutate(width = p_q975 - p_q025)
## Districts without survey clusters: outlined in both panels, with a
## legend entry in panel (a)
nocl_label <- "No survey clusters: estimate from\ncovariates and neighbouring districts"
nocl <- function(show) list(
  geom_sf(data = filter(map_df, n_clusters == 0), aes(colour = nocl_label), fill = NA, linewidth = 0.45,
          inherit.aes = FALSE, show.legend = show),
  scale_colour_manual(values = setNames("#1a1a1a", nocl_label), name = NULL))
## Panel (a): diverging scale centred on the national share (functions.R)
f1 <- (district_map(map_df, "p_mean") + nocl(TRUE) +
         scale_fill_gradientn(colours = div_cols, values = div_vals, limits = c(0, 1),
                               breaks = div_breaks, labels = div_labels, na.value = na_fill, name = NULL) +
         labs(title = "(a) Share in the national poorest 40%,\ncentred on the national share")) |
  (district_map(map_df, "width") + nocl(FALSE) +
     with(capped_scale(map_df$width), seq_fill_unc(name = "Percentage\npoints", limits = limits, breaks = breaks,
                                                      labels = labels, oob = scales::squish)) +
     labs(title = "(b) Uncertainty: width of the 95% interval"))
save_fig(f1 + plot_annotation(caption = cap1, theme = theme_sae()), "fig1_estimate_uncertainty", 10, 6.5)

## ---- Figure 2: precision gain ------------------------------------------
## Uncertainty as the width of the 95% interval (percentage points), not
## the CV, which explodes for shares near 0. Left panel: districts with a
## usable direct estimate (the 01 rule, direct_usable; with the default
## switches these are the districts in the likelihood), direct and M3 for
## the same districts. Right panel: every other district (0 or 1 cluster,
## or 2+ clusters with a direct share of 0 or 1 or no direct variance),
## M3 only. Direct width: 2 x 1.96 x se_direct, clipped to [0, 1].
## Points are unnamed and not linked district by district: a pairing line
## would tie a direct standard error to a named M3 estimate.
panel_in <- "Usable direct estimate:\nsurvey and model, same districts"
panel_out <- "No usable direct estimate:\nmodel only"
gain_df <- m3 %>% left_join(select(d, district_id, p_direct, se_direct, direct_usable), by = "district_id") %>%
  mutate(panel = factor(ifelse(direct_usable, panel_in, panel_out), levels = c(panel_in, panel_out)),
         ## Published cluster groups only: exact counts would identify districts
         clusters = cluster_group(n_clusters))
gain_long <- bind_rows(
  gain_df %>% filter(direct_usable) %>%
    transmute(panel, clusters, width = pmin(pmax(2 * 1.96 * se_direct, 0), 1), estimate = "Direct (survey only)"),
  gain_df %>% transmute(panel, clusters, width = p_q975 - p_q025, estimate = "Spatial Fay-Herriot (M3)")) %>%
  mutate(estimate = factor(estimate, levels = c("Direct (survey only)", "Spatial Fay-Herriot (M3)")))
gain_med <- gain_long %>% group_by(panel, clusters, estimate) %>%
  summarise(districts = n(), median_width = median(width), .groups = "drop")
jit <- position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.75, seed = 1)
dod <- position_dodge(width = 0.75)
f2 <- ggplot(gain_long, aes(clusters, 100 * width, colour = estimate)) +
  geom_point(position = jit, size = 1.5, alpha = 0.6) +
  geom_crossbar(data = gain_med, aes(y = 100 * median_width, ymin = 100 * median_width, ymax = 100 * median_width),
                position = dod, width = 0.6, linewidth = 0.5, middle.linewidth = 0.5, show.legend = FALSE) +
  facet_grid(~ panel, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = c(grey_mark, blue), name = NULL, drop = FALSE) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "DHS clusters in the district", y = "Width of the 95% interval (percentage points)",
       title = "Uncertainty of district estimates, survey alone and with the model",
       caption = paste("Horizontal bars: medians. Direct interval: estimate ± 1.96 standard errors.",
                       "Right panel: districts with 0 or 1 cluster, or 2+ clusters with a direct share of 0 or 1, predicted from covariates and neighbours.",
                       paste0(src_dhs, "."), sep = "\n")) +
  plot_theme() + theme(legend.position = "top", strip.text = element_text(hjust = 0, face = "bold"),
                       panel.spacing.x = unit(1.2, "lines"))
save_fig(f2, "fig2_precision_gain", 8, 5)

## ---- Figure 3: ranked districts ------------------------------------------
rank_df <- m3 %>% arrange(p_mean) %>% mutate(rank = row_number(),
  confident = ifelse(pr_above_50 >= 0.9, "P(share > 50%) at least 0.9", "Other districts"))
n_conf <- sum(rank_df$pr_above_50 >= 0.9)
f3 <- ggplot(rank_df, aes(rank, p_mean, colour = confident)) +
  geom_hline(yintercept = 0.5, colour = "grey55", linewidth = 0.4, linetype = "dashed") +
  geom_linerange(aes(ymin = p_q025, ymax = p_q975), linewidth = 0.5, alpha = 0.7) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = c("Other districts" = blue_light, "P(share > 50%) at least 0.9" = "#0d366b"), name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Districts, ranked by estimated share", y = "Share in national poorest 40% (95% interval)",
       title = paste0(n_conf, " of ", nrow(rank_df), " districts can be identified with confidence as majority-poor"),
       caption = paste0("Spatial Fay-Herriot model (M3). Dashed line: 50%. District names are given in the interactive map.\n",
                        src_dhs, ".")) +
  plot_theme() + theme(legend.position = "top", panel.grid.major.x = element_blank())
save_fig(f3, "fig3_ranked_districts", 8, 4.5)

## ---- Figure 4: province benchmark -------------------------------------
f4 <- benchmark %>%
  mutate(province = reorder(tools::toTitleCase(province), direct)) %>%
  ggplot(aes(y = province)) +
  geom_linerange(aes(xmin = direct_lo, xmax = direct_hi), colour = grey_mark, linewidth = 0.6,
                 position = position_nudge(y = 0.15)) +
  geom_point(aes(x = direct, colour = "DHS direct estimate (95% CI)"), size = 2, position = position_nudge(y = 0.15)) +
  geom_linerange(aes(xmin = model_q025, xmax = model_q975), colour = blue_light, linewidth = 0.6,
                 position = position_nudge(y = -0.15)) +
  geom_point(aes(x = model_mean, colour = "Model, districts aggregated (95% interval)"), size = 2,
             position = position_nudge(y = -0.15)) +
  scale_colour_manual(values = c(grey_mark, blue), name = NULL) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Share in national poorest 40%", y = NULL,
       title = "Province benchmark: model against the survey where the survey is reliable",
       caption = paste0("Model: district posterior draws weighted by Meta 2020 population. ", src_dhs, ".")) +
  plot_theme() + theme(legend.position = "top", panel.grid.major.y = element_blank())
save_fig(f4, "fig4_province_benchmark", 7, 5)

## ---- Figure 5: model check, direct against model (unnamed) ------------
sc <- gain_df %>% filter(direct_usable)
f5 <- ggplot(sc, aes(p_direct, p_mean)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.4, linetype = "dashed") +
  geom_linerange(aes(xmin = pmax(0, p_direct - 1.96 * se_direct), xmax = pmin(1, p_direct + 1.96 * se_direct)),
                 colour = "grey80", linewidth = 0.3) +
  geom_linerange(aes(ymin = p_q025, ymax = p_q975), colour = blue_light, linewidth = 0.3) +
  geom_point(colour = blue, size = 1.6) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  coord_equal() +
  labs(x = "Direct estimate (95% CI)", y = "Model estimate (95% interval)",
       title = "Shrinkage of direct estimates towards the model",
       caption = paste0(nrow(sc), " districts with 2+ clusters. Dashed line: equality.\n", src_dhs, ".")) +
  plot_theme()
save_fig(f5, "fig5_shrinkage", 5.5, 5.5, dir = file.path(dir_work, "figures"))   # not published

## ---- Figure S2: clusters per district and P(share > 50%) ---------------
## Cluster counts mapped in the published groups (0, 1, 2-4, 5+)
s2 <- (district_map(map_df %>% mutate(n_cl = cluster_group(n_clusters)), "n_cl") +
         scale_fill_manual(values = colorRampPalette(c("#f0efec", "#256abf"))(4), name = NULL, drop = FALSE) +
         labs(title = "(a) DHS survey clusters per district")) |
  (district_map(map_df, "pr_above_50") + seq_fill(name = NULL, limits = c(0, 1), labels = scales::percent) +
     labs(title = "(b) P(share > 50%)"))
save_fig(s2 + plot_annotation(caption = paste0(src_dhs, "; geoBoundaries ADM2."), theme = theme_sae()),
         "figS2_clusters_exceedance", 10, 6.5)

## ---- Population in confidently majority-poor districts ---------------
## Districts with P(share > 50%) >= 0.9: their Meta 2020 population, and
## the number of people in the national poorest 40% living there (sum of
## pop x p over those districts, for each M3 posterior draw)
conf <- d$district_id %in% m3$district_id[m3$pr_above_50 >= 0.9]
pop_conf <- sum(d$pop[conf]); pop_share_conf <- pop_conf / sum(d$pop)
poor_conf <- colSums(p_draws[conf, , drop = FALSE] * d$pop[conf])
mn <- \(x) round(x / 1e6, 1)   # millions, to 0.1 million

## ---- README headline numbers -------------------------------------------
gain_m3 <- read_csv(file.path(dir_out, "s03_precision_gain.csv"), show_col_types = FALSE) %>% filter(model == "M3")
write_out(tibble(
  item = c("Districts", "Districts with 2+ clusters and a usable direct estimate", "Districts predicted out of sample",
           "Sampling variances (MOZSAE_VAR)",
           "Median CV, direct (districts in model)", "Median CV, M3 (same districts)",
           "Share of these districts with CV below 20%: direct", "Share with CV below 20%: M3",
           "Districts with P(share > 50%) >= 0.9", "Provinces where the model mean lies in the direct 95% CI"),
  value = c(as.character(c(nrow(d), sum(d$direct_usable), sum(!d$in_model))), var_method,
            as.character(round(c(gain_m3$median_cv_direct, gain_m3$median_cv_model,
                                 gain_m3$share_cv_below_20pct_direct, gain_m3$share_cv_below_20pct_model), 3)),
            as.character(n_conf), paste0(sum(benchmark$model_mean_in_direct_ci), " of ", nrow(benchmark)))) %>%
    ## Median width of the 95% interval (percentage points) by Figure 2 group
    bind_rows(tibble(
      item = c("Population of districts with P(share > 50%) >= 0.9 (million, Meta 2020)",
               "Share of the national population in these districts",
               "People in the national poorest 40% in these districts (million, posterior mean)",
               "People in the national poorest 40% in these districts (million, 95% interval)"),
      value = c(as.character(mn(pop_conf)), as.character(round(pop_share_conf, 3)), as.character(mn(mean(poor_conf))),
                paste(mn(quantile(poor_conf, c(0.025, 0.975))), collapse = " to ")))) %>%
    bind_rows(gain_med %>% transmute(
      item = paste0("Median 95% interval width (pp), ", sub(":.*", "", gsub("\n", " ", panel)), ", ",
                    clusters, " clusters, ", estimate, " (", districts, " districts)"),
      value = as.character(round(100 * median_width, 1)))),
  "s05_readme_numbers")
message("05_outputs.R finished.")
