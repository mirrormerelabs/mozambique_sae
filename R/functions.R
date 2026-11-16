############################################################
## functions.R
## Shared paths, packages and helpers for the Mozambique district
## poverty SAE demonstration. Sourced by every numbered script.
############################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr)
  library(haven); library(sf); library(terra)
})
sf_use_s2(FALSE)  # planar predicates on EPSG:4326 polygons, as in Paper A

stopifnot("Run from the project root (mozambique_sae/)" = file.exists("run_all.R"))

## ---- Input paths ---------------------------------------------------
## Inputs are read in place from the Paper A project and never copied.
## Set MOZSAE_PAPER_A to point elsewhere (e.g. on another machine).
path_paper_a <- Sys.getenv("MOZSAE_PAPER_A",
                           "../Thesis Corrections/Ch 4 reanalysis/analysis")
path_ext  <- file.path(path_paper_a, "data_raw/external")
path_dhs  <- file.path(path_ext, "dhs/MZ_2022-23_DHS_10162026_1741_262354")
file_hr   <- file.path(path_dhs, "MZHR82DT/MZHR82FL.dta")
file_ge   <- file.path(path_dhs, "MZGE81FL/MZGE81FL.shp")
file_adm2 <- file.path(path_ext, "geoboundaries/geoBoundaries-MOZ-ADM2.geojson")
file_adm1_zip <- file.path(path_ext, "gadm41/gadm41_MOZ_1.json.zip")
file_rwi  <- file.path(path_ext, "chi2022/moz_relative_wealth_index.csv")
file_meta_zip <- file.path(path_ext, "meta_population2020/moz_general_2020_csv.zip")
file_tt   <- file.path(path_ext, "nelson2019/travel_time_to_cities_11.tif")  # nearest city >= 50,000

stopifnot("Paper A inputs not found: set MOZSAE_PAPER_A" =
            all(file.exists(c(file_hr, file_ge, file_adm2, file_adm1_zip, file_rwi, file_meta_zip, file_tt))))

## ---- Output paths ----------------------------------------------------
## data_work/: household, cluster and direct-estimate files (DHS terms:
##             never published; listed in .gitignore).
## output/:    district-level model estimates, covariates from open
##             data, summaries and figures. Safe to publish.
## Simulation mode (MOZSAE_SIM=1; see run_sim.R): synthetic direct
## estimates from R/sim_01_synthetic.R replace 01, and every output goes
## to sim/ (excluded from version control), labelled as simulated.
## MOZSAE_SIM_RUN (simulation only) names a subfolder of sim/ for one
## run of a multi-seed study (R/sim_variance_study.R); it is always
## inside sim/.
sim_mode <- identical(Sys.getenv("MOZSAE_SIM"), "1")
sim_run <- Sys.getenv("MOZSAE_SIM_RUN")
stopifnot("MOZSAE_SIM_RUN must be a plain subfolder name" = !grepl("\\.\\.|^/", sim_run))
dir_root <- if (!sim_mode) "." else if (nzchar(sim_run)) file.path("sim", sim_run) else "sim"
dir_work <- file.path(dir_root, "data_work"); dir_out <- file.path(dir_root, "output")
dir_fig <- file.path(dir_out, "figures"); dir_docs <- file.path(dir_root, "docs")
for (d in c(dir_work, dir_out, dir_fig)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
if (sim_mode) message("SIMULATION MODE: reading and writing sim/ only")

## Cluster counts per district are published only in groups (DHS data
## protection): exact counts stay in data_work/.
cluster_levels <- c("0", "1", "2-4", "5+")
cluster_group <- function(n) factor(cut(n, c(-Inf, 0, 1, 4, Inf), labels = cluster_levels), levels = cluster_levels)

write_out <- function(x, name) {
  f <- file.path(dir_out, paste0(name, ".csv"))
  write_csv(x, f); message("Written: ", f, " (", nrow(x), " rows)")
  invisible(f)
}

## ---- Sampling variances --------------------------------------------
## MOZSAE_VAR chooses the sampling variances given to the Fay-Herriot
## models in 03 and 04 (see 03, section 1b, and the README):
##   "smoothed"  (default) design-effect generalised variance function:
##               V_i = deff_pooled * p~_i (1 - p~_i) / n_households_i
##   "direct"    the design-based variance of each direct estimate
##   "smoothed1" as "smoothed", and districts with a single cluster
##               (no direct variance) also enter the likelihood
var_method <- Sys.getenv("MOZSAE_VAR", "smoothed")
stopifnot("MOZSAE_VAR must be direct, smoothed or smoothed1" = var_method %in% c("direct", "smoothed", "smoothed1"))
## MOZSAE_ZERO: districts with 2+ clusters and a direct share of exactly
## 0 or 1 (no direct variance; logit undefined; see 03, section 1):
##   "exclude" (default) predicted out of sample, the rule of the first
##             real run; kept after the simulation study of 16 November
##             2026 (README, Deviations)
##   "include" in the likelihood, with the empirical logit
##             logit((x + 0.5) / (n_eff + 1)) and the smoothed variance
##             (needs MOZSAE_VAR "smoothed" or "smoothed1"); sensitivity
zero_method <- Sys.getenv("MOZSAE_ZERO", "exclude")
stopifnot("MOZSAE_ZERO must be exclude or include" = zero_method %in% c("exclude", "include"))

## ---- Sensitivity runs on the real data -------------------------------
## A real run with MOZSAE_VAR or MOZSAE_ZERO away from its default writes
## everything to data_work/sensitivity/<MOZSAE_VAR>_<MOZSAE_ZERO>/
## (private, never published), never to output/ or docs/. Run it with
## R/sens_real.R, which copies the inputs from 01 and 02 there first.
sens_mode <- !sim_mode && (var_method != "smoothed" || zero_method != "exclude")
if (sens_mode) {
  dir_root <- file.path("data_work", "sensitivity", paste(var_method, zero_method, sep = "_"))
  dir_work <- file.path(dir_root, "data_work"); dir_out <- file.path(dir_root, "output")
  dir_fig <- file.path(dir_out, "figures"); dir_docs <- file.path(dir_root, "docs")
  for (d in c(dir_work, dir_out, dir_fig)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
  message("SENSITIVITY RUN: reading and writing ", dir_root, " only")
}

## ---- Districts -------------------------------------------------------
## geoBoundaries ADM2 has 159 polygons. Two are Likoma and Chizumulu,
## islands in Lake Niassa that belong to Malawi ("Ilha Licom", "Ilha
## Risunodo"; no DHS clusters, no RWI tiles); they are dropped, leaving
## 157 districts. district_id is the row order after dropping, fixed
## here and used everywhere (INLA's spatial index needs 1..n).
not_mozambique <- c("Ilha Licom", "Ilha Risunodo")
n_districts <- 157
read_districts <- function() {
  adm2 <- st_read(file_adm2, quiet = TRUE) %>%
    st_make_valid() %>%
    filter(!shapeName %in% not_mozambique) %>%
    transmute(district_id = row_number(), district = shapeName)
  stopifnot(nrow(adm2) == n_districts, !anyDuplicated(adm2$district_id))
  adm2
}

## Province names: lower case, no accents, Maputo city spelled one way,
## so DHS (hv024) and GADM (NAME_1) labels can be matched.
norm_province <- function(x) {
  x <- tolower(iconv(x, to = "ASCII//TRANSLIT"))
  x <- gsub("[^a-z ]", "", x)
  x <- trimws(gsub("\\s+", " ", x))
  ## GADM 4.1 NAME_1 runs words together and misspells Niassa
  x[x == "cabodelgado"] <- "cabo delgado"; x[x == "nassa"] <- "niassa"
  x[x %in% c("maputo city", "maputocity", "maputo cidade", "cidade de maputo", "maputo cidade capital")] <- "maputo city"
  x[x %in% c("maputo province", "maputo provincia", "provincia de maputo")] <- "maputo"
  x
}

## ---- Plot settings ---------------------------------------------------
## Sequential magnitude: one hue, light to dark (reference blue ramp).
seq_low <- "#cde2fb"; seq_high <- "#0d366b"; na_fill <- "#e8e7e3"
## Diverging scale for the district share (Figure 1a, web map estimate
## layer), centred on the national share of 0.40 (the poorest 40% by
## construction): orange below (less poor than the country as a whole),
## neutral grey at 0.40, blue above (poorer). Each arm spans the full
## colour range, so the arms are balanced although 0.4 to 1 is shorter
## than 0 to 0.4. Orange and blue, with dark ends of similar lightness
## (CIE L 29 and 23) and saturated mid-arm stops, at the middle of each
## arm (L 58 and 55), so that the arms do not turn grey; no red-green
## contrast (checked by simulating the three dichromacies). Interpolated
## in CIE Lab by scales.
## Uncertainty (Figure 1b, web map uncertainty layer): a violet ramp,
## light to dark, sharing no hue with the diverging scale below. A grey
## ramp was tried, but the black no-cluster outlines vanished against
## its darkest districts.
unc_low <- "#e4e0f7"; unc_high <- "#2e2470"
div_mid <- 0.40
div_cols <- c("#7a2a0c", "#d2733e", "#f0efec", "#4a86d0", "#0d366b")
div_vals <- c(0, div_mid / 2, div_mid, (1 + div_mid) / 2, 1)   # positions of div_cols on [0, 1]
div_pal <- scales::gradient_n_pal(div_cols, values = div_vals)   # maps [0, 1] to colours
div_breaks <- seq(0, 1, by = 0.2)
div_labels <- ifelse(abs(div_breaks - div_mid) < 1e-9, "40% (national)", paste0(round(100 * div_breaks), "%"))

## Uncertainty colour scale (Figure 1b, web map), for the width of the
## 95% interval (a proportion; labels in percentage points). If a few
## districts would take most of the range (maximum above 1.25 times the
## 90th percentile), the scale is capped at the 90th percentile rounded
## to 10 points, values above take the end colour (scales::squish), and
## the top label reads "+", as for the covariates.
capped_scale <- function(x) {
  p90 <- quantile(x, 0.9, na.rm = TRUE)[[1]]; top <- max(x, na.rm = TRUE)
  cap <- if (top > 1.25 * p90) max(0.1, round(p90, 1)) else ceiling(top * 10) / 10
  br <- round(unique(c(seq(0, cap, by = if (cap <= 0.5) 0.1 else 0.2), cap)), 2)
  list(limits = c(0, cap), breaks = br,
       labels = paste0(round(100 * br), ifelse(br == round(cap, 2) & top > cap, "+", "")))
}
theme_sae <- function() {
  ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
                   legend.position = "right", plot.title = ggplot2::element_text(size = 10, face = "bold"),
                   plot.caption = ggplot2::element_text(size = 7, colour = "grey35", hjust = 0))
}
