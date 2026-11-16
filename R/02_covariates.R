############################################################
## 02_covariates.R
## District covariates from open data, for all 157 districts:
##   rwi        population-weighted mean Relative Wealth Index
##              (Chi et al., 2022; tile population from Meta 2020,
##              as in Paper A, script 38)
##   log_tt     log(1 + population-weighted mean travel time, minutes,
##              to the nearest city of 50,000+) (Nelson et al., 2019)
##   log_dens   log population density, people per km2 (Meta 2020)
## Population is aggregated onto the travel-time grid (30 arc-seconds,
## about 1 km) and districts are rasterised onto the same grid.
## All inputs are public; outputs are publishable.
############################################################

source("R/functions.R")

adm2 <- read_districts()
adm2_v <- vect(adm2)

## ---- 1. Grid, population and travel time -----------------------------
tt <- rast(file_tt)
tt <- crop(tt, ext(adm2_v) + 0.1)
tt[tt < 0] <- NA
## Values above 48 hours on populated cells are artefacts of the
## friction surface (e.g. Ilha de Mocambique, linked to the mainland by
## a bridge, otherwise gets about 58 hours). Treated as missing; the
## population affected is reported in s02_covariate_checks.
tt_cap <- 48 * 60
tt_raw <- tt
tt[tt > tt_cap] <- NA
names(tt) <- "tt"

meta_tmp <- tempfile("meta_"); dir.create(meta_tmp)
unzip(file_meta_zip, exdir = meta_tmp)
meta <- read_csv(file.path(meta_tmp, "moz_general_2020.csv"), col_types = "ddd",
                 col_names = c("lon", "lat", "pop"), skip = 1)
unlink(meta_tmp, recursive = TRUE)
stopifnot(nrow(meta) > 1e6, !anyNA(meta$pop), all(meta$pop >= 0))

pop_r <- rasterize(as.matrix(meta[, c("lon", "lat")]), tt, values = meta$pop, fun = sum, background = 0)
names(pop_r) <- "pop"
dist_r <- rasterize(adm2_v, tt, field = "district_id")
names(dist_r) <- "district_id"

cells_per_district <- freq(dist_r) %>% as_tibble() %>% transmute(district_id = as.integer(value), cells = count)
stopifnot("Some districts have no grid cell; rasterise with touches = TRUE" =
            nrow(cells_per_district) == n_districts)

z_pop   <- zonal(pop_r, dist_r, fun = "sum", na.rm = TRUE)
pop_tt  <- mask(pop_r, tt)                         # population where travel time is known
z_ttpop <- zonal(c(tt * pop_tt, pop_tt), dist_r, fun = "sum", na.rm = TRUE)
names(z_ttpop) <- c("district_id", "tt_x_pop", "pop_with_tt")

## For the record. The Meta CSV covers the bounding rectangle, so it
## includes neighbouring countries; only district cells are Mozambique.
pop_total_meta <- sum(meta$pop); pop_in_districts <- sum(z_pop$pop)
pop_tt_capped <- global(mask(pop_r, tt_raw > tt_cap, maskvalues = c(0, NA)), "sum", na.rm = TRUE)[[1]]

## ---- 2. Relative Wealth Index, population-weighted --------------------
z14 <- 2^14
tile_x_of <- function(lon) floor((lon + 180) / 360 * z14)
tile_y_of <- function(lat) { r <- lat * pi / 180; floor((1 - log(tan(r) + 1 / cos(r)) / pi) / 2 * z14) }
rwi <- read_csv(file_rwi, col_select = c(latitude, longitude, rwi), show_col_types = FALSE) %>%
  rename(lat = latitude, lon = longitude) %>%
  mutate(tile_key = tile_x_of(lon) * 1e5 + tile_y_of(lat))
stopifnot(!anyDuplicated(rwi$tile_key), !anyNA(rwi$rwi))
tile_pop <- rowsum(meta$pop, tile_x_of(meta$lon) * 1e5 + tile_y_of(meta$lat))
rwi$pop <- as.numeric(tile_pop[match(rwi$tile_key, as.numeric(rownames(tile_pop)))])
rwi$pop[is.na(rwi$pop)] <- 0
rm(meta, tile_pop); invisible(gc())

## District of each tile centre (point in polygon; tiles off the grid
## are matched to the polygon directly)
rwi_pts <- st_as_sf(rwi, coords = c("lon", "lat"), crs = 4326)
rwi$district_id <- adm2$district_id[sapply(st_within(rwi_pts, adm2), \(i) if (length(i)) i[1] else NA_integer_)]
rwi_d <- rwi %>% filter(!is.na(district_id)) %>% group_by(district_id) %>%
  summarise(n_tiles = n(), tile_pop = sum(pop),
            rwi = if (sum(pop) > 0) weighted.mean(rwi, pop) else mean(rwi), .groups = "drop")

## Districts with no RWI tile: value of the nearest tile to a point on
## the district's surface (flagged)
no_rwi <- setdiff(adm2$district_id, rwi_d$district_id)
if (length(no_rwi) > 0) {
  near <- st_nearest_feature(st_point_on_surface(adm2[adm2$district_id %in% no_rwi, ]), rwi_pts)
  rwi_d <- bind_rows(rwi_d, tibble(district_id = no_rwi, n_tiles = 0L, tile_pop = 0, rwi = rwi$rwi[near]))
}

## ---- 3. Assemble ------------------------------------------------------
area_km2 <- expanse(adm2_v, unit = "km")
covariates <- adm2 %>% st_drop_geometry() %>%
  mutate(area_km2 = area_km2) %>%
  left_join(z_pop, by = "district_id") %>%
  left_join(z_ttpop, by = "district_id") %>%
  left_join(rwi_d, by = "district_id") %>%
  mutate(tt_mean = tt_x_pop / pop_with_tt,
         log_tt = log1p(tt_mean),
         log_dens = log(pop / area_km2),
         rwi_nearest_tile = n_tiles == 0) %>%
  select(district_id, district, area_km2, pop, rwi, n_tiles, rwi_nearest_tile, tt_mean, log_tt, log_dens)
stopifnot(nrow(covariates) == n_districts, !anyNA(covariates[c("rwi", "log_tt", "log_dens", "pop")]),
          all(covariates$pop > 0))

saveRDS(covariates, file.path(dir_work, "covariates_district.rds"))
write_out(covariates %>% mutate(across(where(is.double), \(x) round(x, 4))), "s02_covariates_district")
write_out(tibble(item = c("Meta CSV population, whole bounding rectangle incl. neighbouring countries",
                          "Population on Mozambique district cells",
                          "Population on cells with travel time > 48 h (set to missing)",
                          "Districts with no RWI tile (nearest tile used)",
                          "Pearson r, rwi vs log_tt", "Pearson r, rwi vs log_dens", "Pearson r, log_tt vs log_dens"),
                 value = round(c(pop_total_meta, pop_in_districts, pop_tt_capped, length(no_rwi),
                                 cor(covariates$rwi, covariates$log_tt), cor(covariates$rwi, covariates$log_dens),
                                 cor(covariates$log_tt, covariates$log_dens)), 3)),
          "s02_covariate_checks")
## ---- 4. Figure S1: covariate maps (public data) ------------------------
## Colour scales are clamped near the 90th percentile (squished beyond),
## so that rural districts, not the few cities, span the colour range.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("R/figures.R")
cov_map <- adm2 %>% left_join(covariates, by = c("district_id", "district"))
s1 <- (district_map(cov_map, "rwi") +
         seq_fill(name = NULL, limits = c(-0.55, 0.25), oob = scales::squish, breaks = c(-0.5, -0.25, 0, 0.25),
                  labels = c("-0.5", "-0.25", "0", "0.25+")) +
         labs(title = "(a) Relative Wealth Index", subtitle = "Population-weighted district mean")) |
  (district_map(cov_map, "tt_mean") + seq_fill(name = NULL, transform = "log10", limits = c(20, 240), oob = scales::squish,
              breaks = c(20, 60, 120, 240), labels = c("20 or less", "60", "120", "240+")) +
     labs(title = "(b) Travel time to a city of 50,000+", subtitle = "Population-weighted mean, minutes (log scale)")) |
  (district_map(cov_map, "log_dens") +
     seq_fill(name = NULL, limits = log(c(5, 500)), oob = scales::squish, breaks = log(c(5, 20, 100, 500)),
              labels = c("5 or less", "20", "100", "500+")) +
     labs(title = "(c) Population density", subtitle = "People per km\u00b2 (log scale)"))
s1 <- s1 + plot_annotation(caption = src_open, theme = theme_sae())
save_fig(s1, "figS1_covariates", 11, 5.5)

message("02_covariates.R finished.")
