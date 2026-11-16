############################################################
## 01_direct.R
## Direct (design-based) estimates of the district poverty share from
## the Mozambique DHS 2022-23.
##
## Target: share of the population living in households in the two
## poorest national wealth quintiles (hv270 <= 2), i.e. the national
## poorest 40%. Person-weighted: hv005 x hv012 (de jure members), so
## the national estimate should be close to 0.40 by construction
## (DHS quintiles are defined on the de jure household population).
##
## Districts: geoBoundaries ADM2 (157 after dropping two Malawian islands). The GPS README states that
## displacement is restricted to admin2, so a displaced cluster lies in
## the same district as the true one; boundary differences between the
## DHS and geoBoundaries layers are checked and counted below.
##
## DATA PROTECTION: coordinates stay in memory; household, cluster and
## direct-estimate files go to data_work/ only (never published).
############################################################

source("R/functions.R")
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu = "adjust")

## ---- 1. Household recode (design and wealth variables only) ---------
hr <- read_dta(file_hr, col_select = c(hv001, hv002, hv005, hv012, hv022, hv024, hv025, hv270))
stopifnot(nrow(hr) == 14250, !anyDuplicated(hr[c("hv001", "hv002")]), !anyNA(hr))

hh <- hr %>%
  transmute(cluster = as.integer(hv001), household = as.integer(hv002),
            stratum = as.integer(hv022),
            province_dhs = norm_province(as.character(as_factor(hv024))),
            residence = as.character(as_factor(hv025)),
            members = as.integer(hv012),
            w_person = hv005 / 1e6 * hv012,
            poor = as.integer(hv270 <= 2))
n_zero_members <- sum(hh$members == 0)
hh <- filter(hh, members > 0)   # no de jure members: no person weight

## ---- 2. Clusters to districts ----------------------------------------
adm2 <- read_districts()
ge <- st_read(file_ge, quiet = TRUE) %>%
  select(cluster = DHSCLUST, LATNUM, LONGNUM) %>%
  mutate(cluster = as.integer(cluster)) %>%
  filter(!(LATNUM == 0 & LONGNUM == 0))
stopifnot(!anyDuplicated(ge$cluster))

within_idx <- st_within(ge, adm2)
n_hits <- lengths(within_idx)
ge$district_id <- NA_integer_
ge$district_id[n_hits == 1] <- adm2$district_id[unlist(within_idx[n_hits == 1])]
## Points outside every polygon (layer mismatch at borders or coast):
## nearest district. Counted, not dropped.
outside <- which(n_hits == 0)
if (length(outside) > 0) ge$district_id[outside] <- adm2$district_id[st_nearest_feature(ge[outside, ], adm2)]
stopifnot(all(n_hits <= 1), !anyNA(ge$district_id))

cluster_district <- st_drop_geometry(ge) %>% select(cluster, district_id)
hh <- hh %>% inner_join(cluster_district, by = "cluster")

## Province of each district, for the province benchmark in 05. GADM
## 4.1 level 1 places Cidade da Matola (Maputo province) inside
## "MaputoCity" (97% of its area), so GADM alone misassigns it. For
## districts with clusters, the province is therefore the DHS province
## (hv024, the survey's own domains) of most of their clusters; GADM
## (point on the district's surface) is used only for districts without
## clusters. Agreement between hv024 and GADM is counted as a check.
adm1_tmp <- tempfile("adm1_"); dir.create(adm1_tmp)
unzip(file_adm1_zip, exdir = adm1_tmp)
adm1 <- st_read(list.files(adm1_tmp, pattern = "json$", full.names = TRUE), quiet = TRUE) %>%
  st_make_valid() %>% transmute(province = norm_province(NAME_1))
unlink(adm1_tmp, recursive = TRUE)
province_gadm <- st_join(st_point_on_surface(adm2), adm1, join = st_within) %>%
  st_drop_geometry() %>% select(district_id, province_gadm = province)
stopifnot(nrow(province_gadm) == n_districts, !anyNA(province_gadm$province_gadm))
stopifnot("Province names do not match between DHS and GADM; adjust norm_province()" =
            all(unique(hh$province_dhs) %in% province_gadm$province_gadm))
province_dhs_major <- hh %>% distinct(cluster, district_id, province_dhs) %>%
  count(district_id, province_dhs) %>%
  group_by(district_id) %>% slice_max(n, n = 1, with_ties = FALSE) %>% ungroup() %>%
  select(district_id, province_dhs_major = province_dhs)
district_province <- province_gadm %>%
  left_join(province_dhs_major, by = "district_id") %>%
  transmute(district_id, province = coalesce(province_dhs_major, province_gadm),
            from_dhs_differs = !is.na(province_dhs_major) & province_dhs_major != province_gadm)

hh <- hh %>% left_join(district_province %>% select(district_id, province), by = "district_id") %>%
  left_join(province_gadm, by = "district_id")
province_agree <- hh %>% distinct(cluster, province_dhs, province_gadm, province) %>%
  summarise(clusters = n(), agree_gadm = sum(province_dhs == province_gadm), agree_district = sum(province_dhs == province))

## ---- 3. Survey design and national check ------------------------------
des <- svydesign(ids = ~cluster, strata = ~stratum, weights = ~w_person, data = hh, nest = TRUE)
national <- svymean(~poor, des)
stopifnot(abs(coef(national) - 0.40) < 0.03)   # should be about 0.40

## ---- 4. Direct estimates by district and by province ------------------
dir_d <- svyby(~poor, ~district_id, des, svymean, vartype = "se", keep.var = TRUE)
counts <- hh %>% group_by(district_id) %>%
  summarise(n_clusters = n_distinct(cluster), n_households = n(), .groups = "drop")
direct <- adm2 %>% st_drop_geometry() %>%
  left_join(select(district_province, district_id, province), by = "district_id") %>%
  ## pmin/pmax: floating-point error can give a share of 1 + 2e-16
  left_join(tibble(district_id = dir_d$district_id, p_direct = pmin(pmax(dir_d$poor, 0), 1), se_direct = dir_d$se),
            by = "district_id") %>%
  left_join(counts, by = "district_id") %>%
  mutate(n_clusters = coalesce(n_clusters, 0L), n_households = coalesce(n_households, 0L),
         ## Used in the model likelihood only if at least two clusters and
         ## a direct estimate strictly between 0 and 1 with positive SE.
         ## Other districts are predicted (out of sample).
         in_model = n_clusters >= 2 & p_direct > 0 & p_direct < 1 & se_direct > 0,
         in_model = coalesce(in_model, FALSE),
         ## Effective sample size (for the arcsin Fay-Herriot in 04)
         n_eff = ifelse(in_model, p_direct * (1 - p_direct) / se_direct^2, NA_real_))

## By DHS province (hv024): the survey's own estimation domains
dir_p <- svyby(~poor, ~province_dhs, des, svymean, vartype = c("se", "ci"))
province_direct <- tibble(province = dir_p$province_dhs, p_direct = dir_p$poor, se_direct = dir_p$se,
                          lo = dir_p$ci_l, hi = dir_p$ci_u)

## ---- 5. Write ----------------------------------------------------------
saveRDS(hh, file.path(dir_work, "hh.rds"))                  # household level: never published
saveRDS(direct, file.path(dir_work, "direct_district.rds"))  # small-n districts: never published
saveRDS(province_direct, file.path(dir_work, "direct_province.rds"))

summary_01 <- tibble(
  item = c("Households read", "Households with no de jure members (dropped)", "Households analysed",
           "Clusters with coordinates", "Clusters outside every district polygon (nearest district used)",
           "Clusters whose DHS province matches the district's GADM province",
           "Clusters whose DHS province matches the district's assigned province",
           "Districts whose assigned (DHS) province differs from GADM",
           "National estimate, share in poorest 40% (person-weighted)", "National SE",
           "Districts", "Districts with any cluster", "Districts with 2+ clusters",
           "Districts in model likelihood", "Median clusters per sampled district"),
  value = c(nrow(hr), n_zero_members, nrow(hh), nrow(ge), length(outside),
            province_agree$agree_gadm, province_agree$agree_district, sum(district_province$from_dhs_differs),
            coef(national), SE(national),
            nrow(direct), sum(direct$n_clusters > 0), sum(direct$n_clusters >= 2),
            sum(direct$in_model), median(direct$n_clusters[direct$n_clusters > 0])))
summary_01$value <- round(summary_01$value, 4)
write_out(summary_01, "s01_direct_summary")
write_out(province_direct %>% mutate(across(where(is.numeric), \(x) round(x, 4))), "s01_province_direct")
message("01_direct.R finished.")
