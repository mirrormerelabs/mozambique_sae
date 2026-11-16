# Data provenance

## DHS data

- **Survey:** Mozambique Demographic and Health Survey (DHS) 2022-23 (INE and ICF, 2024).
- **Files used:** household recode MZHR82 (`MZHR82FL.dta`) and GPS cluster file MZGE81 (`MZGE81FL.shp`).
- **Access:** a new DHS Program project registration for this small area estimation analysis was approved on 16 November 2026. No DHS-dependent script was run before that date. The pipeline was developed and tested on simulated data only (`run_sim.R`, outputs in `sim/`).
- **Location:** the files are read in place from the Paper A project folder (`MOZSAE_PAPER_A`, by default `../Thesis Corrections/Ch 4 reanalysis/analysis/data_raw/external/dhs/`; see `R/functions.R`). They are never copied into this project.

## Use and protection

- Cluster coordinates are used only to assign each survey cluster to its district (`R/01_direct.R`). They are held in memory and never written to disk.
- No household, cluster or coordinate data are published.
- Household-level data and direct district estimates are kept in `data_work/` only. That folder is excluded from version control and is never published.
- `output/` and `docs/` hold district-level model estimates, province-level direct estimates, covariates from open data, and figures. Direct district estimates appear only as unnamed points in Figure 2, grouped by cluster count; Figure 5 (direct against model estimates) is kept in `data_work/figures/` and not published. Sensitivity runs (`R/sens_real.R`) write only to `data_work/sensitivity/`. Cluster counts per district are published only in groups (0, 1, 2-4, 5+).

## Other inputs (public data)

Relative Wealth Index (Chi et al., 2022), Meta and CIESIN (2020) population, travel time to cities (Nelson et al., 2019), geoBoundaries ADM2 and GADM 4.1 level 1 boundaries. These are also read in place from the Paper A folder.
