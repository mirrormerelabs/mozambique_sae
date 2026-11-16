############################################################
## 06_web.R
## Interactive district map for GitHub Pages (docs/index.html).
##
## Layers (one shown at a time, chosen in the layer control):
##   with model outputs: estimated share in the national poorest 40%
##   (M3), coefficient of variation, P(share > 50%);
##   always: Relative Wealth Index, travel time, population density.
## Hovering over a district shows its name and every value.
## Runs before the models exist (covariate layers only), so the page
## can be tested now; rerun after 03 to add the estimate layers.
##
## Publishing: enable GitHub Pages from the /docs folder of the main
## branch. No cluster location, and no direct estimate, is included.
############################################################

source("R/functions.R")
suppressPackageStartupMessages({ library(leaflet); library(htmlwidgets); library(htmltools) })

adm2 <- read_districts()
covs <- readRDS(file.path(dir_work, "covariates_district.rds"))
est_file <- file.path(dir_out, "s03_district_estimates.csv")
have_models <- file.exists(est_file)

## Simplified outlines, coordinates rounded to 4 decimals (about 10 m),
## keep the page small. Tolerance is in degrees (about 500 m); sf warns
## that it does not simplify longitude/latitude correctly, which is
## acceptable at this scale.
geo <- adm2 %>%
  st_set_precision(1e4) %>%
  st_simplify(dTolerance = 0.005, preserveTopology = TRUE) %>%
  suppressWarnings() %>%
  left_join(select(covs, district_id, pop, rwi, tt_mean, log_dens), by = "district_id") %>%
  mutate(dens = exp(log_dens))
if (have_models) {
  m3 <- read_csv(est_file, show_col_types = FALSE) %>% filter(model == "M3") %>%
    select(district_id, province, clusters, p_mean, p_q025, p_q975, cv, pr_above_50)   # clusters: grouped
  geo <- left_join(geo, m3, by = "district_id")
}

## ---- Hover labels -----------------------------------------------------
pct <- \(x, d = 0) sprintf(paste0("%.", d, "f%%"), 100 * x)
fmt_int <- \(x) formatC(round(x), format = "d", big.mark = ",")
labels <- lapply(seq_len(nrow(geo)), function(i) {
  g <- geo[i, ]
  rows <- c(
    if (have_models) c(
      sprintf("Province: %s", tools::toTitleCase(g$province)),
      sprintf("<b>Share in national poorest 40%%: %s</b> (95%% interval %s to %s)", pct(g$p_mean), pct(g$p_q025), pct(g$p_q975)),
      sprintf("Width of the 95%% interval: %.0f percentage points", 100 * (g$p_q975 - g$p_q025)),
      sprintf("P(share &gt; 50%%): %.2f", g$pr_above_50),
      sprintf("%s: %s", if (sim_mode) "Simulated survey clusters" else "DHS survey clusters", g$clusters),
      if (g$clusters == "0") "<i>No survey clusters: estimate from covariates and neighbouring districts</i>"),
    sprintf("Relative Wealth Index: %.2f", g$rwi),
    sprintf("Travel time to a city of 50,000+: %s minutes", fmt_int(g$tt_mean)),
    sprintf("Population (2020): %s; density %s per km&sup2;", fmt_int(g$pop), fmt_int(g$dens)))
  HTML(paste0("<b>", htmlEscape(g$district), "</b><br>", paste(rows, collapse = "<br>")))
})

## ---- Layers --------------------------------------------------------------
## Covariate colour scales are capped as in Figure S1 (values beyond the
## domain take the end colour), so rural districts span the range; their
## legends show the Figure S1 breaks as swatches.
ramp_blue <- c(seq_low, seq_high); ramp_unc <- c(unc_low, unc_high)
## ramp: colours spread evenly over the domain, or a function on [0, 1]
## (div_pal, functions.R)
layers <- list(
  list(id = "rwi", name = "Relative Wealth Index", v = geo$rwi, ramp = ramp_blue, domain = c(-0.55, 0.25),
       breaks = c(-0.5, -0.25, 0, 0.25), labels = c("-0.5", "-0.25", "0", "0.25+")),
  list(id = "tt", name = "Travel time to a city (minutes)", v = log1p(geo$tt_mean), ramp = ramp_blue,
       domain = log1p(c(20, 240)), breaks = log1p(c(20, 60, 120, 240)), labels = c("20 or less", "60", "120", "240+")),
  list(id = "dens", name = "Population density (per km²)", v = geo$log_dens, ramp = ramp_blue,
       domain = log(c(5, 500)), breaks = log(c(5, 20, 100, 500)), labels = c("5 or less", "20", "100", "500+")))
## Uncertainty: width of the 95% interval, capped as in Figure 1b (the
## CV is unstable for very low shares; it stays in the CSV)
if (have_models) { geo$width <- geo$p_q975 - geo$p_q025; w_sc <- capped_scale(geo$width) }
if (have_models) layers <- c(list(
  ## Diverging scale centred on the national share, as Figure 1a
  list(id = "est", name = "Share in national poorest 40%", v = geo$p_mean, ramp = div_pal, domain = c(0, 1),
       breaks = div_breaks, labels = div_labels),
  list(id = "width", name = "Uncertainty: 95% interval width (percentage points)", v = geo$width, ramp = ramp_unc,
       domain = w_sc$limits, breaks = w_sc$breaks, labels = w_sc$labels),
  list(id = "exc", name = "P(share > 50%)", v = geo$pr_above_50, ramp = ramp_blue, domain = c(0, 1),
       fmt = labelFormat(digits = 1))), layers)

## Districts without survey clusters: outlined in every layer (drawn on
## top, not interactive, so hovering still reaches the district)
nocl <- if (have_models) geo[geo$clusters == "0", ] else geo[0, ]
m <- leaflet(geo, options = leafletOptions(minZoom = 4)) %>%
  addProviderTiles(providers$Esri.WorldGrayCanvas)  # light basemap, no API key needed
for (L in layers) {
  dom <- if (is.null(L$domain)) range(L$v, na.rm = TRUE) else L$domain
  pal <- colorNumeric(L$ramp, domain = dom, na.color = na_fill)
  m <- m %>%
    addPolygons(fillColor = pal(scales::squish(L$v, dom)), fillOpacity = 0.85, color = "white", weight = 0.6,
                label = labels, group = L$name,
                labelOptions = labelOptions(textsize = "12px", direction = "auto"),
                highlightOptions = highlightOptions(weight = 2, color = "#333333", bringToFront = TRUE))
  if (nrow(nocl) > 0)
    m <- addPolygons(m, data = nocl, fill = FALSE, color = "#1a1a1a", weight = 2, opacity = 1, group = L$name,
                     options = pathOptions(interactive = FALSE))
  leg_class <- paste("info legend", paste0("leg-", L$id))
  m <- if (is.null(L$breaks)) {
    addLegend(m, "bottomright", pal = pal, values = dom, title = L$name, labFormat = L$fmt,
              opacity = 0.9, className = leg_class)
  } else {
    addLegend(m, "bottomright", colors = pal(L$breaks), labels = L$labels, title = L$name,
              opacity = 0.9, className = leg_class)
  }
}
group_names <- vapply(layers, `[[`, "", "name")
intro <- paste0(
  "<div style='max-width:300px;font:13px/1.4 sans-serif'>",
  "<b>District poverty in Mozambique</b><br>",
  if (have_models) "Share of each district's population in the national poorest 40% (DHS wealth quintiles 1 and 2), estimated with a spatial Fay-Herriot model in R-INLA. "
  else "Covariate layers; model estimates will be added once the models have run. ",
  "Hover over a district for its values.<br>",
  "<span style='color:#555;font-size:11px'>Sources: ", if (have_models) "Mozambique DHS 2022-23 (INE and ICF, 2024); ", "Chi et al. (2022); ",
  "Meta and CIESIN (2020); Nelson et al. (2019); geoBoundaries.</span></div>")
if (sim_mode) intro <- sub("<b>District poverty in Mozambique</b>",
  "<b style='color:#b03a2e'>SIMULATED DATA, for testing only</b><br><b>District poverty in Mozambique</b>",
  sub("Mozambique DHS 2022-23 (INE and ICF, 2024); ", "simulated direct estimates; ", intro, fixed = TRUE), fixed = TRUE)
m <- m %>%
  addLayersControl(baseGroups = group_names, options = layersControlOptions(collapsed = FALSE)) %>%
  hideGroup(group_names[-1]) %>%
  addControl(HTML(intro), position = "topleft", className = "info intro")  # not "legend", or the switcher hides it
if (nrow(nocl) > 0) m <- addControl(m, HTML(paste0(
  "<div style='font:12px/1.4 sans-serif'><span style='display:inline-block;width:14px;height:10px;",
  "border:2px solid #1a1a1a;vertical-align:middle;margin-right:6px'></span>",
  "No survey clusters: estimate from covariates<br>and neighbouring districts</div>")),
  position = "bottomleft", className = "info nocl")   # shown with every layer

## Show only the legend of the selected layer; outline the estimate
## legend's swatches
legend_ids <- setNames(as.list(paste0("leg-", vapply(layers, `[[`, "", "id"))), group_names)
m <- onRender(m, "
function(el, x, data) {
  var map = this;
  function show(name) {
    el.querySelectorAll('.legend').forEach(function(l) {
      l.style.display = l.classList.contains(data[name]) ? '' : 'none';
    });
  }
  show(Object.keys(data)[0]);
  // outline the swatches of the estimate legend, so the near-white 40% swatch shows
  el.querySelectorAll('.leg-est i').forEach(function(s) { s.style.boxShadow = 'inset 0 0 0 1px #b5b3ae'; });
  map.on('baselayerchange', function(e) { show(e.name); });
}", data = legend_ids)

dir.create(dir_docs, showWarnings = FALSE, recursive = TRUE)
saveWidget(m, file.path(getwd(), dir_docs, "index.html"), selfcontained = FALSE, libdir = "lib",
           title = "District poverty in Mozambique: small area estimates")
invisible(file.create(file.path(dir_docs, ".nojekyll")))
message("Written: ", file.path(dir_docs, "index.html"), " (", if (have_models) "estimates and covariates" else "covariates only", ")")
