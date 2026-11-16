############################################################
## figures.R: shared figure helpers. Sourced after functions.R by
## scripts that draw figures (02, 05, 06).
############################################################

src_open <- "Sources: Chi et al. (2022); Meta and CIESIN (2020); Nelson et al. (2019); geoBoundaries ADM2."
src_dhs  <- if (sim_mode) "SIMULATED DATA, for testing only (not DHS)" else "Mozambique DHS 2022-23 (INE and ICF, 2024)"
blue <- "#2a78d6"; blue_light <- "#9ec5f4"; grey_mark <- "#8c8c8c"
map_theme <- function() theme_sae() + ggplot2::theme(legend.key.height = ggplot2::unit(0.9, "cm"), legend.key.width = ggplot2::unit(0.3, "cm"))
plot_theme <- function() ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
        plot.title = ggplot2::element_text(size = 10, face = "bold"),
        plot.caption = ggplot2::element_text(size = 7, colour = "grey35", hjust = 0))
## dir: output/figures/ (published) by default; data_work/figures/ for
## figures that must not be published
save_fig <- function(p, name, w, h, dir = dir_fig) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
  message("Written: ", file.path(dir, name), ".png/.pdf")
}
seq_fill <- function(...) ggplot2::scale_fill_gradient(low = seq_low, high = seq_high, na.value = na_fill, ...)
seq_fill_unc <- function(...) ggplot2::scale_fill_gradient(low = unc_low, high = unc_high, na.value = na_fill, ...)
district_map <- function(df, var) {
  ggplot2::ggplot(df, ggplot2::aes(fill = .data[[var]])) + ggplot2::geom_sf(colour = "white", linewidth = 0.1) + map_theme()
}

