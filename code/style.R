# =============================================================================
# Shared publication plot style (house rules: knowledge/Academic_plot_style.md)
# Sourced by every figure script. Double-column figures are the default for
# multi-panel displays; PNG 600 dpi, opaque white.
# All paths are relative to the repository root.
# =============================================================================
library(ggplot2)

WIDTH_MM <- c(single_column = 85, double_column = 178)
mm2in <- function(mm) mm / 25.4

okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7")

theme_pub <- function(base_size = 9, axis_title_size = 10) {
  theme_classic(base_size = base_size, base_family = "Helvetica") %+replace%
    theme(
      axis.line         = element_line(colour = "black", linewidth = 0.5),
      axis.ticks        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.10, "cm"),
      axis.title        = element_text(size = axis_title_size),
      axis.text         = element_text(size = base_size, colour = "black"),
      legend.text       = element_text(size = base_size),
      legend.title      = element_blank(),
      legend.position   = "inside",
      legend.position.inside = c(0.98, 0.98),
      legend.justification.inside = c(1, 1),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      panel.grid        = element_blank(),
      plot.title        = element_blank(),
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA)
    )
}

# Panel label "a", "b", ... at top-left of the full panel footprint (tag).
panel_tag <- function(tag) {
  labs(tag = tag) + theme(plot.tag = element_text(size = 10, face = "bold"),
                          plot.tag.position = c(0, 1))
}

save_fig <- function(plot, path, width_mm, height_mm) {
  ggsave(path, plot, width = mm2in(width_mm), height = mm2in(height_mm),
         dpi = 600, bg = "white", device = ragg::agg_png)
}

PROJ    <- "."
FIGDIR  <- file.path("output", "figures")
LOCKDIR <- file.path("output", "data_lock")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOCKDIR, showWarnings = FALSE, recursive = TRUE)
