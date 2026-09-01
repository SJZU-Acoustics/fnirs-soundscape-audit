# =============================================================================
# Figure 1 — Measurement anatomy (178 mm, three panels)
#   a  Mean HbO across the 90-s presentation cycle, first clip of each sequence vs later clips.
#      The fig1a_* intermediates ARE the frozen product of the cycle recompute
#      (which needs the raw delivered .mat series, not deposited); they are
#      read via load_intermediate() and copied through to output/data_lock/
#      byte-identical.
#   b  Amplitude response of the original 0.01-0.10 Hz band-pass, recomputed
#      by code/mod_fig1b_filter_gain.R (0.688 gate passed) into
#      output/data_lock/fig1b_filter_gain.csv — run that module first.
#   c  Re-filter result (valid arm, motion_correction = "none"): cycle-position
#      contrast across four high-pass settings, with the response-window
#      non-detection (min BH q) below it. Shipped intermediate
#      fig1c_refilter.csv.
# Plots exactly what the lock CSVs contain.
# =============================================================================

source("code/load_data.R")
source("code/style.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

POS_COL <- c("opens a run" = okabe_ito[1], "follows a 90-s cycle" = okabe_ito[4])
POS_LAB <- c("First clip of sequence", "Later clips")
BASE <- 8   # three panels in a row: secondary text 8 pt, axis titles 9 pt
ATITLE <- 9
TAG_THEME <- theme(plot.tag = element_text(size = 10, face = "bold",
                                           family = "Helvetica"))

# The read-through intermediates are the lock CSVs: copy them through
# byte-identical so output/data_lock/ holds everything the figure plots.
copy_lock <- function(filename) {
  from <- file.path(INTERMEDIATE_DIR, filename)
  to <- file.path(LOCKDIR, filename)
  file.copy(from, to, overwrite = TRUE)
  stopifnot(tools::md5sum(from) == tools::md5sum(to))
  cat("  locked", to, "\n")
}
fig1a_locks <- c("fig1a_cycle_curve.csv", "fig1a_anchor_check.csv",
                 "fig1a_difference_score_check.csv", "fig1a_event_bins_head.csv")
invisible(lapply(c(fig1a_locks, "fig1c_refilter.csv"), copy_lock))

# ----------------------------------------------------------------- panel a
curve <- load_intermediate("fig1a_cycle_curve.csv") %>%
  mutate(cycle_position = factor(cycle_position,
                                 levels = c("opens a run", "follows a 90-s cycle")))
windows <- tibble(xmin = c(-10, 0, 5), xmax = c(0, 5, 30),
                  shade = c("grey93", "grey87", "grey93"))
p_a <- ggplot(curve, aes(bin_mid, mean, colour = cycle_position,
                         fill = cycle_position)) +
  geom_rect(data = windows,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = windows$shade, colour = NA) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se),
              colour = NA, alpha = 0.25) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = POS_COL, breaks = names(POS_COL),
                      labels = POS_LAB) +
  scale_fill_manual(values = POS_COL, breaks = names(POS_COL),
                    labels = POS_LAB) +
  scale_x_continuous(breaks = seq(0, 90, 30)) +
  coord_cartesian(xlim = c(-10, 90), ylim = c(-0.13, 0.12)) +
  labs(x = "Time from clip onset (s)", y = "HbO (µM)") +
  theme_pub(base_size = BASE, axis_title_size = ATITLE) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.02, 0.99),
        legend.justification.inside = c(0, 1),
        legend.key.spacing.y = unit(0, "mm"),
        legend.key.height = unit(3.2, "mm"),
        legend.margin = margin(0, 0, 0, 0))

# ----------------------------------------------------------------- panel b
gain <- readr::read_csv(file.path(LOCKDIR, "fig1b_filter_gain.csv"),
                        show_col_types = FALSE)
mark <- readr::read_csv(file.path(LOCKDIR, "fig1b_gain_anchor_check.csv"),
                        show_col_types = FALSE)
g90 <- mark[1, ]
p_b <- ggplot(gain, aes(freq_hz, gain)) +
  annotate("rect", xmin = 0.01, xmax = 0.10, ymin = -Inf, ymax = Inf,
           fill = "grey92") +
  geom_hline(yintercept = 1, colour = "grey60", linewidth = 0.3,
             linetype = "dashed") +
  geom_line(linewidth = 0.6, colour = okabe_ito[1]) +
  geom_point(data = g90, colour = okabe_ito[4], size = 1.8) +
  annotate("segment", x = g90$freq_hz, xend = g90$freq_hz,
           y = 0, yend = g90$gain, colour = okabe_ito[4],
           linewidth = 0.4, linetype = "dashed") +
  annotate("text", x = g90$freq_hz * 1.3, y = g90$gain + 0.06,
           label = sprintf("%.3f", g90$gain),
           size = BASE / .pt, colour = "black", hjust = 0) +
  scale_x_log10(breaks = c(0.01, 0.02, 0.05, 0.1, 0.2),
                labels = c("0.01", "0.02", "0.05", "0.1", "0.2")) +
  coord_cartesian(xlim = c(0.005, 0.3), ylim = c(0, 1.2)) +
  labs(x = "Frequency (Hz)", y = "Amplitude gain (|H|²)") +
  theme_pub(base_size = BASE, axis_title_size = ATITLE)

# ----------------------------------------------------------------- panel c
ref <- load_intermediate("fig1c_refilter.csv") %>%
  mutate(setting = factor(sprintf("%g", hpf),
                          levels = c("0.01", "0.005", "0.002", "0"),
                          labels = c("0.01\n(original)", "0.005", "0.002",
                                     "none")))
xlab_c <- "High-pass corner (Hz)"
p_c1 <- ggplot(ref, aes(setting, diff_contrast, group = 1)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(colour = okabe_ito[1], linewidth = 0.6) +
  geom_point(colour = okabe_ito[1], size = 1.8) +
  labs(x = NULL, y = "Cycle-position\ncontrast (µM)") +
  theme_pub(base_size = BASE, axis_title_size = ATITLE) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_c2 <- ggplot(ref, aes(setting, min_response_q_BH, group = 1)) +
  geom_hline(yintercept = 0.05, colour = "grey40", linewidth = 0.3,
             linetype = "dashed") +
  geom_line(colour = okabe_ito[4], linewidth = 0.6) +
  geom_point(colour = okabe_ito[4], size = 1.8) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_discrete(expand = expansion(add = 0.35)) +
  labs(x = xlab_c, y = "min BH q,\nresponse window") +
  theme_pub(base_size = BASE, axis_title_size = ATITLE)

p_c <- p_c1 / p_c2 + plot_layout(heights = c(1, 1))

# tags descend into the nested c stack (c1, c2): give the 4th an empty tag
fig1 <- p_a | p_b | p_c
fig1 <- fig1 + plot_annotation(tag_levels = list(c("a", "b", "c", ""))) &
  TAG_THEME
save_fig(fig1, file.path(FIGDIR, "fig1_anatomy.png"), width_mm = 178, height_mm = 64)
cat("wrote", file.path(FIGDIR, "fig1_anatomy.png"), "\n")
