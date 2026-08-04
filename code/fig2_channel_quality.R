# =============================================================================
# Figure 2 — Channel quality (178 mm, four panels)
#   a  Raw intensity traces for the three failure kinds (shipped intermediates
#      fig2a_*, extracted from the raw SNIRF recordings, which are not
#      deposited): dark channel at 0.01312, live channel with cardiac
#      pulsation (10-s inset), rail-pinned wavelength at 1.000000
#   b  Per-session dropout calendar, detectors 6-10, with the two
#      between-sequence rest intervals and the 14 release onsets
#      (shipped intermediates fig2b_*)
#   c  Mirror-pair montage schematic (hand-built) — the completeness counts in
#      output/data_lock/fig2c_mirror_pairs.csv are RECOMPUTED here from
#      load_channels() + load_dead_channels(): the superseded session-CV mask,
#      against which the counts were declared (27 participants with all four
#      T-pair channels live, 1 with all four TP-pair channels)
#   d  Effective N per ROI under the union mask, RECOMPUTED here from
#      load_channel_mask() into output/data_lock/fig2d_effective_n.csv
# Panels (a)/(b) read the intermediates through; (c)/(d) plot the recomputation.
# =============================================================================

source("code/load_data.R")
source("code/style.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

BASE <- 8   # dense multi-panel figure: secondary text 8 pt, axis titles 9 pt
ATITLE <- 9
TAG_THEME <- theme(plot.tag = element_text(size = 10, face = "bold",
                                           family = "Helvetica"))
WL_COL <- c("760" = okabe_ito[1], "850" = okabe_ito[4])

# The read-through intermediates are the lock CSVs: copy them through
# byte-identical so output/data_lock/ holds everything the figure plots.
copy_lock <- function(filename) {
  from <- file.path(INTERMEDIATE_DIR, filename)
  to <- file.path(LOCKDIR, filename)
  file.copy(from, to, overwrite = TRUE)
  stopifnot(tools::md5sum(from) == tools::md5sum(to))
  cat("  locked", to, "\n")
}
fig2ab_locks <- c("fig2a_dark_trace.csv", "fig2a_live_trace.csv",
                  "fig2a_live_zoom.csv", "fig2a_rail_trace.csv",
                  "fig2b_dark_raster.csv", "fig2b_rest_intervals.csv",
                  "fig2b_release_onsets.csv")
invisible(lapply(fig2ab_locks, copy_lock))

read_trace <- function(name) {
  load_intermediate(name) %>%
    pivot_longer(-time_s, names_to = "wavelength", values_to = "intensity",
                 names_pattern = "intensity_(\\d+)nm")
}
trace_plot <- function(d, ybreaks, ylim, strip_label, show_x = FALSE,
                       show_y_title = FALSE) {
  p <- ggplot(d, aes(time_s, intensity, colour = wavelength)) +
    geom_line(linewidth = 0.35) +
    scale_colour_manual(values = WL_COL, labels = c("760 nm", "850 nm")) +
    scale_y_continuous(breaks = ybreaks, limits = ylim) +
    annotate("text", x = 8, y = Inf, label = strip_label,
             hjust = 0, vjust = 1.4, size = BASE / .pt, colour = "black") +
    labs(x = if (show_x) "Time (s)" else NULL,
         y = if (show_y_title) "Intensity (a.u.)" else NULL) +
    theme_pub(base_size = BASE, axis_title_size = ATITLE) +
    theme(legend.position = "none")
  if (!show_x) p <- p + theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank())
  p
}

# ----------------------------------------------------------------- panel a
dark <- read_trace("fig2a_dark_trace.csv")
live <- read_trace("fig2a_live_trace.csv")
zoom <- read_trace("fig2a_live_zoom.csv")
rail <- read_trace("fig2a_rail_trace.csv")

s1 <- trace_plot(dark, ybreaks = c(0.01310, 0.01312, 0.01314),
                 ylim = c(0.01309, 0.01314), "dark") +
  annotate("text", x = 1250, y = 0.0131320, label = "760 nm", hjust = 1,
           vjust = 0, size = (BASE - 1) / .pt, colour = WL_COL["760"]) +
  annotate("text", x = 1250, y = 0.0131180, label = "850 nm", hjust = 1,
           vjust = 1, size = (BASE - 1) / .pt, colour = WL_COL["850"])

zoom_p <- ggplot(zoom %>% filter(wavelength == "850"),
                 aes(time_s, intensity)) +
  geom_line(linewidth = 0.4, colour = WL_COL["850"]) +
  scale_x_continuous(breaks = c(325, 330, 335)) +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 6.5, axis_title_size = 6.5) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA,
                                    linewidth = 0.4),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0.5, 1, 0.5, 1),
        plot.background = element_rect(fill = "white", colour = NA))
s2 <- trace_plot(live, ybreaks = c(0.4, 0.7, 1.0), ylim = c(0.3, 1.0),
                 "live", show_y_title = TRUE) +
  inset_element(zoom_p, left = 0.38, bottom = 0.14, right = 0.995,
                top = 0.985)

s3 <- trace_plot(rail, ybreaks = c(0.6, 0.8, 1.0), ylim = c(0.5, 1.02),
                 "rail", show_x = TRUE)

p_a <- s1 / s2 / s3 + plot_layout(heights = c(1, 1, 1))

# ----------------------------------------------------------------- panel b
raster <- load_intermediate("fig2b_dark_raster.csv")
rest <- load_intermediate("fig2b_rest_intervals.csv")
releases <- load_intermediate("fig2b_release_onsets.csv")

sess <- raster %>% distinct(subject_id, date) %>%
  arrange(date, subject_id) %>% mutate(sess_idx = row_number())
raster <- raster %>% left_join(sess, by = c("subject_id", "date")) %>%
  mutate(y = -sess_idx + (detector - 8) * 0.18)
releases <- releases %>% left_join(sess, by = "subject_id") %>%
  mutate(y = -sess_idx + (detector - 8) * 0.18)

# date blocks for the y axis
blocks <- sess %>% group_by(date) %>%
  summarise(mid = -mean(sess_idx), n = n(), .groups = "drop")
rest_med <- rest %>% summarise(r1s = median(rest1_start_s),
                               r1e = median(rest1_end_s),
                               r2s = median(rest2_start_s),
                               r2e = median(rest2_end_s))
DET_COL <- c("6" = okabe_ito[1], "7" = okabe_ito[2], "8" = okabe_ito[4],
             "9" = okabe_ito[3], "10" = okabe_ito[6])

p_b <- ggplot(raster %>% filter(frac_dark >= 0.5)) +
  annotate("rect", xmin = rest_med$r1s, xmax = rest_med$r1e,
           ymin = -69.6, ymax = -0.4, fill = "grey88") +
  annotate("rect", xmin = rest_med$r2s, xmax = rest_med$r2e,
           ymin = -69.6, ymax = -0.4, fill = "grey88") +
  geom_tile(aes(x = bin_start_s + 2.5, y = y, fill = factor(detector)),
            height = 0.17, width = 5) +
  geom_point(data = releases, aes(x = onset_s, y = y),
             shape = 4, size = 1.1, stroke = 0.7, colour = "black") +
  scale_fill_manual(values = DET_COL, name = "Detector") +
  scale_x_continuous(breaks = c(0, 300, 600, 900, 1200),
                     expand = expansion(mult = c(0.005, 0.005))) +
  scale_y_continuous(breaks = blocks$mid, labels = blocks$date,
                     expand = expansion(mult = c(0.005, 0.005))) +
  labs(x = "Session time (s)", y = NULL) +
  theme_pub(base_size = BASE, axis_title_size = ATITLE) +
  theme(legend.position = "top",
        legend.justification = "right",
        legend.direction = "horizontal",
        legend.title = element_text(size = BASE),
        legend.key.size = unit(0.32, "cm"),
        axis.text.y = element_text(size = 6.5))

# ----------------------------------------------------------------- panel c
# RECOMPUTED: mirror-pair completeness under the superseded session-CV mask
# (load_dead_channels()). A channel is live for a subject when its pair is
# active in the channels sheet and the cell is not in the CV dead set;
# subjects_pair_complete counts participants with all four channels of the
# pair's mirror group live.
chs <- load_channels()
dead <- load_dead_channels()
dead_set <- dead %>% distinct(subject_id, channel_index)
subjects <- unique(chs$subject_id)

mirror <- tribble(
  ~pair, ~left_ch, ~right_ch,
  "T",  16L, 20L,
  "T",  17L, 21L,
  "TP", 18L, 22L,
  "TP", 19L, 23L)
ch_info <- chs %>% distinct(channel_index, channel_name, detector)
mirror <- mirror %>%
  left_join(ch_info, by = c("left_ch" = "channel_index")) %>%
  rename(left_name = channel_name, left_detector = detector) %>%
  left_join(ch_info, by = c("right_ch" = "channel_index")) %>%
  rename(right_name = channel_name, right_detector = detector)

cell_live <- function(s, ch) {
  act <- chs %>% filter(subject_id == s, channel_index == ch) %>%
    pull(active_pair)
  length(act) == 1 && act == 1 &&
    !any(dead_set$subject_id == s & dead_set$channel_index == ch)
}
cell_dead <- function(s, ch) {
  any(dead_set$subject_id == s & dead_set$channel_index == ch)
}
complete_pair <- map_dfr(c("T", "TP"), function(pp) {
  ch4 <- mirror %>% filter(pair == pp)
  n_ok <- sum(map_lgl(subjects, function(s) {
    all(map_lgl(c(ch4$left_ch, ch4$right_ch), function(ch) cell_live(s, ch)))
  }))
  tibble(pair = pp, subjects_pair_complete = n_ok)
})
mirror <- mirror %>% left_join(complete_pair, by = "pair") %>%
  rowwise() %>%
  mutate(failed_side = if (sum(map_lgl(subjects, cell_dead, ch = right_ch)) >
                             sum(map_lgl(subjects, cell_dead, ch = left_ch))) {
    "right"
  } else {
    "left"
  }) %>%
  ungroup() %>%
  select(pair, left_ch, left_name, left_detector, right_ch, right_name,
         right_detector, failed_side, subjects_pair_complete)

# Gate: the declared CV-mask-era counts.
stopifnot("Fig2c gate failed: mirror-pair completeness is not T=27 / TP=1" =
            identical(mirror$subjects_pair_complete, c(27L, 27L, 1L, 1L)))
readr::write_csv(mirror, file.path(LOCKDIR, "fig2c_mirror_pairs.csv"), eol = "\r\n")
cat("  wrote", file.path(LOCKDIR, "fig2c_mirror_pairs.csv"), "\n")

# Hand-built symmetric temporal montage (verified against the channels sheet:
# ch 16/17 det 7, 18/19 det 9, 20/21 det 8, 22/23 det 10; mirror pairs
# 16/20, 17/21, 18/22, 19/23 — the recomputation above).
det_pos <- tibble(detector = c(7, 9, 8, 10),
                  site = c("T7", "TP7", "T8", "TP8"),
                  x = c(-1, -1, 1, 1), y = c(0.42, -0.42, 0.42, -0.42),
                  failed = c(FALSE, TRUE, TRUE, FALSE))
src_pos <- tibble(source = c(7, 9, 11, 8, 10, 12),
                  x = c(-1.55, -1.55, -1.55, 1.55, 1.55, 1.55),
                  y = c(0.75, 0.05, -0.75, 0.75, 0.05, -0.75))
ch_pos <- tibble(
  channel = 16:23,
  x = c(-1.27, -1.27, -1.27, -1.27, 1.27, 1.27, 1.27, 1.27),
  y = c(0.62, 0.22, -0.22, -0.62, 0.62, 0.22, -0.22, -0.62),
  detector = c(7, 7, 9, 9, 8, 8, 10, 10),
  pair = c("T", "T", "TP", "TP", "T", "T", "TP", "TP")) %>%
  # labels: outer pair above/below the dots, inner pair outside the dots
  mutate(lab_x = x + ifelse(abs(y) < 0.4, sign(x) * 0.22, 0),
         lab_y = y + ifelse(abs(y) > 0.4, sign(y) * 0.19, 0))
pairs <- tibble(x = c(-1.27, -1.27, -1.27, -1.27), xend = 1.27,
                y = c(0.62, 0.22, -0.22, -0.62),
                yend = c(0.62, 0.22, -0.22, -0.62))
sd_links <- tibble(
  x = c(rep(-1.55, 4), rep(1.55, 4)),
  y = c(0.75, 0.05, 0.05, -0.75, 0.75, 0.05, 0.05, -0.75),
  xend = c(rep(-1, 4), rep(1, 4)),
  yend = c(0.42, 0.42, -0.42, -0.42, 0.42, 0.42, -0.42, -0.42))
head_outline <- tibble(t = seq(0, 2 * pi, length.out = 200),
                       x = 1.88 * cos(t), y = 1.88 * sin(t))

p_c <- ggplot() +
  geom_path(data = head_outline, aes(x, y), colour = "grey75",
            linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = 0, y = -1.6, yend = 1.6),
               colour = "grey85", linewidth = 0.4, linetype = "dashed") +
  geom_segment(data = pairs, aes(x = x, xend = xend, y = y, yend = yend),
               colour = "grey70", linewidth = 0.4, linetype = "dashed") +
  geom_segment(data = sd_links, aes(x = x, xend = xend, y = y, yend = yend),
               colour = "grey80", linewidth = 0.4) +
  geom_point(data = src_pos, aes(x, y), shape = 17, size = 2.2,
             colour = "grey45") +
  geom_point(data = det_pos, aes(x, y), shape = 22, size = 3.2, stroke = 1,
             fill = ifelse(det_pos$failed, okabe_ito[4], "white"),
             colour = "black") +
  geom_point(data = ch_pos %>% filter(detector %in% c(7, 10)),
             aes(x, y), shape = 16, size = 2, colour = okabe_ito[1]) +
  geom_point(data = ch_pos %>% filter(detector %in% c(8, 9)),
             aes(x, y), shape = 16, size = 2, colour = "grey55") +
  geom_text(data = ch_pos, aes(lab_x, lab_y, label = channel),
            size = (BASE - 1) / .pt, colour = "black") +
  geom_text(data = det_pos,
            aes(x + ifelse(x < 0, 0.21, -0.21), y,
                label = paste0("d", detector)),
            size = (BASE - 1) / .pt, colour = "grey25",
            hjust = ifelse(det_pos$x < 0, 0, 1)) +
  annotate("text", x = -1.55, y = 1.15, label = "left", size = BASE / .pt,
           colour = "grey30") +
  annotate("text", x = 1.55, y = 1.15, label = "right", size = BASE / .pt,
           colour = "grey30") +
  annotate("text", x = 0, y = -2.02,
           label = "▲ source   ● channel   ■ detector",
           size = (BASE - 0.5) / .pt, colour = "grey30") +
  coord_equal(xlim = c(-2.2, 2.2), ylim = c(-2.2, 2.05), clip = "off") +
  theme_void(base_size = BASE, base_family = "Helvetica") +
  theme(plot.margin = margin(2, 2, 2, 2))

# ----------------------------------------------------------------- panel d
# RECOMPUTED: effective N per ROI under the union mask (load_channel_mask()).
# n_full_complement counts participants whose surviving montage is complete —
# all 13 prefrontal channels, or both channels of the temporal ROI's surviving
# detector pair (detector 7 left / detector 10 right; detectors 8 and 9 fail
# montage-wide) — i.e. n_live >= the complement below.
mask <- load_channel_mask()
live <- mask %>% filter(active_pair == 1, use_in_aggregation == 1) %>%
  dplyr::count(subject_id, roi, name = "n_live") %>%
  complete(subject_id, roi = ROIS, fill = list(n_live = 0L))
FULL_COMPLEMENT <- c(PFC_Frontal = 13L, Left_Temporal = 2L, Right_Temporal = 2L)
effn <- live %>% group_by(roi) %>%
  summarise(n_full_complement = sum(n_live >= FULL_COMPLEMENT[roi]),
            n_any_usable      = sum(n_live >= 1),
            n_none            = sum(n_live == 0),
            n_total           = n(), .groups = "drop") %>%
  mutate(roi = factor(roi, levels = ROIS)) %>% arrange(roi) %>%
  mutate(roi = as.character(roi))

# Gate: the declared union-mask counts.
stopifnot("Fig2d gate failed: effective-N counts differ from the declared ones" =
            identical(effn$n_full_complement, c(61L, 65L, 57L)) &&
            identical(effn$n_any_usable, c(69L, 67L, 58L)) &&
            identical(effn$n_none, c(0L, 2L, 11L)) &&
            identical(effn$n_total, c(69L, 69L, 69L)))
readr::write_csv(effn, file.path(LOCKDIR, "fig2d_effective_n.csv"), eol = "\r\n")
cat("  wrote", file.path(LOCKDIR, "fig2d_effective_n.csv"), "\n")

effn <- effn %>%
  mutate(roi_label = factor(roi, levels = c("Right_Temporal", "Left_Temporal",
                                            "PFC_Frontal"),
                            labels = c("Right temporal", "Left temporal",
                                       "Prefrontal")),
         n_usable = if_else(roi == "PFC_Frontal", n_full_complement,
                            n_any_usable))
p_d <- ggplot(effn, aes(roi_label, n_usable)) +
  geom_col(fill = okabe_ito[1], width = 0.55) +
  geom_text(aes(label = sprintf("%d/69", n_usable)), hjust = -0.15,
            size = BASE / .pt, colour = "black") +
  scale_y_continuous(breaks = c(0, 20, 40, 60),
                     expand = expansion(mult = c(0, 0.16))) +
  coord_flip(clip = "off") +
  labs(x = NULL, y = "Participants (n)") +
  theme_pub(base_size = BASE, axis_title_size = ATITLE) +
  theme(axis.text.y = element_text(size = BASE))

# ----------------------------------------------------------------- assembly
row1 <- p_a | p_c | p_d + plot_layout(widths = c(1.1, 0.85, 1.05))
fig2 <- row1 / p_b + plot_layout(heights = c(1, 1.1))
# tag order: s1, s2, zoom inset, s3 (nested a stack), montage, bars, calendar
fig2 <- fig2 + plot_annotation(tag_levels = list(c("a", "", "", "", "c", "d",
                                                   "b"))) & TAG_THEME
save_fig(fig2, file.path(FIGDIR, "fig2_channel_quality.png"),
         width_mm = 178, height_mm = 150)
cat("wrote", file.path(FIGDIR, "fig2_channel_quality.png"), "\n")
