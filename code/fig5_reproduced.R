# =============================================================================
# Figure 5 — "What did reproduce" (178 mm, two panels)
# (a) Inter-subject time-course reproducibility (A18): mean between-subject
#     correlation of the 0-60 s course for same audio / same position-different
#     audio / neither, for prefrontal and whole head; contrast arrows with the
#     two permutation p values; inset = recomputed permutation null (1,000
#     perms) of the prefrontal position-minus-neither contrast with the
#     observed value marked. Null recomputed by
#     code/recompute_a18_position_null.R (exact RNG-stream replica of the A18
#     module; recomputed p == stored p for all four units).
# (b) Prefrontal encounter trends (A12): per-encounter baseline and response
#     window means. No per-point SE is stored in the A12 outcome files, so none
#     is plotted. Trend slopes (z-scaled encounter) are caption material.
# Source: output/analysis_18_clip_reliability_and_adjectives/neural_clip_reliability.csv
#         output/analysis_12_repetition_structure/run_means.csv
#         output/data_lock/fig5a_position_contrast_null.csv
# =============================================================================

source("code/style.R")
library(dplyr); library(readr); library(tidyr); library(tibble)
library(patchwork)

A18 <- file.path(PROJ, "output/analysis_18_clip_reliability_and_adjectives")
A12 <- file.path(PROJ, "output/analysis_12_repetition_structure")
isc  <- read_csv(file.path(A18, "neural_clip_reliability.csv"), show_col_types = FALSE)
runm <- read_csv(file.path(A12, "run_means.csv"),               show_col_types = FALSE)
null <- read_csv(file.path(LOCKDIR, "fig5a_position_contrast_null.csv"),
                 show_col_types = FALSE)

# ------------------------------------------------------------ panel a -------
# Points carry unit by shape only (black); pair type is the x axis, so no
# colour legend is needed. Grey arrows trace the two contrasts per unit.
a <- isc %>% filter(unit %in% c("Prefrontal", "Whole head")) %>%
  transmute(unit, n_courses,
            `Same audio` = r_same_wave,
            `Same position,\ndifferent audio` = r_pos_ctrl,
            `Neither` = r_diff,
            p_pos = p_position_vs_diff, p_audio = p_audio_vs_position) %>%
  pivot_longer(cols = c(`Same audio`, `Same position,\ndifferent audio`, `Neither`),
               names_to = "pair_type", values_to = "r") %>%
  mutate(pair_type = factor(pair_type,
                            levels = c("Same audio", "Same position,\ndifferent audio",
                                       "Neither")),
         x = as.numeric(pair_type) + if_else(unit == "Prefrontal", -0.13, 0.13))

arrow_df <- tribble(
  ~unit,        ~x1,   ~x2,  ~contrast,
  "Prefrontal", 2 - .13, 1 - .13, "audio",
  "Prefrontal", 3 - .13, 2 - .13, "position",
  "Whole head", 2 + .13, 1 + .13, "audio",
  "Whole head", 3 + .13, 2 + .13, "position"
) %>% left_join(a %>% select(unit, x, r), by = c("unit", "x1" = "x")) %>%
  rename(y1 = r) %>%
  left_join(a %>% select(unit, x, r), by = c("unit", "x2" = "x")) %>%
  rename(y2 = r)

# inset: recomputed permutation null of the prefrontal position contrast
pfc_obs  <- isc %>% filter(unit == "Prefrontal") %>%
  summarise(obs = position_vs_diff, p = p_position_vs_diff)
null_pfc <- null %>% filter(unit == "Prefrontal")
inset <- ggplot(null_pfc, aes(null_position_vs_diff)) +
  geom_histogram(binwidth = 0.002, fill = "grey75", colour = NA) +
  geom_vline(xintercept = pfc_obs$obs, colour = "#D55E00", linewidth = 0.5) +
  scale_x_continuous(breaks = c(-0.01, 0, 0.01)) +
  labs(x = "position-contrast null (r)", y = NULL) +
  theme_void(base_family = "Helvetica") +
  theme(axis.text.x = element_text(size = 6.5, colour = "black",
                                   margin = margin(t = 1)),
        axis.title.x = element_text(size = 6.5, margin = margin(t = 2)),
        axis.ticks.x = element_line(colour = "black", linewidth = 0.3),
        axis.ticks.length.x = unit(0.05, "cm"),
        axis.line.x = element_line(colour = "black", linewidth = 0.4))

pa <- ggplot(a, aes(x, r, shape = unit)) +
  geom_segment(data = arrow_df, aes(x = x1, y = y1, xend = x2, yend = y2),
               inherit.aes = FALSE, colour = "grey40", linewidth = 0.45,
               arrow = arrow(length = unit(0.13, "cm"), type = "open")) +
  geom_point(size = 2.3, stroke = 0.7) +
  geom_text(aes(label = sprintf("%.3f", r),
                hjust = if_else(unit == "Prefrontal", 1.15, -0.15)),
            vjust = 0.5, size = 7 / .pt, family = "Helvetica",
            show.legend = FALSE) +
  annotate("text", x = 0.52, y = c(0.2465, 0.2395),
           label = c("position − neither: p = 0.001",
                     "audio − position: p ≥ 0.85"),
           size = 7 / .pt, family = "Helvetica", hjust = 0) +
  scale_shape_manual(values = c(Prefrontal = 16, `Whole head` = 17)) +
  scale_x_continuous(breaks = 1:3,
                     labels = levels(a$pair_type),
                     limits = c(0.45, 3.55),
                     expand = expansion(add = 0.15)) +
  scale_y_continuous(limits = c(0.183, 0.251), breaks = c(0.19, 0.21, 0.23)) +
  labs(x = NULL, y = "Mean between-subject\ncorrelation (r)") +
  theme_pub() +
  theme(legend.position.inside = c(0.99, 0.99),
        legend.justification = c(1, 1),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.3, "cm")) +
  inset_element(inset, left = 0.05, bottom = 0.03, right = 0.44, top = 0.26,
                align_to = "panel")

# ------------------------------------------------------------ panel b -------
b <- runm %>%
  transmute(encounter,
            Baseline = baseline_PFC_Frontal,
            Response = response_PFC_Frontal) %>%
  pivot_longer(-encounter, names_to = "window", values_to = "mean_uM")

win_cols <- c(Baseline = "#0072B2", Response = "#D55E00")

pb <- ggplot(b, aes(encounter, mean_uM, colour = window, group = window)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_line(linewidth = 0.7, show.legend = FALSE) +
  geom_point(size = 2.2, stroke = 0.7, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%+.3f", mean_uM)),
            size = 7 / .pt, family = "Helvetica", show.legend = FALSE,
            vjust = if_else(b$window == "Baseline", -1.0, 1.8)) +
  annotate("text", x = 3.28, y = b$mean_uM[b$window == "Baseline" & b$encounter == 3],
           label = "Baseline", colour = win_cols["Baseline"],
           size = 8 / .pt, family = "Helvetica", hjust = 0) +
  annotate("text", x = 3.28, y = b$mean_uM[b$window == "Response" & b$encounter == 3],
           label = "Response", colour = win_cols["Response"],
           size = 8 / .pt, family = "Helvetica", hjust = 0) +
  scale_colour_manual(values = win_cols) +
  scale_x_continuous(breaks = 1:3, limits = c(0.8, 4.15),
                     expand = expansion(add = 0.1)) +
  scale_y_continuous(limits = c(-0.058, 0.118), breaks = c(-0.04, 0, 0.04, 0.08)) +
  labs(x = "Encounter (1–3)", y = "Prefrontal window\nmean (µM)") +
  theme_pub()

p <- pa + pb + plot_layout(ncol = 2, widths = c(1.25, 1)) +
     plot_annotation(tag_levels = list(c("a", "", "b"))) &
     theme(plot.tag = element_text(size = 10, face = "bold"),
           plot.tag.position = c(0, 1))

save_fig(p, file.path(FIGDIR, "fig5_reproduced.png"), 178, 68)

# ---- data locks -------------------------------------------------------------
a %>% select(unit, pair_type, r, n_courses, p_pos, p_audio) %>%
  write_csv(file.path(LOCKDIR, "fig5a_clip_reliability.csv"))
b %>% write_csv(file.path(LOCKDIR, "fig5b_encounter_means.csv"))
cat("fig5 done.\n")
