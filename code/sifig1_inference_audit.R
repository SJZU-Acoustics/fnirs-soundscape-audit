# =============================================================================
# SI Figure S1 — "The inference audit" (178 mm, three panels)
# (a) Split-half reliability of the per-person appraisal slope across the 13
#     outcomes tried (B.4), with the prefrontal subject-bootstrap 95% CI
#     (B, 1,000 resamples) overlaid on the HbO prefrontal row.
# (b) Responder-subset selection curves (parts C/D): slope t value against
#     subset size k. 12 out-of-sample curves (grey) vs 3 in-sample curves
#     (ROI colours). The y axis is the t value, not the estimate, because the
#     stored curve-level permutation nulls are max-|t| statistics. Estimates at
#     the anchor k values are in the data lock. Per-k envelopes are not stored;
#     the curve-level in-sample p95 is drawn as the reference.
# (c) Two nulls for the prefrontal response-window comfort variance cell:
#     parametric-bootstrap p95 vs within-subject permutation p95, observed LRT.
#     Full null vectors are not stored, so the p95 values are drawn as markers
#     on a common axis.
# Source: output/analysis_23_responder_heterogeneity/ (code/mod23_responder_audit.R)
# =============================================================================

source("code/style.R")
library(dplyr); library(readr); library(tidyr); library(tibble)
library(ggtext); library(patchwork)

A23 <- file.path("output", "analysis_23_responder_heterogeneity")
b4   <- read_csv(file.path(A23, "B4_reliability_generality.csv"),  show_col_types = FALSE)
bprim<- read_csv(file.path(A23, "B_slope_reliability.csv"),        show_col_types = FALSE)
boot <- read_csv(file.path(A23, "B_reliability_bootstrap.csv"),    show_col_types = FALSE)
curv <- read_csv(file.path(A23, "CD_selection_curves.csv"),        show_col_types = FALSE)
summ <- read_csv(file.path(A23, "CD_selection_summary.csv"),       show_col_types = FALSE)
a2   <- read_csv(file.path(A23, "A2_variance_falsification.csv"),  show_col_types = FALSE)

roi_cols <- c(Prefrontal = "#0072B2", `Left temporal` = "#E69F00",
              `Right temporal` = "#009E73")
roi_lab <- c(PFC_Frontal = "Prefrontal", Left_Temporal = "Left temporal",
             Right_Temporal = "Right temporal", head = "Whole head",
             channel_10 = "Channel 10")

# ------------------------------------------------------------ panel a -------
boot_ci <- quantile(boot$boot_r, c(0.025, 0.975))
b_r     <- bprim$r_enc_pairwise_mean[bprim$roi == "PFC_Frontal"]

unit_order <- c("PFC_Frontal", "Left_Temporal", "Right_Temporal", "head")
mk_rows <- function(chrom, suffix) {
  tibble(chromophore = chrom,
         unit = c(unit_order, if (chrom == "HbO") "channel_10"),
         sfx = suffix)
}
rows <- bind_rows(mk_rows("HbO", ""), mk_rows("HbR", "​"), mk_rows("HbT", "​​"))

a <- b4 %>% mutate(sfx = recode(chromophore, HbO = "", HbR = "​", HbT = "​​")) %>%
  mutate(y_lab = paste0(roi_lab[unit], sfx))

# full y-axis structure: headers + spacer rows (bottom to top)
lev_hbo <- c("Channel 10", "Whole head", "Right temporal", "Left temporal",
             "Prefrontal")
lev_hbr <- paste0(c("Whole head", "Right temporal", "Left temporal",
                    "Prefrontal"), "​")
lev_hbt <- paste0(c("Whole head", "Right temporal", "Left temporal",
                    "Prefrontal"), "​​")
full_levels <- c(lev_hbt, "**HbT**", " ", lev_hbr, "**HbR**", "  ", lev_hbo, "**HbO**")
a <- a %>% mutate(y_key = factor(y_lab, levels = full_levels))

pa <- ggplot(a, aes(y = y_key)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55",
             linewidth = 0.45) +
  geom_linerange(data = tibble(y = "Prefrontal", lo = boot_ci[1], hi = boot_ci[2]),
                 aes(y = y, xmin = lo, xmax = hi), inherit.aes = FALSE,
                 colour = "#0072B2", linewidth = 0.8) +
  geom_point(data = tibble(y = "Prefrontal", x = b_r),
             aes(y = y, x = x), inherit.aes = FALSE,
             colour = "#0072B2", size = 2.2, shape = 18) +
  geom_point(aes(x = split_half_r), size = 1.9, stroke = 0.6) +
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(breaks = c(-0.3, -0.15, 0, 0.15, 0.3),
                     limits = c(-0.37, 0.34)) +
  labs(x = "Split-half reliability (r)", y = NULL) +
  theme_pub(base_size = 8, axis_title_size = 9) +
  theme(axis.text.y = element_markdown(size = 7.5, hjust = 1, lineheight = 0.85),
        axis.ticks.y = element_blank())

# ------------------------------------------------------------ panel b -------
null_p95_insample <- summ %>% filter(arm == "D_in_sample", roi == "PFC_Frontal") %>%
  pull(null_p95)
peak <- summ %>% filter(arm == "D_in_sample", roi == "PFC_Frontal")

b <- curv %>%
  mutate(grp = if_else(arm == "D_in_sample", roi_lab[roi], "Out-of-sample"),
         curve_id = paste(arm, split, roi, sep = "|"),
         grp = factor(grp, levels = c("Prefrontal", "Left temporal",
                                      "Right temporal", "Out-of-sample")))
grp_cols <- c(roi_cols, `Out-of-sample` = "grey60")
grp_lw   <- c(Prefrontal = 1.0, `Left temporal` = 0.55, `Right temporal` = 0.55,
              `Out-of-sample` = 0.45)

pb <- ggplot(b, aes(k, t, group = curve_id, colour = grp, linewidth = grp)) +
  geom_hline(yintercept = c(-1.97, 1.97), linetype = "dashed", colour = "grey55",
             linewidth = 0.4) +
  geom_hline(yintercept = -null_p95_insample, linetype = "dotted",
             colour = "grey30", linewidth = 0.5) +
  geom_line(alpha = 0.85, na.rm = TRUE) +
  geom_point(data = tibble(k = peak$k_at_max, t = -peak$obs_max_abs_t),
             aes(k, t), inherit.aes = FALSE, colour = "#0072B2", size = 2.2) +
  annotate("text", x = 8, y = -6.02,
           label = "peak |t| = 5.71 (k = 27)", size = 7 / .pt,
           family = "Helvetica", hjust = 0, colour = "#0072B2") +
  annotate("text", x = 68.5, y = c(1.97 + 0.75, -null_p95_insample + 0.75),
           label = c("nominal |t| = 1.97",
                     "curve-level null p95 = 6.39"),
           size = 7 / .pt, family = "Helvetica", hjust = 1, colour = "grey30") +
  scale_colour_manual(values = grp_cols, limits = force,
                      labels = c(Prefrontal = "In-sample, prefrontal",
                                 `Left temporal` = "In-sample, left temporal",
                                 `Right temporal` = "In-sample, right temporal",
                                 `Out-of-sample` = "Out-of-sample (12 curves)")) +
  scale_linewidth_manual(values = grp_lw, guide = "none") +
  scale_x_continuous(breaks = c(5, 20, 35, 50, 65), expand = expansion(add = 1.5)) +
  scale_y_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-8.8, 3.2)) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                               override.aes = list(linewidth = c(1.0, 0.55, 0.55, 0.45)))) +
  labs(x = "Subset size (subjects)", y = "Slope t value") +
  theme_pub(base_size = 8, axis_title_size = 9) +
  theme(legend.position.inside = c(0.5, 0.005),
        legend.justification.inside = c(0.5, 0),
        legend.direction = "horizontal",
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.32, "cm"),
        legend.spacing.x = unit(0.15, "cm"),
        legend.spacing.y = unit(0.01, "cm"))

# ------------------------------------------------------------ panel c -------
cell <- a2 %>% filter(window == "response", roi == "PFC_Frontal", axis == "comfort_wi")
c_nulls <- tibble(
  null = c("Parametric<br>bootstrap null", "Within-subject<br>permutation null"),
  y = c(2, 1),
  p95 = c(cell$boot_null_p95, cell$perm_p95)
)

pc <- ggplot(c_nulls, aes(y = y)) +
  geom_segment(aes(x = 0, xend = p95, yend = y), colour = "grey70",
               linewidth = 2.5, lineend = "round") +
  geom_vline(xintercept = cell$lrt_obs, colour = "#D55E00", linewidth = 0.7) +
  geom_point(aes(x = p95), size = 2.4) +
  geom_text(aes(x = p95, label = sprintf("p95 = %.2f", p95)),
            hjust = if_else(c_nulls$p95 < 10, -0.25, 1.1),
            vjust = if_else(c_nulls$p95 < 10, 0.5, -1.1),
            size = 7 / .pt, family = "Helvetica") +
  annotate("text", x = cell$lrt_obs, y = 2.55,
           label = sprintf("observed LRT = %.1f", cell$lrt_obs),
           size = 7 / .pt, family = "Helvetica", hjust = 1.05,
           colour = "#D55E00") +
  scale_y_continuous(breaks = c_nulls$y, labels = c_nulls$null,
                     limits = c(0.45, 2.75), expand = expansion(add = 0.15)) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40), limits = c(0, 43),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "LRT statistic", y = NULL) +
  theme_pub(base_size = 8, axis_title_size = 9) +
  theme(axis.text.y = element_markdown(size = 7.5, hjust = 1),
        axis.ticks.y = element_blank())

p <- pa + pb + pc + plot_layout(ncol = 3, widths = c(0.95, 1.35, 0.85)) +
     plot_annotation(tag_levels = "a") &
     theme(plot.tag = element_text(size = 10, face = "bold"),
           plot.tag.position = c(0, 1))

save_fig(p, file.path(FIGDIR, "sifig1_inference_audit.png"), 178, 85)

# ---- data locks -------------------------------------------------------------
a %>% transmute(chromophore, unit, split_half_r, n_subjects) %>%
  write_csv(file.path(LOCKDIR, "sifig1a_split_half_r.csv"))
tibble(statistic = c("bootstrap_point_r", "bootstrap_ci_lo", "bootstrap_ci_hi"),
       value = c(b_r, boot_ci[1], boot_ci[2])) %>%
  write_csv(file.path(LOCKDIR, "sifig1a_bootstrap_ci.csv"))
b %>% transmute(arm, split, roi, k, n, est, se, t) %>%
  write_csv(file.path(LOCKDIR, "sifig1b_selection_curves.csv"))
tibble(reference = c("nominal_abs_t", "curve_level_null_p95_in_sample_pfc",
                     "peak_abs_t_in_sample_pfc", "peak_at_k"),
       value = c(1.97, null_p95_insample, peak$obs_max_abs_t, peak$k_at_max)) %>%
  write_csv(file.path(LOCKDIR, "sifig1b_references.csv"))
c_nulls %>% mutate(observed_lrt = cell$lrt_obs) %>%
  write_csv(file.path(LOCKDIR, "sifig1c_two_nulls.csv"))

# anchor checks
stopifnot(abs(boot_ci[1] - (-0.272)) < 0.001, abs(boot_ci[2] - 0.232) < 0.001,
          abs(cell$lrt_obs - 37.9) < 0.1, abs(cell$boot_null_p95 - 2.20) < 0.01,
          abs(cell$perm_p95 - 31.30) < 0.01,
          abs(null_p95_insample - 6.39) < 0.01,
          abs(peak$obs_max_abs_t - 5.71) < 0.01, peak$k_at_max == 27)
cat("sifig1 done. nominal-k in-sample PFC:", peak$n_nominal_k, "of", peak$n_k,
    "; out-of-sample curves surviving BH:", sum(summ$q_bh[summ$arm == "C_out_of_sample"] < 0.05),
    "of 12\n")
