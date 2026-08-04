# =============================================================================
# Figure 4 — "The confirmatory family"
# 2 appraisal axes x 3 ROIs = 6 within-subject slopes per analysis window;
# one panel per window (a difference / b response / c baseline), 178 mm.
# Forest style: estimate with analytic 95% CI (estimate +/- t(df)*se), q_BH as
# direct data label per row, bold group headers inside the y axis, zero
# reference line. Marking rule reserved for q < 0.05 (bold q label + filled
# point): NO row enters it in any window.
# Source: output/analysis_04_appraisal_to_haemodynamics/
#   family_by_window.csv  (bootstrap_ci.csv cross-check: prefrontal rows only)
#
# Run from the repository root, after code/mod04_confirmatory_family.R:
#   Rscript code/fig4_confirmatory.R
# =============================================================================

source("code/style.R")
library(dplyr); library(readr); library(tidyr); library(tibble)
library(ggtext); library(patchwork)

A04 <- file.path("output", "analysis_04_appraisal_to_haemodynamics")
fam  <- read_csv(file.path(A04, "family_by_window.csv"), show_col_types = FALSE)
boot <- read_csv(file.path(A04, "bootstrap_ci.csv"),     show_col_types = FALSE)

# Fixed ROI colour mapping (same in every ROI figure)
roi_cols <- c(Prefrontal = "#0072B2", `Left temporal` = "#E69F00",
              `Right temporal` = "#009E73")
roi_lab <- c(PFC_Frontal = "Prefrontal", Left_Temporal = "Left temporal",
             Right_Temporal = "Right temporal")
axis_lab <- c(comfort_wi = "Comfort", richness_wi = "Richness")

d <- fam %>%
  mutate(
    axis  = axis_lab[term],
    roi   = roi_lab[roi],
    roi_c = c("Left temporal" = "Left temporal",
              "Right temporal" = "Right temporal", "Prefrontal" = "Prefrontal")[roi],
    crit  = qt(0.975, df),
    ci_lo = estimate - crit * se,
    ci_hi = estimate + crit * se,
    sig   = q_BH < 0.05
  )

# Y axis: bold header rows + spacer row between the two axis groups.
# Zero-width-space suffixes keep the repeated ROI labels unique for discrete y.
d <- d %>%
  mutate(y_key = paste0(roi, if_else(axis == "Comfort", "", "\u200B")),
         y_key = factor(y_key, levels = c(
           "Right temporal\u200B", "Left temporal\u200B", "Prefrontal\u200B",
           "**Richness**", " ",
           "Right temporal", "Left temporal", "Prefrontal", "**Comfort**")))
full_levels <- c("Right temporal", "Left temporal", "Prefrontal", "**Comfort**",
                 " ", "Right temporal\u200B", "Left temporal\u200B",
                 "Prefrontal\u200B", "**Richness**")
d$y_key <- factor(d$y_key, levels = full_levels)

windows <- c(difference = "a", response = "b", baseline = "c")

panel <- function(win, tag, show_y) {
  dd <- d %>% filter(analysis == win)
  # q-label column sits at a fixed x just right of the widest CI
  x_hi <- max(dd$ci_hi)
  x_lo <- min(dd$ci_lo)
  rng  <- x_hi - x_lo
  pad  <- rng * 0.06
  x_q  <- x_hi + rng * 0.24
  gg <- ggplot(dd, aes(y = y_key)) +
    annotate("richtext", x = x_q, y = length(full_levels) + 0.9,
             label = "*q*<sub>BH</sub>", size = 8 / .pt, family = "Helvetica",
             fill = NA, label.colour = NA) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55",
               linewidth = 0.45) +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi, colour = roi_c),
                  orientation = "y", width = 0, linewidth = 0.6,
                  show.legend = FALSE) +
    geom_point(aes(x = estimate, colour = roi_c, fill = roi_c,
                   size = sig, shape = sig), stroke = 0.7, show.legend = FALSE) +
    geom_text(aes(x = x_q, label = sprintf("%.3f", q_BH), fontface = ifelse(sig, 2, 1)),
              size = 8 / .pt, family = "Helvetica", hjust = 0.5) +
    scale_colour_manual(values = roi_cols) +
    scale_fill_manual(values = roi_cols) +
    scale_size_manual(values = c(`FALSE` = 1.9, `TRUE` = 2.6)) +
    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 16)) +
    scale_y_discrete(drop = FALSE, expand = expansion(add = c(0.7, 1.7))) +
    scale_x_continuous(limits = c(x_lo - pad, x_q + rng * 0.14), n.breaks = 4) +
    labs(x = "Slope (µM per axis unit)", y = NULL) +
    theme_pub(base_size = 8, axis_title_size = 9) +
    theme(axis.text.y = if (show_y) element_markdown(size = 8, hjust = 1)
          else element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y  = if (show_y) element_line(colour = "black", linewidth = 0.5)
          else element_blank())
  gg
}

p <- panel("difference", "a", TRUE) + panel("response", "b", FALSE) +
     panel("baseline", "c", FALSE) +
     plot_layout(ncol = 3, widths = c(1.28, 1, 1)) +
     plot_annotation(tag_levels = "a") &
     theme(plot.tag = element_text(size = 10, face = "bold"),
           plot.tag.position = c(0, 1))

save_fig(p, file.path(FIGDIR, "fig4_confirmatory.png"), 178, 62)

# ---- data locks (one per panel) + bootstrap cross-check ---------------------
for (win in names(windows)) {
  d %>% filter(analysis == win) %>%
    transmute(axis, roi = roi_c, estimate, se, df, ci_lo, ci_hi, q_BH,
              marked_q_lt_0.05 = sig) %>%
    write_csv(file.path(LOCKDIR, paste0("fig4", windows[win], "_", win, ".csv")))
}

chk <- boot %>%
  tidyr::separate(dv, into = c("analysis", "roi"), sep = "_", extra = "merge") %>%
  transmute(analysis, roi = roi_lab[roi], term,
            boot_lo = ci_lo, boot_hi = ci_hi) %>%
  left_join(d %>% mutate(roi = roi_c) %>%
              select(analysis, roi, term, estimate, ci_lo, ci_hi),
            by = c("analysis", "roi", "term"))
cat("\nAnalytic vs bootstrap 95% CI (prefrontal rows):\n")
print(as.data.frame(chk %>% mutate(across(where(is.numeric), ~ round(.x, 4)))),
      row.names = FALSE)
cat("\nq < 0.05 rows marked in any window:", sum(d$sig), "\n")
