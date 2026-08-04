# =============================================================================
# Figure 6 — "The chain, quantified" (85 mm single column: 4 configs x 3 ROIs
# is low density; near-square panel)
# Lambda (smallest detectable peak response, 80% power) under the nested
# chain C0 ideal -> C1 +high-pass -> C2 +channel loss -> C3 +dependence,
# one line per ROI (corrected-ROI primary rows of lambda_chain.csv).
# Reference: 0.073 µM peak implied by the observed prefrontal difference-window
# slope (-0.052 µM window-metric, /0.71 -> -0.073; see mod22_chain.R report) —
# magnitude plotted, must fall below prefrontal lambda at C0 and C3.
# Source: output/analysis_22_chain/lambda_chain.csv (code/mod22_chain.R)
# =============================================================================

source("code/style.R")
library(dplyr); library(readr); library(tidyr); library(tibble)

A22 <- file.path("output", "analysis_22_chain")
lam <- read_csv(file.path(A22, "lambda_chain.csv"), show_col_types = FALSE)

roi_cols <- c(Prefrontal = "#0072B2", `Left temporal` = "#E69F00",
              `Right temporal` = "#009E73")
roi_lab <- c(PFC_Frontal = "Prefrontal", Left_Temporal = "Left temporal",
             Right_Temporal = "Right temporal")
cfg_lab <- c("Ideal", "+ High-pass", "+ Channel\nloss", "+ Dependence")

d <- lam %>% filter(se_source == "corrected ROI (primary)") %>%
  mutate(roi = roi_lab[roi]) %>%
  pivot_longer(cols = starts_with("C"), names_to = "config", values_to = "lambda") %>%
  mutate(x = as.integer(substr(config, 2, 2)) + 1L)

REF <- 0.073  # |peak implied by the observed prefrontal slope| (mod22 report)

p <- ggplot(d, aes(x, lambda, colour = roi, shape = roi, group = roi)) +
  geom_hline(yintercept = REF, linetype = "dashed", colour = "grey40",
             linewidth = 0.45) +
  annotate("text", x = 2.3, y = 0.057,
           label = "implied by observed prefrontal slope (0.073)",
           size = 7 / .pt, family = "Helvetica", hjust = 0.5,
           colour = "grey25") +
  geom_line(linewidth = 0.7, show.legend = FALSE) +
  geom_point(size = 2.2, stroke = 0.7) +
  geom_text(data = d %>% filter(x == 4),
            aes(label = sprintf("%.3f", lambda)),
            hjust = -0.35, size = 7 / .pt, family = "Helvetica",
            show.legend = FALSE) +
  scale_colour_manual(values = roi_cols, limits = names(roi_cols)) +
  scale_shape_manual(values = c(Prefrontal = 16, `Left temporal` = 17,
                                `Right temporal` = 15),
                     limits = names(roi_cols)) +
  scale_x_continuous(breaks = 1:4, labels = cfg_lab,
                     limits = c(0.7, 4.6), expand = expansion(add = 0.05)) +
  scale_y_continuous(breaks = c(0.05, 0.10, 0.15, 0.20),
                     limits = c(0.04, 0.255),
                     expand = expansion(mult = c(0, 0.02))) +
  guides(colour = guide_legend(override.aes = list(linewidth = 0.7))) +
  labs(x = NULL, y = "Detectable peak response\nλ (µM per axis unit)") +
  theme_pub() +
  theme(axis.text.x = element_text(size = 7.5),
        legend.position.inside = c(0.03, 0.99),
        legend.justification.inside = c(0, 1),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.3, "cm"))

save_fig(p, file.path(FIGDIR, "fig6_chain.png"), 85, 85)

# ---- data lock --------------------------------------------------------------
d %>% transmute(roi, config, x, lambda, reference_line = REF) %>%
  write_csv(file.path(LOCKDIR, "fig6a_lambda_chain.csv"))

# anchor check: reference must sit below prefrontal lambda at every config
pfc <- d %>% filter(roi == "Prefrontal")
stopifnot(all(pfc$lambda > REF))
cat("fig6 done. Prefrontal lambda min:", round(min(pfc$lambda), 4),
    "> reference", REF, "\n")
