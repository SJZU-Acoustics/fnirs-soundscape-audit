# =============================================================================
# Module 05 — Alternative response geometry.
#
# Module 04 tested linear within-subject Comfort and Richness slopes. This pass
# asks whether that representation hides either:
#   (1) curvature / an axis interaction in the measured appraisal plane; or
#   (2) a non-monotonic four-level nominal-quadrant response.
#
# The family was fixed before this analysis ran: within each window,
# 3 ROIs x 2 alternative representations = 6 omnibus tests, Benjamini-Hochberg
# corrected together. Individual quadratic terms are descriptive unless their
# parent omnibus test survives.
#
# Outputs (output/analysis_05_alternative_response_geometry/):
#   roi_rebuild_gate.csv             aggregation gate (hard check, not a result)
#   alternative_geometry_omnibus.csv the 18 omnibus tests, BH within each window
#   quadratic_terms.csv              fixed effects of the quadratic extension
#   nominal_quadrant_terms.csv       fixed effects of the nominal-quadrant fit
#
# Run from the repository root: Rscript code/mod05_geometry.R
# =============================================================================

source("code/load_data.R")

OUT <- file.path("output", "analysis_05_alternative_response_geometry")
con <- open_log(file.path(OUT, "mod05_geometry_report.txt"))

WINDOWS <- c(
  difference = "mean_response_minus_baseline_m5_30",
  response   = "mean_response_m5_30",
  baseline   = "mean_baseline_m10_0"
)

# Subject x condition x ROI window means joined to the canonical design. The
# aggregation includes every active channel in the delivered ROI definition
# (the module-04 corrected ROI is a separate, explicitly declared variant).
load_condition_roi_analysis <- function(chromophore = "HbO") {
  map <- load_channels() %>%
    filter(active_pair == 1) %>%
    distinct(subject_id, channel_index, roi)
  load_condition_metrics(active_only = TRUE) %>%
    filter(.data$chromophore == !!chromophore) %>%
    left_join(map, by = c("subject_id", "channel_index")) %>%
    filter(!is.na(roi)) %>%
    group_by(subject_id, marker, roi) %>%
    summarise(
      across(c(mean_baseline_m10_0, mean_stim_m0_30, mean_response_m5_30,
               mean_response_minus_baseline_m5_30), mean),
      n_channels = n(),
      .groups = "drop"
    ) %>%
    inner_join(.analysis_design("full"), by = c("subject_id", "marker"))
}

d <- load_condition_roi_analysis(chromophore = "HbO") %>%
  mutate(quadrant_f = factor(quadrant_design,
                             levels = c("Q1", "Q2", "Q3", "Q4"))) %>%
  group_by(subject_id, roi) %>%
  mutate(
    comfort2_wi = comfort_wi^2 - mean(comfort_wi^2),
    richness2_wi = richness_wi^2 - mean(richness_wi^2),
    comfort_x_richness_wi =
      comfort_wi * richness_wi - mean(comfort_wi * richness_wi)
  ) %>%
  ungroup()

rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n",
                        strrep("=", 74), "\n", sep = "")

cat("Module 05 — alternative response geometry\n")
cat("Haemoglobin unit:", attr(load_p24(), "hb_unit"), "\n")
cat("Family per window: 3 ROIs x 2 representations = 6 omnibus tests, BH.\n")

# Verify that the condition-channel aggregation reproduces the delivered ROI
# difference scores. This is a hard gate, not a result.
gate <- d %>%
  select(subject_id, marker, roi, mean_response_minus_baseline_m5_30) %>%
  pivot_wider(names_from = roi,
              values_from = mean_response_minus_baseline_m5_30) %>%
  inner_join(load_p24() %>%
               select(subject_id, marker, starts_with("HbO_")),
             by = c("subject_id", "marker"))
gate_check <- map_dfr(ROIS, function(rr) {
  tibble(
    roi = rr,
    n_compared = sum(!is.na(gate[[rr]]) &
                       !is.na(gate[[paste0("HbO_", rr)]])),
    max_abs_diff = max(abs(gate[[rr]] - gate[[paste0("HbO_", rr)]]),
                       na.rm = TRUE)
  )
})
stopifnot(max(gate_check$max_abs_diff) < 1e-12)
write_outcome(gate_check, file.path(OUT, "roi_rebuild_gate.csv"))

omnibus <- list()
quad_terms <- list()
quadratic_terms <- list()

for (win in names(WINDOWS)) {
  y <- unname(WINDOWS[[win]])
  rule(paste("WINDOW:", win))

  for (rr in ROIS) {
    dd <- d %>% filter(roi == rr, !is.na(.data[[y]]))

    # Quadratic extension of the same measured-axis model used by module 04.
    m_linear <- p24_lmm(
      dd, "comfort_wi + richness_wi", y, reml = FALSE
    )
    m_quad <- p24_lmm(
      dd,
      paste(
        "comfort_wi + richness_wi + comfort2_wi + richness2_wi +",
        "comfort_x_richness_wi"
      ),
      y,
      reml = FALSE
    )
    lr <- anova(m_linear, m_quad)
    omnibus[[length(omnibus) + 1]] <- tibble(
      window = win,
      roi = rr,
      representation = "quadratic measured axes",
      test = "likelihood-ratio extension",
      df_num = lr$Df[2],
      statistic = lr$Chisq[2],
      p = lr$`Pr(>Chisq)`[2],
      AIC_reference = AIC(m_linear),
      AIC_alternative = AIC(m_quad)
    )
    quadratic_terms[[length(quadratic_terms) + 1]] <-
      fixef_table(
        m_quad,
        keep = c("comfort_wi", "richness_wi", "comfort2_wi",
                 "richness2_wi", "comfort_x_richness_wi")
      ) %>%
      mutate(window = win, roi = rr, .before = 1)

    # Nominal quadrant makes no equal-spacing or monotonicity assumption.
    m_quad_nominal <- p24_lmm(dd, "quadrant_f", y)
    aq <- anova(m_quad_nominal)
    omnibus[[length(omnibus) + 1]] <- tibble(
      window = win,
      roi = rr,
      representation = "nominal quadrant",
      test = "three-df omnibus F",
      df_num = aq$NumDF[1],
      statistic = aq$`F value`[1],
      p = aq$`Pr(>F)`[1],
      AIC_reference = NA_real_,
      AIC_alternative = AIC(m_quad_nominal)
    )
    quad_terms[[length(quad_terms) + 1]] <-
      fixef_table(m_quad_nominal) %>%
      filter(startsWith(term, "quadrant_f")) %>%
      mutate(window = win, roi = rr, .before = 1)
  }
}

omnibus <- bind_rows(omnibus) %>%
  group_by(window) %>%
  mutate(q_BH = bh(p), family_size = n()) %>%
  ungroup() %>%
  arrange(factor(window, levels = names(WINDOWS)), p)
quadratic_terms <- bind_rows(quadratic_terms) %>%
  left_join(
    omnibus %>%
      filter(representation == "quadratic measured axes") %>%
      select(window, roi, parent_omnibus_p = p, parent_omnibus_q = q_BH),
    by = c("window", "roi")
  )
quad_terms <- bind_rows(quad_terms) %>%
  left_join(
    omnibus %>%
      filter(representation == "nominal quadrant") %>%
      select(window, roi, parent_omnibus_p = p, parent_omnibus_q = q_BH),
    by = c("window", "roi")
  )

rule("OMNIBUS RESULTS")
print(as.data.frame(omnibus), digits = 4)
cat("\nSurviving q < 0.05:", sum(omnibus$q_BH < 0.05), "of",
    nrow(omnibus), "across the three separate six-test families.\n")

rule("QUADRATIC TERMS (DESCRIPTIVE UNLESS THE PARENT OMNIBUS SURVIVES)")
print(as.data.frame(
  quadratic_terms %>%
    select(window, roi, term, estimate, se, p,
           parent_omnibus_p, parent_omnibus_q)
), digits = 4)

write_outcome(omnibus, file.path(OUT, "alternative_geometry_omnibus.csv"))
write_outcome(quadratic_terms, file.path(OUT, "quadratic_terms.csv"))
write_outcome(quad_terms, file.path(OUT, "nominal_quadrant_terms.csv"))

close_log(con)
cat("Module 05 done.\n")
