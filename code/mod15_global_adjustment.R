# Module 15 — Spatial specificity of the appraisal association under global adjustment
#
# Port of exploration/analysis_15_spatial_specificity (batch-5 declaration):
#
# A11 found the [20,30) comfort peak "appears head-wide" and every downstream
# report carries that as a caveat. It has never been tested: across A01-A14,
# HbO_head is only ever an OUTCOME, never a covariate. A03's systemic exclusion
# does not reach this finding — it asked whether the MEAN OFFSET was graded by
# source-detector distance, before A11/A14 existed, across active channels
# spanning only 21.9-34.2 mm. This probe has NO short-separation channels
# (shortest active pair 21.9 mm; the two excluded pairs, 6 and 13, are the
# LONGEST at 45.1/45.4 mm), so the standard systemic regressor does not exist.
#
# The substitute is a LEAVE-THAT-UNIT-OUT global regressor: for target unit u,
# g_-u is the mean over that subject x condition's active channels NOT in u.
# Leave-out is mandatory, not fastidious: 13 of 21 active channels are
# prefrontal, so adjusting prefrontal for the whole-head mean would regress it
# partly on itself and bias the test against the very unit under study.
#
# Declared families (fixed before the original analysis ran):
#   PRIMARY   2 axes x 3 ROIs = 6 adjusted slopes on the response window [5,30), BH
#   SECONDARY 6 on the frozen difference score; 6 on the pre-stimulus baseline (BH each)
#   CHANNEL   21 channels x 2 axes = 42 on the response window, BH (mirrors A06)
#   Descriptive: the global component's own slope (declared target, A13 precedent,
#   still not a claim carrier); spatial-uniformity omnibus; corrected LT.
#
# Run from the repository root (Rscript code/mod15_global_adjustment.R).
# Reads only the Mendeley workbook via load_data.R; writes
# output/analysis_15_spatial_specificity/.

source("code/load_data.R")

OUT <- file.path("output", "analysis_15_spatial_specificity")
con <- open_log(file.path(OUT, "spatial_specificity_report.txt"))
rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")

AXES    <- c("comfort_wi", "richness_wi")
WINDOWS <- c(response   = "mean_response_m5_30",
             difference = "mean_response_minus_baseline_m5_30",
             baseline   = "mean_baseline_m10_0")

# ------------------------------------------------------------------ 0. the data
rule("0. DATA + PROBE GEOMETRY")

cd <- load_condition_analysis("HbO")           # subject x condition x channel, active only
chan_map <- load_channels() %>% filter(active_pair == 1) %>%
  distinct(subject_id, channel_index, roi, sd_distance_mm)
cd <- cd %>% left_join(chan_map, by = c("subject_id", "channel_index")) %>% filter(!is.na(roi))

cat("rows:", nrow(cd), " subjects:", n_distinct(cd$subject_id),
    " channels/cell (median):", median(table(paste(cd$subject_id, cd$marker))), "\n")

sd_summary <- chan_map %>% distinct(channel_index, sd_distance_mm) %>% arrange(sd_distance_mm)
cat("\nActive source-detector distances (mm): min", min(sd_summary$sd_distance_mm),
    " max", max(sd_summary$sd_distance_mm), "\n")
cat("No channel below 15 mm => no short-separation reference exists on this probe:",
    sum(sd_summary$sd_distance_mm < 15) == 0, "\n")

# helper: within/between decomposition of any cell-level regressor
decomp2 <- function(d, v) {
  d %>% group_by(subject_id) %>%
    mutate(!!paste0(v, "_bw") := mean(.data[[v]], na.rm = TRUE)) %>% ungroup() %>%
    mutate(!!paste0(v, "_wi") := .data[[v]] - .data[[paste0(v, "_bw")]])
}

# --------------------------------------------- 1. ROI families, global-adjusted
rule("1. ROI FAMILIES — leave-that-ROI-out global adjustment")

roi_rows <- list()
for (wn in names(WINDOWS)) {
  ycol <- WINDOWS[[wn]]
  for (rr in ROIS) {
    # target ROI mean and the leave-that-ROI-out global mean, per cell
    cell <- cd %>%
      group_by(subject_id, marker) %>%
      summarise(y   = mean(.data[[ycol]][roi == rr]),
                g   = mean(.data[[ycol]][roi != rr]),
                n_u = sum(roi == rr), n_g = sum(roi != rr), .groups = "drop") %>%
      filter(is.finite(y), is.finite(g), n_u > 0, n_g > 0) %>%
      inner_join(.analysis_design("full"), by = c("subject_id", "marker")) %>%
      decomp2("g")

    m_adj <- p24_lmm(cell, paste(c(AXES, "g_wi", "g_bw"), collapse = " + "), "y")
    m_raw <- p24_lmm(cell, paste(AXES, collapse = " + "), "y")
    fa <- fixef_table(m_adj, keep = AXES); fr <- fixef_table(m_raw, keep = AXES)
    gterm <- fixef_table(m_adj, keep = "g_wi")

    for (ax in AXES) {
      roi_rows[[length(roi_rows) + 1]] <- tibble(
        window = wn, unit = rr, axis = ax,
        n_cells = nrow(cell), n_chan_unit = round(mean(cell$n_u), 1),
        est_unadjusted = fr$estimate[fr$term == ax], p_unadjusted = fr$p[fr$term == ax],
        est_adjusted   = fa$estimate[fa$term == ax], se_adjusted = fa$se[fa$term == ax],
        p_adjusted     = fa$p[fa$term == ax],
        g_slope = gterm$estimate, g_p = gterm$p,
        shrinkage = 1 - fa$estimate[fa$term == ax] / fr$estimate[fr$term == ax])
    }
  }
}
roi_tab <- bind_rows(roi_rows) %>% group_by(window) %>%
  mutate(q_adjusted = bh(p_adjusted), q_unadjusted = bh(p_unadjusted)) %>% ungroup()

for (wn in names(WINDOWS)) {
  cat("\n--", wn, "window — family of 6, BH within window --\n")
  print(as.data.frame(roi_tab %>% filter(window == wn) %>%
    select(unit, axis, est_unadjusted, p_unadjusted, q_unadjusted,
           est_adjusted, p_adjusted, q_adjusted, shrinkage)), digits = 3, row.names = FALSE)
}
write_outcome(roi_tab, file.path(OUT, "roi_global_adjusted.csv"))

# ------------------------------------------- 2. the global component's own slope
rule("2. THE GLOBAL COMPONENT ITSELF (declared target; not a claim carrier)")

glob_rows <- list()
for (wn in names(WINDOWS)) {
  ycol <- WINDOWS[[wn]]
  cell <- cd %>% group_by(subject_id, marker) %>%
    summarise(y = mean(.data[[ycol]]), n_chan = n(), .groups = "drop") %>%
    inner_join(.analysis_design("full"), by = c("subject_id", "marker"))
  m <- p24_lmm(cell, paste(AXES, collapse = " + "), "y")
  ft <- fixef_table(m, keep = AXES)
  for (ax in AXES) glob_rows[[length(glob_rows) + 1]] <- tibble(
    window = wn, axis = ax, estimate = ft$estimate[ft$term == ax],
    se = ft$se[ft$term == ax], p = ft$p[ft$term == ax], n_cells = nrow(cell))
}
glob_tab <- bind_rows(glob_rows)
print(as.data.frame(glob_tab), digits = 3, row.names = FALSE)
write_outcome(glob_tab, file.path(OUT, "global_component_slopes.csv"))

# ------------------------------------------------- 3. channel family (response)
rule("3. CHANNEL FAMILY — 21 x 2 = 42, leave-that-channel-out adjustment, response window")

ycol <- WINDOWS[["response"]]
chans <- sort(unique(cd$channel_index))
ch_rows <- list()
for (ci in chans) {
  cell <- cd %>% group_by(subject_id, marker) %>%
    summarise(y = mean(.data[[ycol]][channel_index == ci]),
              g = mean(.data[[ycol]][channel_index != ci]),
              has = any(channel_index == ci), .groups = "drop") %>%
    filter(has, is.finite(y), is.finite(g)) %>%
    inner_join(.analysis_design("full"), by = c("subject_id", "marker")) %>%
    decomp2("g")
  m_adj <- p24_lmm(cell, paste(c(AXES, "g_wi", "g_bw"), collapse = " + "), "y")
  m_raw <- p24_lmm(cell, paste(AXES, collapse = " + "), "y")
  fa <- fixef_table(m_adj, keep = AXES); fr <- fixef_table(m_raw, keep = AXES)
  for (ax in AXES) ch_rows[[length(ch_rows) + 1]] <- tibble(
    channel_index = ci, axis = ax, n_cells = nrow(cell),
    est_unadjusted = fr$estimate[fr$term == ax], p_unadjusted = fr$p[fr$term == ax],
    est_adjusted = fa$estimate[fa$term == ax], se_adjusted = fa$se[fa$term == ax],
    p_adjusted = fa$p[fa$term == ax])
}
ch_tab <- bind_rows(ch_rows) %>%
  mutate(q_adjusted = bh(p_adjusted), q_unadjusted = bh(p_unadjusted)) %>%
  left_join(chan_map %>% distinct(channel_index, roi, sd_distance_mm), by = "channel_index")

cat("\nTop 10 by adjusted p (family of 42, BH):\n")
print(as.data.frame(ch_tab %>% arrange(p_adjusted) %>% head(10) %>%
  select(channel_index, roi, axis, est_unadjusted, p_unadjusted, q_unadjusted,
         est_adjusted, p_adjusted, q_adjusted)), digits = 3, row.names = FALSE)
cat("\nSurvivors at q<0.05, adjusted:", sum(ch_tab$q_adjusted < 0.05),
    " | unadjusted:", sum(ch_tab$q_unadjusted < 0.05), "\n")
cat("\nChannel 10 (the A06/A07/A08 candidate), both axes:\n")
print(as.data.frame(ch_tab %>% filter(channel_index == 10) %>%
  select(channel_index, roi, axis, est_unadjusted, p_unadjusted,
         est_adjusted, p_adjusted, q_adjusted)), digits = 3, row.names = FALSE)
write_outcome(ch_tab, file.path(OUT, "channel_global_adjusted.csv"))

# ------------------------------------------------ 4. spatial-uniformity omnibus
rule("4. SPATIAL-UNIFORMITY OMNIBUS — does comfort modulate the spatial PROFILE?")

# With a cell-level random intercept (1|subject:marker) the cell mean — i.e. the
# global component — is absorbed. The comfort x channel interaction is then a
# purely WITHIN-CELL contrast: it asks whether any channel departs from the
# global response. A flat profile is the signature of a global/systemic source.
long <- cd %>% mutate(y = .data[[WINDOWS[["response"]]]],
                      cell = paste(subject_id, marker, sep = "_"),
                      channel_f = factor(channel_index)) %>%
  filter(is.finite(y))
cat("rows:", nrow(long), " cells:", n_distinct(long$cell), " channels:", nlevels(long$channel_f), "\n")

f0 <- y ~ comfort_wi + richness_wi + channel_f + (1 | subject_id) + (1 | cell)
f1 <- y ~ comfort_wi * channel_f + richness_wi * channel_f + (1 | subject_id) + (1 | cell)
m0 <- lmer(f0, data = long, REML = FALSE,
           control = lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4)))
m1 <- lmer(f1, data = long, REML = FALSE,
           control = lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4)))
lrt <- anova(m0, m1)
print(lrt)
unif <- tibble(test = "comfort/richness x channel interaction (within-cell)",
               df = lrt$Df[2], chisq = lrt$Chisq[2], p = lrt$`Pr(>Chisq)`[2])
cat("\nReading: p >= 0.05 => the comfort effect is spatially FLAT once the cell\n",
    "(global) level is absorbed — the signature of a global/systemic source.\n", sep = "")
write_outcome(unif, file.path(OUT, "spatial_uniformity_omnibus.csv"))

# --------------------------------------- 5. corrected left-temporal, adjusted
rule("5. CORRECTED LEFT-TEMPORAL ROI (channels 18/19 dropped) UNDER THE SAME ADJUSTMENT")

corr_rows <- list()
for (wn in names(WINDOWS)) {
  ycol <- WINDOWS[[wn]]
  cell <- cd %>% filter(!channel_index %in% FLAT_CHANNELS) %>%
    group_by(subject_id, marker) %>%
    summarise(y = mean(.data[[ycol]][roi == "Left_Temporal"]),
              g = mean(.data[[ycol]][roi != "Left_Temporal"]),
              n_u = sum(roi == "Left_Temporal"), .groups = "drop") %>%
    filter(is.finite(y), is.finite(g), n_u > 0) %>%
    inner_join(.analysis_design("full"), by = c("subject_id", "marker")) %>%
    decomp2("g")
  m_adj <- p24_lmm(cell, paste(c(AXES, "g_wi", "g_bw"), collapse = " + "), "y")
  m_raw <- p24_lmm(cell, paste(AXES, collapse = " + "), "y")
  fa <- fixef_table(m_adj, keep = AXES); fr <- fixef_table(m_raw, keep = AXES)
  for (ax in AXES) corr_rows[[length(corr_rows) + 1]] <- tibble(
    window = wn, unit = "Left_Temporal_corrected", axis = ax,
    n_chan_unit = round(mean(cell$n_u), 1),
    est_unadjusted = fr$estimate[fr$term == ax], p_unadjusted = fr$p[fr$term == ax],
    est_adjusted = fa$estimate[fa$term == ax], p_adjusted = fa$p[fa$term == ax])
}
corr_tab <- bind_rows(corr_rows)
print(as.data.frame(corr_tab), digits = 3, row.names = FALSE)
write_outcome(corr_tab, file.path(OUT, "left_temporal_corrected_adjusted.csv"))

# ------------------------------- 6. CR2 audit of the ADJUSTED channel family
rule("6. CR2 SUBJECT-CLUSTER-ROBUST AUDIT OF THE ADJUSTED CHANNEL FAMILY (A08's test)")

# A08 established that CR2 is the binding correction for this design: it raised
# ch10's SE by 35% and took the unchanged 42-test family from one survivor to
# none. The adjusted family must face the same audit, or A15's channel result
# would be graded on an easier scale than A08's.
suppressPackageStartupMessages(library(clubSandwich))

# Bound explicitly: the section-5 loop leaves the shared `ycol` on the LAST
# window it visited (baseline), so re-using it here would silently audit the
# wrong window. Named separately so no later section can shadow it either.
YCOL_RESPONSE <- unname(WINDOWS[["response"]])
stopifnot(YCOL_RESPONSE == "mean_response_m5_30")

cr2_rows <- list()
for (ci in chans) {
  cell <- cd %>% group_by(subject_id, marker) %>%
    summarise(y = mean(.data[[YCOL_RESPONSE]][channel_index == ci]),
              g = mean(.data[[YCOL_RESPONSE]][channel_index != ci]),
              has = any(channel_index == ci), .groups = "drop") %>%
    filter(has, is.finite(y), is.finite(g)) %>%
    inner_join(.analysis_design("full"), by = c("subject_id", "marker")) %>%
    decomp2("g")
  # CR2 is applied to the SAME random-intercept fit as section 3 (A08's method):
  # it recomputes the standard error under arbitrary within-subject covariance
  # and leaves the point estimate untouched. Fitting an OLS here instead would
  # silently change the estimand from the within-subject slope to a pooled one.
  fit <- p24_lmm(cell, paste(c(AXES, "g_wi", "g_bw"), collapse = " + "), "y")
  ct <- as.data.frame(clubSandwich::coef_test(
    fit, vcov = "CR2", cluster = cell$subject_id, test = "Satterthwaite"))
  ct$term <- rownames(ct)
  for (ax in AXES) {
    i <- which(ct$term == ax)
    cr2_rows[[length(cr2_rows) + 1]] <- tibble(
      channel_index = ci, axis = ax,
      est = ct$beta[i], se_cr2 = ct$SE[i], df = ct$df_Satt[i], p_cr2 = ct$p_Satt[i])
  }
}
cr2_tab <- bind_rows(cr2_rows) %>% mutate(q_cr2 = bh(p_cr2)) %>%
  left_join(ch_tab %>% select(channel_index, axis, se_adjusted, p_adjusted),
            by = c("channel_index", "axis")) %>%
  mutate(se_inflation = se_cr2 / se_adjusted - 1)

cat("\nTop 6 by CR2 p (family of 42, BH):\n")
print(as.data.frame(cr2_tab %>% arrange(p_cr2) %>% head(6) %>%
  select(channel_index, axis, est, se_adjusted, se_cr2, se_inflation, p_adjusted,
         p_cr2, q_cr2)), digits = 3, row.names = FALSE)
cat("\nCR2 survivors at q<0.05:", sum(cr2_tab$q_cr2 < 0.05), "\n")
cat("\nChannel 10 under CR2, global-adjusted:\n")
print(as.data.frame(cr2_tab %>% filter(channel_index == 10)), digits = 3, row.names = FALSE)
write_outcome(cr2_tab, file.path(OUT, "channel_adjusted_cr2.csv"))

rule("DONE")
close_log(con)
