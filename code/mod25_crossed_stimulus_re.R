# =============================================================================
# Module 25 — The stimulus sample as a dependence structure.
#
# POST-REVIEW SENSITIVITY (2026-08-05). This is NOT one of the campaign's
# declared estimator families and does not appear in main-text Table 2. It was
# run after the analysis campaign closed, in response to a review point that
# the primary model carries only a participant random intercept while the same
# 16 clips are heard by many participants.
#
# The declared region-of-interest family is refitted with crossed participant
# and clip random intercepts, (1 | subject_id) + (1 | clip_id), in place of the
# participant intercept alone, on all three analysis windows. Each of the 16
# clips was heard by 16-18 participants and each participant heard four, with
# clips nested within stimulus group, so participant and clip are crossed
# within group.
#
# Supports main-text Section 3.3 (the fourth link, dependence structure) and
# Supplementary Section "The stimulus sample as a dependence structure" with
# its Table S4.
#
# Outputs (output/analysis_25_crossed_stimulus_re/):
#   crossed_vs_subject_only.csv  18 tests, both random-effect structures
#   tableS4_crossed_re.csv       the SI Table S4 body
#
# Run from the repository root: Rscript code/mod25_crossed_stimulus_re.R
# =============================================================================

source("code/load_data.R")

OUT <- file.path("output", "analysis_25_crossed_stimulus_re")
con <- open_log(file.path(OUT, "mod25_crossed_stimulus_re_report.txt"))

AXES    <- c("comfort_wi", "richness_wi")
ROI_LAB <- c(PFC_Frontal = "Prefrontal", Left_Temporal = "Left temporal",
             Right_Temporal = "Right temporal")
WIN_LAB <- c(difference = "Difference", response = "Response", baseline = "Baseline")

rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")
cat("Module 25 - crossed participant + clip random intercepts.\n")
cat("A POST-REVIEW SENSITIVITY, not a declared family. Slopes in uM per axis unit.\n")

d  <- load_p24()
bm <- load_block_metrics(active_only = TRUE)

# ROI window means per subject x condition, by the same path as module 04.
roi_blocks <- bm %>% filter(chromophore == "HbO") %>%
  left_join(load_channels() %>% filter(active_pair == 1) %>%
              distinct(subject_id, channel_index, roi),
            by = c("subject_id", "channel_index")) %>%
  filter(!is.na(roi)) %>%
  group_by(subject_id, event_index, marker, roi) %>%
  summarise(baseline = mean(baseline_m10_0), response = mean(response_m5_30),
            difference = mean(response_minus_baseline_m5_30), .groups = "drop")

cells <- roi_blocks %>% group_by(subject_id, marker, roi) %>%
  summarise(across(c(baseline, response, difference), mean), .groups = "drop") %>%
  left_join(d %>% select(subject_id, marker, clip_id, all_of(AXES)),
            by = c("subject_id", "marker"))

rule("0. CROSSING STRUCTURE")
cat("clips:", dplyr::n_distinct(cells$clip_id),
    "| participants per clip:", paste(range(table(cells$clip_id[cells$roi == "PFC_Frontal"])), collapse = "-"),
    "| clips per participant:", paste(range(table(cells$subject_id[cells$roi == "PFC_Frontal"])), collapse = "-"), "\n")

#' One window, one random-effect structure: the declared family of 6.
fam <- function(win, keep_clip) {
  rows <- map_dfr(names(ROI_LAB), function(rr) {
    dd <- cells %>% filter(roi == rr)
    m  <- p24_lmm(dd, paste(AXES, collapse = " + "), win, keep_clip = keep_clip)
    vc <- as.data.frame(lme4::VarCorr(m))
    fixef_table(m, keep = AXES) %>%
      mutate(roi = rr,
             clip_sd = if (keep_clip) vc$sdcor[vc$grp == "clip_id"] else NA_real_,
             .before = 1)
  })
  rows %>% mutate(q_BH = bh(p), window = win, family_size = nrow(rows))
}

cmp <- map_dfr(names(WIN_LAB), function(w) {
  a <- fam(w, FALSE); b <- fam(w, TRUE)
  tibble(window = w, roi = a$roi, axis = a$term,
         estimate_subject = a$estimate, estimate_crossed = b$estimate,
         se_subject = a$se, se_crossed = b$se,
         se_inflation_pct = 100 * (b$se / a$se - 1),
         p_subject = a$p, p_crossed = b$p,
         q_subject = a$q_BH, q_crossed = b$q_BH,
         clip_sd = b$clip_sd)
})

rule("1. THE DECLARED FAMILY UNDER BOTH RANDOM-EFFECT STRUCTURES")
print(as.data.frame(cmp), digits = 3)
write_outcome(cmp, file.path(OUT, "crossed_vs_subject_only.csv"))

rule("2. WHAT CHANGES")
summ <- cmp %>% group_by(window) %>%
  summarise(survivors_subject = sum(q_subject < 0.05),
            survivors_crossed = sum(q_crossed < 0.05),
            max_se_inflation_pct = max(se_inflation_pct),
            max_clip_sd = max(clip_sd),
            max_abs_estimate_shift = max(abs(estimate_crossed - estimate_subject)),
            .groups = "drop")
print(as.data.frame(summ), digits = 3)

cat("\nSurvivors at q < 0.05, all 18 tests: subject-only",
    sum(cmp$q_subject < 0.05), "| crossed", sum(cmp$q_crossed < 0.05), "\n")
cat("Clip variance is at the boundary (0) wherever clip_sd rounds to 0.0000;\n")
cat("the crossed model then reduces to the reported one.\n")

# SI Table S4 body, in display order and rounding.
tableS4 <- cmp %>%
  transmute(window = unname(WIN_LAB[window]), region = unname(ROI_LAB[roi]),
            axis = ifelse(axis == "comfort_wi", "Comfort", "Richness"),
            estimate_subject = round(estimate_subject, 4),
            estimate_crossed = round(estimate_crossed, 4),
            se_subject = round(se_subject, 4), se_crossed = round(se_crossed, 4),
            q_subject = round(q_subject, 3), q_crossed = round(q_crossed, 3),
            clip_sd = round(clip_sd, 4))
write_outcome(tableS4, file.path(OUT, "tableS4_crossed_re.csv"))

close_log(con)
