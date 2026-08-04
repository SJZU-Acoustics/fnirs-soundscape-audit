# M24 — the temporal ROIs under the dark-FRACTION mask.
#
# The session-CV criterion the earlier correction pass (M20) used cannot see
# INTERMITTENT dropout. A channel alternating between tissue light and the
# instrument's dark floor has MORE variation than a clean one, so any screen
# built from variation keeps it — the same failure mode as the collector's own
# SNR = mean/SD, which scored a permanently dead channel at ~6,700 against a
# threshold of 2. The two screens fail in opposite directions and neither sees
# the middle. The replacement is absolute: a sample is dark when both
# wavelengths sit in the dark floor AND agree to 1e-4, and a cell is excluded
# at >= 5% dark (plus the single rail-saturated cell a dark-floor test cannot
# see). The mask is recomputed from the workbook's channel_mask sheet by
# load_channel_mask(), which verifies the recomputation against the deposited
# flags.
#
# This is a CORRECTION pass, not a discovery pass — the second one, and it
# mirrors M20 exactly so the three mask states (none / CV / dark-fraction) are
# directly comparable on the same declared families. Nothing new is targeted.
#
# Reading rule, fixed in advance and deliberately not symmetric, carried over
# from M20: the declared families were null, so a stricter mask can only remove
# an artefactual result or leave the null standing. A survivor appearing under
# the new mask would be a NEW finding on a corrected measure, declared as such,
# never inherited.

source("code/load_data.R")

OUT <- file.path("output", "analysis_24_dark_fraction_mask")
con <- open_log(file.path(OUT, "dark_fraction_mask_report.txt"))

rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")
AXES <- c("comfort_wi", "richness_wi")

cat("M24 — temporal ROIs under the dark-fraction mask.\n")
cat("All concentrations in uM; slopes in uM per unit of the axis.\n")

# ---------------------------------------------------------------- local masks

#' Subject x channel cells the dark-fraction criterion excludes (the
#' mask_excluded component of the union mask), with the sheet's frac_dark.
dark_fraction_cells <- function() {
  load_channel_mask() %>%
    filter(mask_excluded == 1) %>%
    distinct(subject_id, channel_index, roi, detector, frac_dark)
}

#' Rebuild subject x condition ROI means under the superseded session-CV mask,
#' from the condition x channel layer — needed only for section 4's
#' mask-state comparison of the prefrontal line. Mirrors the M20 helper.
load_p24_roi_masked_cv <- function(chromophore = "HbO") {
  dead <- load_dead_channels() %>% distinct(subject_id, channel_index)
  map <- load_channels() %>% filter(active_pair == 1) %>%
    distinct(subject_id, channel_index, roi)
  agg <- load_condition_metrics(active_only = TRUE) %>%
    filter(.data$chromophore == !!chromophore) %>%
    anti_join(dead, by = c("subject_id", "channel_index")) %>%
    left_join(map, by = c("subject_id", "channel_index")) %>%
    filter(!is.na(roi)) %>%
    group_by(subject_id, marker, roi) %>%
    summarise(value = mean(mean_response_minus_baseline_m5_30),
              n_live = n(), .groups = "drop")
  wide_v <- agg %>% mutate(roi = paste0(chromophore, "m_", roi)) %>%
    select(subject_id, marker, roi, value) %>%
    pivot_wider(names_from = roi, values_from = value)
  wide_n <- agg %>% mutate(roi = paste0("n_live_", roi)) %>%
    select(subject_id, marker, roi, n_live) %>%
    pivot_wider(names_from = roi, values_from = n_live)
  load_p24() %>%
    left_join(wide_v, by = c("subject_id", "marker")) %>%
    left_join(wide_n, by = c("subject_id", "marker"))
}

# ================================================================== 0. gates
rule("0. GATES — the rebuild must be the builder's, and the mask must be real")

cat("(a) does the rebuild reproduce the delivered ROI columns when nothing is\n")
cat("    excluded? If not, the aggregation differs from the builder's and\n")
cat("    nothing below is interpretable.\n\n")
gate_a <- check_roi_rebuild()
print(as.data.frame(gate_a), digits = 3)
stopifnot(max(gate_a$max_abs_diff) < 1e-9)
cat("\n    GATE A PASSED (worst discrepancy", signif(max(gate_a$max_abs_diff), 3), "uM).\n")

cat("\n(b) is the criterion specific? The dark floor must be absent from\n")
cat("    detectors that never failed, or the cut is measuring something else.\n\n")
all_cells <- load_channel_mask()
spec <- all_cells %>%
  group_by(detector) %>%
  summarise(n_cells = n(),
            median_frac_dark = median(frac_dark),
            n_ge_5pct = sum(frac_dark >= DARK_FRACTION_CUT),
            n_ge_99pct = sum(frac_dark >= 0.99), .groups = "drop") %>%
  arrange(detector)
print(as.data.frame(spec), digits = 3)
cat("\n    Detectors 1-5 are the never-failed set. If they carry zero cells at\n")
cat("    the cut, the criterion is not firing on ordinary signal.\n")
clean_hits <- spec %>% filter(detector %in% 1:5) %>% pull(n_ge_5pct) %>% sum()
cat("    cells >= 5% dark on detectors 1-5:", clean_hits, "\n")
stopifnot(clean_hits == 0)
cat("    GATE B PASSED.\n")

cat("\n(c) what the new mask adds over the CV mask it superseded:\n\n")
new  <- dark_fraction_cells()
old  <- load_dead_channels() %>% distinct(subject_id, channel_index)
added <- new %>% anti_join(old, by = c("subject_id", "channel_index"))
lost  <- old %>% anti_join(new %>% distinct(subject_id, channel_index),
                           by = c("subject_id", "channel_index"))
cat("    CV mask cells        :", nrow(old), "\n")
cat("    dark-fraction cells  :", nrow(new), "\n")
cat("    newly excluded       :", nrow(added), "\n")
cat("    excluded before, not now:", nrow(lost), "\n\n")
print(as.data.frame(added %>% group_by(roi, detector) %>%
                      summarise(n_cells = n(), .groups = "drop")), digits = 3)
write_outcome(added %>% arrange(desc(frac_dark)), file.path(OUT, "cells_newly_excluded.csv"))

cat("\n    the cells the CV mask kept that are almost wholly dark:\n\n")
print(as.data.frame(added %>% filter(frac_dark >= 0.5) %>% arrange(desc(frac_dark))), digits = 3)

# ================================================= 1. what the ROIs now contain
rule("1. ROI COMPOSITION — how many live channels each subject actually has")

map <- load_channels() %>% filter(active_pair == 1) %>%
  distinct(subject_id, channel_index, roi)
# the union mask: dark-fraction cells plus the wholesale detector-8 exclusion
masked_keys <- bind_rows(
  new %>% distinct(subject_id, channel_index),
  load_channels() %>% filter(detector %in% UNUSABLE_DETECTORS) %>%
    distinct(subject_id, channel_index)) %>% distinct()

live <- map %>% anti_join(masked_keys, by = c("subject_id", "channel_index")) %>%
  dplyr::count(subject_id, roi, name = "n_live") %>%
  complete(subject_id, roi, fill = list(n_live = 0L))

comp <- live %>% dplyr::count(roi, n_live) %>%
  pivot_wider(names_from = n_live, values_from = n, values_fill = 0) %>%
  arrange(roi)
print(as.data.frame(comp))
write_outcome(live, file.path(OUT, "live_channels_per_subject.csv"))

cat("\nsubjects with NO live channel in an ROI:\n\n")
empty <- live %>% filter(n_live == 0) %>% arrange(roi, subject_id)
print(as.data.frame(empty %>% dplyr::count(roi, name = "n_subjects")))
cat("\n")
print(as.data.frame(empty), max = 200)
cat("\nThe CV-mask pass recorded FOUR subjects with no live right-temporal\n")
cat("channel. The count above is what the corrected criterion gives.\n")

# ============================================ 2. the declared families, re-run
rule("2. DECLARED FAMILIES UNDER THE NEW MASK — 2 axes x 3 ROIs, three windows")

bm <- load_block_metrics(active_only = TRUE) %>% filter(chromophore == "HbO")
mk <- masked_keys %>% mutate(.dead = TRUE)

roi_blocks <- bm %>%
  left_join(mk, by = c("subject_id", "channel_index")) %>%
  filter(is.na(.dead)) %>%
  left_join(map, by = c("subject_id", "channel_index")) %>%
  filter(!is.na(roi)) %>%
  group_by(subject_id, event_index, sequence, marker, onset_sec, roi) %>%
  summarise(baseline = mean(baseline_m10_0), response = mean(response_m5_30),
            difference = mean(response_minus_baseline_m5_30), .groups = "drop")

cells <- roi_blocks %>% group_by(subject_id, marker, roi) %>%
  summarise(across(c(baseline, response, difference), mean), .groups = "drop") %>%
  pivot_wider(names_from = roi, values_from = c(baseline, response, difference),
              names_sep = "_") %>%
  left_join(load_p24() %>% select(subject_id, marker, all_of(AXES)),
            by = c("subject_id", "marker"))

family_scan <- function(data, dv_prefix, label) {
  rows <- map_dfr(ROIS, function(rr) {
    y <- paste0(dv_prefix, rr)
    if (!y %in% names(data)) return(tibble())
    m <- p24_lmm(data, paste(AXES, collapse = " + "), y)
    fixef_table(m, keep = AXES) %>% mutate(roi = rr, .before = 1)
  })
  rows %>% mutate(q_BH = bh(p), family_size = nrow(rows), analysis = label)
}

fam <- map_dfr(c("difference_", "response_", "baseline_"), function(pref) {
  family_scan(cells, pref, sub("_$", "", pref))
})
print(as.data.frame(fam %>% select(analysis, roi, term, estimate, se, p, q_BH)), digits = 3)
write_outcome(fam, file.path(OUT, "families_dark_fraction_mask.csv"))

cat("\nsurvivors at q < 0.05, by window:\n")
print(as.data.frame(fam %>% group_by(analysis) %>%
        summarise(n_tests = n(), n_survivors = sum(q_BH < 0.05), min_q = min(q_BH))),
      digits = 3)
cat("\nTOTAL survivors across the 18 declared tests:", sum(fam$q_BH < 0.05), "\n")

# ================================== 3. side by side with the two earlier states
rule("3. THE THREE MASK STATES SIDE BY SIDE")

prev <- file.path("output", "analysis_20_detector_dropout", "families_masked.csv")
if (file.exists(prev)) {
  a20 <- read_csv(prev, show_col_types = FALSE) %>%
    select(analysis, roi, term, est_cv = estimate, se_cv = se, p_cv = p, q_cv = q_BH)
  cmp <- fam %>%
    select(analysis, roi, term, est_new = estimate, se_new = se,
           p_new = p, q_new = q_BH) %>%
    left_join(a20, by = c("analysis", "roi", "term")) %>%
    mutate(d_est = est_new - est_cv)
  print(as.data.frame(cmp %>%
    select(analysis, roi, term, est_cv, est_new, d_est, p_cv, p_new, q_new)),
    digits = 3)
  write_outcome(cmp, file.path(OUT, "families_cv_vs_dark_fraction.csv"))
  cat("\nlargest absolute movement in an estimate:",
      signif(max(abs(cmp$d_est), na.rm = TRUE), 3), "uM\n")
  cat("tests changing side of q = 0.05:",
      sum((cmp$q_cv < 0.05) != (cmp$q_new < 0.05), na.rm = TRUE), "\n")
} else {
  cat("M20 outcome not found (run code/mod20_detector_dropout.R first);",
      "comparison skipped.\n")
}

# =============================== 4. does the prefrontal line move at all?
rule("4. THE PREFRONTAL LINE — the ROI every surviving conclusion rests on")

pf_frozen <- load_p24_roi_masked_cv("HbO") %>%
  select(subject_id, marker, HbOm_PFC_Frontal, n_live_PFC_Frontal)
pf_new <- load_p24_roi_masked("HbO") %>%
  select(subject_id, marker, new_PFC = HbOm_PFC_Frontal,
         n_live_new = n_live_PFC_Frontal)
pfc <- pf_frozen %>% left_join(pf_new, by = c("subject_id", "marker"))
cat("prefrontal cell values, CV mask versus dark-fraction mask:\n\n")
print(as.data.frame(pfc %>% summarise(
  n_cells = n(),
  r = cor(HbOm_PFC_Frontal, new_PFC, use = "complete.obs"),
  max_abs_diff = max(abs(HbOm_PFC_Frontal - new_PFC), na.rm = TRUE),
  median_abs_diff = median(abs(HbOm_PFC_Frontal - new_PFC), na.rm = TRUE),
  n_subjects_13ch = n_distinct(subject_id[n_live_new == 13]),
  n_subjects_fewer = n_distinct(subject_id[n_live_new < 13]))), digits = 4)

cat("\nchannel 10 is a single channel, so no ROI mask touches it. Confirming it\n")
cat("is not itself a newly-flagged cell in any subject:\n\n")
ch10 <- all_cells %>% filter(channel_index == 10) %>%
  summarise(n = n(), max_frac_dark = max(frac_dark), n_ge_cut = sum(frac_dark >= DARK_FRACTION_CUT))
print(as.data.frame(ch10), digits = 3)

rule("VERDICT")
cat("The union mask (dark-fraction cells plus wholesale detector 8) is the\n")
cat("release default: load_p24_roi_masked() applies it, and M20's CV-mask\n")
cat("outputs are regenerated alongside so the three mask states remain\n")
cat("directly comparable. Survivors across the 18 declared tests under the\n")
cat("new mask:", sum(fam$q_BH < 0.05), ".\n")

close_log(con)
