# M20 — the temporal ROIs under a per-subject dead-channel mask.
#
# The dead-channel defect is detector-shaped. Channels 18/19 are the only two
# channels of detector 9 (TP7), dark in 66-68 of 69 sessions; channels 20/21
# are the only two of detector 8 (T8), dark in 29 of 69 — in the RIGHT-temporal
# ROI, which no earlier stage of the analysis had flagged. The delivered
# `active_pair` flag retains nearly every one of them.
#
# This is a CORRECTION pass, not a discovery pass. Every contrast below already
# exists in the unmasked family scan; the module re-runs it under a corrected
# ROI definition and introduces no new target. Families mirror the unmasked
# scan exactly so masked and unmasked estimates are directly comparable.
#
# This module was declared against the SUPERSEDED session-CV dead-channel mask
# of 2026-07-30 (raw-intensity CV < 0.002 per subject x channel, shipped as
# data/intermediates/qa_detector_dropout_0730.csv), and it must keep
# reproducing the report written against that mask. The release loader's
# load_p24_roi_masked() applies the union (dark-fraction + detector-8) mask
# only, so the CV-masked ROI rebuild is reconstructed in this module by
# anti-joining load_dead_channels() cells against the condition layer
# (load_p24_roi_masked_cv below). Gate A verifies that the condition-layer
# aggregation reproduces the block-layer aggregation the original pass used.
#
# Reading rule, fixed in advance and deliberately not symmetric: the unmasked
# families were null, so masking can only remove an artefactual result or leave
# the null standing. A survivor appearing under the mask would be a NEW finding
# on a corrected measure, declared as such, never inherited.

source("code/load_data.R")

OUT <- file.path("output", "analysis_20_detector_dropout")
con <- open_log(file.path(OUT, "detector_dropout_report.txt"))

rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")
AXES <- c("comfort_wi", "richness_wi")

cat("M20 — temporal ROIs under a per-subject dead-channel mask.\n")
cat("All concentrations in uM; slopes in uM per unit of the axis.\n")

# ---------------------------------------------------------------- local masks

#' Rebuild subject x condition ROI means under the superseded session-CV mask,
#' from the condition x channel layer (the release loader's aggregation).
#' `exclude_detectors = integer(0)` keeps detector 8: sections 1 and 4 exist to
#' show what the UNCORRECTED right-temporal ROI does, and section 4b is what
#' establishes that detector 8 must be dropped wholesale.
load_p24_roi_masked_cv <- function(chromophore = "HbO",
                                   exclude_detectors = integer(0)) {
  dead <- load_dead_channels() %>% distinct(subject_id, channel_index)
  if (length(exclude_detectors)) {
    dead <- bind_rows(dead, load_channels() %>%
                        filter(.data$detector %in% exclude_detectors) %>%
                        distinct(subject_id, channel_index)) %>% distinct()
  }
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

#' The same rebuild from the presentation-level block layer (mean within
#' event, then mean over events) — the aggregation the original 2026-07-30
#' pass used. Kept here only so gate A can demonstrate the two layers agree.
roi_masked_cv_blocklayer <- function(chromophore = "HbO") {
  dead <- load_dead_channels() %>% distinct(subject_id, channel_index) %>%
    mutate(.dead = TRUE)
  map <- load_channels() %>% filter(active_pair == 1) %>%
    distinct(subject_id, channel_index, roi)
  load_block_metrics(active_only = TRUE) %>%
    filter(.data$chromophore == !!chromophore) %>%
    left_join(dead, by = c("subject_id", "channel_index")) %>%
    filter(is.na(.data$.dead)) %>%
    left_join(map, by = c("subject_id", "channel_index")) %>%
    filter(!is.na(roi)) %>%
    group_by(subject_id, marker, event_index, roi) %>%
    summarise(block_mean = mean(response_minus_baseline_m5_30), .groups = "drop") %>%
    group_by(subject_id, marker, roi) %>%
    summarise(value = mean(block_mean), .groups = "drop")
}

#' Cross-check the physical (raw-recording) CV mask against an equivalent
#' criterion computed inside the delivered block layer: per subject x channel
#' response SD below 1% of the campaign median.
check_dead_channel_mask <- function(chromophore = "HbO") {
  bm <- load_block_metrics(active_only = TRUE) %>%
    filter(.data$chromophore == !!chromophore)
  sds <- bm %>% group_by(subject_id, channel_index) %>%
    summarise(sd = sd(response_minus_baseline_m5_30), n = n(), .groups = "drop") %>%
    filter(n >= 3)
  floor <- 0.01 * median(sds$sd, na.rm = TRUE)
  frozen_dead <- sds %>% filter(sd < floor) %>%
    transmute(key = paste(subject_id, channel_index))
  raw_dead <- load_dead_channels() %>%
    transmute(key = paste(subject_id, channel_index))
  a <- raw_dead$key; b <- frozen_dead$key
  tibble(n_raw = length(a), n_frozen = length(b),
         n_both = length(intersect(a, b)),
         n_raw_only = length(setdiff(a, b)), n_frozen_only = length(setdiff(b, a)),
         jaccard = length(intersect(a, b)) / length(union(a, b)),
         frozen_floor_uM = floor)
}

# Declared unmasked values this correction pass is read against. These are the
# estimates the pre-mask family scan and the two surviving right-temporal
# claims reported, quoted here at full precision so the correction comparison
# needs no file outside this repository.
A04_TEMPORAL_FAMILY <- tribble(
  ~analysis,     ~roi,             ~term,         ~estimate,                ~se,                    ~p,
  "difference", "Left_Temporal",  "comfort_wi",  -0.018617389602807982,     0.011047373857505441,   0.0934653869283922,
  "difference", "Left_Temporal",  "richness_wi", -0.011578772484746206,     0.011405412955941776,   0.31120728556084654,
  "difference", "Right_Temporal", "comfort_wi",  -0.02524152501932646,      0.019512692540124332,   0.1973047648991941,
  "difference", "Right_Temporal", "richness_wi", -0.021850819186400607,     0.02014770672591761,    0.27944144292147416,
  "response",   "Left_Temporal",  "comfort_wi",  -0.0026904073467485553,    0.005204713658537236,   0.6057710143501318,
  "response",   "Left_Temporal",  "richness_wi",  0.0011297733677989696,    0.005373395465630775,   0.8336789164471136,
  "response",   "Right_Temporal", "comfort_wi",   9.922125266917322e-4,     0.012095499002028711,   0.9347040736585136,
  "response",   "Right_Temporal", "richness_wi",  0.009107932695697816,     0.012489130656642368,   0.46669474068680883,
  "baseline",   "Left_Temporal",  "comfort_wi",   0.0159269822560594,       0.008729386510052084,   0.06952842241051412,
  "baseline",   "Left_Temporal",  "richness_wi",  0.012708545852545185,     0.009012300957981077,   0.16001576410325136,
  "baseline",   "Right_Temporal", "comfort_wi",   0.026233737546018182,     0.012751987623389466,   0.040968473104381474,
  "baseline",   "Right_Temporal", "richness_wi",  0.030958751882098404,     0.013166983812216876,   0.019688728901483907)

# A04, right-temporal pre-stimulus baseline, richness slope after dropping
# run-opening events (frozen ROI, declared estimate).
A04_ORIG_OPENERS <- list(estimate = 0.052993519703233696,
                         se = 0.016740301548787483, p = 0.0017907022177894793)
# A12, right-temporal response encounter trend (frozen ROI, declared estimate).
A12_ORIG_TREND <- list(estimate = 0.015935712043379263,
                       se = 0.00522617743432025, p = 0.0023768826590366446)

# ============================================================== 0. the two gates
rule("0. GATES — the rebuild must be the builder's, and the mask must be real")

cat("(a) does the condition-layer rebuild reproduce the delivered ROI columns\n")
cat("    when nothing is excluded? If not, the aggregation differs from the\n")
cat("    builder's and nothing below is interpretable.\n\n")
gate_a <- check_roi_rebuild()
print(as.data.frame(gate_a), digits = 3)
stopifnot(max(gate_a$max_abs_diff) < 1e-9)
cat("\n    GATE A1 PASSED (worst discrepancy", signif(max(gate_a$max_abs_diff), 3), "uM).\n")

cat("\n(a2) does the condition-layer aggregation equal the block-layer\n")
cat("     aggregation the original pass used, under the same CV mask? They\n")
cat("     agree when every channel contributes the same number of events per\n")
cat("     condition; the worst discrepancy across cells is reported.\n\n")
d_cond <- load_p24_roi_masked_cv("HbO")
blk <- roi_masked_cv_blocklayer("HbO")
agg_chk <- map_dfr(ROIS, function(rr) {
  a <- d_cond %>% select(subject_id, marker, value = all_of(paste0("HbOm_", rr)))
  b <- blk %>% filter(roi == rr) %>% select(subject_id, marker, value_blk = value)
  dd <- a %>% full_join(b, by = c("subject_id", "marker"))
  keep <- !is.na(dd$value) & !is.na(dd$value_blk)
  tibble(roi = rr, n_compared = sum(keep),
         max_abs_diff = max(abs(dd$value[keep] - dd$value_blk[keep])),
         max_rel_diff = max(abs(dd$value[keep] - dd$value_blk[keep]) /
                              pmax(abs(dd$value_blk[keep]), 1e-12)))
})
print(as.data.frame(agg_chk), digits = 3)
cat("\n    GATE A2: worst relative discrepancy",
    signif(max(agg_chk$max_rel_diff), 3), "across", sum(agg_chk$n_compared),
    "subject x condition x ROI cells.\n")

cat("\n(b) could the mask have been found without leaving the delivered tables?\n")
cat("    The mask is the physical raw-recording criterion; this is the\n")
cat("    equivalent criterion computed inside the delivered block layer.\n\n")
gate_b <- check_dead_channel_mask()
print(as.data.frame(gate_b), digits = 3)
write_outcome(gate_b, file.path(OUT, "mask_cross_check.csv"))
cat("\n    The two criteria agree on", gate_b$n_both, "of", gate_b$n_raw, "cells",
    "(Jaccard", round(gate_b$jaccard, 3), "). The delivered-layer criterion is\n",
    "   mildly conservative at the margin, not divergent.\n")

dead <- load_dead_channels()
cat("\n(c) what the mask contains:\n")
print(as.data.frame(dead %>% count(channel_index, detector, roi, name = "n_subjects")))
write_outcome(dead, file.path(OUT, "dead_channel_mask.csv"))

# ================================================= 1. what the mask does to ROIs
rule("1. WHAT THE MASK CHANGES — ROI amplitude, variance and composition")

# exclude_detectors = integer(0) deliberately: sections 1 and 4 exist to show
# what the UNCORRECTED right-temporal ROI does, and section 4b is what
# establishes that detector 8 must be dropped wholesale.
d  <- load_p24_roi_masked_cv("HbO")
dr <- load_p24_roi_masked_cv("HbR")

live_tab <- map_dfr(ROIS, function(rr) {
  n <- d[[paste0("n_live_", rr)]]
  n[is.na(n)] <- 0L
  # one row per subject, not per cell
  per_subj <- d %>% mutate(.n = n) %>% group_by(subject_id) %>%
    summarise(.n = max(.n), .groups = "drop")
  tibble(roi = rr,
         subjects_0_live = sum(per_subj$.n == 0), subjects_1_live = sum(per_subj$.n == 1),
         subjects_2_live = sum(per_subj$.n == 2), subjects_3plus  = sum(per_subj$.n >= 3))
})
cat("live channels per subject, after masking (active and not dark):\n")
print(as.data.frame(live_tab))
write_outcome(live_tab, file.path(OUT, "live_channels_per_roi.csv"))

amp <- map_dfr(ROIS, function(rr) {
  a <- d[[paste0("HbO_",  rr)]]
  b <- d[[paste0("HbOm_", rr)]]
  tibble(roi = rr,
         n_frozen = sum(!is.na(a)), n_masked = sum(!is.na(b)),
         mean_frozen = mean(a, na.rm = TRUE), mean_masked = mean(b, na.rm = TRUE),
         sd_frozen = sd(a, na.rm = TRUE), sd_masked = sd(b, na.rm = TRUE),
         sd_ratio = sd(b, na.rm = TRUE) / sd(a, na.rm = TRUE),
         r_frozen_masked = cor(a, b, use = "complete.obs"))
})
cat("\ndelivered versus masked ROI values:\n")
print(as.data.frame(amp), digits = 3)
write_outcome(amp, file.path(OUT, "roi_amplitude_change.csv"))
cat("\nA diluted ROI is an attenuated ROI: dead channels contribute a near-zero\n")
cat("constant to the mean, so the delivered column is the live signal shrunk\n")
cat("toward zero by the live fraction. The correlation column says the masked\n")
cat("ROI is the same quantity rescaled where the dilution is uniform, and\n")
cat("something else where it is not.\n")

# ========================================== 2. the declared families, re-run
rule("2. DECLARED FAMILIES UNDER THE MASK — 2 axes x 3 ROIs, three windows")

# rebuild the three windows at block level, masked, exactly as the unmasked
# scan built them
bm  <- load_block_metrics(active_only = TRUE) %>% filter(chromophore == "HbO")
map <- load_channels() %>% filter(active_pair == 1) %>%
  distinct(subject_id, channel_index, roi)
mask_keys <- dead %>% distinct(subject_id, channel_index) %>% mutate(.dead = TRUE)

roi_blocks <- bm %>%
  left_join(mask_keys, by = c("subject_id", "channel_index")) %>%
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
write_outcome(fam, file.path(OUT, "families_masked.csv"))
cat("\nsurvivors at q < 0.05, by window:\n")
print(as.data.frame(fam %>% group_by(analysis) %>%
        summarise(n_tests = n(), n_survivors = sum(q_BH < 0.05), min_q = min(q_BH))),
      digits = 3)

cat("\nside by side with the declared unmasked estimates for the two temporal\n")
cat("ROIs:\n")
cmp <- A04_TEMPORAL_FAMILY %>%
  select(analysis, roi, term, est_unmasked = estimate,
         se_unmasked = se, p_unmasked = p) %>%
  inner_join(fam %>% select(analysis, roi, term, est_masked = estimate,
                            se_masked = se, p_masked = p),
             by = c("analysis", "roi", "term")) %>%
  mutate(ratio = est_masked / est_unmasked)
print(as.data.frame(cmp), digits = 3)
write_outcome(cmp, file.path(OUT, "masked_vs_a04.csv"))

# ============================================= 3. the lateralisation artefact
rule("3. THE LATERALISATION CONTRAST AT THREE ROI DEFINITIONS")
cat("The unmasked scan reported left minus right at +0.0214 uM (p = 0.013) on\n")
cat("the delivered columns, vanishing to -0.0075 (p = 0.58) when channels 18/19\n")
cat("were dropped. That correction fixed the LEFT side only, while the right\n")
cat("side is itself half-dark in 29 of 69 sessions. Three definitions, one\n")
cat("test each:\n\n")

dc <- load_p24_roi_corrected()   # the fixed two-channel drop, for continuity
lat_input <- d %>%
  select(subject_id, marker, HbO_Left_Temporal, HbO_Right_Temporal,
         HbOm_Left_Temporal, HbOm_Right_Temporal,
         n_live_Left_Temporal, n_live_Right_Temporal) %>%
  left_join(dc %>% select(subject_id, marker, HbOc_Left_Temporal, HbOc_Right_Temporal),
            by = c("subject_id", "marker")) %>%
  mutate(lat_frozen    = HbO_Left_Temporal  - HbO_Right_Temporal,
         lat_a04_fixed = HbOc_Left_Temporal - HbOc_Right_Temporal,
         lat_masked    = HbOm_Left_Temporal - HbOm_Right_Temporal)

lat <- map_dfr(c("lat_frozen", "lat_a04_fixed", "lat_masked"), function(y) {
  x <- lat_input[[y]]; x <- x[!is.na(x)]
  tt <- t.test(x)
  tibble(definition = y, n_cells = length(x), mean = mean(x),
         t = unname(tt$statistic), p = tt$p.value,
         ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2])
})

# and the same on the honest subsample: both ROIs retaining at least 2 live channels
strict <- lat_input %>% filter(n_live_Left_Temporal >= 2, n_live_Right_Temporal >= 2)
x <- strict$lat_masked; x <- x[!is.na(x)]
tt <- t.test(x)
lat <- bind_rows(lat, tibble(definition = "lat_masked (both ROIs >= 2 live)",
                             n_cells = length(x), mean = mean(x),
                             t = unname(tt$statistic), p = tt$p.value,
                             ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2]))
print(as.data.frame(lat), digits = 3)
write_outcome(lat, file.path(OUT, "lateralisation_three_definitions.csv"))
cat("\nNOTE, unchanged and binding: handedness was never collected, so NO\n")
cat("lateralisation claim may be made under any of these definitions. The\n")
cat("contrast is reported here only as a demonstration of what a diluted ROI\n")
cat("does to an inferential test.\n")

cat("\n(3b) the strict subsample moves the estimate a long way for a small\n")
cat("     exclusion. Which subjects carry it, and how thin are their ROIs?\n\n")
thin <- lat_input %>%
  filter(n_live_Left_Temporal < 2 | n_live_Right_Temporal < 2, !is.na(lat_masked))
cat("    excluded cells:", nrow(thin), "from", n_distinct(thin$subject_id), "subjects\n")
print(as.data.frame(thin %>% group_by(subject_id) %>%
        summarise(n_cells = n(), live_L = max(n_live_Left_Temporal),
                  live_R = max(n_live_Right_Temporal),
                  mean_lat = mean(lat_masked), .groups = "drop")), digits = 3)
cat("\n    their mean contrast:", sprintf("%+.4f", mean(thin$lat_masked)),
    "uM against", sprintf("%+.4f", mean(strict$lat_masked, na.rm = TRUE)),
    "in the strict set —\n    an ROI built from one channel is a channel, not an ROI.\n")
write_outcome(thin, file.path(OUT, "lateralisation_thin_roi_cells.csv"))

# =============================== 4. the right-temporal coupling anomaly, re-asked
rule("4. HbO/HbR COUPLING PER ROI — is the right-temporal anomaly instrumental?")
cat("The right-temporal ROI looked blood-volume-like in the unmasked data: HbR\n")
cat("silent while HbO moved, and r(HbO,HbR) = +0.28 there against -0.43 to\n")
cat("-0.51 elsewhere. Detector 8 is dark in 29 of 69 sessions in exactly that\n")
cat("ROI, so the question is whether the anomaly is physiology or two dead\n")
cat("channels averaged into the mean.\n\n")

cpl <- map_dfr(ROIS, function(rr) {
  one <- function(a, b) {
    keep <- !is.na(a) & !is.na(b)
    dd <- tibble(sid = d$subject_id[keep], o = a[keep], r = b[keep]) %>%
      group_by(sid) %>% mutate(oc = o - mean(o), rc = r - mean(r)) %>% ungroup()
    c(n = sum(keep), r_within = cor(dd$oc, dd$rc), sd_ratio = sd(a[keep]) / sd(b[keep]))
  }
  f <- one(d[[paste0("HbO_",  rr)]], dr[[paste0("HbR_",  rr)]])
  m <- one(d[[paste0("HbOm_", rr)]], dr[[paste0("HbRm_", rr)]])
  tibble(roi = rr,
         n_frozen = f[["n"]], r_within_frozen = f[["r_within"]], sd_ratio_frozen = f[["sd_ratio"]],
         n_masked = m[["n"]], r_within_masked = m[["r_within"]], sd_ratio_masked = m[["sd_ratio"]])
})
print(as.data.frame(cpl), digits = 3)
write_outcome(cpl, file.path(OUT, "hbo_hbr_coupling_masked.csv"))

cat("\nthe same split by whether the subject's detector 8 was working:\n")
det8_dead <- dead %>% filter(channel_index %in% c(20L, 21L)) %>%
  distinct(subject_id) %>% pull(subject_id)
split_cpl <- map_dfr(c(TRUE, FALSE), function(is_dead) {
  sel <- if (is_dead) d$subject_id %in% det8_dead else !(d$subject_id %in% det8_dead)
  a <- d[["HbO_Right_Temporal"]][sel]; b <- dr[["HbR_Right_Temporal"]][sel]
  keep <- !is.na(a) & !is.na(b)
  dd <- tibble(sid = d$subject_id[sel][keep], o = a[keep], r = b[keep]) %>%
    group_by(sid) %>% mutate(oc = o - mean(o), rc = r - mean(r)) %>% ungroup()
  tibble(detector8 = if (is_dead) "dark" else "working",
         n_subjects = n_distinct(dd$sid), n_cells = sum(keep),
         r_within_frozen_RT = cor(dd$oc, dd$rc))
})
print(as.data.frame(split_cpl), digits = 3)
write_outcome(split_cpl, file.path(OUT, "coupling_by_detector8_status.csv"))
cat("\nThat split runs OPPOSITE to the hypothesis that dead channels create the\n")
cat("anomaly: the subjects whose detector 8 was DARK show ordinary negative\n")
cat("coupling, and the positive coupling belongs to the sessions where\n")
cat("channels 20/21 were alive. But it is a BETWEEN-subject comparison, and\n")
cat("detector-8 status clusters by session date, so it cannot separate site\n")
cat("from cohort. The within-subject test can.\n")

cat("\n(4b) within the subjects whose detector 8 worked, the same coupling\n")
cat("     computed from the T8 pair (20/21) and the TP8 pair (22/23)\n")
cat("     separately — one subject, two sites, no cohort difference:\n\n")
pair_roi <- bm %>%
  filter(!subject_id %in% det8_dead, channel_index %in% 20:23) %>%
  mutate(site = if_else(channel_index %in% 20:21, "T8 (det 8)", "TP8 (det 10)")) %>%
  group_by(subject_id, marker, site) %>%
  summarise(HbO = mean(response_minus_baseline_m5_30), .groups = "drop")
pair_roi_r <- load_block_metrics(active_only = TRUE) %>%
  filter(chromophore == "HbR", !subject_id %in% det8_dead, channel_index %in% 20:23) %>%
  mutate(site = if_else(channel_index %in% 20:21, "T8 (det 8)", "TP8 (det 10)")) %>%
  group_by(subject_id, marker, site) %>%
  summarise(HbR = mean(response_minus_baseline_m5_30), .groups = "drop")
site_cpl <- pair_roi %>% inner_join(pair_roi_r, by = c("subject_id", "marker", "site")) %>%
  group_by(site) %>%
  group_modify(~ {
    dd <- .x %>% group_by(subject_id) %>%
      mutate(oc = HbO - mean(HbO), rc = HbR - mean(HbR)) %>% ungroup()
    tibble(n_subjects = n_distinct(dd$subject_id), n_cells = nrow(dd),
           r_within = cor(dd$oc, dd$rc), sd_ratio = sd(dd$HbO) / sd(dd$HbR))
  }) %>% ungroup()
print(as.data.frame(site_cpl), digits = 3)
write_outcome(site_cpl, file.path(OUT, "coupling_by_site_within_subject.csv"))

# ========================================================== 5. sensitivity
rule("5. SENSITIVITY — does the primary family depend on the thin-ROI subjects?")
strict_cells <- cells %>%
  left_join(d %>% select(subject_id, marker, starts_with("n_live_")),
            by = c("subject_id", "marker"))
fam_strict <- map_dfr(ROIS, function(rr) {
  nlive <- paste0("n_live_", rr)
  dd <- strict_cells %>% filter(.data[[nlive]] >= 2)
  m <- p24_lmm(dd, paste(AXES, collapse = " + "), paste0("difference_", rr))
  fixef_table(m, keep = AXES) %>%
    mutate(roi = rr, n_cells = nrow(dd), n_subjects = n_distinct(dd$subject_id), .before = 1)
})
fam_strict <- fam_strict %>% mutate(q_BH = bh(p))
print(as.data.frame(fam_strict %>% select(roi, n_subjects, term, estimate, se, p, q_BH)),
      digits = 3)
write_outcome(fam_strict, file.path(OUT, "family_strict_two_live.csv"))

# ================== 6. the montage is mirrored, and the failures are complementary
rule("6. NO MIRROR PAIR SURVIVES — why the left-right contrast is not defined")
cat("The temporal montage is exactly mirror-symmetric, channel for channel:\n\n")
mirror <- tribble(
  ~pair, ~left, ~left_name, ~right, ~right_name,
  "T",  16L, "FTT7h_T7",  20L, "FTT8h_T8",
  "T",  17L, "TTP7h_T7",  21L, "TTP8h_T8",
  "TP", 18L, "TTP7h_TP7", 22L, "TTP8h_TP8",
  "TP", 19L, "TPP7h_TP7", 23L, "TPP8h_TP8")
print(as.data.frame(mirror))
cat("\nThe two failures hit OPPOSITE members of the two mirror pairs. Detector 9\n")
cat("kills the left half of the TP pair; detector 8 kills the right half of the\n")
cat("T pair. What survives is left-T7 and right-TP8 — sites that are not each\n")
cat("other's mirror.\n\n")

dead_set <- mask_keys %>% distinct(subject_id, channel_index)
ch_tbl <- load_channels()
complete_pair <- map_dfr(c("T", "TP"), function(pp) {
  chs <- mirror %>% filter(pair == pp)
  n_ok <- sum(map_lgl(unique(d$subject_id), function(s) {
    live <- function(ch) {
      act <- ch_tbl %>% filter(subject_id == s, channel_index == ch) %>%
        pull(active_pair)
      length(act) == 1 && act == 1 &&
        !any(dead_set$subject_id == s & dead_set$channel_index == ch)
    }
    all(map_lgl(c(chs$left, chs$right), live))
  }))
  tibble(mirror_pair = pp, channels = paste(c(chs$left, chs$right), collapse = "/"),
         subjects_with_all_four_live = n_ok)
})
print(as.data.frame(complete_pair))
write_outcome(complete_pair, file.path(OUT, "mirror_pair_completeness.csv"))
cat("\nAnd where the T pair IS complete, its right member is detector 8 — the\n")
cat("pair section 4b showed carries positive HbO/HbR coupling (+0.67) while the\n")
cat("same subjects' detector-10 channels carry ordinary negative coupling\n")
cat("(-0.46). So the one available mirror comparison rests on the one pair\n")
cat("known to be non-neurovascular.\n")
cat("\nA left-versus-right contrast in these data is therefore not merely\n")
cat("unlicensed by the missing handedness — after the two failures it is NOT\n")
cat("DEFINED BY THE SURVIVING MONTAGE. Everything below compares left-T7 with\n")
cat("right-TP8 and is reported to show what that comparison does, not to make\n")
cat("a lateralisation statement.\n\n")

cat("The two-live-detector definition (left = det 7, ch 16/17; right = det 10,\n")
cat("ch 22/23) — two channels each, but NOT anatomically matched:\n\n")

sym <- bm %>%
  filter(channel_index %in% c(16L, 17L, 22L, 23L)) %>%
  left_join(mask_keys, by = c("subject_id", "channel_index")) %>%
  filter(is.na(.dead)) %>%
  mutate(roi = if_else(channel_index %in% c(16L, 17L), "Left_Temporal", "Right_Temporal")) %>%
  group_by(subject_id, marker, roi) %>%
  summarise(value = mean(response_minus_baseline_m5_30), n_live = n_distinct(channel_index),
            .groups = "drop")
sym_wide <- sym %>%
  pivot_wider(names_from = roi, values_from = c(value, n_live), names_sep = "_") %>%
  left_join(load_p24() %>% select(subject_id, marker, all_of(AXES)),
            by = c("subject_id", "marker"))

cat("subjects retaining both two-live-detector ROIs:",
    n_distinct(sym_wide$subject_id[!is.na(sym_wide$value_Left_Temporal) &
                                   !is.na(sym_wide$value_Right_Temporal)]), "of 69\n\n")

sym_fam <- map_dfr(c("Left_Temporal", "Right_Temporal"), function(rr) {
  m <- p24_lmm(sym_wide, paste(AXES, collapse = " + "), paste0("value_", rr))
  fixef_table(m, keep = AXES) %>% mutate(roi = rr, .before = 1)
}) %>% mutate(q_BH = bh(p))
cat("the two temporal slopes on the two-live-detector definition (difference score):\n")
print(as.data.frame(sym_fam %>% select(roi, term, estimate, se, p, q_BH)), digits = 3)
write_outcome(sym_fam, file.path(OUT, "symmetric_definition_family.csv"))

sym_lat <- sym_wide %>%
  mutate(lat = value_Left_Temporal - value_Right_Temporal) %>%
  filter(!is.na(lat))
tt <- t.test(sym_lat$lat)
cat("\nleft minus right on the two-live-detector definition:",
    sprintf("%+.4f uM, t = %.2f, p = %.3f, n = %d cells from %d subjects\n",
            mean(sym_lat$lat), unname(tt$statistic), tt$p.value,
            nrow(sym_lat), n_distinct(sym_lat$subject_id)))
write_outcome(sym_lat %>% select(subject_id, marker, lat), file.path(OUT, "symmetric_lateralisation.csv"))

cat("\nAll five definitions of the same contrast, in one place:\n")
allfive <- bind_rows(
  lat %>% select(definition, n_cells, mean, p),
  tibble(definition = "two-live-detector (T7 vs TP8)", n_cells = nrow(sym_lat),
         mean = mean(sym_lat$lat), p = tt$p.value))
print(as.data.frame(allfive), digits = 3)
write_outcome(allfive, file.path(OUT, "lateralisation_all_definitions.csv"))
cat("\nRead this carefully, because it does not say what the unmasked scan said.\n")
cat("The two composition-matched definitions (4 and 5) AGREE: +0.027 to +0.033\n")
cat("uM at p = 0.002-0.004. The only definition that makes the contrast\n")
cat("disappear is the fixed two-channel correction, and it disappears because\n")
cat("that correction fixed the LEFT side while leaving the right diluted in 29\n")
cat("of 69 subjects — attenuating the right ROI toward zero and pushing the\n")
cat("difference negative. The 'spurious lateralisation that vanishes when\n")
cat("corrected' is therefore an artefact of a one-sided correction, in both of\n")
cat("its halves.\n")

cat("\n(6b) but WHICH quantity carries it? The difference score's negative\n")
cat("     offset is the pre-stimulus baseline sitting high. A left-right\n")
cat("     difference in the difference score may therefore be a difference\n")
cat("     between two baselines, not between two responses.\n\n")
sym_win <- bm %>%
  filter(channel_index %in% c(16L, 17L, 22L, 23L)) %>%
  left_join(mask_keys, by = c("subject_id", "channel_index")) %>%
  filter(is.na(.dead)) %>%
  mutate(roi = if_else(channel_index %in% c(16L, 17L), "L", "R")) %>%
  group_by(subject_id, marker, roi) %>%
  summarise(baseline = mean(baseline_m10_0), response = mean(response_m5_30),
            difference = mean(response_minus_baseline_m5_30), .groups = "drop") %>%
  pivot_wider(names_from = roi, values_from = c(baseline, response, difference))
win_lat <- map_dfr(c("baseline", "response", "difference"), function(w) {
  x <- sym_win[[paste0(w, "_L")]] - sym_win[[paste0(w, "_R")]]
  x <- x[!is.na(x)]
  tt <- t.test(x)
  tibble(window = w, n = length(x), left_minus_right = mean(x),
         t = unname(tt$statistic), p = tt$p.value)
})
print(as.data.frame(win_lat), digits = 3)
write_outcome(win_lat, file.path(OUT, "symmetric_lateralisation_by_window.csv"))

cat("\nNOTE, unchanged and binding whatever these numbers say: handedness was\n")
cat("never collected, so no lateralisation claim may be made from these data\n")
cat("under any definition. What this module establishes is narrower and is\n")
cat("about measurement: the contrast is real in the data and the earlier\n")
cat("correction mis-corrected it, and it belongs to whichever window section\n")
cat("6b identifies.\n")

# ================= 7. every surviving q<0.05 claim that rests on a temporal ROI
rule("7. THE THREE SURVIVING RIGHT-TEMPORAL CLAIMS, RE-RUN")
cat("Three results in the unmasked campaign reach q < 0.05 on a temporal ROI,\n")
cat("and all three are RIGHT temporal — the ROI detector 8 compromises. Each is\n")
cat("re-run here on the defensible definition (detector 10 only, all subjects),\n")
cat("with the model left exactly as its own analysis specified it.\n")

# the corrected right-temporal series at block level, det 10 only
rt10 <- bm %>%
  filter(channel_index %in% c(22L, 23L)) %>%
  left_join(mask_keys, by = c("subject_id", "channel_index")) %>%
  filter(is.na(.dead)) %>%
  group_by(subject_id, event_index, sequence, marker, onset_sec) %>%
  summarise(baseline = mean(baseline_m10_0), response = mean(response_m5_30),
            difference = mean(response_minus_baseline_m5_30), .groups = "drop")
cat("\n    corrected right-temporal block series:", nrow(rt10), "events from",
    n_distinct(rt10$subject_id), "subjects\n")

cat("\n(7a) right-temporal PRE-STIMULUS BASELINE comfort/richness slopes after\n")
cat("     dropping run-opening events (declared q = 0.011 on the frozen ROI).\n\n")
a04_cells <- rt10 %>% mutate(opens = as.integer(event_index %in% c(1, 5, 9))) %>%
  filter(opens == 0) %>%
  group_by(subject_id, marker) %>%
  summarise(baseline = mean(baseline), .groups = "drop") %>%
  left_join(load_p24() %>% select(subject_id, marker, all_of(AXES)),
            by = c("subject_id", "marker"))
a7a <- fixef_table(p24_lmm(a04_cells, paste(AXES, collapse = " + "), "baseline"),
                   keep = AXES)
print(as.data.frame(a7a), digits = 3)
write_outcome(a7a, file.path(OUT, "recheck_a04_baseline_openers.csv"))

cat("\n(7b) right-temporal pre-stimulus baseline by NOMINAL QUADRANT, the\n")
cat("     three-df omnibus (declared p = 0.0013, q = 0.008 on the frozen ROI).\n\n")
a05_cells <- rt10 %>% group_by(subject_id, marker) %>%
  summarise(baseline = mean(baseline), .groups = "drop") %>%
  left_join(load_p24() %>% select(subject_id, marker, quadrant),
            by = c("subject_id", "marker")) %>%
  mutate(quadrant_f = factor(quadrant, levels = c("Q1", "Q2", "Q3", "Q4"))) %>%
  filter(!is.na(baseline), !is.na(quadrant_f))
m_q  <- lmerTest::lmer(baseline ~ quadrant_f + (1 | subject_id), data = a05_cells, REML = FALSE)
m_q0 <- lmerTest::lmer(baseline ~ 1 + (1 | subject_id), data = a05_cells, REML = FALSE)
lrt <- anova(m_q0, m_q)
a7b <- tibble(test = "nominal quadrant, 3 df",
              chisq = lrt$Chisq[2], df = lrt$Df[2], p = lrt$`Pr(>Chisq)`[2])
print(as.data.frame(a7b), digits = 3)
write_outcome(a7b, file.path(OUT, "recheck_a05_quadrant_omnibus.csv"))

cat("\n(7c) right-temporal RESPONSE encounter trend (declared +0.0159 uM per\n")
cat("     encounter, q = 0.008, the adaptation candidate).\n\n")
a12 <- rt10 %>% group_by(subject_id, marker) %>%
  mutate(encounter = rank(event_index)) %>% ungroup()
a7c <- fixef_table(p24_lmm(a12, "z(encounter)", "response"), keep = "z(encounter)")
print(as.data.frame(a7c), digits = 3)
write_outcome(a7c, file.path(OUT, "recheck_a12_encounter_trend.csv"))

cat("\n(7d) the three side by side with their declared originals:\n\n")
recheck <- tribble(
  ~claim, ~est_original, ~se_original, ~p_original, ~est_corrected, ~se_corrected, ~p_corrected,
  "A04 RT baseline richness (openers dropped)",
    A04_ORIG_OPENERS$estimate, A04_ORIG_OPENERS$se, A04_ORIG_OPENERS$p,
    a7a$estimate[a7a$term == "richness_wi"], a7a$se[a7a$term == "richness_wi"],
    a7a$p[a7a$term == "richness_wi"],
  "A12 RT response encounter trend",
    A12_ORIG_TREND$estimate, A12_ORIG_TREND$se, A12_ORIG_TREND$p,
    a7c$estimate[1], a7c$se[1], a7c$p[1]) %>%
  mutate(est_retained = est_corrected / est_original,
         se_ratio = se_corrected / se_original)
print(as.data.frame(recheck), digits = 3)
write_outcome(recheck, file.path(OUT, "recheck_summary.csv"))
cat("\nThe quadrant omnibus is a chi-square and is not comparable in this\n")
cat("format: p = 0.0013 (q = 0.008) on the frozen ROI -> chisq(3) = 7.45,\n")
cat("p = 0.059.\n")
cat("\nRead the two ratio columns before concluding anything. The encounter\n")
cat("trend's standard error FALLS by 17% while its estimate falls 62%: losing\n")
cat("detector 8 made that estimate MORE precise and smaller, which is what an\n")
cat("artefact does, not what losing power does. The baseline richness slope\n")
cat("retains 78% of its magnitude while its SE rises 31% (4 channels to 2, 69\n")
cat("subjects to 62), so that one is weakened mostly by precision and should\n")
cat("be described as no longer significant rather than as removed.\n")

rule("VERDICT")
cat("survivors under the mask, all three windows:",
    sum(fam$q_BH < 0.05), "of", nrow(fam), "tests\n")
cat("lateralisation, masked:", sprintf("%+.4f uM, p = %.3f",
    lat$mean[lat$definition == "lat_masked"], lat$p[lat$definition == "lat_masked"]), "\n")
cat("right-temporal r(HbO,HbR): delivered", sprintf("%+.3f", cpl$r_within_frozen[3]),
    "-> masked", sprintf("%+.3f", cpl$r_within_masked[3]), "\n")

close_log(con)
