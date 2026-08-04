# mod14_detrended.R — Baseline-coupling anatomy: fast event alignment or slow session state?
#
# The strangest finding in the campaign is still unexplained: the PRE-stimulus
# baseline — the negative control — carries the strongest, most consistent
# positive comfort slope (the current clip's later-rated comfort predicts the
# baseline BEFORE it plays; not carry-over, not slot, not anticipation). Two
# mechanisms remain:
#   (a) SLOW: within-person session-state co-fluctuation — mood/arousal drifts
#       across sequences and moves both the haemodynamic level and how clips
#       are later rated.
#   (b) FAST: the association is aligned to the event itself (expectation, or
#       an artefact of how post-session ratings were reconstructed).
# The two are separable because the ratings have no sequence-level variance
# (every sequence contains all four markers): detrending each person's 12-event
# physiology at sequence level leaves only the fast, event-aligned arm.
# Declared reading rule: a fast-residual slope at the unde-trended magnitude
# means event-aligned coupling; a collapse means slow state co-fluctuation.
#
# Declared family: comfort_wi slope on the sequence-detrended series, 2 windows
# (baseline; response) x 3 ROIs = 6 tests, BH. Richness, unde-trended
# block-level slopes, the lead-lag scan and the slow between-person arm are
# descriptive (the person-side module owns the formal between-person families;
# the unde-trended inferential anchors live with the condition-level families).
#
# Outputs (output/analysis_14_baseline_coupling_anatomy/):
#   family_detrended_comfort_slopes.csv, undetrended_slopes_descriptive.csv,
#   fast_vs_raw_comparison_descriptive.csv, lead_lag_joint_descriptive.csv,
#   lead_lag_scan_descriptive.csv, slow_between_person_descriptive.csv,
#   richness_detrended_descriptive.csv
# Display anchor: main-text Table 2 and SI Table S1 quote the family minimum
# (min q 0.061 unmasked; the SI text additionally quotes 0.059 under the
# channel mask, which is a different variant not regenerated here).

if (!exists("load_p24")) source("code/load_data.R")

OUT <- file.path("output", "analysis_14_baseline_coupling_anatomy")
con <- open_log(file.path(OUT, "baseline_coupling_anatomy_report.txt"))

rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")

# ------------------------------------------------- 1. block-level ROI window series
# active channels -> ROI mean per event
rule("1. BLOCK-LEVEL ROI WINDOW SERIES")
map <- load_channels() %>% filter(active_pair == 1) %>%
  distinct(subject_id, channel_index, roi)
ev <- load_block_metrics(active_only = TRUE) %>%
  filter(chromophore == "HbO") %>%
  left_join(map, by = c("subject_id", "channel_index")) %>%
  filter(!is.na(roi)) %>%
  summarise(.by = c(subject_id, event_index, sequence, marker, roi),
            baseline = mean(baseline_m10_0), response = mean(response_m5_30))

design <- load_p24_block() %>%
  select(subject_id, event_index, comfort_wi, richness_wi, comfort_bw, richness_bw)
d <- ev %>% left_join(design, by = c("subject_id", "event_index"))
stopifnot(!any(is.na(d$comfort_wi)))
cat("events:", nrow(d), "; subjects:", n_distinct(d$subject_id), "\n")

# ------------------------------------------------- 2. sequence-level detrend
rule("2. SEQUENCE-LEVEL DETREND (person x sequence means removed)")
d <- d %>%
  group_by(subject_id, sequence, roi) %>%
  mutate(baseline_fast = baseline - mean(baseline),
         response_fast = response - mean(response)) %>%
  ungroup()
cat("fast residuals: SD baseline", round(sd(d$baseline_fast), 4),
    "uM vs raw", round(sd(d$baseline), 4),
    "; response", round(sd(d$response_fast), 4), "vs", round(sd(d$response), 4), "\n")

# ------------------------------------------------- 3. declared family: 2 x 3 = 6, BH
rule("3. DECLARED FAMILY — DETRENDED COMFORT SLOPES, 2 WINDOWS x 3 ROIS = 6")
family <- map_dfr(ROIS, function(rr) {
  dd <- d %>% filter(roi == rr)
  map_dfr(c("baseline_fast", "response_fast"), function(w) {
    m <- p24_lmm(dd, "comfort_wi", w)
    fixef_table(m, keep = "comfort_wi") %>%
      mutate(roi = rr, window = str_remove(w, "_fast"))
  })
}) %>% mutate(q_BH = bh(p), family_size = n())
print(as.data.frame(family %>%
      select(window, roi, estimate, se, t, p, q_BH)), digits = 3)
write_outcome(family, file.path(OUT, "family_detrended_comfort_slopes.csv"))
cat("\nsurviving q < 0.05:", sum(family$q_BH < 0.05), "of", nrow(family), "\n")
cat("family minimum:", min(family$q_BH), "\n")

# ------------------------------------------------- 4. descriptive: unde-trended block
# slopes (same unit as the family; the condition-level values are the
# inferential anchor and are quoted in the report, not re-tested here)
rule("4. DESCRIPTIVE — UNDE-TRENDED BLOCK-LEVEL SLOPES FOR LIKE-FOR-LIKE COMPARISON")
raw_slopes <- map_dfr(ROIS, function(rr) {
  dd <- d %>% filter(roi == rr)
  map_dfr(c("baseline", "response"), function(w) {
    m <- p24_lmm(dd, "comfort_wi", w)
    fixef_table(m, keep = "comfort_wi") %>%
      mutate(roi = rr, window = w)
  })
}) %>% mutate(note = "descriptive anchor comparison, no correction")
print(as.data.frame(raw_slopes %>% select(window, roi, estimate, p)), digits = 3)
write_outcome(raw_slopes, file.path(OUT, "undetrended_slopes_descriptive.csv"))

comp <- family %>%
  select(window, roi, fast_estimate = estimate, fast_p = p, fast_q = q_BH) %>%
  left_join(raw_slopes %>% select(window, roi, raw_estimate = estimate, raw_p = p),
            by = c("window", "roi")) %>%
  mutate(fast_share_of_raw = fast_estimate / raw_estimate)
cat("\nfast-residual slope as share of the unde-trended slope:\n")
print(as.data.frame(comp), digits = 3)
write_outcome(comp, file.path(OUT, "fast_vs_raw_comparison_descriptive.csv"))

# ------------------------------------------------- 5. descriptive: richness on detrended
rule("5. DESCRIPTIVE — RICHNESS ON THE DETRENDED SERIES")
rich <- map_dfr(ROIS, function(rr) {
  dd <- d %>% filter(roi == rr)
  map_dfr(c("baseline_fast", "response_fast"), function(w) {
    m <- p24_lmm(dd, "richness_wi", w)
    fixef_table(m, keep = "richness_wi") %>%
      mutate(roi = rr, window = str_remove(w, "_fast"))
  })
}) %>% mutate(note = "descriptive, no correction")
print(as.data.frame(rich %>% select(window, roi, estimate, p)), digits = 3)
write_outcome(rich, file.path(OUT, "richness_detrended_descriptive.csv"))

# ------------------------------------------------- 6. descriptive: lead-lag scan
rule("6. DESCRIPTIVE — LEAD-LAG SCAN (detrended baseline vs comfort of k-1, k, k+1)")
lag_scan <- map_dfr(ROIS, function(rr) {
  dd <- d %>% filter(roi == rr) %>% arrange(subject_id, event_index)
  map_dfr(c(-1, 0, 1), function(lag) {
    dd2 <- dd %>%
      group_by(subject_id) %>%
      mutate(cw_lag = if (lag == -1) dplyr::lag(comfort_wi)
             else if (lag == 1) dplyr::lead(comfort_wi) else comfort_wi) %>%
      ungroup() %>% filter(!is.na(cw_lag))
    m <- p24_lmm(dd2, "cw_lag", "baseline_fast")
    fixef_table(m, keep = "cw_lag") %>%
      mutate(roi = rr, lag = lag,
             definition = c("-1" = "previous event's rating",
                            "0" = "current event's rating",
                            "1" = "next event's rating")[as.character(lag)])
  })
}) %>% mutate(note = "descriptive, no correction")
print(as.data.frame(lag_scan %>% select(roi, lag, definition, estimate, p)), digits = 3)
write_outcome(lag_scan, file.path(OUT, "lead_lag_scan_descriptive.csv"))

# adjacent-lag predictors correlate, so the scan alone cannot apportion: joint model
rule("6b. DESCRIPTIVE — JOINT LAG MODEL (k-1, k, k+1 together)")
lag_joint <- map_dfr(ROIS, function(rr) {
  dd <- d %>% filter(roi == rr) %>% arrange(subject_id, event_index) %>%
    group_by(subject_id) %>%
    mutate(cw_prev = dplyr::lag(comfort_wi), cw_next = dplyr::lead(comfort_wi)) %>%
    ungroup() %>% filter(!is.na(cw_prev), !is.na(cw_next))
  m <- p24_lmm(dd, "cw_prev + comfort_wi + cw_next", "baseline_fast")
  fixef_table(m, keep = c("cw_prev", "comfort_wi", "cw_next")) %>% mutate(roi = rr)
}) %>% mutate(note = "descriptive, no correction")
print(as.data.frame(lag_joint %>% select(roi, term, estimate, p)), digits = 3)
write_outcome(lag_joint, file.path(OUT, "lead_lag_joint_descriptive.csv"))

# ------------------------------------------------- 7. descriptive: slow between-person arm
rule("7. DESCRIPTIVE — SLOW ARM: person baseline level/drift vs comfort_bw")
slow <- d %>%
  summarise(.by = c(subject_id, roi),
            level = mean(baseline),
            drift = coef(lm(baseline ~ event_index))[["event_index"]],
            comfort_bw = first(comfort_bw))
slow_cors <- map_dfr(ROIS, function(rr) {
  ss <- slow %>% filter(roi == rr)
  tibble(roi = rr,
         r_level_comfortbw = cor(ss$level, ss$comfort_bw),
         r_drift_comfortbw = cor(ss$drift, ss$comfort_bw))
}) %>% mutate(note = "descriptive, no correction")
print(as.data.frame(slow_cors), digits = 3)
write_outcome(slow_cors, file.path(OUT, "slow_between_person_descriptive.csv"))

close_log(con)
cat("\nmod14 done.\n")
