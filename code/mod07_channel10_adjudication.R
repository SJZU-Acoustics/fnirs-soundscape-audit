# mod07_channel10_adjudication.R — Adjudication of the spatial scan's sole
# survivor.
#
# Release port of the stage-2 exploration pass A07
# (exploration/analysis_07_focal_signal_adjudication).
#
# Target selection is finished before this module runs: module 06 found exactly
# one q < 0.05 result in the 42-test response-window family, HbO channel 10 x
# within-subject Comfort. This module cannot select a substitute channel or
# window; it asks whether the survivor is:
#   - temporally plausible (absent before/at onset, present after 5 s);
#   - physiologically convergent in HbR rather than HbO/HbT alone;
#   - stable across repetitions, stimulus sets and form;
#   - robust to the declared quality/sample switches and single subjects; and
#   - locally coherent across channels sharing either optode.
#
# The bootstrap resamples subjects 1,000 times. All other follow-ups are
# coefficient-stability or convergence diagnostics, not second discovery
# screens.
#
# Run from the repository root:  Rscript code/mod07_channel10_adjudication.R

source("code/load_data.R")

OUT <- file.path("output", "analysis_07_focal_signal_adjudication")
TARGET_CHANNEL <- 10L
TARGET_DV <- "mean_response_m5_30"
AXES <- c("comfort_wi", "richness_wi")

con <- open_log(file.path(OUT, "focal_signal_adjudication_report.txt"))
rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n",
                        strrep("=", 74), "\n", sep = "")

# Hard dependency/guard: the target must be module 06's one and only
# response-window survivor. This prevents a later rerun from silently
# following a new winner.
a06_file <- file.path("output", "analysis_06_channelwise_spatial_scan",
                      "channel_scan_response.csv")
if (!file.exists(a06_file)) {
  stop("Module 06 output not found at ", a06_file, "\n",
       "Run code/mod06_channel_scan.R first.", call. = FALSE)
}
a06 <- readr::read_csv(a06_file, show_col_types = FALSE)
a06_survivors <- a06 %>% filter(q_BH < 0.05)
stopifnot(
  nrow(a06_survivors) == 1L,
  a06_survivors$channel_index[[1]] == TARGET_CHANNEL,
  a06_survivors$term[[1]] == "comfort_wi"
)

cat("Module 07 — focal signal adjudication\n")
cat("Fixed target: channel 10, HbO response [5,30), within-subject Comfort.\n")
cat("Target inherited from module 06; no replacement selection is permitted.\n")

# ------------------------------------------------------ 1. temporal profile
rule("1. TEMPORAL PROFILE AT THE FIXED CHANNEL")

cc <- load_condition_analysis("HbO") %>%
  filter(channel_index == TARGET_CHANNEL) %>%
  mutate(first5_m0_5 = 6 * mean_stim_m0_30 - 5 * mean_response_m5_30)

WINDOWS <- c(
  baseline_m10_0 = "mean_baseline_m10_0",
  onset_m0_5 = "first5_m0_5",
  whole_stimulus_m0_30 = "mean_stim_m0_30",
  response_m5_30 = "mean_response_m5_30",
  difference_m5_30_minus_baseline =
    "mean_response_minus_baseline_m5_30"
)

window_profile <- map_dfr(names(WINDOWS), function(win) {
  y <- unname(WINDOWS[[win]])
  m <- p24_lmm(cc, paste(AXES, collapse = " + "), y)
  fixef_table(m, keep = AXES) %>%
    mutate(
      window = win,
      ci_lo = estimate - qnorm(0.975) * se,
      ci_hi = estimate + qnorm(0.975) * se,
      .before = 1
    )
})
print(as.data.frame(window_profile), digits = 4)
write_outcome(window_profile, file.path(OUT, "temporal_window_profile.csv"))

# ------------------------------------------ 2. cluster bootstrap, fixed target
rule("2. SUBJECT-CLUSTER BOOTSTRAP (1,000 RESAMPLES)")
boot <- boot_lmm(
  cc,
  paste(AXES, collapse = " + "),
  TARGET_DV,
  terms = AXES,
  n = BOOT_N,
  seed = 2407
)
print(as.data.frame(boot), digits = 4)
write_outcome(boot, file.path(OUT, "bootstrap_ci.csv"))

# ------------------------------------------------------ 3. sample sensitivity
rule("3. DECLARED SAMPLE / QUALITY SENSITIVITY")

sample_masks <- list(
  full = rep(TRUE, nrow(cc)),
  observed_assignment = cc$event_assignment == "positional",
  drop_collector_qc = cc$flag_qc_check == 0,
  drop_onset_displaced = cc$flag_onset_displaced == 0,
  qc_clean = cc$flag_qc_check == 0 & cc$flag_onset_displaced == 0,
  on_quadrant_only = cc$flag_off_quadrant_stimulus == 0,
  unique_q4_waveforms = cc$flag_shared_q4_waveform == 0
)

sensitivity <- imap_dfr(sample_masks, function(keep, label) {
  dd <- cc[keep, ]
  m <- p24_lmm(dd, paste(AXES, collapse = " + "), TARGET_DV)
  fixef_table(m, keep = "comfort_wi") %>%
    mutate(
      sample = label,
      n_cells = nrow(dd),
      n_subjects = n_distinct(dd$subject_id),
      .before = 1
    )
})
full_est <- sensitivity$estimate[sensitivity$sample == "full"]
sensitivity <- sensitivity %>%
  mutate(effect_retained_pct = 100 * estimate / full_est)
print(as.data.frame(sensitivity), digits = 4)
write_outcome(sensitivity, file.path(OUT, "sample_sensitivity.csv"))

# ----------------------------------------------------- 4. single-subject influence
rule("4. LEAVE-ONE-SUBJECT-OUT INFLUENCE")

subject_flags <- cc %>%
  distinct(subject_id, flag_qc_check, flag_onset_displaced, event_assignment)
jackknife <- map_dfr(unique(cc$subject_id), function(ss) {
  dd <- cc %>% filter(subject_id != ss)
  m <- p24_lmm(dd, paste(AXES, collapse = " + "), TARGET_DV)
  x <- fixef_table(m, keep = "comfort_wi")
  tibble(
    dropped_subject = ss,
    estimate = x$estimate,
    se = x$se,
    p = x$p
  )
}) %>%
  left_join(subject_flags, by = c("dropped_subject" = "subject_id")) %>%
  mutate(delta_from_full = estimate - full_est) %>%
  arrange(desc(abs(delta_from_full)))

jackknife_summary <- jackknife %>%
  summarise(
    n_deletions = n(),
    min_estimate = min(estimate),
    max_estimate = max(estimate),
    all_negative = all(estimate < 0),
    max_abs_change = max(abs(delta_from_full)),
    min_p = min(p),
    max_p = max(p)
  )
print(as.data.frame(jackknife_summary), digits = 4)
cat("\nMost influential deletions:\n")
print(as.data.frame(head(jackknife, 10)), digits = 4)
write_outcome(jackknife, file.path(OUT, "leave_one_subject_out.csv"))
write_outcome(jackknife_summary, file.path(OUT, "leave_one_subject_out_summary.csv"))

# ------------------------------------------ 5. repetition / design heterogeneity
rule("5. REPETITION, STIMULUS-SET AND FORM STABILITY")

bb <- load_block_analysis("HbO") %>%
  filter(channel_index == TARGET_CHANNEL)

run_slopes <- map_dfr(1:3, function(rr) {
  dd <- bb %>% filter(run_index == rr)
  m <- p24_lmm(dd, "comfort_wi + richness_wi + position_f",
               "response_m5_30")
  fixef_table(m, keep = AXES) %>%
    mutate(run = rr, n_cells = nrow(dd),
           n_subjects = n_distinct(dd$subject_id), .before = 1)
})

subgroup_slopes <- bind_rows(
  map_dfr(sort(unique(cc$group)), function(gg) {
    dd <- cc %>% filter(group == gg)
    m <- p24_lmm(dd, paste(AXES, collapse = " + "), TARGET_DV)
    fixef_table(m, keep = AXES) %>%
      mutate(split = "stimulus_group", level = as.character(gg),
             n_subjects = n_distinct(dd$subject_id), .before = 1)
  }),
  map_dfr(sort(unique(cc$form)), function(ff) {
    dd <- cc %>% filter(form == ff)
    m <- p24_lmm(dd, paste(AXES, collapse = " + "), TARGET_DV)
    fixef_table(m, keep = AXES) %>%
      mutate(split = "sequence_form", level = as.character(ff),
             n_subjects = n_distinct(dd$subject_id), .before = 1)
  })
)

# Targeted omnibus tests for heterogeneity in the Comfort coefficient while
# allowing Richness to vary across the same grouping factor.
heterogeneity_lrt <- function(data, y, factor_name, label, extra = NULL) {
  f <- paste0("factor(", factor_name, ")")
  reduced <- paste(c(
    "comfort_wi",
    paste0("richness_wi * ", f),
    extra
  ), collapse = " + ")
  full <- paste(c(
    paste0("comfort_wi * ", f),
    paste0("richness_wi * ", f),
    extra
  ), collapse = " + ")
  m0 <- p24_lmm(data, reduced, y, reml = FALSE)
  m1 <- p24_lmm(data, full, y, reml = FALSE)
  a <- anova(m0, m1)
  tibble(
    dimension = label,
    df = a$Df[2],
    statistic = a$Chisq[2],
    p = a$`Pr(>Chisq)`[2]
  )
}

heterogeneity <- bind_rows(
  heterogeneity_lrt(
    bb, "response_m5_30", "run_index", "repetition/run",
    extra = "position_f"
  ),
  heterogeneity_lrt(
    cc, TARGET_DV, "group", "stimulus group"
  ),
  heterogeneity_lrt(
    cc, TARGET_DV, "form", "sequence form"
  )
) %>%
  mutate(q_BH = bh(p), family_size = n())

cat("Per-run slopes (slot controlled):\n")
print(as.data.frame(run_slopes), digits = 4)
cat("\nComfort-slope heterogeneity tests:\n")
print(as.data.frame(heterogeneity), digits = 4)
write_outcome(run_slopes, file.path(OUT, "run_stability.csv"))
write_outcome(subgroup_slopes, file.path(OUT, "group_form_slopes.csv"))
write_outcome(heterogeneity, file.path(OUT, "heterogeneity_tests.csv"))

# -------------------------------------- 6. stimulus versus idiosyncratic rating
rule("6. STIMULUS-CONSENSUAL VERSUS IDIOSYNCRATIC COMFORT")

si_model <- p24_lmm(
  cc,
  "comfort_stim + comfort_idio + richness_stim + richness_idio",
  TARGET_DV
)
stim_idio <- fixef_table(
  si_model,
  keep = c("comfort_stim", "comfort_idio",
           "richness_stim", "richness_idio")
)
print(as.data.frame(stim_idio), digits = 4)
write_outcome(stim_idio, file.path(OUT, "stimulus_vs_idiosyncratic.csv"))

# --------------------------------------------------- 7. chromophore convergence
rule("7. CHROMOPHORE CONVERGENCE AT CHANNEL 10")

by_chromophore <- map_dfr(CHROMOPHORES, function(chr) {
  dd <- load_condition_analysis(chr) %>%
    filter(channel_index == TARGET_CHANNEL)
  m <- p24_lmm(dd, paste(AXES, collapse = " + "), TARGET_DV)
  fixef_table(m, keep = AXES) %>%
    mutate(chromophore = chr, .before = 1)
})

coupling <- map_dfr(CHROMOPHORES, load_condition_analysis) %>%
  filter(channel_index == TARGET_CHANNEL) %>%
  select(subject_id, marker, chromophore, all_of(TARGET_DV)) %>%
  pivot_wider(names_from = chromophore, values_from = all_of(TARGET_DV)) %>%
  group_by(subject_id) %>%
  mutate(across(c(HbO, HbR, HbT), ~ .x - mean(.x), .names = "{.col}_c")) %>%
  ungroup() %>%
  summarise(
    n_cells = sum(complete.cases(HbO_c, HbR_c)),
    r_HbO_HbR_within = cor(HbO_c, HbR_c, use = "complete.obs"),
    sd_ratio_HbO_HbR = sd(HbO_c, na.rm = TRUE) / sd(HbR_c, na.rm = TRUE)
  )

print(as.data.frame(by_chromophore), digits = 4)
cat("\nWithin-subject channel-10 response-window coupling:\n")
print(as.data.frame(coupling), digits = 4)
write_outcome(by_chromophore, file.path(OUT, "chromophore_slopes.csv"))
write_outcome(coupling, file.path(OUT, "chromophore_coupling.csv"))

# --------------------------------------------- 8. optode neighbourhood context
rule("8. OPTODE-NEIGHBOURHOOD CONTEXT")

meta <- load_channels() %>%
  filter(active_pair == 1) %>%
  group_by(channel_index) %>%
  summarise(
    roi = first(roi),
    source = first(source),
    detector = first(detector),
    n_active_subjects = n(),
    .groups = "drop"
  )
target_meta <- meta %>% filter(channel_index == TARGET_CHANNEL)
neighbour_channels <- meta %>%
  filter(
    source == target_meta$source[[1]] |
      detector == target_meta$detector[[1]]
  ) %>%
  pull(channel_index) %>%
  sort()

neighbour_individual <- a06 %>%
  filter(channel_index %in% neighbour_channels, term == "comfort_wi") %>%
  select(channel_index, roi, source, detector, n_subjects,
         estimate, se, p, q_BH) %>%
  arrange(channel_index)

neighbour_cells <- load_condition_analysis("HbO") %>%
  filter(channel_index %in% neighbour_channels) %>%
  group_by(subject_id, marker) %>%
  summarise(
    across(c(mean_baseline_m10_0, mean_stim_m0_30, mean_response_m5_30,
             mean_response_minus_baseline_m5_30),
           ~ mean(.x, na.rm = TRUE)),
    comfort_wi = first(comfort_wi),
    richness_wi = first(richness_wi),
    n_channels = n(),
    .groups = "drop"
  )
neighbour_aggregate <- map_dfr(
  c(
    response = "mean_response_m5_30",
    baseline = "mean_baseline_m10_0",
    difference = "mean_response_minus_baseline_m5_30"
  ),
  function(y) {
    m <- p24_lmm(neighbour_cells, paste(AXES, collapse = " + "), y)
    fixef_table(m, keep = AXES)
  },
  .id = "window"
)

cat("Channel 10 shares either source 6 or detector 4 with channels:",
    paste(neighbour_channels, collapse = ", "), "\n")
cat("\nIndividual response-window Comfort coefficients:\n")
print(as.data.frame(neighbour_individual), digits = 4)
cat("\nPost-selection neighbourhood aggregates (descriptive only):\n")
print(as.data.frame(neighbour_aggregate), digits = 4)
write_outcome(neighbour_individual,
              file.path(OUT, "neighbour_channel_slopes.csv"))
write_outcome(neighbour_aggregate,
              file.path(OUT, "neighbourhood_aggregate.csv"))

close_log(con)
cat("Module 07 done.\n")
