# mod06_channel_scan.R — Channel-wise spatial scan.
#
# Release port of the stage-2 exploration pass A06
# (exploration/analysis_06_channelwise_spatial_scan).
#
# The ROI-level modules average channels into three ROIs. This module asks
# whether that aggregation hides a spatially focal association. It scans every
# active channel without excluding the two near-flat channels (18/19).
#
# Family fixed before running: within each window, 21 active channels x 2
# measured axes = 42 tests, BH corrected together. The three windows remain
# separate families; baseline is the negative-control family.
#
# Run from the repository root:  Rscript code/mod06_channel_scan.R

source("code/load_data.R")

OUT <- file.path("output", "analysis_06_channelwise_spatial_scan")
con <- open_log(file.path(OUT, "channelwise_spatial_scan_report.txt"))

WINDOWS <- c(
  difference = "mean_response_minus_baseline_m5_30",
  response   = "mean_response_m5_30",
  baseline   = "mean_baseline_m10_0"
)
AXES <- c("comfort_wi", "richness_wi")

# Subject x condition x ROI window means rebuilt from the condition x channel
# layer, joined to the canonical design. The aggregation includes every active
# channel in the ROI — the contextual estimates below must describe the same
# aggregation the scan challenges.
condition_roi_analysis <- function(chromophore = "HbO") {
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
    inner_join(
      load_p24() %>% select(subject_id, marker, comfort_wi, richness_wi),
      by = c("subject_id", "marker")
    )
}

d <- load_condition_analysis(chromophore = "HbO")
channels <- sort(unique(d$channel_index))
stopifnot(length(channels) == 21L,
          setequal(as.integer(channels), setdiff(1:23, c(6, 13))))

channel_meta <- load_channels() %>%
  group_by(channel_index, roi) %>%
  summarise(
    source = first(source),
    detector = first(detector),
    n_active_subjects = sum(active_pair == 1),
    mean_sd_distance_mm = mean(sd_distance_mm[active_pair == 1]),
    .groups = "drop"
  ) %>%
  mutate(dead_channel_A03 = as.integer(channel_index %in% FLAT_CHANNELS))

rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n",
                        strrep("=", 74), "\n", sep = "")

cat("Module 06 — channel-wise spatial screen\n")
cat("21 active channels; both axes enter the same model.\n")
cat("Family per window: 21 channels x 2 axes = 42 tests, BH.\n")
cat("Channels 18/19 remain in the family and are flagged, not removed.\n")

scan <- map_dfr(names(WINDOWS), function(win) {
  y <- unname(WINDOWS[[win]])
  map_dfr(channels, function(cc) {
    dd <- d %>% filter(channel_index == cc, !is.na(.data[[y]]))
    m <- p24_lmm(dd, paste(AXES, collapse = " + "), y)
    fixef_table(m, keep = AXES) %>%
      mutate(
        window = win,
        channel_index = cc,
        n_cells = nrow(dd),
        n_subjects = n_distinct(dd$subject_id),
        .before = 1
      )
  })
}) %>%
  left_join(channel_meta, by = "channel_index") %>%
  group_by(window) %>%
  mutate(
    q_BH = bh(p),
    family_size = n(),
    rank_in_42 = min_rank(p),
    ci_lo = estimate - qnorm(0.975) * se,
    ci_hi = estimate + qnorm(0.975) * se
  ) %>%
  ungroup() %>%
  arrange(factor(window, levels = names(WINDOWS)), p)

stopifnot(all(scan %>% count(window) %>% pull(n) == 42L))

family_summary <- scan %>%
  group_by(window) %>%
  summarise(
    family_size = first(family_size),
    min_p = min(p),
    min_q = min(q_BH),
    n_q_lt_05 = sum(q_BH < 0.05),
    winning_channel = channel_index[which.min(q_BH)],
    winning_axis = term[which.min(q_BH)],
    .groups = "drop"
  )

spatial_summary <- scan %>%
  group_by(window, term, roi) %>%
  summarise(
    n_channels = n(),
    median_estimate = median(estimate),
    n_negative = sum(estimate < 0),
    n_positive = sum(estimate > 0),
    best_channel = channel_index[which.min(p)],
    best_estimate = estimate[which.min(p)],
    best_p = min(p),
    best_q_family42 = q_BH[which.min(p)],
    .groups = "drop"
  )

# Directly quantify the aggregation choice being challenged: the same response
# model at the three ROIs, on the table reconstructed from the channel layer.
# These are contextual estimates, not extra family members.
roi_d <- condition_roi_analysis(chromophore = "HbO")
roi_response <- map_dfr(ROIS, function(rr) {
  dd <- roi_d %>% filter(roi == rr)
  m <- p24_lmm(dd, paste(AXES, collapse = " + "),
               "mean_response_m5_30")
  fixef_table(m, keep = AXES) %>%
    mutate(roi = rr, n_cells = nrow(dd), .before = 1)
})

rule("FAMILY SUMMARY")
print(as.data.frame(family_summary), digits = 4)

rule("ALL FDR SURVIVORS")
survivors <- scan %>% filter(q_BH < 0.05)
if (nrow(survivors) == 0) {
  cat("None.\n")
} else {
  print(as.data.frame(
    survivors %>%
      select(window, channel_index, roi, source, detector, term,
             estimate, se, ci_lo, ci_hi, p, q_BH, n_subjects)
  ), digits = 4)
}

rule("TOP FIVE TESTS IN EACH WINDOW")
print(as.data.frame(
  scan %>%
    group_by(window) %>%
    slice_min(p, n = 5, with_ties = FALSE) %>%
    ungroup() %>%
    select(window, channel_index, roi, source, detector, term,
           estimate, se, p, q_BH, dead_channel_A03)
), digits = 4)

rule("ROI RESPONSE-WINDOW CONTEXT")
print(as.data.frame(roi_response), digits = 4)

write_outcome(scan, file.path(OUT, "channel_scan_all.csv"))
write_outcome(
  scan %>% filter(window == "response"),
  file.path(OUT, "channel_scan_response.csv")
)
write_outcome(family_summary, file.path(OUT, "family_summary.csv"))
write_outcome(spatial_summary, file.path(OUT, "spatial_summary.csv"))
write_outcome(roi_response, file.path(OUT, "roi_response_context.csv"))

close_log(con)
cat("Module 06 done.\n")
