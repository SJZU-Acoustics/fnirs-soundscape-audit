# =============================================================================
# mod22_chain.R — demonstrating the chain: does each validity link degrade the
# same quantity?
#
# Four links form a chain rather than a list. This module tests it on one
# scalar with a unit:
#
#   lambda = the smallest true within-subject appraisal slope this study could
#            have detected, at alpha = 0.05 and 80% power, in uM per axis unit
#          = 2.80 * SE / a
#
# where `a` is the fraction of a true effect the measurement chain delivers to
# the estimate, and SE is the standard error the real data actually produce.
# Every link either shrinks `a` or inflates SE, so every link raises lambda.
#
# WHY `a` IS COMPUTED, NOT SIMULATED. In the no-motion-correction arm the
# pipeline from optical density onward is linear (OD -> band-pass -> modified
# Beer-Lambert), so an injected response superposes exactly and its recovery
# can be read off analytically. This does NOT hold for the wavelet arm, which
# is nonlinear — and that is the arm the motion-correction audit rejected.
# Stated in the report.
#
# This module cannot rescue the null and does not try: lambda is a property of
# the design and the measurement chain, computed without reference to whether
# any effect exists in the data.
#
# Outputs (output/analysis_22_chain/): recovery_per_subject.csv,
# channel_dilution.csv, se_by_roi.csv, se_by_roi_corrected.csv,
# lambda_chain.csv, per_link_multiplier.csv, estimability.csv, chain_report.txt
# =============================================================================

source("code/load_data.R")
suppressPackageStartupMessages(library(clubSandwich))
# alias rather than attach: signal::filter masks dplyr::filter
butter   <- signal::butter
filtfilt <- signal::filtfilt

OUT <- file.path("output", "analysis_22_chain")
con <- open_log(file.path(OUT, "chain_report.txt"))

rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")
AXES <- c("comfort_wi", "richness_wi")

# ---- the collector's own machinery -------------------------------------------
FS  <- 11.0                                   # sampling rate (delivered series)
NYQ <- FS / 2
BF  <- butter(3, c(0.01, 0.10) / NYQ, type = "pass")

canonical_hrf <- function(peak = 5, dt = 1 / FS, len = 60) {
  t <- seq(0, len, by = dt)
  a1 <- peak + 1; a2 <- 16
  h <- dgamma(t, shape = a1, rate = 1) - dgamma(t, shape = a2, rate = 1) / 6
  h / max(h)
}
HRF <- canonical_hrf()

# Peak of an isolated unit-weight event, used to normalise the injection so that
# a weight of 1 means "a stimulus-locked response whose PEAK is 1 uM". Without
# this the injected amplitude is an arbitrary convolution sum and lambda has no
# physical unit.
UNIT_PEAK <- local({
  dur <- round(30 * FS)
  x <- c(numeric(round(10 * FS)), rep(1, dur), numeric(round(60 * FS)))
  max(stats::convolve(x, rev(HRF), type = "open")[seq_along(x)])
})

#' A 30-s boxcar at each onset, convolved with the HRF, weighted per event, and
#' scaled so a unit weight gives a 1 uM PEAK response.
inject <- function(onsets, weights, n_samp) {
  x <- numeric(n_samp)
  dur <- round(30 * FS)
  for (i in seq_along(onsets)) {
    s <- round(onsets[i] * FS) + 1
    if (s < 1 || s > n_samp) next
    e <- min(s + dur - 1, n_samp)
    x[s:e] <- x[s:e] + weights[i]
  }
  stats::convolve(x, rev(HRF), type = "open")[seq_len(n_samp)] / UNIT_PEAK
}

#' Window metric on a series: mean[5,30) minus mean[-10,0), per event.
window_metric <- function(series, onsets) {
  vapply(onsets, function(o) {
    idx <- function(a, b) {
      i <- (round((o + a) * FS) + 1):(round((o + b) * FS))
      i <- i[i >= 1 & i <= length(series)]
      if (!length(i)) return(NA_real_)
      mean(series[i])
    }
    idx(5, 30) - idx(-10, 0)
  }, numeric(1))
}

#' Rebuild subject x condition ROI means from the presentation x channel layer
#' under the session-CV dead-channel mask (load_dead_channels()) plus the
#' detector exclusion — the mask the chain and the responder-heterogeneity
#' audit were declared against. Mirrors the default load_p24_roi_masked()
#' aggregation (channel mean within presentation, then mean over the cell's
#' presentations) but applies the declared mask in place of the union mask.
roi_masked_cv <- function(chromophore = "HbO", exclude_detectors = UNUSABLE_DETECTORS) {
  dead <- load_dead_channels() %>% distinct(subject_id, channel_index)
  if (length(exclude_detectors)) {
    dropped <- load_channels() %>% filter(.data$detector %in% exclude_detectors) %>%
      distinct(subject_id, channel_index)
    dead <- bind_rows(dead, dropped) %>% distinct(subject_id, channel_index)
  }
  dead <- dead %>% mutate(.dead = TRUE)
  bm <- load_block_metrics(active_only = TRUE) %>%
    filter(.data$chromophore == !!chromophore)
  map <- load_channels() %>% filter(active_pair == 1) %>%
    distinct(subject_id, channel_index, roi)
  agg <- bm %>%
    left_join(dead, by = c("subject_id", "channel_index")) %>%
    filter(is.na(.data$.dead)) %>%
    left_join(map, by = c("subject_id", "channel_index")) %>%
    filter(!is.na(roi)) %>%
    group_by(subject_id, marker, event_index, roi) %>%
    summarise(block_mean = mean(response_minus_baseline_m5_30),
              n_live = n(), .groups = "drop") %>%
    group_by(subject_id, marker, roi) %>%
    summarise(value = mean(block_mean), n_live = max(n_live), .groups = "drop")
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

cat("A22 — the chain. lambda = smallest detectable PEAK response amplitude,\n")
cat("in uM per unit of the appraisal axis; smaller is better.\n")
cat(sprintf("unit-peak normalisation constant: %.3f\n", UNIT_PEAK))
cat("Sampling rate", FS, "Hz; band-pass 0.01-0.10 Hz, 3rd-order Butterworth,\n")
cat("zero-phase — the collector's own setting, reproduced from its definition.\n")

# ================================================ 1. recovery of a known effect
rule("1. RECOVERY `a` — how much of a true effect reaches the estimate?")
cat("Per subject: a 30-s boxcar at each frozen onset, convolved with the\n")
cat("canonical HRF, weighted by that event's within-subject-centred comfort,\n")
cat("unit amplitude. `a` is the slope of the resulting window metric on\n")
cat("comfort_wi — exactly the declared estimand, on a signal we control.\n\n")

ev <- load_events() %>% filter(!is.na(onset_sec)) %>%
  left_join(load_p24() %>% select(subject_id, marker, comfort_wi),
            by = c("subject_id", "marker")) %>%
  filter(!is.na(comfort_wi)) %>%
  arrange(subject_id, event_index)

rec <- ev %>% group_by(subject_id) %>% group_modify(function(d, key) {
  n_samp <- ceiling((max(d$onset_sec) + 90) * FS)
  w <- d$comfort_wi
  clean    <- inject(d$onset_sec, w, n_samp)
  filtered <- as.numeric(filtfilt(BF, clean))
  m_clean <- window_metric(clean, d$onset_sec)
  m_filt  <- window_metric(filtered, d$onset_sec)
  ok <- !is.na(m_clean) & !is.na(m_filt)
  tibble(a_nofilter = unname(coef(lm(m_clean[ok] ~ w[ok]))[2]),
         a_filter   = unname(coef(lm(m_filt[ok]  ~ w[ok]))[2]),
         n_events   = sum(ok))
}) %>% ungroup()

a_nofilter <- mean(rec$a_nofilter)
a_filter   <- mean(rec$a_filter)
cat(sprintf("no filter : a = %.4f  (SD across subjects %.4f)\n",
            a_nofilter, sd(rec$a_nofilter)))
cat(sprintf("0.01 Hz HP: a = %.4f  (SD %.4f)\n", a_filter, sd(rec$a_filter)))
cat(sprintf("\nthe filter delivers %.1f%% of what an unfiltered chain would\n",
            100 * a_filter / a_nofilter))
write_outcome(rec, file.path(OUT, "recovery_per_subject.csv"))

# ================================================ 2. channel loss dilutes `a`
rule("2. CHANNEL LOSS — a dark detector cannot record a cortical response")
cat("An injected cortical response appears only in channels that received\n")
cat("light. An ROI mean over n_total channels of which n_live are alive\n")
cat("therefore delivers n_live/n_total of it. Measured per ROI per subject.\n\n")

# NB `.d` as a column name silently partial-matches mutate()'s own `.data`
# argument; use a plain name.
dead <- load_dead_channels() %>% distinct(subject_id, channel_index) %>%
  mutate(is_dark = TRUE)
chan <- load_channels() %>% filter(active_pair == 1) %>%
  select(subject_id, channel_index, roi) %>%
  left_join(dead, by = c("subject_id", "channel_index")) %>%
  mutate(dark = !is.na(is_dark)) %>% select(-is_dark)

cat("Dilution is computed on the FROZEN ROI definition — the configuration whose\n")
cat("lambda this is — counting only channels that were DARK. Detector 8's live\n")
cat("sessions are left in, even though the detector audit shows their signal is\n")
cat("not neurovascular, because they do carry the injected response. That makes\n")
cat("this link's contribution a CONSERVATIVE floor, not a ceiling.\n\n")

dilution <- chan %>% group_by(subject_id, roi) %>%
  summarise(n_total = n(), n_live = sum(!dark), .groups = "drop") %>%
  mutate(frac = n_live / n_total) %>%
  group_by(roi) %>%
  summarise(mean_frac = mean(frac), n_subjects_full = sum(frac == 1),
            n_subjects_half_or_less = sum(frac <= 0.5), .groups = "drop")
print(as.data.frame(dilution), digits = 3)
write_outcome(dilution, file.path(OUT, "channel_dilution.csv"))

# ================================================= 3. SE from the real data
rule("3. SE — measured on the real data, model-based and CR2")
cat("The dependence link does not change the point estimate; it changes what\n")
cat("the standard error should be. CR2 is the honest choice on the unadjusted\n")
cat("family, where it inflates SEs by a median 105%.\n\n")

d <- load_p24()
se_tab <- map_dfr(ROIS, function(rr) {
  y <- paste0("HbO_", rr)
  dd <- d %>% filter(!is.na(.data[[y]]))
  m <- p24_lmm(dd, paste(AXES, collapse = " + "), y)
  ft <- fixef_table(m, keep = "comfort_wi")
  cr <- coef_test(m, vcov = "CR2", cluster = dd$subject_id, test = "Satterthwaite")
  cr <- cr[rownames(cr) == "comfort_wi" | cr$Coef == "comfort_wi", ]
  tibble(roi = rr, se_model = ft$se[1], se_cr2 = cr$SE[1],
         inflation = cr$SE[1] / ft$se[1])
})
print(as.data.frame(se_tab), digits = 3)
write_outcome(se_tab, file.path(OUT, "se_by_roi.csv"))

cat("\nThe right-temporal CR2 SE is SMALLER than the model-based one, which is\n")
cat("not what dependence normally does. That ROI is the compromised one, so the\n")
cat("same comparison is repeated on the corrected definition (detector 10 only):\n\n")
# The declared session-CV mask plus the detector-8 exclusion (see roi_masked_cv).
dm <- roi_masked_cv("HbO")
se_corrected <- map_dfr(ROIS, function(rr) {
  y <- paste0("HbOm_", rr)
  dd <- dm %>% filter(!is.na(.data[[y]]))
  m <- p24_lmm(dd, paste(AXES, collapse = " + "), y)
  ft <- fixef_table(m, keep = "comfort_wi")
  cr <- coef_test(m, vcov = "CR2", cluster = dd$subject_id, test = "Satterthwaite")
  cr <- cr[cr$Coef == "comfort_wi", ]
  tibble(roi = rr, se_model = ft$se[1], se_cr2 = cr$SE[1],
         inflation = cr$SE[1] / ft$se[1])
})
print(as.data.frame(se_corrected), digits = 3)
write_outcome(se_corrected, file.path(OUT, "se_by_roi_corrected.csv"))
cat("\nIf the reversal disappears on the corrected ROI, the dependence link does\n")
cat("behave as a link everywhere and the exception was another detector-8\n")
cat("artefact. If it persists, the link is genuinely ROI-specific.\n")

# ==================================================== 4. the chain, in lambda
rule("4. THE CHAIN — lambda under nested configurations")
cat("lambda = 2.80 * SE / a, in uM per axis unit. C0 is the study this design\n")
cat("could have been; each row adds exactly one link.\n\n")

cat("PRIMARY TABLE USES THE CORRECTED-ROI SEs, because the frozen right-temporal\n")
cat("SE is itself a detector-8 artefact (its CR2 inflation is 0.74, and 1.01 once\n")
cat("the detector is removed). The dilution factor still comes from the FROZEN\n")
cat("definition, because dilution is what the study actually did. The frozen-SE\n")
cat("variant is kept in the outcome file as a sensitivity.\n\n")

Z <- qnorm(0.975) + qnorm(0.80)          # 2.80
lam <- map_dfr(ROIS, function(rr) {
  se  <- se_corrected$se_model[se_corrected$roi == rr]
  se2 <- se_corrected$se_cr2[se_corrected$roi == rr]
  fr  <- dilution$mean_frac[dilution$roi == rr]
  tibble(
    roi = rr,
    C0_ideal            = Z * se  / a_nofilter,
    C1_plus_filter      = Z * se  / a_filter,
    C2_plus_channelloss = Z * se  / (a_filter * fr),
    C3_plus_dependence  = Z * se2 / (a_filter * fr))
}) %>% mutate(total_inflation = C3_plus_dependence / C0_ideal)
print(as.data.frame(lam), digits = 3)
lam_frozen <- map_dfr(ROIS, function(rr) {
  se  <- se_tab$se_model[se_tab$roi == rr]
  se2 <- se_tab$se_cr2[se_tab$roi == rr]
  fr  <- dilution$mean_frac[dilution$roi == rr]
  tibble(roi = rr, se_source = "frozen ROI (sensitivity)",
         C0_ideal = Z * se / a_nofilter, C1_plus_filter = Z * se / a_filter,
         C2_plus_channelloss = Z * se / (a_filter * fr),
         C3_plus_dependence = Z * se2 / (a_filter * fr))
})
write_outcome(bind_rows(lam %>% mutate(se_source = "corrected ROI (primary)"),
                        lam_frozen), file.path(OUT, "lambda_chain.csv"))

cat("\nper-link multiplier (what each link alone does to lambda):\n")
mult <- lam %>% transmute(roi,
  filter      = C1_plus_filter / C0_ideal,
  channelloss = C2_plus_channelloss / C1_plus_filter,
  dependence  = C3_plus_dependence / C2_plus_channelloss,
  compounded  = total_inflation)
print(as.data.frame(mult), digits = 3)
write_outcome(mult, file.path(OUT, "per_link_multiplier.csv"))

cat("\nreading rule, fixed before the run: the chain claim is supported only if\n")
cat("EACH link raises lambda on its own and the compounded value is materially\n")
cat("larger than any single link's.\n")
each_ok <- all(mult$filter > 1) && all(mult$channelloss >= 1) && all(mult$dependence > 1)
cat("  every link raises lambda in every ROI:", each_ok, "\n")
stopifnot(each_ok)
cat("  largest single-link multiplier:", signif(max(c(mult$filter, mult$channelloss,
                                                      mult$dependence)), 3),
    " vs compounded:", signif(max(mult$compounded), 3), "\n")

# ============================================ 5. the fourth link, different
rule("5. THE FOURTH LINK — estimability, which does not live on this scale")
cat("Channel loss also removed the planned left-right contrast outright: the\n")
cat("montage is mirror-symmetric and the two failures took opposite members of\n")
cat("the two mirror pairs. For that contrast lambda is not large, it is\n")
cat("undefined — no matched comparison exists in any usable subject. Reporting\n")
cat("it as one more multiplier would misrepresent it.\n\n")
mirror <- tibble(
  mirror_pair = c("T (16/17 vs 20/21)", "TP (18/19 vs 22/23)"),
  subjects_all_four_live = c(27L, 1L),
  note = c("right member is detector 8 — non-neurovascular when it registers",
           "left member is detector 9 — dark in 66-68 of 69 sessions"))
print(as.data.frame(mirror))
write_outcome(mirror, file.path(OUT, "estimability.csv"))

rule("VERDICT")
for (rr in ROIS) {
  r <- lam[lam$roi == rr, ]
  cat(sprintf("%-15s lambda %.3f -> %.3f uM  (x%.2f)\n", rr,
              r$C0_ideal, r$C3_plus_dependence, r$total_inflation))
}
cat("\nFOR SCALE, on the same footing. The prefrontal difference-score comfort\n")
cat("slope is -0.052 uM in WINDOW-METRIC units. Dividing by the chain's own\n")
cat("recovery converts it to the underlying peak amplitude it would imply:\n")
implied_peak <- -0.052 / a_filter
cat(sprintf("  observed window slope -0.052 uM  ->  underlying peak %.3f uM\n",
            implied_peak))
ratio_full  <- abs(implied_peak) / lam$C3_plus_dependence[lam$roi == "PFC_Frontal"]
ratio_ideal <- abs(implied_peak) / lam$C0_ideal[lam$roi == "PFC_Frontal"]
cat(sprintf("\nThat implied effect is %.2f x lambda under the full chain, and %.2f x\n",
            ratio_full, ratio_ideal))
cat("lambda under the ideal configuration. BOTH ARE BELOW 1, and that is the\n")
cat("most consequential number in this analysis:\n\n")
cat("  * An effect of the size actually observed at prefrontal was BELOW this\n")
cat("    study's detection threshold. Non-detection of it is what an\n")
cat("    underpowered design produces, not evidence that it is absent.\n")
cat("  * It was below threshold even in the IDEAL configuration, so this is not\n")
cat("    only a story about the measurement chain: the design was underpowered\n")
cat("    for effects of this magnitude before any link degraded it, and the\n")
cat("    chain then made it about 1.7x worse at prefrontal.\n")
cat("  * The synthesis sentence \"the study is not underpowered\" therefore\n")
cat("    cannot stand as written and must be narrowed — see the report.\n")

# anchor: the implied peak must fall below prefrontal lambda at every config
pfc_lam <- lam %>% filter(roi == "PFC_Frontal")
stopifnot(abs(implied_peak) < pfc_lam$C0_ideal,
          abs(implied_peak) < pfc_lam$C1_plus_filter,
          abs(implied_peak) < pfc_lam$C2_plus_channelloss,
          abs(implied_peak) < pfc_lam$C3_plus_dependence)

close_log(con)
