# Module 17 — Design-matched canonical-HRF GLM (second level, from shipped betas)
#
# Port of exploration/analysis_17_hrf_glm (batch-5 declaration).
#
# The FIRST level of this analysis fits one GLM per subject x channel on the
# delivered continuous concentration series (the stage-1 builder's own input).
# That series, and the subject crosswalk needed to align it with the design
# table, are NOT part of the Mendeley deposit, so the first level cannot be
# re-run from the workbook alone. What CAN be deposited is its entire output:
# the per-subject GLM coefficients. This module therefore reads the shipped
# intermediate
#   data/intermediates/a17_glm_subject_betas.csv
# — one row per subject_id x channel_index x chromophore x HRF peak, carrying
#   beta_main      unmodulated boxcar beta (response-existence check)
#   beta_comfort   comfort-modulated boxcar beta
#   beta_richness  richness-modulated boxcar beta
# at the canonical peak (5 s; HbO and HbR) and at the four swept peaks
# (8/11/14/17 s; HbO only) — and recomputes every group-level statistic of the
# original analysis from those betas alone. All betas are in uM per rating
# point (modulators) or uM (main boxcar).
#
# What it changes relative to A01-A16 is the ESTIMATOR. Every earlier analysis
# reads a difference between two window means on a series whose 90-s stimulus
# cycle (0.0111 Hz) sits just above the 0.01 Hz high-pass corner. A GLM with a
# canonical haemodynamic response asks the same question with a model of the
# expected response shape instead of two arbitrary windows.
#
# The design matrix is passed through the SAME 0.01-0.10 Hz band-pass as the
# data (3rd-order Butterworth, zero-phase filtfilt, verified stable at
# max|pole| = 0.9973). Without this the model is misspecified against a filtered
# series: the filter's gain at the stimulus-cycle frequency is 0.688, so an
# unfiltered regressor would be fitted to data whose response has been attenuated
# to 69% of its amplitude. That filter characterisation is a deterministic fact
# of the filter definition and is recomputed below (section 2).
#
# Declared families (fixed before the original analysis ran):
#   PRIMARY  2 axes x 3 ROIs = 6, BH   (mirrors A04)
#   CHANNEL  21 channels x 2 axes = 42, BH   (mirrors A06)
#   Descriptive: unmodulated boxcar beta (response-existence check); ch10;
#   HbR; GLM-vs-window-difference agreement.
#
# Run from the repository root (Rscript code/mod17_hrf_glm.R). Writes
# output/analysis_17_hrf_glm/.

source("code/load_data.R")
# `signal` is NOT attached: signal::filter would mask dplyr::filter and silently
# break every pipeline in this script. Called namespaced instead.
butter   <- signal::butter
filtfilt <- signal::filtfilt

OUT <- file.path("output", "analysis_17_hrf_glm")
con <- open_log(file.path(OUT, "hrf_glm_report.txt"))
rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")

FS <- 10; NYQ <- FS / 2
BF <- butter(3, c(0.01, 0.10) / NYQ, type = "pass")   # the collector's band
CLIP_SEC <- 30
AXES <- c("comfort_wi", "richness_wi")
CANONICAL_PEAK <- 5
SWEEP_PEAKS <- c(5, 8, 11, 14, 17)

# ------------------------------------------------------- canonical HRF (SPM form)
# double gamma: peak ~6 s, undershoot ~16 s, ratio 1/6 — the standard used for
# fNIRS GLMs as well as fMRI. `peak` is the time-to-peak in seconds
# (dgamma(shape=a) peaks at a-1), so the canonical peak of 5 s is shape 6.
canonical_hrf <- function(peak = 5, dt = 1 / FS, len = 60) {
  t <- seq(0, len, by = dt)
  a1 <- peak + 1; a2 <- a1 + 10
  h <- dgamma(t, shape = a1, rate = 1) - dgamma(t, shape = a2, rate = 1) / 6
  h / max(h)
}

# ------------------------------------------------------------------- 0. the betas
rule("0. SHIPPED FIRST-LEVEL BETAS")

betas <- load_intermediate("a17_glm_subject_betas.csv")
chan_map <- load_channels() %>% filter(active_pair == 1) %>%
  distinct(subject_id, channel_index, roi)
betas <- betas %>% left_join(chan_map, by = c("subject_id", "channel_index"))
cat("beta rows:", nrow(betas), " subjects:", n_distinct(betas$subject_id),
    " peaks:", paste(sort(unique(betas$hrf_peak_s)), collapse = "/"), "\n")

# --------------------------------------------- 2. filter characterisation (fact)
rule("2. WHAT THE COLLECTOR'S FILTER DOES TO THE STIMULUS CYCLE")
tprobe <- seq(0, 1231, by = 1 / FS)
for (per in c(90, 60, 30, 20)) {
  s <- sin(2 * pi * tprobe / per); yv <- filtfilt(BF, s)
  cat(sprintf("  period %3d s (%.4f Hz): amplitude gain %.3f\n", per, 1 / per,
              max(abs(yv[500:11000]))))
}
cat("\nThe 90-s presentation cycle is passed at 0.688 of its amplitude — the\n",
    "quantified form of A03's filter-corner observation.\n", sep = "")

# --------------------------------------------------------- 3. PRIMARY ROI family
rule("3. PRIMARY FAMILY — 2 axes x 3 ROIs = 6, BH (mirrors A04)")

roi_beta <- betas %>% filter(chromophore == "HbO", hrf_peak_s == CANONICAL_PEAK,
                             !is.na(roi)) %>%
  group_by(subject_id, roi) %>%
  summarise(comfort = mean(beta_comfort), richness = mean(beta_richness),
            main = mean(beta_main), n_ch = n(), .groups = "drop")

roi_tab <- map_dfr(ROIS, function(rr) {
  d <- roi_beta %>% filter(roi == rr)
  map_dfr(c("comfort", "richness"), function(ax) {
    tt <- t.test(d[[ax]])
    tibble(roi = rr, axis = paste0(ax, "_wi"), n_subjects = nrow(d),
           mean_beta = unname(tt$estimate), se = tt$stderr,
           t = unname(tt$statistic), p = tt$p.value,
           ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2])
  })
}) %>% mutate(q_BH = bh(p))
print(as.data.frame(roi_tab), digits = 3, row.names = FALSE)
cat("\nSurvivors at q<0.05:", sum(roi_tab$q_BH < 0.05), "\n")
write_outcome(roi_tab, file.path(OUT, "roi_family_glm.csv"))

# ------------------------------------- 4. response-existence check (descriptive)
rule("4. IS THERE ANY MODELLED RESPONSE AT ALL? (unmodulated boxcar, descriptive)")
main_tab <- map_dfr(ROIS, function(rr) {
  d <- roi_beta %>% filter(roi == rr); tt <- t.test(d$main)
  tibble(roi = rr, n = nrow(d), mean_beta_main = unname(tt$estimate),
         t = unname(tt$statistic), p = tt$p.value)
})
print(as.data.frame(main_tab), digits = 3, row.names = FALSE)
write_outcome(main_tab, file.path(OUT, "main_effect_boxcar.csv"))

# ------------------------------------------------------- 5. CHANNEL family (42)
rule("5. CHANNEL FAMILY — 21 x 2 = 42, BH (mirrors A06)")

ch_tab <- betas %>% filter(chromophore == "HbO", hrf_peak_s == CANONICAL_PEAK) %>%
  group_by(channel_index) %>%
  group_modify(~ {
    map_dfr(c("comfort", "richness"), function(ax) {
      v <- .x[[paste0("beta_", ax)]]
      tt <- t.test(v)
      tibble(axis = paste0(ax, "_wi"), n_subjects = length(v),
             mean_beta = unname(tt$estimate), t = unname(tt$statistic), p = tt$p.value)
    })
  }) %>% ungroup() %>% mutate(q_BH = bh(p)) %>%
  left_join(chan_map %>% distinct(channel_index, roi), by = "channel_index")

cat("\nTop 8 by p (family of 42, BH):\n")
print(as.data.frame(ch_tab %>% arrange(p) %>% head(8)), digits = 3, row.names = FALSE)
cat("\nSurvivors at q<0.05:", sum(ch_tab$q_BH < 0.05), "\n")
cat("\nChannel 10:\n")
print(as.data.frame(ch_tab %>% filter(channel_index == 10)), digits = 3, row.names = FALSE)
write_outcome(ch_tab, file.path(OUT, "channel_family_glm.csv"))

# ------------------------------------------------------- 6. HbR (descriptive)
rule("6. HbR — chromophore convergence (descriptive)")
hbr_tab <- betas %>% filter(chromophore == "HbR", hrf_peak_s == CANONICAL_PEAK,
                            !is.na(roi)) %>%
  group_by(subject_id, roi) %>%
  summarise(comfort = mean(beta_comfort), .groups = "drop") %>%
  group_by(roi) %>% summarise(mean_beta = mean(comfort),
                              t = t.test(comfort)$statistic,
                              p = t.test(comfort)$p.value, .groups = "drop")
print(as.data.frame(hbr_tab), digits = 3, row.names = FALSE)
write_outcome(hbr_tab, file.path(OUT, "hbr_descriptive.csv"))

# ------------------------------- 7. agreement with the window-difference measure
rule("7. DOES THE GLM AGREE WITH THE WINDOW-DIFFERENCE MEASURE? (descriptive)")
win_slope <- load_condition_analysis("HbO") %>%
  group_by(channel_index) %>%
  group_modify(~ {
    m <- p24_lmm(.x, paste(AXES, collapse = " + "), "mean_response_m5_30")
    ft <- fixef_table(m, keep = "comfort_wi")
    tibble(window_slope = ft$estimate)
  }) %>% ungroup()
agree <- ch_tab %>% filter(axis == "comfort_wi") %>%
  select(channel_index, glm_beta = mean_beta) %>%
  inner_join(win_slope, by = "channel_index")
cat("Pearson r across", nrow(agree), "channels:",
    round(cor(agree$glm_beta, agree$window_slope), 3), "\n")
cat("Spearman:", round(cor(agree$glm_beta, agree$window_slope, method = "spearman"), 3), "\n")
write_outcome(agree, file.path(OUT, "glm_vs_window_agreement.csv"))

# ------------------------------------------------- 8. HRF-delay sweep (descriptive)
rule("8. HRF-DELAY SWEEP — is the canonical shape simply the wrong one?")

# A11 localised the comfort association to [20,30) s. A canonical HRF convolved
# with a 30-s boxcar reaches plateau within ~10-15 s, so a canonical-only GLM
# could miss a genuinely later response by shape mismatch alone. This sweeps the
# time-to-peak and reports the whole curve. DESCRIPTIVE and uncorrected: it is a
# search over shapes, declared as such, and may not carry a claim.

cat("plateau timing of the convolved 30-s boxcar, by HRF peak:\n")
for (pk in SWEEP_PEAKS) {
  h <- canonical_hrf(pk)
  bc <- stats::convolve(c(rep(1, CLIP_SEC * FS), rep(0, 600)), rev(h), type = "open")
  cat(sprintf("  peak %2d s -> convolved max at %.1f s\n", pk, (which.max(bc) - 1) / FS))
}

sweep_rows <- list()
for (pk in SWEEP_PEAKS) {
  sbd <- betas %>% filter(chromophore == "HbO", hrf_peak_s == pk)
  # prefrontal ROI and ch10, the two a-priori targets (A04/A06)
  pfc <- sbd %>% filter(roi == "PFC_Frontal") %>%
    group_by(subject_id) %>% summarise(b = mean(beta_comfort), .groups = "drop")
  c10 <- sbd %>% filter(channel_index == 10)
  sweep_rows[[length(sweep_rows) + 1]] <- tibble(
    hrf_peak_s = pk,
    pfc_beta = mean(pfc$b), pfc_t = t.test(pfc$b)$statistic, pfc_p = t.test(pfc$b)$p.value,
    ch10_beta = mean(c10$beta_comfort), ch10_t = t.test(c10$beta_comfort)$statistic,
    ch10_p = t.test(c10$beta_comfort)$p.value)
}
sweep_tab <- bind_rows(sweep_rows)
print(as.data.frame(sweep_tab), digits = 3, row.names = FALSE)
cat("\nSmallest prefrontal p over the sweep:", signif(min(sweep_tab$pfc_p), 3),
    "at peak", sweep_tab$hrf_peak_s[which.min(sweep_tab$pfc_p)], "s\n")
cat("Smallest ch10 p over the sweep:", signif(min(sweep_tab$ch10_p), 3),
    "at peak", sweep_tab$hrf_peak_s[which.min(sweep_tab$ch10_p)], "s\n")
write_outcome(sweep_tab, file.path(OUT, "hrf_delay_sweep.csv"))

rule("DONE")
close_log(con)
