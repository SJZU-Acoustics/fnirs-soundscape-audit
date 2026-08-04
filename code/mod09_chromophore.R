# mod09_chromophore.R — Chromophore triangulation of the ROI families
#
# The declared 2 axes x 3 ROIs family ran on HbO only. The HbO/HbR MEANS are
# physiologically incoherent (HbO at p ~ 1e-19, HbR flat) while the trial-level
# FLUCTUATIONS are properly coupled (within-subject r ~ -0.5). What has never
# been tested is whether the condition-level appraisal slopes exist on HbR or
# HbT at all, and whether the HbO:HbR slope signature, cell by cell, looks
# neurovascular (opposite signs, |HbO| : |HbR| ~ 2.3 : 1) or blood-volume-like
# (same sign, HbT tracking HbO).
#
# Declared families:
#   per window (frozen difference score; response [5,30) alone):
#     2 axes x 3 ROIs x 2 chromophores (HbR, HbT) = 12 tests, BH, q < 0.05.
#   HbO is quoted from the HbO family, not re-tested.
#   The HbO:HbR slope-ratio diagnostics are subject-bootstrap coefficient
#   diagnostics, not significance tests, and carry no correction.
#
# Outputs (output/analysis_09_chromophore_triangulation/):
#   family_difference_hbr_hbt.csv, family_response_hbr_hbt.csv,
#   hbt_identity_check.csv, slope_ratio_coupling.csv
# Display anchors: main-text Table 2 and SI Table S1 quote the two family
# minima (min q 0.364 difference / 0.290 response); Table 1 quotes the
# HbT = HbO + HbR identity to 1e-16.

if (!exists("load_p24")) source("code/load_data.R")

OUT <- file.path("output", "analysis_09_chromophore_triangulation")
con <- open_log(file.path(OUT, "chromophore_triangulation_report.txt"))

AXES <- c("comfort_wi", "richness_wi")
rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")

# subject x condition x ROI means for one chromophore, joined to the canonical
# analysis design — built here from the public loaders.
condition_roi_long <- function(chromophore) {
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

d <- load_p24()                       # 276 rows; difference-score columns per ROI
resp <- map_dfr(CHROMOPHORES, ~ condition_roi_long(.x) %>%
                  mutate(chromophore = .x))   # response-window ROI means, all three

# ------------------------------------------------ 1. the two declared families
# Difference score: wide columns {HbR,HbT}_<ROI> in the analysis table.
fam_diff <- map_dfr(c("HbR", "HbT"), function(cc) {
  map_dfr(ROIS, function(rr) {
    y <- paste0(cc, "_", rr)
    dd <- d %>% filter(!is.na(.data[[y]]))
    m <- p24_lmm(dd, paste(AXES, collapse = " + "), y)
    fixef_table(m, keep = AXES) %>%
      mutate(chromophore = cc, roi = rr, n_rows = nrow(dd), .before = 1)
  })
}) %>% mutate(q_BH = bh(p), family_size = n(), window = "difference")

rule("1a. DIFFERENCE SCORE (frozen DV) — HbR and HbT, family of 12")
print(as.data.frame(fam_diff %>%
      select(chromophore, roi, term, estimate, se, t, p, q_BH)), digits = 3)
write_outcome(fam_diff, file.path(OUT, "family_difference_hbr_hbt.csv"))

# Response window: `resp` is long over roi, one mean_response_m5_30 per
# subject x marker x roi — subset to the ROI before fitting.
rule("1b. RESPONSE WINDOW [5,30) — HbR and HbT, family of 12")
fam_resp <- map_dfr(c("HbR", "HbT"), function(cc) {
  map_dfr(ROIS, function(rr) {
    dd <- resp %>% filter(chromophore == cc, roi == rr,
                          !is.na(mean_response_m5_30))
    m <- p24_lmm(dd, paste(AXES, collapse = " + "), "mean_response_m5_30")
    fixef_table(m, keep = AXES) %>%
      mutate(chromophore = cc, roi = rr, n_rows = nrow(dd), .before = 1)
  })
}) %>% mutate(q_BH = bh(p), family_size = n(), window = "response")
print(as.data.frame(fam_resp %>%
      select(chromophore, roi, term, estimate, se, t, p, q_BH)), digits = 3)
write_outcome(fam_resp, file.path(OUT, "family_response_hbr_hbt.csv"))

cat("\nsurviving q < 0.05: difference", sum(fam_diff$q_BH < 0.05),
    "; response", sum(fam_resp$q_BH < 0.05), "\n")
cat("family minima: difference", min(fam_diff$q_BH),
    "; response", min(fam_resp$q_BH), "\n")

# ---------------------------------- 2. identity check: HbT slope == HbO + HbR
# HbT = HbO + HbR holds exactly in the data (builder guarantee); with subject
# fixed effects the slope identity then holds exactly by linearity of OLS.
rule("2. INTERNAL IDENTITY — HbT = HbO + HbR in data and in slopes")
row_gap <- map_dfr(ROIS, function(rr) {
  tibble(roi = rr,
         max_abs_row_gap = max(abs(d[[paste0("HbT_", rr)]] -
                                 d[[paste0("HbO_", rr)]] -
                                 d[[paste0("HbR_", rr)]]), na.rm = TRUE))
})
print(as.data.frame(row_gap), digits = 3)

ident_rows <- map_dfr(ROIS, function(rr) {
  dd <- d %>% filter(!is.na(.data[[paste0("HbT_", rr)]]))
  co <- map(c("HbO", "HbR", "HbT"),
            ~ coef(lm(as.formula(paste0(.x, "_", rr,
                      " ~ comfort_wi + richness_wi + subject_id")), data = dd)))
  ax <- AXES
  tibble(roi = rr, term = ax,
         hbt = vapply(ax, function(a) co[[3]][[a]], numeric(1)),
         hbo_plus_hbr = vapply(ax, function(a) co[[1]][[a]] + co[[2]][[a]],
                               numeric(1)))
}) %>% mutate(abs_diff = abs(hbt - hbo_plus_hbr))
print(as.data.frame(ident_rows), digits = 6)
cat("\nmax |slope identity gap|:", max(ident_rows$abs_diff), "uM\n")
write_outcome(ident_rows %>% left_join(row_gap, by = "roi"),
              file.path(OUT, "hbt_identity_check.csv"))

# ----------------- 3. slope-ratio coupling diagnostic (subject bootstrap)
# Per ROI, comfort axis: refit HbO and HbR on 1,000 subject-cluster resamples and
# look at the JOINT distribution of (beta_HbO, beta_HbR):
#   neurovascular deactivation  -> opposite signs, beta_HbR / beta_HbO ~ -1/2.3
#   blood-volume decrease       -> same sign,      beta_HbR / beta_HbO ~ +1/2.3
#   HbR silent                  -> beta_HbR ~ 0 regardless of beta_HbO
# Subject fixed effects (lm) inside the bootstrap: the point estimates coincide
# with the LMM's, and speed matters at 1,000 x 2 windows x 3 ROIs x 2 fits.
rule("3. HbO:HbR SLOPE-RATIO COUPLING — comfort axis, 1,000 subject bootstraps")

boot_pair <- function(data, y_hbo, y_hbr, n = BOOT_N, seed = 24) {
  dd <- data %>% filter(!is.na(.data[[y_hbo]]), !is.na(.data[[y_hbr]]))
  subs <- unique(dd$subject_id)
  set.seed(seed)
  out <- matrix(NA_real_, nrow = n, ncol = 2, dimnames = list(NULL, c("hbo", "hbr")))
  for (b in seq_len(n)) {
    pick <- sample(subs, length(subs), replace = TRUE)
    db <- map_dfr(seq_along(pick), function(i) dd %>% filter(subject_id == pick[i]))
    f1 <- try(lm(as.formula(paste(y_hbo, "~ comfort_wi + richness_wi + subject_id")),
                 data = db), silent = TRUE)
    f2 <- try(lm(as.formula(paste(y_hbr, "~ comfort_wi + richness_wi + subject_id")),
                 data = db), silent = TRUE)
    if (inherits(f1, "try-error") || inherits(f2, "try-error")) next
    out[b, ] <- c(coef(f1)["comfort_wi"], coef(f2)["comfort_wi"])
  }
  as_tibble(out) %>% filter(!is.na(hbo))
}

diag_rows <- list()
for (rr in ROIS) {
  bs <- boot_pair(d, paste0("HbO_", rr), paste0("HbR_", rr))
  ratio <- bs$hbr / bs$hbo
  diag_rows[[length(diag_rows) + 1]] <- tibble(
    window = "difference", roi = rr,
    n_boot = nrow(bs),
    hbo_median = median(bs$hbo), hbr_median = median(bs$hbr),
    p_hbo_negative = mean(bs$hbo < 0), p_hbr_positive = mean(bs$hbr > 0),
    p_opposite_signs = mean(bs$hbo * bs$hbr < 0),
    ratio_median = median(ratio[is.finite(ratio)]),
    ratio_q05 = quantile(ratio[is.finite(ratio)], 0.05),
    ratio_q95 = quantile(ratio[is.finite(ratio)], 0.95))
  cat(sprintf("\n%s / %s: HbO median %.4f (P<0 %.2f), HbR median %.4f (P>0 %.2f), P(opposite) %.2f, ratio median %.2f [%.2f, %.2f]\n",
              "difference", rr, median(bs$hbo), mean(bs$hbo < 0), median(bs$hbr),
              mean(bs$hbr > 0), mean(bs$hbo * bs$hbr < 0),
              median(ratio[is.finite(ratio)]),
              quantile(ratio[is.finite(ratio)], 0.05),
              quantile(ratio[is.finite(ratio)], 0.95)))
}

# response-window version: wide per-ROI HbO/HbR response means
resp_wide <- resp %>%
  select(subject_id, marker, chromophore, roi, mean_response_m5_30,
         comfort_wi, richness_wi) %>%
  mutate(key = paste0(chromophore, "_", roi)) %>%
  select(-chromophore, -roi) %>%
  pivot_wider(names_from = key, values_from = mean_response_m5_30)
for (rr in ROIS) {
  y_hbo <- paste0("HbO_", rr); y_hbr <- paste0("HbR_", rr)
  if (!all(c(y_hbo, y_hbr) %in% names(resp_wide))) next
  bs <- boot_pair(resp_wide, y_hbo, y_hbr)
  ratio <- bs$hbr / bs$hbo
  diag_rows[[length(diag_rows) + 1]] <- tibble(
    window = "response", roi = rr,
    n_boot = nrow(bs),
    hbo_median = median(bs$hbo), hbr_median = median(bs$hbr),
    p_hbo_negative = mean(bs$hbo < 0), p_hbr_positive = mean(bs$hbr > 0),
    p_opposite_signs = mean(bs$hbo * bs$hbr < 0),
    ratio_median = median(ratio[is.finite(ratio)]),
    ratio_q05 = quantile(ratio[is.finite(ratio)], 0.05),
    ratio_q95 = quantile(ratio[is.finite(ratio)], 0.95))
  cat(sprintf("\n%s / %s: HbO median %.4f (P<0 %.2f), HbR median %.4f (P>0 %.2f), P(opposite) %.2f, ratio median %.2f [%.2f, %.2f]\n",
              "response", rr, median(bs$hbo), mean(bs$hbo < 0), median(bs$hbr),
              mean(bs$hbr > 0), mean(bs$hbo * bs$hbr < 0),
              median(ratio[is.finite(ratio)]),
              quantile(ratio[is.finite(ratio)], 0.05),
              quantile(ratio[is.finite(ratio)], 0.95)))
}
diag <- bind_rows(diag_rows)
cat("\nphysiological benchmarks: neurovascular ~ -1/2.3 = -0.43; blood-volume ~ +0.43\n")
write_outcome(diag, file.path(OUT, "slope_ratio_coupling.csv"))

close_log(con)
cat("\nmod09 done.\n")
