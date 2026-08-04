# A12 — Repetition structure of the appraisal-haemodynamic relation
#
# Every condition is heard three times (once per sequence), and the ratings were
# collected once, AFTER the session, on re-listening. Two questions:
#   (i)  does the haemodynamic response habituate across the three hearings
#        (repetition suppression / session progression)? Because every sequence
#        contains all four markers, encounter number == run number, so repetition
#        and time-on-task are one contrast here and are reported as such.
#   (ii) is the appraisal->HbO coupling a first-hearing phenomenon? The primary
#        family (2 axes x 3 ROIs) is fitted separately at encounters 1, 2 and 3.
#   (iii) a single targeted contrast: Comfort slope at encounter 1 minus
#        encounters 2-3.
#
# Declared families: (i) 2 windows x 3 ROIs = 6;
# (ii) 3 encounters x 2 axes x 3 ROIs = 18, BH together; (iii) 3 ROIs = 3, BH.
# Single presentations are noisy (ICC 0.135), so (ii)/(iii) are read as
# coefficient patterns with subject-bootstrap support, never as three separate
# discoveries.

source("code/load_data.R")

OUT  <- file.path("output", "analysis_12_repetition_structure")
con  <- open_log(file.path(OUT, "repetition_structure_report.txt"))

rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")
AXES <- c("comfort_wi", "richness_wi")

# event-level ROI table: difference score wide, plus baseline/response from channels
blk <- load_p24_block() %>%          # 826 rows, difference-score ROI columns
  mutate(run_index = ceiling(event_index / 4)) %>%
  group_by(subject_id, marker) %>%
  mutate(encounter = rank(event_index)) %>%   # chronological hearing number
  ungroup()

bm_roi <- load_block_analysis("HbO") %>%
  left_join(load_channels() %>% filter(active_pair == 1) %>%
              distinct(subject_id, channel_index, roi),
            by = c("subject_id", "channel_index")) %>%
  filter(!is.na(roi)) %>%
  group_by(subject_id, event_index, roi) %>%
  summarise(baseline = mean(baseline_m10_0),
            response = mean(response_m5_30), .groups = "drop") %>%
  pivot_wider(names_from = roi,
              values_from = c(baseline, response), names_glue = "{.value}_{roi}")

blk <- blk %>% left_join(bm_roi, by = c("subject_id", "event_index"))

# gate: my per-event aggregation must satisfy difference == response - baseline
# against the delivered wide columns exactly (same channels, same mean)
gate <- max(abs((blk$response_PFC_Frontal - blk$baseline_PFC_Frontal) -
                blk$HbO_PFC_Frontal), na.rm = TRUE)
stopifnot("per-event aggregation gate failed" = gate < 1e-9)
cat("aggregation gate (response - baseline == delivered difference): max gap",
    format(gate, digits = 2), "\n")

# encounter == run for every event but the two lost-slot cells (P12/P14 lose M4's
# first slot, so that marker's later hearings rank one lower); analyse by ENCOUNTER
# (the declared variable) and record the divergence count.
n_div <- sum(blk$encounter != blk$run_index)
cat("events where encounter != run (lost-slot descendants):", n_div, "of",
    nrow(blk), "\n")
cat("events per encounter:", paste(table(blk$encounter), collapse = " / "), "\n")

# ============================================================ (i) habituation
rule("(i) ENCOUNTER TREND (repetition == session progression here) — family of 6")
hab <- map_dfr(ROIS, function(rr) {
  imap_dfr(c(response = paste0("response_", rr),
             baseline = paste0("baseline_", rr)), function(y, win) {
    dd <- blk %>% filter(!is.na(.data[[y]]))
    m <- p24_lmm(dd, "z(encounter)", y)
    fixef_table(m, keep = "z(encounter)") %>%
      mutate(window = win, roi = rr, n_events = nrow(dd), .before = 1)
  })
}) %>% mutate(q_BH = bh(p), family_size = n())
print(as.data.frame(hab %>% select(window, roi, estimate, se, t, p, q_BH)), digits = 3)
write_outcome(hab, file.path(OUT, "run_trend_by_window_roi.csv"))

# encounter means for the record
run_means <- blk %>%
  group_by(encounter) %>%
  summarise(across(c(HbO_PFC_Frontal, HbO_Left_Temporal, HbO_Right_Temporal,
                     response_PFC_Frontal, baseline_PFC_Frontal),
                   ~ mean(.x, na.rm = TRUE)), n = n(), .groups = "drop")
print(as.data.frame(run_means), digits = 3)
write_outcome(run_means, file.path(OUT, "run_means.csv"))

# ==================================================== (ii) per-encounter families
rule("(ii) APPRAISAL FAMILY PER ENCOUNTER — 3 x (2 axes x 3 ROIs) = 18, BH")
per_run <- map_dfr(1:3, function(k) {
  dd <- blk %>% filter(encounter == k)
  map_dfr(ROIS, function(rr) {
    y <- paste0("HbO_", rr)
    ddd <- dd %>% filter(!is.na(.data[[y]]))
    m <- p24_lmm(ddd, paste(AXES, collapse = " + "), y)
    fixef_table(m, keep = AXES) %>%
      mutate(encounter = k, roi = rr, n_events = nrow(ddd),
             n_subjects = n_distinct(ddd$subject_id), .before = 1)
  })
}) %>% mutate(q_BH = bh(p), family_size = n())
print(as.data.frame(per_run %>%
      filter(term == "comfort_wi") %>%
      select(encounter, roi, estimate, se, t, p, q_BH)), digits = 3)
cat("\nrichness:\n")
print(as.data.frame(per_run %>%
      filter(term == "richness_wi") %>%
      select(encounter, roi, estimate, se, t, p, q_BH)), digits = 3)
write_outcome(per_run, file.path(OUT, "per_encounter_family.csv"))

# ================================================== (iii) first-vs-later contrast
# per-encounter comfort slopes in one model (emtrends), then the linear contrast
# b1 - (b2 + b3)/2 with Satterthwaite df — no rank deficiency.
rule("(iii) COMFORT SLOPE: ENCOUNTER 1 vs 2-3 — family of 3")
suppressPackageStartupMessages(library(emmeans))
emm_options(lmer.df = "satterthwaite")
contrast <- map_dfr(ROIS, function(rr) {
  y <- paste0("HbO_", rr)
  dd <- blk %>% filter(!is.na(.data[[y]])) %>%
    mutate(enc_f = factor(encounter))
  m <- p24_lmm(dd, "enc_f + richness_wi + comfort_wi:enc_f", y)
  et <- emtrends(m, ~ enc_f, var = "comfort_wi")
  ct <- summary(contrast(et, list(first_vs_later = c(1, -0.5, -0.5))),
                infer = TRUE)
  tibble(roi = rr, estimate = ct$estimate, se = ct$SE, df = ct$df,
         t = ct$t.ratio, p = ct$p.value)
}) %>% mutate(q_BH = bh(p), family_size = n())
print(as.data.frame(contrast %>% select(roi, estimate, se, t, p, q_BH)), digits = 3)
write_outcome(contrast, file.path(OUT, "first_vs_later_contrast.csv"))

# descriptive post-hoc: the per-encounter pattern in (ii) is enc1 ~= enc2, enc3 -> 0,
# so also quote the third-vs-earlier contrast. NOT part of the declared family.
third_vs_earlier <- map_dfr(ROIS, function(rr) {
  y <- paste0("HbO_", rr)
  dd <- blk %>% filter(!is.na(.data[[y]])) %>%
    mutate(enc_f = factor(encounter))
  m <- p24_lmm(dd, "enc_f + richness_wi + comfort_wi:enc_f", y)
  et <- emtrends(m, ~ enc_f, var = "comfort_wi")
  ct <- summary(contrast(et, list(third_vs_earlier = c(-0.5, -0.5, 1))),
                infer = TRUE)
  tibble(roi = rr, estimate = ct$estimate, se = ct$SE, df = ct$df,
         t = ct$t.ratio, p = ct$p.value)
}) %>% mutate(note = "descriptive post-hoc contrast; no correction")
cat("\nDESCRIPTIVE (post-hoc, uncorrected): encounter 3 vs encounters 1-2\n")
print(as.data.frame(third_vs_earlier), digits = 3)
write_outcome(third_vs_earlier, file.path(OUT, "third_vs_earlier_descriptive.csv"))

# ------------------------------------------ bootstrap support for (ii) and (iii)
rule("BOOTSTRAP — comfort slope per encounter and the first-vs-later contrast")
boot_rows <- map_dfr(ROIS, function(rr) {
  y <- paste0("HbO_", rr)
  dd <- blk %>% filter(!is.na(.data[[y]])) %>%
    select(subject_id, encounter, all_of(y), comfort_wi, richness_wi)
  subs <- unique(dd$subject_id)
  set.seed(24)
  est <- matrix(NA_real_, BOOT_N, 5,
                dimnames = list(NULL, c("enc1", "enc2", "enc3", "contrast", "tve")))
  for (b in seq_len(BOOT_N)) {
    pick <- sample(subs, length(subs), replace = TRUE)
    db <- map_dfr(seq_along(pick), function(i) dd %>% filter(subject_id == pick[i]))
    slopes <- map_dbl(1:3, function(k) {
      f <- try(lm(as.formula(paste(y, "~ comfort_wi + richness_wi + subject_id")),
                  data = db %>% filter(encounter == k)), silent = TRUE)
      if (inherits(f, "try-error")) return(NA_real_)
      coef(f)["comfort_wi"]
    })
    est[b, 1:3] <- slopes
    est[b, 4] <- slopes[1] - mean(slopes[2:3])
    est[b, 5] <- slopes[3] - mean(slopes[1:2])
  }
  as_tibble(est) %>%
    summarise(across(everything(), list(
      median = ~ median(.x, na.rm = TRUE),
      ci_lo = ~ quantile(.x, 0.025, na.rm = TRUE),
      ci_hi = ~ quantile(.x, 0.975, na.rm = TRUE),
      p_neg = ~ mean(.x < 0, na.rm = TRUE)), .names = "{.col}_{.fn}")) %>%
    mutate(tve_p_pos = 1 - tve_p_neg, roi = rr, .before = 1)
})
print(as.data.frame(boot_rows), digits = 3)
write_outcome(boot_rows, file.path(OUT, "bootstrap_run_slopes.csv"))

close_log(con)
cat("\nA12 done.\n")
