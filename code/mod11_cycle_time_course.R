# A11 — Fine-grained cycle time course: the declared families on the 5-s bin
#        table
#
# Where in the 90-s presentation cycle does the appraisal association live, and
# does the prefrontal / channel-10 Comfort association have an HRF-plausible
# time course (absent at 0-5 s, rising into [10,30), decaying after clip end)
# or a step/constant one (artefact-like)? Declared family: 2 targets
# (prefrontal ROI; channel 10) x 14 bins = 28 within-subject Comfort-slope
# tests, BH, q < 0.05. Event-referenced bins are primary; as-recorded bins, the
# richness axis, and the other spatial summaries are descriptive. One growth
# contrast per target (slope in [15,30) minus [0,15)), family of 2, BH.
#
# PROVENANCE. This module reads the SHIPPED INTERMEDIATE
# data/intermediates/condition_bin_values_hbo.csv: the subject x condition x
# 5-s-bin HbO table (event-referenced and as-recorded values for the two
# primary targets and three descriptive summaries, with the appraisal
# decompositions attached). It is derived from the delivered continuous
# concentration series epoched on the frozen onsets; that series is not part of
# the deposit, so the series-extraction and validation-gate sections of the
# original pipeline cannot be re-run here and are not ported. The bin table
# passed the original pipeline's mandatory gate (the three frozen window means
# recomputed from the same series read reproduced the frozen block metrics to
# 2.46e-11), and mod18 re-asserts the inheritance independently before use.
#
# The HbR bin table is not shipped (the intermediate is HbO-only, as deposited),
# so the HbR descriptive family of the original pipeline cannot be regenerated
# here; every HbO family is computed below.

source("code/load_data.R")

OUT  <- file.path("output", "analysis_11_cycle_time_course")
con  <- open_log(file.path(OUT, "cycle_time_course_report.txt"))

rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")
AXES <- c("comfort_wi", "richness_wi")
PRIMARY <- c("PFC_Frontal", "ch10")

hbo <- load_intermediate("condition_bin_values_hbo.csv")
stopifnot(all(hbo$chromophore == "HbO"))

# bin geometry, recovered from the shipped table itself
bin_lut <- hbo %>% distinct(bin, bin_lo) %>% arrange(bin_lo)
BIN_LABELS <- bin_lut$bin
EDGES <- c(bin_lut$bin_lo, max(bin_lut$bin_lo) + 5)   # 14 half-open bins

design <- hbo %>%
  distinct(subject_id, marker, comfort, richness, comfort_wi, richness_wi,
           comfort_bw, richness_bw)

cat("shipped bin table:", n_distinct(hbo$subject_id), "subjects;",
    n_distinct(paste(hbo$subject_id, hbo$marker)), "condition cells;",
    n_distinct(hbo$bin), "bins;", n_distinct(hbo$target), "targets\n")

# --------------------------------------------- 2. the declared 28-test family
rule("2. COMFORT SLOPE PER BIN — event-referenced, family 2 targets x 14 bins = 28")
fam <- map_dfr(PRIMARY, function(tg) {
  map_dfr(BIN_LABELS, function(bb) {
    dd <- hbo %>% filter(target == tg, bin == bb)
    m <- p24_lmm(dd, paste(AXES, collapse = " + "), "referenced")
    fixef_table(m, keep = "comfort_wi") %>%
      mutate(target = tg, bin = bb, bin_lo = EDGES[which(BIN_LABELS == bb)],
             n_rows = nrow(dd), n_subjects = n_distinct(dd$subject_id), .before = 1)
  })
}) %>% mutate(q_BH = bh(p), family_size = n())
print(as.data.frame(fam %>%
      select(target, bin, estimate, se, t, p, q_BH, n_subjects)), digits = 3)
write_outcome(fam, file.path(OUT, "bin_slopes_comfort.csv"))
cat("\nsurviving q < 0.05:", sum(fam$q_BH < 0.05), "of", nrow(fam), "\n")

# --------------------------------------------- descriptive: richness + raw + others
rule("3. DESCRIPTIVE — richness axis, as-recorded levels, other targets")
desc_rich <- map_dfr(PRIMARY, function(tg) {
  map_dfr(BIN_LABELS, function(bb) {
    dd <- hbo %>% filter(target == tg, bin == bb)
    m <- p24_lmm(dd, paste(AXES, collapse = " + "), "referenced")
    fixef_table(m, keep = "richness_wi") %>%
      mutate(target = tg, bin = bb, bin_lo = EDGES[which(BIN_LABELS == bb)], .before = 1)
  })
})
write_outcome(desc_rich, file.path(OUT, "bin_slopes_richness_descriptive.csv"))

desc_all <- map_dfr(c("Left_Temporal (descriptive)", "Right_Temporal (descriptive)",
                      "head (descriptive)"), function(tg) {
  map_dfr(BIN_LABELS, function(bb) {
    dd <- hbo %>% filter(target == tg, bin == bb)
    m <- p24_lmm(dd, paste(AXES, collapse = " + "), "referenced")
    fixef_table(m, keep = AXES) %>%
      mutate(target = tg, bin = bb, bin_lo = EDGES[which(BIN_LABELS == bb)], .before = 1)
  })
})
write_outcome(desc_all, file.path(OUT, "bin_slopes_other_targets_descriptive.csv"))

raw_desc <- map_dfr(PRIMARY, function(tg) {
  map_dfr(BIN_LABELS, function(bb) {
    dd <- hbo %>% filter(target == tg, bin == bb)
    m <- p24_lmm(dd, paste(AXES, collapse = " + "), "raw")
    fixef_table(m, keep = "comfort_wi") %>%
      mutate(target = tg, bin = bb, bin_lo = EDGES[which(BIN_LABELS == bb)], .before = 1)
  })
})
write_outcome(raw_desc, file.path(OUT, "bin_slopes_asrecorded_descriptive.csv"))

# --------------------------------------------- 4. growth contrast, family of 2
rule("4. GROWTH CONTRAST — slope in [15,30) minus [0,15), 2 targets, BH")
macro <- hbo %>%
  filter(target %in% PRIMARY, bin_lo >= 0, bin_lo < 30) %>%
  mutate(phase = ifelse(bin_lo < 15, "early", "late")) %>%
  group_by(subject_id, marker, target, phase) %>%
  summarise(referenced = mean(referenced), .groups = "drop") %>%
  left_join(design, by = c("subject_id", "marker"))

growth <- map_dfr(PRIMARY, function(tg) {
  dd <- macro %>% filter(target == tg)
  slopes <- map_dfr(c("early", "late"), function(ph) {
    m <- p24_lmm(dd %>% filter(phase == ph), paste(AXES, collapse = " + "), "referenced")
    fixef_table(m, keep = "comfort_wi") %>% mutate(phase = ph)
  })
  # bootstrap the late-minus-early difference (subject cluster, subject FE lm)
  wide <- dd %>% select(subject_id, marker, phase, referenced,
                        comfort_wi, richness_wi) %>%
    pivot_wider(names_from = phase, values_from = referenced)
  subs <- unique(wide$subject_id)
  set.seed(24)
  diffs <- map_dbl(seq_len(BOOT_N), function(b) {
    pick <- sample(subs, length(subs), replace = TRUE)
    db <- map_dfr(seq_along(pick), function(i) wide %>% filter(subject_id == pick[i]))
    s <- map_dbl(c("early", "late"), function(ph) {
      f <- try(lm(as.formula(paste(ph, "~ comfort_wi + richness_wi + subject_id")),
                  data = db), silent = TRUE)
      if (inherits(f, "try-error")) return(NA_real_)
      coef(f)["comfort_wi"]
    })
    s[2] - s[1]
  })
  diffs <- diffs[!is.na(diffs)]
  tibble(target = tg,
         slope_early = slopes$estimate[slopes$phase == "early"],
         slope_late = slopes$estimate[slopes$phase == "late"],
         diff_median = median(diffs),
         ci_lo = quantile(diffs, 0.025), ci_hi = quantile(diffs, 0.975),
         p_boot = 2 * min(mean(diffs <= 0), mean(diffs >= 0)))
}) %>% mutate(q_BH = bh(p_boot), family_size = n())
print(as.data.frame(growth), digits = 3)
write_outcome(growth, file.path(OUT, "growth_contrast.csv"))

close_log(con)
cat("\nA11 done.\n")
