# =============================================================================
# Module 04 — Does subjective appraisal predict cortical haemodynamics?
#
# The campaign's primary question, run under the constraints the earlier audit
# established:
#
#   - the primary DV's comfort slope is 68% pre-stimulus baseline, so the
#     difference score cannot carry a stimulus-evoked reading on its own;
#   - M6 (nominal Q4) never opens a run, so cycle position is confounded with
#     the comfortable low-richness condition and pushes the comfort slope
#     negative for a purely temporal reason;
#   - channels 18/19 are dead and halve the left-temporal ROI.
#
# Declared family: 2 axes x 3 ROIs = 6 within-subject slope tests per analysis
# window, Benjamini-Hochberg, significance at q < 0.05.
#
# Outputs (output/analysis_04_appraisal_to_haemodynamics/):
#   family_frozen_dv.csv          declared primary family on the frozen DV
#   family_by_window.csv          the same family on each window separately
#   cycle_control_block_level.csv cycle position as a block-level covariate
#   cycle_control_drop_m6.csv     family refit without M6
#   cycle_control_drop_openers.csv family refit without run-opening events
#   left_temporal_corrected.csv   frozen versus corrected left-temporal ROI
#   lateralisation.csv            left-minus-right t-tests
#   bootstrap_ci.csv              cluster bootstrap on the key slopes
#
# Purely exploratory sections of the original pass (carry-over lags,
# stimulus-versus-idiosyncratic decomposition, head-mean summaries) are not
# part of the confirmatory displays and are not ported.
#
# Run from the repository root: Rscript code/mod04_confirmatory_family.R
# =============================================================================

source("code/load_data.R")

OUT <- file.path("output", "analysis_04_appraisal_to_haemodynamics")
con <- open_log(file.path(OUT, "mod04_confirmatory_family_report.txt"))

d  <- load_p24()
bm <- load_block_metrics(active_only = TRUE)

rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")
cat("Module 04 — appraisal to haemodynamics. All slopes in uM per unit of the axis\n")
cat("(axis range is about +/-0.93, SD 0.38), unless a z-scored DV is stated.\n")

AXES <- c("comfort_wi", "richness_wi")

#' The declared family: one model per ROI, both axes in it, BH across the 6.
family_scan <- function(data, dv_prefix, label, extra = NULL) {
  rows <- map_dfr(ROIS, function(rr) {
    y <- paste0(dv_prefix, rr)
    if (!y %in% names(data)) return(tibble())
    fixed <- paste(c(AXES, extra), collapse = " + ")
    m <- p24_lmm(data, fixed, y)
    fixef_table(m, keep = AXES) %>% mutate(roi = rr, .before = 1)
  })
  rows <- rows %>% mutate(q_BH = bh(p), family_size = nrow(rows), analysis = label)
  rows
}

# ================================================ 1. the declared primary family
rule("1. PRIMARY FAMILY — frozen DV (response minus baseline), 2 axes x 3 ROIs")
fam1 <- family_scan(d, "HbO_", "frozen DV")
print(as.data.frame(fam1 %>% select(roi, term, estimate, se, t, p, q_BH)), digits = 3)
write_outcome(fam1, file.path(OUT, "family_frozen_dv.csv"))
cat("\nfamily size", unique(fam1$family_size), "; surviving q < 0.05:",
    sum(fam1$q_BH < 0.05), "\n")

# ============================================ 2. which window carries the slope
rule("2. WHICH WINDOW CARRIES IT? (the audit showed 68% is pre-stimulus)")
roi_blocks <- bm %>% filter(chromophore == "HbO") %>%
  left_join(load_channels() %>% filter(active_pair == 1) %>%
              distinct(subject_id, channel_index, roi),
            by = c("subject_id", "channel_index")) %>%
  filter(!is.na(roi)) %>%
  group_by(subject_id, event_index, sequence, marker, onset_sec, roi) %>%
  summarise(baseline = mean(baseline_m10_0), response = mean(response_m5_30),
            difference = mean(response_minus_baseline_m5_30), .groups = "drop")

cells <- roi_blocks %>% group_by(subject_id, marker, roi) %>%
  summarise(across(c(baseline, response, difference), mean), .groups = "drop") %>%
  pivot_wider(names_from = roi, values_from = c(baseline, response, difference),
              names_sep = "_") %>%
  left_join(d %>% select(subject_id, marker, all_of(AXES), comfort, richness,
                         comfort_stim, comfort_idio, richness_stim, richness_idio,
                         group, form, quadrant),
            by = c("subject_id", "marker"))

win_fam <- map_dfr(c("difference_", "response_", "baseline_"), function(pref) {
  family_scan(cells, pref, sub("_$", "", pref))
})
print(as.data.frame(win_fam %>% select(analysis, roi, term, estimate, se, p, q_BH)),
      digits = 3)
write_outcome(win_fam, file.path(OUT, "family_by_window.csv"))
cat("\nEach window is its own family of 6. 'baseline' is the control: the clip has\n")
cat("not played yet, so a slope there is not a stimulus-evoked response.\n")

# ============================== 3. the cycle-position confound (M6 never opens)
rule("3. THE CYCLE-POSITION CONFOUND — three ways of removing it")
roi_blocks <- roi_blocks %>%
  mutate(opens = as.integer(event_index %in% c(1, 5, 9)))
cat("run-opening events by condition (M6 can never open a run):\n")
print(as.data.frame(roi_blocks %>% filter(roi == "PFC_Frontal") %>%
        count(marker, opens) %>%
        pivot_wider(names_from = opens, values_from = n, names_prefix = "opens_")))

cat("\n(a) block-level model with cycle position as a covariate\n")
blk <- roi_blocks %>%
  left_join(d %>% select(subject_id, marker, all_of(AXES)), by = c("subject_id", "marker"))
a3a <- map_dfr(ROIS, function(rr) {
  dd <- blk %>% filter(roi == rr)
  m <- p24_lmm(dd, paste(c(AXES, "opens"), collapse = " + "), "difference")
  fixef_table(m, keep = c(AXES, "opens")) %>% mutate(roi = rr, .before = 1)
})
a3a <- a3a %>% mutate(q_BH = ifelse(term %in% AXES, bh(p[term %in% AXES])[
  match(paste(roi, term), paste(roi, term)[term %in% AXES])], NA_real_))
print(as.data.frame(a3a %>% select(roi, term, estimate, se, p)), digits = 3)
write_outcome(a3a, file.path(OUT, "cycle_control_block_level.csv"))

cat("\n(b) drop M6 — the three remaining conditions are balanced on cycle position\n")
cells_noM6 <- cells %>% filter(marker != "M6") %>%
  group_by(subject_id) %>%
  mutate(comfort_wi = comfort - mean(comfort), richness_wi = richness - mean(richness)) %>%
  ungroup()
cat("   n =", nrow(cells_noM6), "cells,", n_distinct(cells_noM6$subject_id), "subjects",
    "(axes re-centred within the retained conditions)\n")
b_fam <- map_dfr(c("difference_", "response_", "baseline_"), function(pref) {
  family_scan(cells_noM6, pref, paste0("M6 dropped: ", sub("_$", "", pref)))
})
print(as.data.frame(b_fam %>% select(analysis, roi, term, estimate, se, p, q_BH)), digits = 3)
write_outcome(b_fam, file.path(OUT, "cycle_control_drop_m6.csv"))

cat("\n(c) drop the run-opening presentations, keeping all four conditions\n")
cells_hom <- roi_blocks %>% filter(opens == 0) %>%
  group_by(subject_id, marker, roi) %>%
  summarise(across(c(baseline, response, difference), mean),
            n_blocks = n(), .groups = "drop") %>%
  pivot_wider(names_from = roi, values_from = c(baseline, response, difference, n_blocks),
              names_sep = "_") %>%
  left_join(d %>% select(subject_id, marker, comfort, richness),
            by = c("subject_id", "marker")) %>%
  group_by(subject_id) %>%
  mutate(comfort_wi = comfort - mean(comfort), richness_wi = richness - mean(richness)) %>%
  ungroup()
cat("   n =", nrow(cells_hom), "cells; presentations per cell:",
    paste(names(table(cells_hom$n_blocks_PFC_Frontal)), collapse = "/"), "\n")
c_fam <- map_dfr(c("difference_", "response_", "baseline_"), function(pref) {
  family_scan(cells_hom, pref, paste0("openers dropped: ", sub("_$", "", pref)))
})
print(as.data.frame(c_fam %>% select(analysis, roi, term, estimate, se, p, q_BH)), digits = 3)
write_outcome(c_fam, file.path(OUT, "cycle_control_drop_openers.csv"))

# ============================================ 4. the corrected left-temporal ROI
rule("4. CORRECTED LEFT-TEMPORAL ROI (channels 18/19 dropped) — PROVISIONAL")
cat("First, does the rebuild reproduce the delivered ROI columns when nothing is\n")
cat("excluded? If not, the aggregation differs from the builder's and nothing\n")
cat("below is interpretable.\n\n")
print(as.data.frame(check_roi_rebuild()), digits = 3)

dc <- load_p24_roi_corrected()
cat("\nfrozen versus corrected left-temporal, same model:\n")
cmp <- map_dfr(c("HbO_Left_Temporal", "HbOc_Left_Temporal"), function(y) {
  m <- p24_lmm(dc, paste(AXES, collapse = " + "), y)
  fixef_table(m, keep = AXES) %>% mutate(dv = y, .before = 1)
})
print(as.data.frame(cmp), digits = 3)
write_outcome(cmp, file.path(OUT, "left_temporal_corrected.csv"))

cat("\nlateralisation: left minus right, frozen versus corrected\n")
dc <- dc %>% mutate(lat_frozen = HbO_Left_Temporal - HbO_Right_Temporal,
                    lat_corrected = HbOc_Left_Temporal - HbOc_Right_Temporal)
lat <- map_dfr(c("lat_frozen", "lat_corrected"), function(y) {
  x <- dc[[y]]; x <- x[!is.na(x)]
  tt <- t.test(x)
  tibble(dv = y, n = length(x), mean = mean(x), t = unname(tt$statistic), p = tt$p.value)
})
print(as.data.frame(lat), digits = 3)
write_outcome(lat, file.path(OUT, "lateralisation.csv"))
cat("\nNOTE handedness was never collected, so no lateralisation claim can be made\n")
cat("from these data regardless of which ROI definition is used.\n")

# ===================================================== 5. bootstrap the survivors
rule("5. CLUSTER BOOTSTRAP (1,000 resamples of subjects) ON THE KEY SLOPES")
boot_targets <- tribble(
  ~data,   ~dv,                     ~label,
  "cells", "difference_PFC_Frontal", "frozen DV, prefrontal",
  "cells", "response_PFC_Frontal",   "response window only, prefrontal",
  "cells", "baseline_PFC_Frontal",   "pre-stimulus baseline, prefrontal",
)
bt <- pmap_dfr(boot_targets, function(data, dv, label) {
  boot_lmm(cells, paste(AXES, collapse = " + "), dv, terms = AXES) %>%
    mutate(dv = dv, label = label, .before = 1)
})
print(as.data.frame(bt), digits = 3)
write_outcome(bt, file.path(OUT, "bootstrap_ci.csv"))

close_log(con)
cat("Module 04 done.\n")
