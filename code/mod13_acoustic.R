# mod13_acoustic.R — Stimulus-acoustic drivers of the late-window response
#
# Every association elsewhere in the campaign runs through the POST-SESSION
# appraisal, whose causal direction is ambiguous (sound -> appraisal -> brain,
# or state -> both). This analysis bypasses appraisal: do the PHYSICAL features
# of the 16 condition clips predict the [20,30) event-referenced response?
# Sound -> brain has an unambiguous direction, so a predictive feature is the
# first strictly evoked evidence available without rebuilding the pipeline.
#
# Features (computed on the mono mixdown of each 30-s WAV): RMS level in dB
# (relative; no absolute level exists or is reported), spectral centroid (Hz),
# normalised spectral flux, low-frequency power ratio (< 200 Hz). Declared
# family: 4 features x 2 targets (prefrontal ROI; whole head — a declared
# target here, not a rule-2 claim carrier) = 8 z(feature) slope tests, BH.
# Descriptive: ch10 / LT / RT, [5,30) robustness, feature correlation matrix,
# clip-level correlations, joint models with comfort_wi.
#
# DATA NOTE. The stimulus WAVs are not part of the deposit, so the features are
# not recomputed here: the per-clip feature table is shipped as the intermediate
# data/intermediates/clip_acoustic_features.csv and copied through to the output
# (clip_acoustic_features.csv) so every downstream table regenerates from the
# same numbers. The [20,30) outcome comes from the validated condition bin
# table, shipped as data/intermediates/condition_bin_values_hbo.csv.
#
# Outputs (output/analysis_13_stimulus_acoustic_drivers/):
#   clip_acoustic_features.csv (copy-through of the shipped intermediate),
#   family_acoustic_slopes_20_30.csv, joint_feature_comfort_descriptive.csv,
#   clip_level_correlations_descriptive.csv, feature_correlations.csv,
#   acoustic_slopes_descriptive.csv
# Display anchor: main-text Table 2 and SI Table S1 quote the family minimum
# (min q 0.289).

if (!exists("load_p24")) source("code/load_data.R")

OUT <- file.path("output", "analysis_13_stimulus_acoustic_drivers")
con <- open_log(file.path(OUT, "stimulus_acoustic_drivers_report.txt"))

rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")

# ================================================== 1. features: shipped intermediate
rule("1. CLIP ACOUSTIC FEATURES (shipped intermediate — WAVs not deposited)")
stim <- load_stimuli()
feats <- load_intermediate("clip_acoustic_features.csv")
stopifnot(nrow(feats) == 16, all(feats$clip_id %in% stim$clip_id))
cat("features for", nrow(feats), "design cells;",
    n_distinct(feats$clip_id), "clips\n")
print(as.data.frame(feats %>% select(-sr, -bits)), digits = 4)
write_outcome(feats, file.path(OUT, "clip_acoustic_features.csv"))

fmat <- feats %>% select(rms_db, spec_centroid_hz, spec_flux, lowfreq_ratio)
cat("\nfeature correlation matrix (16 cells, 15 unique waveforms):\n")
print(round(cor(fmat), 3))
write_outcome(as.data.frame(cor(fmat)) %>% rownames_to_column("feature"),
              file.path(OUT, "feature_correlations.csv"))

# ================================================== 2. outcome: [20,30) event-referenced
rule("2. OUTCOME WINDOWS FROM THE CONDITION BIN TABLE")
bins <- load_intermediate("condition_bin_values_hbo.csv")
resp2030 <- bins %>%
  filter(bin_lo %in% c(20, 25)) %>%
  summarise(.by = c(subject_id, marker, target),
            resp2030 = mean(referenced),
            comfort_wi = first(comfort_wi), richness_wi = first(richness_wi),
            comfort = first(comfort), richness = first(richness))
resp530 <- bins %>%
  filter(bin_lo %in% c(5, 10, 15, 20, 25)) %>%
  summarise(.by = c(subject_id, marker, target), resp530 = mean(referenced))
cells <- resp2030 %>% left_join(resp530, by = c("subject_id", "marker", "target"))
cat("cells:", nrow(cells), " (", n_distinct(cells$subject_id), "subjects x",
    n_distinct(cells$marker), "markers x", n_distinct(cells$target), "targets )\n")

design <- load_p24() %>% distinct(subject_id, group)
cells <- cells %>%
  left_join(design, by = "subject_id") %>%
  left_join(stim %>% select(group, marker, clip_id), by = c("group", "marker")) %>%
  left_join(feats %>% select(clip_id, rms_db, spec_centroid_hz, spec_flux, lowfreq_ratio),
            by = "clip_id")
stopifnot(!any(is.na(cells$rms_db)))

FEATURES <- c("rms_db", "spec_centroid_hz", "spec_flux", "lowfreq_ratio")
TARGETS_PRIMARY <- c("PFC_Frontal", "head (descriptive)")
target_label <- c("PFC_Frontal" = "prefrontal ROI", "head (descriptive)" = "whole head")

# ================================================== 3. declared family: 4 x 2 = 8, BH
rule("3. DECLARED FAMILY — 4 FEATURES x 2 TARGETS = 8 TESTS, BH")
family <- map_dfr(TARGETS_PRIMARY, function(tg) {
  map_dfr(FEATURES, function(f) {
    dd <- cells %>% filter(target == tg)
    m <- p24_lmm(dd %>% mutate(fv = as.numeric(scale(.data[[f]]))), "fv", "resp2030")
    fixef_table(m, keep = "fv") %>%
      mutate(target = target_label[[tg]], feature = f, window = "[20,30) referenced")
  })
}) %>% mutate(q_BH = bh(p), family_size = n())
print(as.data.frame(family %>%
      select(target, feature, estimate, se, t, p, q_BH)), digits = 3)
write_outcome(family, file.path(OUT, "family_acoustic_slopes_20_30.csv"))
cat("\nsurviving q < 0.05:", sum(family$q_BH < 0.05), "of", nrow(family), "\n")
cat("family minimum:", min(family$q_BH), "\n")

# ================================================== 4. descriptive layers
rule("4. DESCRIPTIVE — OTHER TARGETS AND THE [5,30) RESPONSE WINDOW")
descr_targets <- map_dfr(c("ch10", "Left_Temporal (descriptive)",
                           "Right_Temporal (descriptive)"), function(tg) {
  map_dfr(FEATURES, function(f) {
    dd <- cells %>% filter(target == tg)
    m <- p24_lmm(dd %>% mutate(fv = as.numeric(scale(.data[[f]]))), "fv", "resp2030")
    fixef_table(m, keep = "fv") %>% mutate(target = tg, feature = f, window = "[20,30)")
  })
})
descr_530 <- map_dfr(TARGETS_PRIMARY, function(tg) {
  map_dfr(FEATURES, function(f) {
    dd <- cells %>% filter(target == tg)
    m <- p24_lmm(dd %>% mutate(fv = as.numeric(scale(.data[[f]]))), "fv", "resp530")
    fixef_table(m, keep = "fv") %>% mutate(target = target_label[[tg]], feature = f,
                                           window = "[5,30) referenced")
  })
})
descr <- bind_rows(descr_targets, descr_530) %>% mutate(note = "descriptive, no correction")
print(as.data.frame(descr %>% select(target, feature, window, estimate, p)), digits = 3)
write_outcome(descr, file.path(OUT, "acoustic_slopes_descriptive.csv"))

rule("5. DESCRIPTIVE — JOINT MODELS: FEATURE AND COMFORT TOGETHER")
joint <- map_dfr(TARGETS_PRIMARY, function(tg) {
  map_dfr(FEATURES, function(f) {
    dd <- cells %>% filter(target == tg) %>%
      mutate(fv = as.numeric(scale(.data[[f]])))
    m <- p24_lmm(dd, "fv + comfort_wi", "resp2030")
    fixef_table(m, keep = c("fv", "comfort_wi")) %>%
      mutate(target = target_label[[tg]], feature = f)
  })
}) %>% mutate(note = "descriptive, no correction")
print(as.data.frame(joint %>% select(target, feature, term, estimate, p)), digits = 3)
write_outcome(joint, file.path(OUT, "joint_feature_comfort_descriptive.csv"))

rule("6. DESCRIPTIVE — CLIP-LEVEL CORRELATIONS (n = 16 CELLS)")
clip_lvl <- cells %>%
  filter(target == "PFC_Frontal") %>%
  summarise(.by = c(clip_id, group, marker),
            resp = mean(resp2030), comfort = mean(comfort), richness = mean(richness)) %>%
  left_join(feats %>% select(clip_id, rms_db, spec_centroid_hz, spec_flux, lowfreq_ratio),
            by = "clip_id")
clip_cors <- map_dfr(FEATURES, function(f) {
  tibble(feature = f,
         r_response = cor(clip_lvl[[f]], clip_lvl$resp),
         r_comfort = cor(clip_lvl[[f]], clip_lvl$comfort),
         r_richness = cor(clip_lvl[[f]], clip_lvl$richness))
})
print(as.data.frame(clip_cors), digits = 3)
write_outcome(clip_cors, file.path(OUT, "clip_level_correlations_descriptive.csv"))

close_log(con)
cat("\nmod13 done.\n")
