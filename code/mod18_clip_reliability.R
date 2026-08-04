# A18 — Neural clip reliability, and the adjective screen
#
# Two closeouts.
#
# (i) NEURAL CLIP RELIABILITY. Listeners agree about clips at ICC(2,k)
#     0.96/0.92 on the subjective side, and clip explains 61% of comfort
#     variance. The neural counterpart here is model-free — no HRF assumed, no
#     appraisal involved: is the response SHAPE reproducible across people
#     hearing the same recording?
#
#     Subjects in the same `group` heard the same four clips, so marker <-> clip
#     is 1:1 within a group. Statistic: mean between-subject correlation of the
#     event-locked 0-60 s course for the SAME clip, minus the mean for DIFFERENT
#     clips, over same-group subject pairs. Null: permute each subject's four clip
#     labels (1,000). Family: 3 ROIs + whole head = 4, BH.
#
# (ii) ADJECTIVE SCREEN. An axis is a composite, and a single adjective could
#     couple more strongly than either axis. 8 adjectives x 3 ROIs = 24
#     within-subject slope tests on the response window, BH.
#
# Provenance: part (i) inherits the shipped subject x condition x 5-s-bin HbO
# table (data/intermediates/condition_bin_values_hbo.csv; see mod11's
# provenance note). Section 0 re-asserts that inheritance with an independent
# check against the delivered difference-score columns before anything is read.

source("code/load_data.R")

OUT  <- file.path("output", "analysis_18_clip_reliability_and_adjectives")
con  <- open_log(file.path(OUT, "clip_reliability_and_adjectives_report.txt"))
rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")

set.seed(24)
N_PERM <- BOOT_N
BINS_POST <- c("[0,5)", "[5,10)", "[10,15)", "[15,20)", "[20,25)", "[25,30)",
               "[30,35)", "[35,40)", "[40,45)", "[45,50)", "[50,55)", "[55,60)")
UNITS <- c("PFC_Frontal", "Left_Temporal (descriptive)",
           "Right_Temporal (descriptive)", "head (descriptive)")
UNIT_LABEL <- c("PFC_Frontal" = "Prefrontal", "Left_Temporal (descriptive)" = "Left temporal",
                "Right_Temporal (descriptive)" = "Right temporal",
                "head (descriptive)" = "Whole head")

bins <- load_intermediate("condition_bin_values_hbo.csv")

# Subject x condition x ROI window means joined to the canonical design, from
# the condition x channel sheet (mean over each ROI's active channels).
load_condition_roi_analysis <- function(chromophore = "HbO") {
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

# ------------------------------------------------------------------ 0. the gate
rule("0. INHERITANCE CHECK — the bins must reproduce the delivered difference score")

# mean of the referenced bins over [5,30) == response_m5_30 - baseline_m10_0,
# which is exactly the delivered HbO_PFC_Frontal column.
chk <- bins %>%
  filter(target == "PFC_Frontal",
         bin %in% c("[5,10)", "[10,15)", "[15,20)", "[20,25)", "[25,30)")) %>%
  group_by(subject_id, marker) %>% summarise(rebuilt = mean(referenced), .groups = "drop") %>%
  inner_join(load_p24() %>% select(subject_id, marker, frozen = HbO_PFC_Frontal),
             by = c("subject_id", "marker")) %>%
  mutate(absdiff = abs(rebuilt - frozen))
cat("cells compared:", nrow(chk), "  worst absolute difference:",
    signif(max(chk$absdiff), 3), "uM\n")
gate_ok <- max(chk$absdiff) < 1e-9
cat("INHERITANCE GATE PASSED:", gate_ok, "\n")
write_outcome(chk %>% summarise(n = n(), max_abs_diff = max(absdiff)),
              file.path(OUT, "gate_bin_inheritance.csv"))
if (!gate_ok) stop("Inheritance gate failed — the bin table does not reproduce the delivered columns.")

# ---------------------------------------------------- 1. neural clip reliability
rule("1. NEURAL CLIP RELIABILITY — does the same recording produce the same course?")

subj_group <- load_subjects() %>% select(subject_id, group)

# WAVEFORM IDENTITY, not group x marker. Groups 1 and 3 share their Q4 source
# recording, so a group-based definition would score that genuinely-same audio
# as "different".
wave_id <- load_p24() %>%
  select(subject_id, marker, clip_id, flag_shared_q4_waveform) %>%
  mutate(waveform = if_else(flag_shared_q4_waveform == 1, "shared_Q4_wdt16", clip_id))
cat("design cells:", n_distinct(paste(wave_id$clip_id)),
    " distinct waveforms:", n_distinct(wave_id$waveform), "\n")

isc_one <- function(unit) {
  w <- bins %>% filter(target == unit, bin %in% BINS_POST) %>%
    select(subject_id, marker, bin, referenced) %>%
    mutate(bin = factor(bin, levels = BINS_POST)) %>%
    arrange(subject_id, marker, bin) %>%
    pivot_wider(names_from = bin, values_from = referenced) %>%
    inner_join(subj_group, by = "subject_id") %>%
    inner_join(wave_id %>% select(subject_id, marker, waveform),
               by = c("subject_id", "marker"))
  M <- as.matrix(w[, BINS_POST])
  ok <- stats::complete.cases(M) & apply(M, 1, function(v) sd(v) > 0)
  w <- w[ok, ]; M <- M[ok, , drop = FALSE]

  # row-standardise once: correlation between any two courses is then a dot product
  Z <- t(scale(t(M)))
  C <- (Z %*% t(Z)) / (ncol(Z) - 1)

  n <- nrow(w)
  same_sub <- outer(w$subject_id, w$subject_id, "==")
  pairs    <- !same_sub & upper.tri(C)          # any two courses from different people

  # THE CONFOUND, and the control for it.
  # Marker is partly confounded with within-cycle position for EVERY subject
  # identically: with S1 = M3-M4-M5-M6, S2 = M5-M3-M6-M4, S3 = M4-M6-M5-M3, each
  # marker occupies a FIXED set of within-sequence positions (M3 {1,2,4}, M4
  # {2,4,1}, M5 {3,1,3}, M6 {4,3,2}) whatever the form, since every subject hears
  # all three sequences. So a same-marker pair shares cycle position as well as
  # audio, and a within-subject label permutation breaks both at once.
  #
  # The discriminator is the between-subjects stimulus set. Same marker in a
  # DIFFERENT group means the SAME cycle position but a DIFFERENT recording.
  #   same_wave  : same audio  + same position   (the effect of interest)
  #   pos_ctrl   : different audio + same position (position alone)
  #   diff_marker: different audio + different position
  # Audio-driven reliability requires same_wave > pos_ctrl, not merely
  # same_wave > diff_marker.
  stat_for <- function(lab_wave, lab_mark) {
    same_wave <- outer(lab_wave, lab_wave, "==")
    same_mark <- outer(lab_mark, lab_mark, "==")
    a  <- mean(C[pairs & same_wave])                       # same audio
    pc <- mean(C[pairs & same_mark & !same_wave])          # same position only
    b  <- mean(C[pairs & !same_mark])                      # neither
    c(audio_vs_position = a - pc, audio_vs_diff = a - b, position_vs_diff = pc - b,
      r_same_wave = a, r_pos_ctrl = pc, r_diff = b)
  }
  obs <- stat_for(w$waveform, w$marker)

  # Null: permute each subject's four condition labels, carrying waveform and
  # marker together (they are one label — the condition the subject heard).
  idx_by_sub <- split(seq_len(n), w$subject_id)
  null <- matrix(NA_real_, N_PERM, 3,
                 dimnames = list(NULL, c("audio_vs_position", "audio_vs_diff",
                                         "position_vs_diff")))
  for (j in seq_len(N_PERM)) {
    lw <- w$waveform; lm_ <- w$marker
    for (i in idx_by_sub) {
      o <- sample.int(length(i)); lw[i] <- w$waveform[i][o]; lm_[i] <- w$marker[i][o]
    }
    null[j, ] <- stat_for(lw, lm_)[1:3]
  }
  tibble(unit = unname(UNIT_LABEL[unit]), n_courses = n,
         r_same_wave = obs[["r_same_wave"]], r_pos_ctrl = obs[["r_pos_ctrl"]],
         r_diff = obs[["r_diff"]],
         audio_vs_position = obs[["audio_vs_position"]],
         p_audio_vs_position = (1 + sum(null[, "audio_vs_position"] >=
                                          obs[["audio_vs_position"]])) / (N_PERM + 1),
         audio_vs_diff = obs[["audio_vs_diff"]],
         p_audio_vs_diff = (1 + sum(null[, "audio_vs_diff"] >=
                                      obs[["audio_vs_diff"]])) / (N_PERM + 1),
         position_vs_diff = obs[["position_vs_diff"]],
         p_position_vs_diff = (1 + sum(null[, "position_vs_diff"] >=
                                         obs[["position_vs_diff"]])) / (N_PERM + 1))
}

isc <- map_dfr(UNITS, isc_one) %>%
  mutate(q_audio_vs_position = bh(p_audio_vs_position),
         q_audio_vs_diff     = bh(p_audio_vs_diff))

cat("\nMean between-subject correlation of the 0-60 s course, by pair type:\n")
print(as.data.frame(isc %>% select(unit, n_courses, r_same_wave, r_pos_ctrl, r_diff)),
      digits = 3, row.names = FALSE)

cat("\nDECLARED TEST (same audio vs different audio, family of 4, BH) —",
    "confounded with cycle position:\n")
print(as.data.frame(isc %>% select(unit, audio_vs_diff, p_audio_vs_diff, q_audio_vs_diff)),
      digits = 3, row.names = FALSE)
cat("  survivors at q<0.05:", sum(isc$q_audio_vs_diff < 0.05), "\n")

cat("\nCONTROLLED TEST (same audio vs SAME POSITION, different audio; family of 4, BH):\n")
print(as.data.frame(isc %>% select(unit, audio_vs_position, p_audio_vs_position,
                                   q_audio_vs_position)), digits = 3, row.names = FALSE)
cat("  survivors at q<0.05:", sum(isc$q_audio_vs_position < 0.05), "\n")

cat("\nHOW MUCH IS POSITION ALONE (same position, different audio vs neither):\n")
print(as.data.frame(isc %>% select(unit, position_vs_diff, p_position_vs_diff)),
      digits = 3, row.names = FALSE)
cat("\nSubjective-side counterpart: clip-mean positions reliable at",
    "ICC(2,k) 0.96 (comfort) / 0.92 (richness);\nclip explains 61% of comfort variance.\n")
write_outcome(isc, file.path(OUT, "neural_clip_reliability.csv"))

# ------------------------------------------------------------ 2. adjective screen
rule("2. ADJECTIVE SCREEN — 8 adjectives x 3 ROIs = 24, BH, response window")

adj_cols <- paste0("adj_", ADJECTIVES$zh)
ratings <- load_ratings() %>% select(subject_id, marker, all_of(adj_cols))
names(ratings)[-(1:2)] <- ADJECTIVES$en          # ASCII names for the formulae

# within-subject decomposition of each adjective, same convention as the axes
for (a in ADJECTIVES$en) {
  ratings <- ratings %>% group_by(subject_id) %>%
    mutate(!!paste0(a, "_wi") := .data[[a]] - mean(.data[[a]], na.rm = TRUE)) %>% ungroup()
}

roi_resp <- load_condition_roi_analysis("HbO") %>%
  select(subject_id, marker, roi, y = mean_response_m5_30) %>%
  inner_join(ratings, by = c("subject_id", "marker"))

adj_tab <- expand_grid(adjective = ADJECTIVES$en, roi = ROIS) %>%
  pmap_dfr(function(adjective, roi) {
    d <- roi_resp %>% filter(.data$roi == !!roi)
    m <- p24_lmm(d, paste0(adjective, "_wi"), "y")
    ft <- fixef_table(m, keep = paste0(adjective, "_wi"))
    tibble(adjective = adjective,
           zh = ADJECTIVES$zh[ADJECTIVES$en == adjective],
           axis_role = ADJECTIVES$axis_role[ADJECTIVES$en == adjective],
           roi = roi, n = nrow(d), estimate = ft$estimate, se = ft$se, p = ft$p)
  }) %>% mutate(q_BH = bh(p))

cat("\nAll 24, sorted by p:\n")
print(as.data.frame(adj_tab %>% arrange(p) %>%
  select(adjective, zh, roi, estimate, se, p, q_BH)), digits = 3, row.names = FALSE)
cat("\nSurvivors at q<0.05:", sum(adj_tab$q_BH < 0.05), "\n")

cat("\nFor comparison, the composite axes on the same window:\n")
axis_cmp <- map_dfr(ROIS, function(rr) {
  d <- load_condition_roi_analysis("HbO") %>% filter(roi == rr)
  m <- p24_lmm(d, "comfort_wi + richness_wi", "mean_response_m5_30")
  fixef_table(m, keep = c("comfort_wi", "richness_wi")) %>% mutate(roi = rr)
})
print(as.data.frame(axis_cmp %>% select(roi, term, estimate, se, p)),
      digits = 3, row.names = FALSE)
write_outcome(adj_tab, file.path(OUT, "adjective_screen.csv"))
write_outcome(axis_cmp, file.path(OUT, "axis_comparison.csv"))

rule("DONE")
close_log(con)
