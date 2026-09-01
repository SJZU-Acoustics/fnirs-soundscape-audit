# =============================================================================
# make_tables.R — regenerate the NUMBERS behind every manuscript table.
#
# Writes one CSV per manuscript table to output/tables/:
#   table1_parameters.csv          (main Table 1 — measurement parameters)
#   table2_multiplicity.csv        (main Table 2 — the multiplicity ledger)
#   table3_nondetection.csv        (main Table 3 — seven reasons)
#   tableS1_estimator_sweep.csv    (SI Table S1 — per-window family sweep)
#   tableS2_channel10_adjudication.csv (SI Table S2 — channel-10 adjudication)
#   tableS3_refilter_gates.csv     (SI Table S3 — gate-C frozen vs re-filter)
#
# Every value is pulled programmatically from the regenerated result CSVs in
# output/ (or, for two raw-QA-derived quantities, from the shipped
# intermediates via load_intermediate()/load_channels()). Where the quoted
# value is a Benjamini-Hochberg family minimum the stored q column of the
# module output is used (never re-adjusted). Row order mirrors the .tex
# tables. A missing source CSV (partial output tree before a full run_all.R
# pass) yields NA values plus a note, never a crash.
#
# Run from the repository root (run_all.R sources this after the modules).
# =============================================================================

source("code/load_data.R")

TAB_DIR <- file.path("output", "tables")
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------- safe reading
.missing <- character(0)

rd <- function(path, ...) {
  if (!file.exists(path)) {
    .missing <<- unique(c(.missing, path))
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE, ...)
}

an <- function(...) file.path("output", ...)
dl <- function(...) file.path("output", "data_lock", ...)

MISSING_NOTE <- "source not yet regenerated (partial output tree)"

# Minimum of a column, optionally over filtered rows; NA-safe.
min_of <- function(d, col) if (is.null(d) || !nrow(d)) NA_real_ else min(d[[col]])

# NULL-safe filter (a missing source stays NULL instead of crashing).
flt <- function(d, ...) if (is.null(d)) NULL else dplyr::filter(d, ...)

# Minimum over several candidates, NA if all are missing.
min_all <- function(...) {
  v <- unlist(list(...)); v <- v[!is.na(v)]
  if (!length(v)) NA_real_ else min(v)
}

# ------------------------------------------------------------ Table 1
# Measurement parameters (table1_parameters.tex), plus the quantities behind
# the "Outcome" column. Long format: one row per quoted quantity.

t1_rows <- list()
t1 <- function(row_label, quantity, value, source_output, note = "") {
  t1_rows[[length(t1_rows) + 1L]] <<-
    tibble(row_label, quantity, value = as.numeric(value), source_output, note)
}

# Row 1 — filter gain at the 0.0111 Hz fundamental and the 60/30/20-s
# components (mod_fig1b_filter_gain.R; DISPLAY_LOCK deviation 1: only the
# 0.688 gate is binding, components read "approximately 1.0").
g <- rd(dl("fig1b_gain_anchor_check.csv"))
if (!is.null(g)) {
  pick <- function(pat) g$gain[grepl(pat, g$component)][1]
  t1("filter amplitude response", "gain at 0.0111 Hz fundamental",
     pick("90-s cycle"), "output/data_lock/fig1b_gain_anchor_check.csv")
  t1("filter amplitude response", "gain at 60-s component",
     pick("60-s"), "output/data_lock/fig1b_gain_anchor_check.csv",
     "least-squares quadrature readout; manuscript says 'essentially unity' (DISPLAY_LOCK deviation 1)")
  t1("filter amplitude response", "gain at 30-s component",
     pick("30-s"), "output/data_lock/fig1b_gain_anchor_check.csv")
  t1("filter amplitude response", "gain at 20-s component",
     pick("20-s"), "output/data_lock/fig1b_gain_anchor_check.csv")
} else {
  for (q in c("gain at 0.0111 Hz fundamental", "gain at 60-s component",
              "gain at 30-s component", "gain at 20-s component"))
    t1("filter amplitude response", q, NA, "output/data_lock/fig1b_gain_anchor_check.csv", MISSING_NOTE)
}

# Rows 2-3 — dark-channel QA from the shipped dark trace (raw SNIRF recordings
# are not deposited; the trace is the shipped intermediate behind Fig 2a).
dk <- rd(dl("fig2a_dark_trace.csv"))
if (!is.null(dk)) {
  snr760 <- mean(dk$intensity_760nm) / sd(dk$intensity_760nm)
  snr850 <- mean(dk$intensity_850nm) / sd(dk$intensity_850nm)
  t1("dark channel under the delivered screen", "dark-floor mean, 760 nm",
     mean(dk$intensity_760nm), "output/data_lock/fig2a_dark_trace.csv")
  t1("dark channel under the delivered screen", "dark-floor mean, 850 nm",
     mean(dk$intensity_850nm), "output/data_lock/fig2a_dark_trace.csv")
  t1("dark channel under the delivered screen", "SNR (mean/SD), 760 nm", snr760,
     "output/data_lock/fig2a_dark_trace.csv",
     "order-of-magnitude anchor only; the quoted ~6,700 is the raw-recording QA statistic")
  t1("dark channel under the delivered screen", "SNR (mean/SD), 850 nm", snr850,
     "output/data_lock/fig2a_dark_trace.csv",
     "order-of-magnitude anchor only; the quoted ~6,700 is the raw-recording QA statistic")
  t1("dRange window", "dark floor / lower dRange bound (factor)",
     mean(c(dk$intensity_760nm, dk$intensity_850nm)) / 0.001,
     "output/data_lock/fig2a_dark_trace.csv")
} else {
  for (q in c("dark-floor mean, 760 nm", "dark-floor mean, 850 nm",
              "SNR (mean/SD), 760 nm", "SNR (mean/SD), 850 nm",
              "dark floor / lower dRange bound (factor)"))
    t1(if (grepl("dRange", q)) "dRange window" else "dark channel under the delivered screen",
       q, NA, "output/data_lock/fig2a_dark_trace.csv", MISSING_NOTE)
}

# Row 4 — union mask: 289 excluded cells; declared families null (0/18).
ms <- rd(an("mask_summary.csv"))
mval <- function(metric) if (is.null(ms)) NA_real_ else ms$value[ms$metric == metric][1]
t1("union channel mask", "excluded subject x channel cells",
   mval("excluded_cells_dark_fraction_criterion"), "output/mask_summary.csv",
   if (is.null(ms)) MISSING_NOTE else "")
f24 <- rd(an("analysis_24_dark_fraction_mask", "families_dark_fraction_mask.csv"))
t1("union channel mask", "declared families surviving q < 0.05 under the mask (of 18)",
   if (is.null(f24)) NA_real_ else sum(f24$q_BH < 0.05),
   "output/analysis_24_dark_fraction_mask/families_dark_fraction_mask.csv",
   if (is.null(f24)) MISSING_NOTE else "")
t1("union channel mask", "declared-family tests under the mask (total)",
   if (is.null(f24)) NA_real_ else nrow(f24),
   "output/analysis_24_dark_fraction_mask/families_dark_fraction_mask.csv",
   if (is.null(f24)) MISSING_NOTE else "")

# Row 5 — re-filter gate: waveform agreement of the valid (no motion
# correction) arm, median r over HbO channels (A19 gate B, from the shipped
# gate_b_reproduction intermediate).
gb <- tryCatch(load_intermediate("gate_b_reproduction.csv"),
               error = function(e) NULL)
t1("wavelet motion-correction gate", "median waveform r, no-MC arm (HbO)",
   if (is.null(gb)) NA_real_
   else median(gb$r[gb$motion_correction == "none" & gb$chromophore == "HbO"]),
   "data/intermediates/gate_b_reproduction.csv",
   if (is.null(gb)) "shipped intermediate missing" else "")

# Row 6 — shortest active source-detector separation.
ch <- tryCatch(load_channels(), error = function(e) NULL)
t1("montage", "shortest active pair (mm)",
   if (is.null(ch)) NA_real_ else min(ch$sd_distance_mm[ch$active_pair == 1], na.rm = TRUE),
   "workbook channels sheet (load_channels())",
   if (is.null(ch)) "workbook unavailable" else "")

# Row 7 — mirror-pair survival (A20's CV-mask-era counts, per DISPLAY_LOCK
# deviation 3).
mp <- rd(an("analysis_20_detector_dropout", "mirror_pair_completeness.csv"))
if (is.null(mp)) mp <- rd(dl("fig2c_mirror_pairs.csv"))
if (!is.null(mp)) {
  if ("subjects_with_all_four_live" %in% names(mp)) {
    t1("mirror-symmetric montage", "subjects with all four T-pair channels live",
       mp$subjects_with_all_four_live[mp$mirror_pair == "T"][1],
       "output/analysis_20_detector_dropout/mirror_pair_completeness.csv")
    t1("mirror-symmetric montage", "subjects with all four TP-pair channels live",
       mp$subjects_with_all_four_live[mp$mirror_pair == "TP"][1],
       "output/analysis_20_detector_dropout/mirror_pair_completeness.csv")
  } else {
    comp <- tapply(mp$subjects_pair_complete, mp$pair, max)
    t1("mirror-symmetric montage", "subjects with all four T-pair channels live",
       comp[["T"]], "output/data_lock/fig2c_mirror_pairs.csv")
    t1("mirror-symmetric montage", "subjects with all four TP-pair channels live",
       comp[["TP"]], "output/data_lock/fig2c_mirror_pairs.csv")
  }
} else {
  t1("mirror-symmetric montage", "subjects with all four T-pair channels live",
     NA, "output/analysis_20_detector_dropout/mirror_pair_completeness.csv", MISSING_NOTE)
  t1("mirror-symmetric montage", "subjects with all four TP-pair channels live",
     NA, "output/analysis_20_detector_dropout/mirror_pair_completeness.csv", MISSING_NOTE)
}

# Row 8 — HbT arithmetic identity.
hb <- rd(an("analysis_09_chromophore_triangulation", "hbt_identity_check.csv"))
t1("HbT identity", "max |HbT - (HbO + HbR)| (uM)",
   if (is.null(hb)) NA_real_ else max(hb$max_abs_row_gap),
   "output/analysis_09_chromophore_triangulation/hbt_identity_check.csv",
   if (is.null(hb)) MISSING_NOTE else "")

table1 <- bind_rows(t1_rows)
write_outcome(table1, file.path(TAB_DIR, "table1_parameters.csv"))

# ------------------------------------------------------------ Table 2
# The multiplicity ledger (table2_multiplicity.tex). One row per printed row,
# in print order; min_q pulled from the stored q column of the module output.

t2_rows <- list()
t2 <- function(row_label, tests, correction, min_q, survivors, source_output, note = "") {
  t2_rows[[length(t2_rows) + 1L]] <<-
    tibble(row_label, tests, correction, min_q = as.numeric(min_q),
           survivors, source_output, note)
}

# Confirmatory family (A04): family_by_window.csv, min q per analysis window.
fw <- rd(an("analysis_04_appraisal_to_haemodynamics", "family_by_window.csv"))
SRC04 <- "output/analysis_04_appraisal_to_haemodynamics/family_by_window.csv"
nt04 <- if (is.null(fw)) MISSING_NOTE else ""
t2("ROI within-subject slopes, difference score (A04)", "6", "BH",
   min_of(fw %>% filter(analysis == "difference"), "q_BH"), "0", SRC04, nt04)
t2("ROI within-subject slopes, response window (A04)", "6", "BH",
   min_of(fw %>% filter(analysis == "response"), "q_BH"), "0", SRC04, nt04)
t2("ROI within-subject slopes, pre-stimulus baseline (A04)", "6", "BH",
   min_of(fw %>% filter(analysis == "baseline"), "q_BH"), "0", SRC04, nt04)

# A05 alternative response geometry: omnibus q per window, minimum quoted in
# Table 2 (the baseline-window survivor is withdrawn — footnote a).
a5 <- rd(an("analysis_05_alternative_response_geometry", "alternative_geometry_omnibus.csv"))
t2("Alternative response geometry, 3 windows (A05)", "6 x 3", "BH per window",
   min_of(a5, "q_BH"), "0 (footnote a)",
   "output/analysis_05_alternative_response_geometry/alternative_geometry_omnibus.csv",
   if (is.null(a5)) MISSING_NOTE else
     "Table 2 quotes the cross-window minimum; the apparent baseline survivor was carried by detector 8 and is withdrawn")

# A06 channel-wise scan: family_summary.csv min over windows.
a6 <- rd(an("analysis_06_channelwise_spatial_scan", "family_summary.csv"))
t2("Channel-wise spatial scan, 3 windows (A06)", "42 x 3", "BH per window",
   min_of(a6, "min_q"), "0",
   "output/analysis_06_channelwise_spatial_scan/family_summary.csv",
   if (is.null(a6)) MISSING_NOTE else "")

# A08 CR2 robust audit: response-window CR2 min q.
a8 <- rd(an("analysis_08_dependence_structure_audit", "family_summary_CR2.csv"))
t2("Cluster-robust (CR2) channel audit (A08)", "42", "BH",
   min_of(a8 %>% filter(window == "response"), "CR2_min_q"), "0",
   "output/analysis_08_dependence_structure_audit/family_summary_CR2.csv",
   if (is.null(a8)) MISSING_NOTE else "")

# A09 HbR/HbT triangulation: minimum over the two window families.
a9d <- rd(an("analysis_09_chromophore_triangulation", "family_difference_hbr_hbt.csv"))
a9r <- rd(an("analysis_09_chromophore_triangulation", "family_response_hbr_hbt.csv"))
t2("HbR and HbT slope triangulation, 2 windows (A09)", "12 x 2", "BH per window",
   min(c(min_of(a9d, "q_BH"), min_of(a9r, "q_BH")), na.rm = TRUE), "0",
   "output/analysis_09_chromophore_triangulation/family_{difference,response}_hbr_hbt.csv",
   if (is.null(a9d) && is.null(a9r)) MISSING_NOTE else "")
if (is.null(a9d) && is.null(a9r)) t2_rows[[length(t2_rows)]]$min_q <- NA_real_

# A10 trait reactivity — DOCUMENTED DEVIATION: the deposit anonymises exact
# age into five-year bands, so this family regenerates to min q = 0.576 where
# the manuscript quotes 0.273. Emit the regenerated value; do not force 0.273.
a10r <- rd(an("analysis_10_person_side", "between_person_reactivity.csv"))
t2("Trait reactivity (A10)", "12", "BH",
   min_of(a10r, "q_BH"), "0",
   "output/analysis_10_person_side/between_person_reactivity.csv",
   if (is.null(a10r)) MISSING_NOTE else
     "deposit anonymisation: exact age not in the workbook (age bands) — family regenerates to 0.576 where the manuscript quotes 0.273; see README")

a10m <- rd(an("analysis_10_person_side", "trait_moderation_response.csv"))
t2("Trait x appraisal moderation (A10)", "24", "BH",
   min_of(a10m, "q_BH"), "0",
   "output/analysis_10_person_side/trait_moderation_response.csv",
   if (is.null(a10m)) MISSING_NOTE else "")

# A11 cycle time course.
a11b <- rd(an("analysis_11_cycle_time_course", "bin_slopes_comfort.csv"))
t2("Cycle time course, 5-s bins (A11)", "28", "BH",
   min_of(a11b, "q_BH"), "0",
   "output/analysis_11_cycle_time_course/bin_slopes_comfort.csv",
   if (is.null(a11b)) MISSING_NOTE else "")
a11g <- rd(an("analysis_11_cycle_time_course", "growth_contrast.csv"))
t2("Time-course growth contrasts (A11)", "2", "BH",
   min_of(a11g, "q_BH"), "0",
   "output/analysis_11_cycle_time_course/growth_contrast.csv",
   if (is.null(a11g)) MISSING_NOTE else "")

# A12 repetition structure.
a12t <- rd(an("analysis_12_repetition_structure", "run_trend_by_window_roi.csv"))
t2("Encounter (repetition) trends (A12)", "6", "BH",
   min_of(a12t, "q_BH"), "2 (footnote a)",
   "output/analysis_12_repetition_structure/run_trend_by_window_roi.csv",
   if (is.null(a12t)) MISSING_NOTE else
     "third survivor (right-temporal trend) withdrawn as carried by detector 8")
a12e <- rd(an("analysis_12_repetition_structure", "per_encounter_family.csv"))
t2("Per-encounter appraisal slopes (A12)", "18", "BH",
   min_of(a12e, "q_BH"), "0",
   "output/analysis_12_repetition_structure/per_encounter_family.csv",
   if (is.null(a12e)) MISSING_NOTE else "")
a12f <- rd(an("analysis_12_repetition_structure", "first_vs_later_contrast.csv"))
t2("First-versus-later coupling contrast (A12)", "3", "BH",
   min_of(a12f, "q_BH"), "0",
   "output/analysis_12_repetition_structure/first_vs_later_contrast.csv",
   if (is.null(a12f)) MISSING_NOTE else "")

# A13 stimulus acoustic features.
a13 <- rd(an("analysis_13_stimulus_acoustic_drivers", "family_acoustic_slopes_20_30.csv"))
t2("Stimulus acoustic features (A13)", "8", "BH",
   min_of(a13, "q_BH"), "0",
   "output/analysis_13_stimulus_acoustic_drivers/family_acoustic_slopes_20_30.csv",
   if (is.null(a13)) MISSING_NOTE else "")

# A14 (Table 2): the masked variant from the A21 correction sweep.
a14m <- rd(an("analysis_21_correction_sweep", "a14_family_masked.csv"))
t2("Sequence-detrended baseline coupling (A14, masked)", "6", "BH",
   min_of(a14m, "q_BH"), "0",
   "output/analysis_21_correction_sweep/a14_family_masked.csv",
   if (is.null(a14m)) MISSING_NOTE else "")

# A15 leave-one-out global adjustment: Table 2 quotes the smallest q across
# the three windows (per-window minima in Table S1).
a15r <- rd(an("analysis_15_spatial_specificity", "roi_global_adjusted.csv"))
t2("Global-adjusted ROI slopes, 3 windows (A15)", "6 x 3", "BH per window",
   min_of(a15r, "q_adjusted"), "0",
   "output/analysis_15_spatial_specificity/roi_global_adjusted.csv",
   if (is.null(a15r)) MISSING_NOTE else
     "minimum across the three windows (baseline, right temporal richness)")
a15c <- rd(an("analysis_15_spatial_specificity", "channel_adjusted_cr2.csv"))
t2("Global-adjusted channel scan (A15)", "42", "BH",
   min_of(a15c, "q_cr2"), "0",
   "output/analysis_15_spatial_specificity/channel_adjusted_cr2.csv",
   if (is.null(a15c)) MISSING_NOTE else "")

# A16 ridge: the declared primary family is the four comfort tests.
a16 <- rd(an("analysis_16_multivariate_pattern", "primary_comfort_pfc13.csv"))
t2("Multivariate ridge, leave-one-subject-out (A16)", "4", "BH",
   min_of(a16, "q_BH"), "0",
   "output/analysis_16_multivariate_pattern/primary_comfort_pfc13.csv",
   if (is.null(a16)) MISSING_NOTE else "")

# A17 canonical-HRF GLM: minimum over the ROI and channel families.
a17r <- rd(an("analysis_17_hrf_glm", "roi_family_glm.csv"))
a17c <- rd(an("analysis_17_hrf_glm", "channel_family_glm.csv"))
a17_min <- suppressWarnings(min(c(min_of(a17r, "q_BH"), min_of(a17c, "q_BH")), na.rm = TRUE))
t2("Canonical-HRF GLM, ROI and channel families (A17)", "6 + 42", "BH per family",
   if (is.infinite(a17_min)) NA_real_ else a17_min, "0",
   "output/analysis_17_hrf_glm/{roi,channel}_family_glm.csv",
   if (is.null(a17r) && is.null(a17c)) MISSING_NOTE else "")

# A18 clip reliability: permutation family, minimum audio-vs-position p.
a18 <- rd(an("analysis_18_clip_reliability_and_adjectives", "neural_clip_reliability.csv"))
t2("Clip time-course reliability, audio vs position control (A18)", "4", "permutation",
   min_of(a18, "p_audio_vs_position"), "0",
   "output/analysis_18_clip_reliability_and_adjectives/neural_clip_reliability.csv",
   if (is.null(a18)) MISSING_NOTE else "minimum permutation p (Table 2 prints p >= 0.29)")
a18a <- rd(an("analysis_18_clip_reliability_and_adjectives", "adjective_screen.csv"))
t2("Adjective-level slopes (A18)", "24", "BH",
   min_of(a18a, "q_BH"), "0",
   "output/analysis_18_clip_reliability_and_adjectives/adjective_screen.csv",
   if (is.null(a18a)) MISSING_NOTE else "")

# A24 masked ROI families: Table 2 quotes the smallest q across the three
# windows (per-window minima in Table S1).
t2("Masked ROI families, 3 windows (A24)", "6 x 3", "BH per window",
   min_of(f24, "q_BH"), "0",
   "output/analysis_24_dark_fraction_mask/families_dark_fraction_mask.csv",
   if (is.null(f24)) MISSING_NOTE else
     "minimum across the three windows (baseline window)")

# A23 responder heterogeneity.
a23v <- rd(an("analysis_23_responder_heterogeneity", "A2_variance_falsification.csv"))
t2("Slope variance components, 3 windows (A23)", "6 x 3", "BH + permutation gate",
   min_of(a23v, "q_perm_F1"), "0",
   "output/analysis_23_responder_heterogeneity/A2_variance_falsification.csv",
   if (is.null(a23v)) MISSING_NOTE else
     "permutation-gated q (five nominal bootstrap cells overturned; footnote c)")
a23s <- rd(an("analysis_23_responder_heterogeneity", "CD_selection_summary.csv"))
t2("Out-of-sample responder selection (A23)", "12", "max-|t| permutation",
   min_of(a23s %>% filter(arm == "C_out_of_sample"), "q_bh"), "0",
   "output/analysis_23_responder_heterogeneity/CD_selection_summary.csv",
   if (is.null(a23s)) MISSING_NOTE else "")

table2 <- bind_rows(t2_rows)
write_outcome(table2, file.path(TAB_DIR, "table2_multiplicity.csv"))

# ------------------------------------------------------------ Table 3
# Seven reasons non-detection is not absence (table3_nondetection.tex).
# Long format; rows 1 and 4 carry no quantity in the manuscript.

t3_rows <- list()
t3 <- function(row_label, quantity, value, source_output, note = "") {
  t3_rows[[length(t3_rows) + 1L]] <<-
    tibble(row_label, quantity, value = as.numeric(value), source_output, note)
}

# Reason 2 — waveform agreement, median r (same gate-B quantity as Table 1).
t3("2. rebuilt-vs-delivered waveform agreement", "median r, no-MC arm (HbO)",
   if (is.null(gb)) NA_real_
   else median(gb$r[gb$motion_correction == "none" & gb$chromophore == "HbO"]),
   "data/intermediates/gate_b_reproduction.csv",
   if (is.null(gb)) "shipped intermediate missing" else "")

# Reason 3 — cross-arm estimate agreement (A21 gate-C decomposition).
gcd <- rd(an("analysis_21_correction_sweep", "a19_gate_c_decomposed.csv"))
SRC21G <- "output/analysis_21_correction_sweep/a19_gate_c_decomposed.csv"
t3("3. rebuilt-vs-delivered estimate agreement", "r over 12 ROI x axis cells, as published",
   if (is.null(gcd)) NA_real_ else gcd$r_both_arms_12[gcd$roi_definition == "as published"][1],
   SRC21G, if (is.null(gcd)) MISSING_NOTE else "")
t3("3. rebuilt-vs-delivered estimate agreement", "r, corrected region (valid arm)",
   if (is.null(gcd)) NA_real_ else gcd$r_valid_arm_6[gcd$roi_definition == "corrected"][1],
   SRC21G, if (is.null(gcd)) MISSING_NOTE else "")

# Reason 5 — no short-separation channel.
t3("5. no short-separation channel", "shortest active pair (mm)",
   if (is.null(ch)) NA_real_ else min(ch$sd_distance_mm[ch$active_pair == 1], na.rm = TRUE),
   "workbook channels sheet (load_channels())",
   if (is.null(ch)) "workbook unavailable" else "")

# Reason 6 — right-temporal region empty in 11 of 69 participants.
t3("6. temporal channel loss", "participants with no usable right-temporal channel",
   mval("right_temporal_subjects_no_channel_usable"), "output/mask_summary.csv",
   if (is.null(ms)) MISSING_NOTE else "")
t3("6. temporal channel loss", "total participants", 69, "design constant", "")

# Reason 7 — underpowering: implied peak against lambda (full chain / ideal).
lc <- rd(an("analysis_22_chain", "lambda_chain.csv"))
lam <- if (is.null(lc)) NULL else
  lc %>% filter(se_source == "corrected ROI (primary)", roi == "PFC_Frontal")
SRC22 <- "output/analysis_22_chain/lambda_chain.csv"
t3("7. underpowering", "lambda, full chain (uM per axis unit)",
   if (is.null(lam)) NA_real_ else lam$C3_plus_dependence[1], SRC22,
   if (is.null(lam)) MISSING_NOTE else "")
t3("7. underpowering", "lambda, ideal chain (uM per axis unit)",
   if (is.null(lam)) NA_real_ else lam$C0_ideal[1], SRC22,
   if (is.null(lam)) MISSING_NOTE else "")
# The implied peak is printed (not CSV-stored) by mod22_chain.R; parse it from
# the regenerated chain report. Fall back to the Fig 6 data-lock reference.
imp <- NA_real_
cr <- an("analysis_22_chain", "chain_report.txt")
if (file.exists(cr)) {
  txt <- readLines(cr, warn = FALSE)
  hit <- grep("underlying peak\\s+-?[0-9.]+", txt, value = TRUE)
  if (length(hit)) imp <- as.numeric(sub(".*underlying peak\\s+(-?[0-9.]+).*", "\\1", hit[1]))
}
imp_src <- "output/analysis_22_chain/chain_report.txt"
if (is.na(imp)) {
  f6 <- rd(dl("fig6a_lambda_chain.csv"))
  if (!is.null(f6)) { imp <- -abs(f6$reference_line[1]); imp_src <- "output/data_lock/fig6a_lambda_chain.csv" }
}
t3("7. underpowering", "implied peak response (uM per axis unit)", imp, imp_src,
   if (is.na(imp)) MISSING_NOTE else "")

table3 <- bind_rows(t3_rows)
write_outcome(table3, file.path(TAB_DIR, "table3_nondetection.csv"))

# ------------------------------------------------------------ Table S1
# Estimator sweep (tableS1_estimator_sweep.tex). Per-window minima as
# difference / response / baseline columns where the family runs per window.

s1_rows <- list()
s1 <- function(row_label, tests, correction,
               min_q = NA_real_, min_q_difference = NA_real_,
               min_q_response = NA_real_, min_q_baseline = NA_real_,
               source_output, note = "") {
  s1_rows[[length(s1_rows) + 1L]] <<-
    tibble(row_label, tests, correction,
           min_q = as.numeric(min_q),
           min_q_difference = as.numeric(min_q_difference),
           min_q_response = as.numeric(min_q_response),
           min_q_baseline = as.numeric(min_q_baseline),
           source_output, note)
}

s1("A04 difference score", "6", "BH",
   min_q = min_of(fw %>% filter(analysis == "difference"), "q_BH"),
   source_output = SRC04, note = nt04)
s1("A04 response window", "6", "BH",
   min_q = min_of(fw %>% filter(analysis == "response"), "q_BH"),
   source_output = SRC04, note = nt04)
s1("A04 pre-stimulus baseline (negative control)", "6", "BH",
   min_q = min_of(fw %>% filter(analysis == "baseline"), "q_BH"),
   source_output = SRC04, note = nt04)

wmin <- function(d, col) {
  if (is.null(d)) return(c(NA_real_, NA_real_, NA_real_))
  by <- if ("window" %in% names(d)) d$window else d$analysis
  v <- tapply(d[[col]], by, min)
  unname(v[c("difference", "response", "baseline")])
}
a5w <- wmin(a5, "q_BH")
s1("A05 quadratic axes; nominal quadrant, 3 windows", "6 x 3", "BH / window",
   min_q_difference = a5w[1], min_q_response = a5w[2], min_q_baseline = a5w[3],
   source_output = "output/analysis_05_alternative_response_geometry/alternative_geometry_omnibus.csv",
   note = if (is.null(a5)) MISSING_NOTE else "baseline-window survivor carried by detector 8, withdrawn")

a9dw <- min_of(a9d, "q_BH"); a9rw <- min_of(a9r, "q_BH")
s1("A09 HbR and HbT ROI slopes, 2 windows", "12 x 2", "BH / window",
   min_q_difference = a9dw, min_q_response = a9rw,
   source_output = "output/analysis_09_chromophore_triangulation/family_{difference,response}_hbr_hbt.csv",
   note = if (is.null(a9d) && is.null(a9r)) MISSING_NOTE else "")

a6w <- if (is.null(a6)) c(NA_real_, NA_real_, NA_real_) else {
  v <- tapply(a6$min_q, a6$window, min); unname(v[c("difference", "response", "baseline")])
}
s1("A06 channel-wise scan, 21 channels x 2 axes, 3 windows", "42 x 3", "BH / window",
   min_q_difference = a6w[1], min_q_response = a6w[2], min_q_baseline = a6w[3],
   source_output = "output/analysis_06_channelwise_spatial_scan/family_summary.csv",
   note = if (is.null(a6)) MISSING_NOTE else "")
s1("A08 CR2 cluster-robust re-run, response", "42", "BH",
   min_q = min_of(a8 %>% filter(window == "response"), "CR2_min_q"),
   source_output = "output/analysis_08_dependence_structure_audit/family_summary_CR2.csv",
   note = if (is.null(a8)) MISSING_NOTE else "")

s1("A10 trait reactivity, 4 traits x 3 ROIs", "12", "BH",
   min_q = min_of(a10r, "q_BH"),
   source_output = "output/analysis_10_person_side/between_person_reactivity.csv",
   note = if (is.null(a10r)) MISSING_NOTE else
     "deposit anonymisation: exact age not in the workbook — regenerates to 0.576 where S1 quotes 0.273; see README")
s1("A10 trait x appraisal moderation", "24", "BH",
   min_q = min_of(a10m, "q_BH"),
   source_output = "output/analysis_10_person_side/trait_moderation_response.csv",
   note = if (is.null(a10m)) MISSING_NOTE else "")

s1("A11 comfort slopes, 2 targets x 14 bins", "28", "BH",
   min_q = min_of(a11b, "q_BH"),
   source_output = "output/analysis_11_cycle_time_course/bin_slopes_comfort.csv",
   note = if (is.null(a11b)) MISSING_NOTE else "")
s1("A11 growth contrasts", "2", "BH",
   min_q = min_of(a11g, "q_BH"),
   source_output = "output/analysis_11_cycle_time_course/growth_contrast.csv",
   note = if (is.null(a11g)) MISSING_NOTE else "")

s1("A12 encounter trends, 2 windows x 3 ROIs", "6", "BH",
   min_q = min_of(a12t, "q_BH"),
   source_output = "output/analysis_12_repetition_structure/run_trend_by_window_roi.csv",
   note = if (is.null(a12t)) MISSING_NOTE else "")
s1("A12 per-encounter appraisal slopes", "18", "BH",
   min_q = min_of(a12e, "q_BH"),
   source_output = "output/analysis_12_repetition_structure/per_encounter_family.csv",
   note = if (is.null(a12e)) MISSING_NOTE else "")
s1("A12 first-versus-later coupling contrast", "3", "BH",
   min_q = min_of(a12f, "q_BH"),
   source_output = "output/analysis_12_repetition_structure/first_vs_later_contrast.csv",
   note = if (is.null(a12f)) MISSING_NOTE else "")

s1("A13 feature slopes, 4 features x 2 targets", "8", "BH",
   min_q = min_of(a13, "q_BH"),
   source_output = "output/analysis_13_stimulus_acoustic_drivers/family_acoustic_slopes_20_30.csv",
   note = if (is.null(a13)) MISSING_NOTE else "")

# S1 A14: the unmasked detrended family from analysis_14 (Table 2 uses the
# masked variant from analysis_21).
a14 <- rd(an("analysis_14_baseline_coupling_anatomy", "family_detrended_comfort_slopes.csv"))
s1("A14 detrended comfort slopes, 2 x 3", "6", "BH",
   min_q = min_of(a14, "q_BH"),
   source_output = "output/analysis_14_baseline_coupling_anatomy/family_detrended_comfort_slopes.csv",
   note = if (is.null(a14)) MISSING_NOTE else "")

a15w <- if (is.null(a15r)) rep(NA_real_, 3) else
  c(min_of(a15r %>% filter(window == "difference"), "q_adjusted"),
    min_of(a15r %>% filter(window == "response"), "q_adjusted"),
    min_of(a15r %>% filter(window == "baseline"), "q_adjusted"))
s1("A15 adjusted ROI slopes, 3 windows", "6 x 3", "BH / window",
   min_q_difference = a15w[1], min_q_response = a15w[2], min_q_baseline = a15w[3],
   source_output = "output/analysis_15_spatial_specificity/roi_global_adjusted.csv",
   note = if (is.null(a15r)) MISSING_NOTE else "")
s1("A15 adjusted channel scan, response", "42", "BH (CR2)",
   min_q = min_of(a15c, "q_cr2"),
   source_output = "output/analysis_15_spatial_specificity/channel_adjusted_cr2.csv",
   note = if (is.null(a15c)) MISSING_NOTE else "")

s1("A16 ridge, nested LOSO, 2 x 2", "4", "BH",
   min_q = min_of(a16, "q_BH"),
   source_output = "output/analysis_16_multivariate_pattern/primary_comfort_pfc13.csv",
   note = if (is.null(a16)) MISSING_NOTE else "")

s1("A17 ROI parametric-modulator family", "6", "BH",
   min_q = min_of(a17r, "q_BH"),
   source_output = "output/analysis_17_hrf_glm/roi_family_glm.csv",
   note = if (is.null(a17r)) MISSING_NOTE else "")
s1("A17 channel parametric-modulator family", "42", "BH",
   min_q = min_of(a17c, "q_BH"),
   source_output = "output/analysis_17_hrf_glm/channel_family_glm.csv",
   note = if (is.null(a17c)) MISSING_NOTE else "")

s1("A18 clip reliability, audio vs position", "4", "permutation",
   min_q = min_of(a18, "p_audio_vs_position"),
   source_output = "output/analysis_18_clip_reliability_and_adjectives/neural_clip_reliability.csv",
   note = if (is.null(a18)) MISSING_NOTE else
     "minimum permutation p; S1 prints the range audio p = 0.29-0.91")
s1("A18 adjective slopes, 8 x 3 ROIs", "24", "BH",
   min_q = min_of(a18a, "q_BH"),
   source_output = "output/analysis_18_clip_reliability_and_adjectives/adjective_screen.csv",
   note = if (is.null(a18a)) MISSING_NOTE else "")

f24w <- wmin(f24, "q_BH")
s1("A24 ROI slopes under the mask, 3 windows", "6 x 3", "BH / window",
   min_q_difference = f24w[1], min_q_response = f24w[2], min_q_baseline = f24w[3],
   source_output = "output/analysis_24_dark_fraction_mask/families_dark_fraction_mask.csv",
   note = if (is.null(f24)) MISSING_NOTE else "")

s1("A23 slope variance components, 3 windows", "6 x 3", "BH + perm. gate",
   min_q = min_of(a23v, "q_perm_F1"),
   source_output = "output/analysis_23_responder_heterogeneity/A2_variance_falsification.csv",
   note = if (is.null(a23v)) MISSING_NOTE else "")
s1("A23 out-of-sample responder selection", "12", "max-|t| perm.",
   min_q = min_of(a23s %>% filter(arm == "C_out_of_sample"), "q_bh"),
   source_output = "output/analysis_23_responder_heterogeneity/CD_selection_summary.csv",
   note = if (is.null(a23s)) MISSING_NOTE else "")

tableS1 <- bind_rows(s1_rows)
write_outcome(tableS1, file.path(TAB_DIR, "tableS1_estimator_sweep.csv"))

# ------------------------------------------------------------ Table S2
# Channel-10 adjudication (tableS2_channel10_adjudication.tex). Long format:
# one row per quoted quantity.

s2_rows <- list()
s2 <- function(analysis, quantity, value, source_output, note = "") {
  s2_rows[[length(s2_rows) + 1L]] <<-
    tibble(analysis, quantity, value = as.numeric(value), source_output, note)
}

# A06 — the nominal-era survivor.
c10 <- rd(an("analysis_06_channelwise_spatial_scan", "channel_scan_response.csv"))
c10r <- if (is.null(c10)) NULL else
  c10 %>% filter(window == "response", channel_index == 10, term == "comfort_wi")
SRC06 <- "output/analysis_06_channelwise_spatial_scan/channel_scan_response.csv"
nt06 <- if (is.null(c10r)) MISSING_NOTE else ""
s2("A06 model-based 42-test response family", "estimate (uM)",
   if (is.null(c10r)) NA_real_ else c10r$estimate[1], SRC06, nt06)
s2("A06 model-based 42-test response family", "p",
   if (is.null(c10r)) NA_real_ else c10r$p[1], SRC06, nt06)
s2("A06 model-based 42-test response family", "q",
   if (is.null(c10r)) NA_real_ else c10r$q_BH[1], SRC06, nt06)

# A07 — adjudication.
b7 <- rd(an("analysis_07_focal_signal_adjudication", "bootstrap_ci.csv"))
b7r <- if (is.null(b7)) NULL else b7 %>% filter(term == "comfort_wi")
SRC07B <- "output/analysis_07_focal_signal_adjudication/bootstrap_ci.csv"
s2("A07 adjudication", "bootstrap CI lo", if (is.null(b7r)) NA_real_ else b7r$ci_lo[1],
   SRC07B, if (is.null(b7r)) MISSING_NOTE else "")
s2("A07 adjudication", "bootstrap CI hi", if (is.null(b7r)) NA_real_ else b7r$ci_hi[1],
   SRC07B, if (is.null(b7r)) MISSING_NOTE else "")
l7 <- rd(an("analysis_07_focal_signal_adjudication", "leave_one_subject_out_summary.csv"))
SRC07L <- "output/analysis_07_focal_signal_adjudication/leave_one_subject_out_summary.csv"
s2("A07 adjudication", "leave-one-subject-out deletions",
   if (is.null(l7)) NA_real_ else l7$n_deletions[1], SRC07L,
   if (is.null(l7)) MISSING_NOTE else "")
s2("A07 adjudication", "all LOSO estimates negative (1/0)",
   if (is.null(l7)) NA_real_ else as.numeric(l7$all_negative[1]), SRC07L,
   if (is.null(l7)) MISSING_NOTE else "")
s7 <- rd(an("analysis_07_focal_signal_adjudication", "sample_sensitivity.csv"))
s7r <- if (is.null(s7)) NULL else s7 %>% filter(sample == "qc_clean", term == "comfort_wi")
SRC07S <- "output/analysis_07_focal_signal_adjudication/sample_sensitivity.csv"
s2("A07 adjudication", "quality-clean estimate (uM)",
   if (is.null(s7r)) NA_real_ else s7r$estimate[1], SRC07S,
   if (is.null(s7r)) MISSING_NOTE else "")
s2("A07 adjudication", "quality-clean p",
   if (is.null(s7r)) NA_real_ else s7r$p[1], SRC07S,
   if (is.null(s7r)) MISSING_NOTE else "")
h7 <- rd(an("analysis_07_focal_signal_adjudication", "chromophore_slopes.csv"))
h7r <- if (is.null(h7)) NULL else h7 %>% filter(chromophore == "HbR", term == "comfort_wi")
s2("A07 adjudication", "HbR p", if (is.null(h7r)) NA_real_ else h7r$p[1],
   "output/analysis_07_focal_signal_adjudication/chromophore_slopes.csv",
   if (is.null(h7r)) MISSING_NOTE else "")

# A08 — CR2 shift.
s8 <- rd(an("analysis_08_dependence_structure_audit", "channel10_inference_shift.csv"))
SRC08 <- "output/analysis_08_dependence_structure_audit/channel10_inference_shift.csv"
nt08 <- if (is.null(s8)) MISSING_NOTE else ""
s2("A08 CR2 robust, unchanged family", "estimate (uM)",
   if (is.null(s8)) NA_real_ else s8$estimate[1], SRC08, nt08)
s2("A08 CR2 robust, unchanged family", "CR2 p",
   if (is.null(s8)) NA_real_ else s8$CR2_p[1], SRC08, nt08)
s2("A08 CR2 robust, unchanged family", "CR2 q",
   if (is.null(s8)) NA_real_ else s8$CR2_q[1], SRC08, nt08)
s2("A08 CR2 robust, unchanged family", "SE increase (%)",
   if (is.null(s8)) NA_real_ else s8$SE_increase_pct[1], SRC08, nt08)

# A15 — family leader under global adjustment.
a15c10 <- if (is.null(a15c)) NULL else
  a15c %>% filter(channel_index == 10, axis == "comfort_wi")
SRC15 <- "output/analysis_15_spatial_specificity/channel_adjusted_cr2.csv"
nt15 <- if (is.null(a15c10)) MISSING_NOTE else ""
s2("A15 leave-one-out global adjustment, CR2", "estimate (uM)",
   if (is.null(a15c10)) NA_real_ else a15c10$est[1], SRC15, nt15)
s2("A15 leave-one-out global adjustment, CR2", "p",
   if (is.null(a15c10)) NA_real_ else a15c10$p_cr2[1], SRC15, nt15)
s2("A15 leave-one-out global adjustment, CR2", "q",
   if (is.null(a15c10)) NA_real_ else a15c10$q_cr2[1], SRC15, nt15)

# A16 — cross-validated single-channel benchmark.
b16 <- rd(an("analysis_16_multivariate_pattern", "single_channel_benchmark.csv"))
b16r <- if (is.null(b16)) NULL else b16 %>% filter(window == "response", channel_index == 10)
SRC16 <- "output/analysis_16_multivariate_pattern/single_channel_benchmark.csv"
nt16 <- if (is.null(b16r)) MISSING_NOTE else ""
s2("A16 cross-validated single-channel benchmark", "r",
   if (is.null(b16r)) NA_real_ else b16r$r_obs[1], SRC16, nt16)
s2("A16 cross-validated single-channel benchmark", "p (perm.)",
   if (is.null(b16r)) NA_real_ else b16r$p_perm[1], SRC16, nt16)

# A17 — canonical-HRF GLM at five delays.
a17c10 <- if (is.null(a17c)) NULL else a17c %>% filter(channel_index == 10, axis == "comfort_wi")
d17 <- rd(an("analysis_17_hrf_glm", "hrf_delay_sweep.csv"))
SRC17 <- "output/analysis_17_hrf_glm/{channel_family_glm,hrf_delay_sweep}.csv"
s2("A17 canonical-HRF GLM, five delays", "estimate at 5-s delay (uM)",
   if (is.null(d17)) NA_real_ else d17$ch10_beta[d17$hrf_peak_s == 5][1], SRC17,
   if (is.null(d17)) MISSING_NOTE else "")
s2("A17 canonical-HRF GLM, five delays", "p at 5-s delay",
   if (is.null(d17)) NA_real_ else d17$ch10_p[d17$hrf_peak_s == 5][1], SRC17,
   if (is.null(d17)) MISSING_NOTE else "")
s2("A17 canonical-HRF GLM, five delays", "q (channel family)",
   if (is.null(a17c10)) NA_real_ else a17c10$q_BH[1], SRC17,
   if (is.null(a17c10)) MISSING_NOTE else "")
s2("A17 canonical-HRF GLM, five delays", "minimum p over delays",
   if (is.null(d17)) NA_real_ else min(d17$ch10_p), SRC17,
   if (is.null(d17)) MISSING_NOTE else "")
s2("A17 canonical-HRF GLM, five delays", "maximum p over delays",
   if (is.null(d17)) NA_real_ else max(d17$ch10_p), SRC17,
   if (is.null(d17)) MISSING_NOTE else "")

# A19 — re-filter corners (valid, no-motion-correction arm).
c19 <- rd(an("analysis_19_refilter", "channel10_by_setting.csv"))
SRC19 <- "output/analysis_19_refilter/channel10_by_setting.csv"
for (hp in c("0.01", "0.005", "0.002", "0")) {
  rr <- if (is.null(c19)) NULL else
    c19 %>% filter(motion_correction == "none", hpf == as.numeric(hp))
  lab <- paste0("hpf ", if (hp == "0") "none" else paste0(hp, " Hz"))
  s2(paste0("A19 re-filter, ", lab), "response estimate (uM)",
     if (is.null(rr) || !nrow(rr)) NA_real_ else rr$response_est[1], SRC19,
     if (is.null(rr) || !nrow(rr)) MISSING_NOTE else "")
  s2(paste0("A19 re-filter, ", lab), "response p",
     if (is.null(rr) || !nrow(rr)) NA_real_ else rr$response_p[1], SRC19,
     if (is.null(rr) || !nrow(rr)) MISSING_NOTE else "")
}

# A21 — void-test bound, 34 live-channel tests.
v21 <- rd(an("analysis_21_correction_sweep", "channel_family_sensitivity.csv"))
v21r <- if (is.null(v21)) NULL else v21 %>% filter(family_size == 34)
SRC21 <- "output/analysis_21_correction_sweep/channel_family_sensitivity.csv"
nt21 <- if (is.null(v21r) || !nrow(v21r)) MISSING_NOTE else ""
s2("A21 void-test bound, 34 live-channel tests", "estimate (uM)",
   if (is.null(v21r) || !nrow(v21r)) NA_real_ else v21r$estimate[1], SRC21, nt21)
s2("A21 void-test bound, 34 live-channel tests", "p",
   if (is.null(v21r) || !nrow(v21r)) NA_real_ else v21r$p[1], SRC21, nt21)
s2("A21 void-test bound, 34 live-channel tests", "q",
   if (is.null(v21r) || !nrow(v21r)) NA_real_ else v21r$q_recomputed[1], SRC21, nt21)

tableS2 <- bind_rows(s2_rows)
write_outcome(tableS2, file.path(TAB_DIR, "tableS2_channel10_adjudication.csv"))

# ------------------------------------------------------------ Table S3
# Gate-C frozen vs re-filter (tableS3_refilter_gates.tex). The frozen side is
# stored in gate_c_frozen_vs_refilter.csv for the difference-score row; the
# other frozen quantities come from the A03/A04/A06 campaign outputs.

gc <- rd(an("analysis_19_refilter", "gate_c_frozen_vs_refilter.csv"))
bs <- rd(an("analysis_19_refilter", "baseline_share_by_setting.csv"))
cp <- rd(an("analysis_19_refilter", "cycle_position_by_setting.csv"))
df19 <- rd(an("analysis_19_refilter", "declared_family_by_setting.csv"))
lmm03 <- rd(an("analysis_03_signal_anatomy", "cycle_position_lmm.csv"))

gc_row <- gc %>% filter(mc == "none", roi == "PFC_Frontal", axis == "comfort_wi")
fw_base <- if (is.null(fw)) NULL else
  fw %>% filter(analysis == "baseline", roi == "PFC_Frontal", term == "comfort_wi")
df_base <- if (is.null(df19)) NULL else
  df19 %>% filter(mc == "none", hpf == 0.01, window == "baseline",
                  roi == "PFC_Frontal", axis == "comfort_wi")
bs_pfc <- if (is.null(bs)) NULL else bs %>% filter(mc == "none", hpf == 0.01, roi == "PFC_Frontal")
cp_none <- if (is.null(cp)) NULL else cp %>% filter(motion_correction == "none", hpf == 0.01)
lmm_opens <- if (is.null(lmm03)) NULL else lmm03 %>% filter(term == "opens")
c19_001 <- if (is.null(c19)) NULL else
  c19 %>% filter(motion_correction == "none", hpf == 0.01)

# Frozen baseline share of the difference-score comfort slope: the A03-style
# exact decomposition (raw-uM windows, response-minus-baseline = response -
# baseline), recomputed from the shipped presentation-level intermediate.
frozen_share <- NA_real_
bm <- tryCatch(load_block_metrics(active_only = TRUE) %>%
                 filter(chromophore == "HbO"), error = function(e) NULL)
if (!is.null(bm)) {
  alt <- bm %>%
    group_by(subject_id, sequence, sequence_position, marker) %>%
    summarise(rmb  = mean(response_minus_baseline_m5_30),
              base = mean(baseline_m10_0), .groups = "drop") %>%
    group_by(subject_id, marker) %>%
    summarise(rmb = mean(rmb), base = mean(base), .groups = "drop") %>%
    left_join(load_p24() %>% select(subject_id, marker, comfort_wi, richness_wi),
              by = c("subject_id", "marker"))
  e_rmb  <- fixef(p24_lmm(alt, "comfort_wi + richness_wi", "rmb"))["comfort_wi"]
  e_base <- fixef(p24_lmm(alt, "comfort_wi + richness_wi", "base"))["comfort_wi"]
  frozen_share <- unname(-e_base / e_rmb)
}

s3_rows <- list()
s3 <- function(campaign_result, frozen_estimate, frozen_p,
               refilter_estimate, refilter_p, source_output, note = "") {
  s3_rows[[length(s3_rows) + 1L]] <<-
    tibble(campaign_result,
           frozen_estimate = as.numeric(frozen_estimate),
           frozen_p = as.numeric(frozen_p),
           refilter_estimate = as.numeric(refilter_estimate),
           refilter_p = as.numeric(refilter_p),
           source_output, note)
}

s3("Prefrontal difference-score comfort slope (uM)",
   if (is.null(gc) || !nrow(gc_row)) NA_real_ else gc_row$frozen_estimate[1],
   if (is.null(gc) || !nrow(gc_row)) NA_real_ else gc_row$frozen_p[1],
   if (is.null(gc) || !nrow(gc_row)) NA_real_ else gc_row$refilter_estimate[1],
   if (is.null(gc) || !nrow(gc_row)) NA_real_ else gc_row$refilter_p[1],
   "output/analysis_19_refilter/gate_c_frozen_vs_refilter.csv",
   if (is.null(gc)) MISSING_NOTE else "")
s3("Prefrontal baseline comfort slope (uM)",
   if (is.null(fw_base)) NA_real_ else fw_base$estimate[1],
   if (is.null(fw_base)) NA_real_ else fw_base$p[1],
   if (is.null(df_base) || !nrow(df_base)) NA_real_ else df_base$estimate[1],
   if (is.null(df_base) || !nrow(df_base)) NA_real_ else df_base$p[1],
   "output/analysis_04_appraisal_to_haemodynamics/family_by_window.csv; output/analysis_19_refilter/declared_family_by_setting.csv",
   if (is.null(fw_base) || is.null(df_base)) MISSING_NOTE else "")
s3("Baseline share of the difference slope (proportion)",
   frozen_share, NA_real_,
   if (is.null(bs_pfc) || !nrow(bs_pfc)) NA_real_ else bs_pfc$baseline_share[1],
   NA_real_,
   "computed from data/intermediates/tidy_block_metrics.csv (frozen side); output/analysis_19_refilter/baseline_share_by_setting.csv (re-filter side)",
   paste0(if (is.na(frozen_share)) "frozen side not regenerable; " else "",
          if (is.null(bs_pfc)) MISSING_NOTE else ""))
s3("Cycle-position contrast (uM)",
   if (is.null(lmm_opens) || !nrow(lmm_opens)) NA_real_ else lmm_opens$estimate[1],
   if (is.null(lmm_opens) || !nrow(lmm_opens)) NA_real_ else lmm_opens$p[1],
   if (is.null(cp_none) || !nrow(cp_none)) NA_real_ else cp_none$diff_contrast[1],
   if (is.null(cp_none) || !nrow(cp_none)) NA_real_ else cp_none$diff_p[1],
   "output/analysis_03_signal_anatomy/cycle_position_lmm.csv; output/analysis_19_refilter/cycle_position_by_setting.csv",
   paste0(if (is.null(lmm_opens)) paste(MISSING_NOTE, "(frozen side); ") else "",
          if (is.null(cp_none)) MISSING_NOTE else ""))
s3("Channel-10 response-window comfort slope (uM)",
   if (is.null(c10r)) NA_real_ else c10r$estimate[1],
   if (is.null(c10r)) NA_real_ else c10r$p[1],
   if (is.null(c19_001) || !nrow(c19_001)) NA_real_ else c19_001$response_est[1],
   if (is.null(c19_001) || !nrow(c19_001)) NA_real_ else c19_001$response_p[1],
   "output/analysis_06_channelwise_spatial_scan/channel_scan_response.csv; output/analysis_19_refilter/channel10_by_setting.csv",
   paste0(if (is.null(c10r)) paste(MISSING_NOTE, "(frozen side); ") else "",
          if (is.null(c19_001)) MISSING_NOTE else ""))

tableS3 <- bind_rows(s3_rows)
write_outcome(tableS3, file.path(TAB_DIR, "tableS3_refilter_gates.csv"))

# ------------------------------------------------------------ summary
cat("\nmake_tables.R: wrote", length(list.files(TAB_DIR, pattern = "\\.csv$")),
    "table CSVs to", TAB_DIR, "\n")
if (length(.missing)) {
  cat("sources not yet regenerated (values NA, flagged in the note column):\n")
  for (m in .missing) cat("  -", m, "\n")
}
cat("make_tables.R done.\n")
