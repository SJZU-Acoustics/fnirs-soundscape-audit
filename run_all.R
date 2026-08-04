#!/usr/bin/env Rscript
# =============================================================================
# run_all.R — master script for the fnirs-soundscape-audit release.
#
# Reproduces every figure, table and headline statistic of
#   "When design timing, preprocessing and channel quality undermine
#    environmental-appraisal fNIRS"
# from the Mendeley Data workbook (place it in data/ — see README.md) plus the
# small shipped intermediates in data/intermediates/ (see data/README.md).
#
# Usage, from the repository root:
#   Rscript run_all.R
# Outputs are written to output/ (figures/, data_lock/, tables/, and one
# folder per analysis module). The final verification pass checks the
# regenerated display values against the manuscript's quoted numbers and the
# ROI-rebuild gate below.
# =============================================================================

t0 <- Sys.time()
source("code/load_data.R")

# Synthesis lock: the condition-layer ROI rebuild must reproduce the delivered
# ROI columns when nothing is excluded (micromol/L, max abs diff < 1e-9).
gate <- check_roi_rebuild()
print(gate)
stopifnot(all(gate$max_abs_diff < 1e-9))
cat("ROI rebuild gate passed (max abs diff", max(gate$max_abs_diff), "uM).\n\n")

run_module <- function(path) {
  cat(strrep("=", 74), "\n", path, "\n", strrep("=", 74), "\n", sep = "")
  t1 <- Sys.time()
  source(path)
  cat(sprintf("  [%.1f s]\n\n", as.numeric(difftime(Sys.time(), t1, units = "secs"))))
}

# ---- analysis modules (workbook-recomputed) --------------------------------
analysis_modules <- c(
  "code/mod_mask_summary.R",          # union channel mask from channel_mask
  "code/mod02_appraisal_structure.R", # appraisal space (Fig 3)
  "code/mod03_signal_anatomy.R",      # run-opener / cycle-position anatomy
  "code/mod04_confirmatory_family.R", # declared ROI families (Fig 4, Table 2)
  "code/mod05_geometry.R",            # alternative response geometry
  "code/mod06_channel_scan.R",        # channel-wise spatial scan
  "code/mod07_channel10_adjudication.R",
  "code/mod08_dependence_cr2.R",      # CR2 cluster-robust audit
  "code/mod09_chromophore.R",         # HbR/HbT triangulation, HbT identity
  "code/mod10_person_side.R",         # trait reactivity and moderation
  "code/mod11_cycle_time_course.R",   # 5-s bin families (from shipped bins)
  "code/mod12_repetition.R",          # encounter trends (Fig 5b)
  "code/mod13_acoustic.R",            # stimulus acoustic features
  "code/mod14_detrended.R",           # sequence-detrended coupling
  "code/mod15_global_adjustment.R",   # leave-one-out global adjustment
  "code/mod16_ridge.R",               # cross-validated multivariate pattern
  "code/mod17_hrf_glm.R",             # canonical-HRF GLM (from shipped betas)
  "code/mod18_clip_reliability.R",    # clip reliability (Fig 5a)
  "code/mod19_refilter.R",            # re-filter gates (from shipped means)
  "code/mod20_detector_dropout.R",    # corrected-ROI rechecks (CV-mask era)
  "code/mod21_correction_sweep.R",    # void-test bound, gate-C decomposition
  "code/mod22_chain.R",               # the lambda chain (Fig 6)
  "code/mod23_responder_audit.R",     # responder heterogeneity (SI Fig S1)
  "code/mod24_dark_fraction_mask.R"   # union-mask declared families
)

# ---- figure scripts ----------------------------------------------------------
figure_scripts <- c(
  "code/mod_fig1b_filter_gain.R",     # filter amplitude response (pure recompute)
  "code/fig1_anatomy.R",
  "code/fig2_channel_quality.R",
  "code/fig3_appraisal_space.R",
  "code/fig4_confirmatory.R",
  "code/recompute_a18_position_null.R",
  "code/fig5_reproduced.R",
  "code/fig6_chain.R",
  "code/sifig1_inference_audit.R"
)

# ---- tables + verification ---------------------------------------------------
closing_scripts <- c(
  "code/make_tables.R",
  "code/verify_displays.R"
)

for (m in c(analysis_modules, figure_scripts, closing_scripts)) {
  if (file.exists(m)) run_module(m) else cat("SKIPPED (not present):", m, "\n\n")
}

# ---- optional raw-recording QA (documentation-grade) ------------------------
# The raw SNIRF recordings are not deposited (see README.md). The Python QA
# scripts behind the channel mask and the re-filter intermediates ship in
# code/python/ for inspection; they run only if the raw recordings are
# supplied locally, and are not part of the standard reproduction path.
if (dir.exists("data/raw_snirf")) {
  message("data/raw_snirf/ found — running raw-recording QA scripts.")
  system("python3 code/python/qa_detector_dropout.py")
  system("python3 code/python/qa_dropout_onset_0803.py")
  system("python3 code/python/refilter_pipeline.py")
} else {
  message("data/raw_snirf/ not present — raw-recording QA scripts skipped ",
          "(standard path; see README.md).")
}

cat(sprintf("\nrun_all.R finished in %.1f min.\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
