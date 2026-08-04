# =============================================================================
# load_data.R — single data entry for the fnirs-soundscape-audit release.
#
# Every module sources this file and takes its data from the loaders below.
# Data source: the Mendeley Data workbook
#   data/Soundscape_appraisal_fNIRS_prefrontal_temporal_data.xlsx
# (DOI 10.17632/cr2tpv999j, CC BY 4.0). Sheets are read with col_types = "text"
# and types are re-inferred through a CSV round-trip, so "NA" and "" both map
# to NA and the XLSX round-trip is loss-free at the 15-significant-figure
# precision Excel stores.
#
# UNITS. The workbook stores haemoglobin in mol/L (order 1e-7). Every loader
# returns it in micromol/L (uM) — multiplied by 1e6 once, at the entry point.
# Never scale again downstream.
#
# All paths are relative to the repository root; run everything via run_all.R.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(lmerTest)
  library(broom)
  library(broom.mixed)
})

XLSX_PATH <- file.path("data", "Soundscape_appraisal_fNIRS_prefrontal_temporal_data.xlsx")
INTERMEDIATE_DIR <- file.path("data", "intermediates")

BOOT_N <- 1000   # house default

ROIS <- c("PFC_Frontal", "Left_Temporal", "Right_Temporal")
ROI_LABELS <- c(PFC_Frontal    = "Prefrontal/frontal",
                Left_Temporal  = "Left temporal",
                Right_Temporal = "Right temporal")
CHROMOPHORES <- c("HbO", "HbR", "HbT")

# The eight rating adjectives, in the order the circumplex formulae use them.
ADJECTIVES <- tribble(
  ~zh,       ~en,           ~axis_role,
  "舒适的", "comfortable", "comfort_pole",
  "混乱的", "chaotic",     "comfort_pole_neg",
  "丰富的", "rich",        "richness_pole",
  "单调的", "monotonous",  "richness_pole_neg",
  "有趣的", "interesting", "diagonal",
  "无聊的", "boring",      "diagonal",
  "热闹的", "lively",      "diagonal",
  "平淡的", "bland",       "diagonal",
)

# Exact circumplex constants.
CIRCUMPLEX_K <- 1 / (4 + 4 * sqrt(2))
CIRCUMPLEX_D <- CIRCUMPLEX_K * sqrt(2) / 2

# ------------------------------------------------------------- sheet reading

# The workbook uses English adjective column names; the analysis code uses the
# Chinese originals the rating form carried. Renamed here, once.
# dplyr::rename() reads a named vector as new_name = "old_name".
.ADJ_RENAME <- c("adj_舒适的" = "adj_comfortable", "adj_热闹的" = "adj_lively",
                 "adj_有趣的" = "adj_interesting", "adj_单调的" = "adj_monotonous",
                 "adj_平淡的" = "adj_bland", "adj_混乱的" = "adj_chaotic",
                 "adj_丰富的" = "adj_varied", "adj_无聊的" = "adj_boring")

#' Read one workbook sheet. col_types = "text" at the XLSX boundary, then a
#' CSV round-trip lets readr re-infer column types (readr maps both "NA" and
#' "" to NA, so the round-trip is loss-free).
read_sheet <- function(name) {
  if (!file.exists(XLSX_PATH)) {
    stop("Workbook not found at ", XLSX_PATH, "\n",
         "Download it from Mendeley Data (DOI 10.17632/cr2tpv999j) and place ",
         "it in data/ — see README.", call. = FALSE)
  }
  raw <- readxl::read_excel(XLSX_PATH, sheet = name, col_types = "text")
  tf <- tempfile(fileext = ".csv")
  on.exit(unlink(tf), add = TRUE)
  readr::write_csv(raw, tf)
  d <- readr::read_csv(tf, show_col_types = FALSE)
  renamed <- .ADJ_RENAME[.ADJ_RENAME %in% names(d)]
  if (length(renamed)) d <- d %>% rename(all_of(renamed))
  d
}

#' Read a shipped intermediate CSV (data/intermediates/). These are the small
#  derived tables whose source (the raw SNIRF recordings and the delivered
#  continuous concentration series) is not part of the Mendeley deposit;
#  README documents each one.
load_intermediate <- function(filename, ...) {
  readr::read_csv(file.path(INTERMEDIATE_DIR, filename), show_col_types = FALSE, ...)
}

# ------------------------------------------------------------- unit conversion

# Concentration columns follow two naming conventions:
# WIDE  — one column per chromophore x spatial unit, e.g. HbO_head
# LONG  — window-statistic value columns, e.g. baseline_m10_0
.CONC_PATTERN <- "^(HbO|HbR|HbT)_|^(mean_)?(baseline|stim|response)_m[0-9]|^(mean_)?response_minus_baseline"

.to_micromolar <- function(d) {
  cc <- grep(.CONC_PATTERN, names(d), value = TRUE)
  if (length(cc) == 0) {
    stop("No concentration column recognised — check .CONC_PATTERN against: ",
         paste(names(d), collapse = ", "), call. = FALSE)
  }
  d <- d %>% mutate(across(all_of(cc), ~ .x * 1e6))
  attr(d, "hb_unit") <- "umol/L"
  attr(d, "hb_columns_converted") <- cc
  d
}

# The workbook reports age in five-year bands ("25-29"); the analysis uses the
# band's lower bound as the numeric age (the frozen pipeline held exact integer
# ages, which the deposit anonymises — see README). The open band "<=19" has no
# lower bound and maps to NA.
.add_age <- function(d) {
  if ("age_band" %in% names(d) && !"age" %in% names(d)) {
    d <- d %>% mutate(age = suppressWarnings(as.integer(substr(age_band, 1, 2))),
                      .after = age_band)
  }
  d
}

# ------------------------------------------------------- design decompositions

# Within/between decomposition of a within-subject predictor, plus the split of
# the within part into the component shared by everyone who heard the same clip
# (stimulus-driven) and the residual (idiosyncratic appraisal).
.decompose <- function(d, var) {
  bw   <- paste0(var, "_bw");   wi   <- paste0(var, "_wi")
  stim <- paste0(var, "_stim"); idio <- paste0(var, "_idio")
  clip_mean <- paste0(var, "_clip_mean")
  d %>%
    group_by(subject_id) %>% mutate(!!bw := mean(.data[[var]], na.rm = TRUE)) %>% ungroup() %>%
    mutate(!!wi := .data[[var]] - .data[[bw]]) %>%
    group_by(clip_id) %>% mutate(!!clip_mean := mean(.data[[var]], na.rm = TRUE)) %>% ungroup() %>%
    group_by(group) %>% mutate(!!stim := .data[[clip_mean]] - mean(.data[[clip_mean]], na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(!!idio := .data[[wi]] - .data[[stim]])
}

.add_design <- function(d) {
  d %>%
    mutate(
      subject   = factor(subject_id),
      clip      = factor(clip_id),
      group_f   = factor(group),
      form_f    = factor(form),
      marker    = factor(marker, levels = c("M3", "M4", "M5", "M6")),
      quadrant  = factor(quadrant, levels = c("Q1", "Q2", "Q3", "Q4")),
      sex       = factor(sex_code, levels = c(1, 2), labels = c("male", "female")),
      inferred_assignment = as.integer(event_assignment == "inferred"),
    )
}

# ----------------------------------------------------------------- the loaders

#' The primary analysis table: one row per subject x condition, 276 rows.
#' Haemoglobin in uM. `sample` is the only sanctioned way to change the row set:
#'   "full" (default) / "observed" / "on_quadrant" / "qc_clean" / "covariates".
load_p24 <- function(sample = c("full", "observed", "on_quadrant", "qc_clean", "covariates")) {
  sample <- match.arg(sample)
  d <- read_sheet("analysis_table") %>% .add_age() %>% .to_micromolar() %>%
    .add_design() %>% .decompose("comfort") %>% .decompose("richness")

  n0 <- nrow(d)
  d <- switch(sample,
    full        = d,
    observed    = d %>% filter(event_assignment == "positional"),
    on_quadrant = d %>% filter(flag_off_quadrant_stimulus == 0),
    qc_clean    = d %>% filter(flag_qc_check == 0, flag_onset_displaced == 0),
    covariates  = d %>% filter(has_covariates == 1)
  )
  attr(d, "hb_unit") <- "umol/L"
  attr(d, "sample")  <- sample
  attr(d, "n_dropped") <- n0 - nrow(d)
  d
}

#' Trial-level table: one row per subject x presentation, 826 rows, HbO only.
load_p24_block <- function() {
  read_sheet("block_level") %>% .add_age() %>% .to_micromolar() %>% .add_design() %>%
    .decompose("comfort") %>% .decompose("richness") %>%
    mutate(sequence = factor(sequence, levels = c("S1", "S2", "S3")))
}

#' subject x condition x channel x chromophore, 19,044 rows (sheet
#' channel_condition) — the layer for channel-wise spatial scans and for the
#' masked ROI aggregation below.
load_condition_metrics <- function(active_only = TRUE) {
  d <- read_sheet("channel_condition") %>% .to_micromolar()
  if (active_only) d <- d %>% filter(active_pair == 1)
  d
}

#' subject x presentation x channel x chromophore, 56,994 rows, HbO/HbR/HbT.
#' NOT in the Mendeley workbook: this is the presentation-level layer the
#' responder-heterogeneity audit works on, shipped as an intermediate derived
#' from the delivered concentration series (see README).
load_block_metrics <- function(active_only = TRUE) {
  d <- load_intermediate("tidy_block_metrics.csv") %>% .to_micromolar()
  if (active_only) d <- d %>% filter(active_pair == 1)
  d %>% mutate(sequence = factor(sequence, levels = c("S1", "S2", "S3")),
               marker = factor(marker, levels = c("M3", "M4", "M5", "M6")))
}

# The canonical design/predictor columns used when a lower-level haemodynamic
# table is joined to the subject x condition analysis table.
.analysis_design <- function(sample = c("full", "observed", "on_quadrant",
                                        "qc_clean", "covariates")) {
  sample <- match.arg(sample)
  load_p24(sample) %>%
    transmute(subject_id, marker, group, form,
           quadrant_design = as.character(quadrant),
           clip_id, comfort, richness,
           comfort_bw, comfort_wi, comfort_stim, comfort_idio,
           richness_bw, richness_wi, richness_stim, richness_idio,
           sex_code, age, who5_total_0_25, nss_total_keyed, has_covariates,
           event_assignment, flag_onset_displaced, flag_qc_check,
           flag_off_quadrant_stimulus, flag_shared_q4_waveform)
}

#' Active subject x condition x channel rows joined to the canonical design.
load_condition_analysis <- function(chromophore = "HbO",
                                    sample = c("full", "observed", "on_quadrant",
                                               "qc_clean", "covariates")) {
  sample <- match.arg(sample)
  load_condition_metrics(active_only = TRUE) %>%
    filter(.data$chromophore == !!chromophore) %>%
    inner_join(.analysis_design(sample), by = c("subject_id", "marker"))
}

#' Active subject x presentation x channel rows joined to the canonical design.
load_block_analysis <- function(chromophore = "HbO",
                                sample = c("full", "observed", "on_quadrant",
                                           "qc_clean", "covariates")) {
  sample <- match.arg(sample)
  load_block_metrics(active_only = TRUE) %>%
    filter(.data$chromophore == !!chromophore) %>%
    inner_join(.analysis_design(sample), by = c("subject_id", "marker")) %>%
    mutate(
      run_index = ceiling(event_index / 4),
      position_in_run = event_index - 4 * (run_index - 1),
      run_f = factor(run_index),
      position_f = factor(position_in_run),
      first5_m0_5 = 6 * stim_m0_30 - 5 * response_m5_30
    )
}

#' subject x clip: the 8 raw adjectives plus recomputed and delivered axes.
load_ratings <- function() read_sheet("subjective_ratings")

#' One row per subject: design cell, covariates, channel count, flags.
load_subjects <- function() read_sheet("subjects") %>% .add_age()

#' subject x channel: ROI membership, source-detector distance, active flag.
load_channels <- function() read_sheet("channels")

#' group x quadrant: clip id, recording site and sources, off-quadrant flags.
load_stimuli <- function() read_sheet("stimuli")

#' subject x scheduled slot: the acquisition record.
load_events <- function() read_sheet("events")

# ------------------------------------------------------- dead-channel masking
#
# THE MASK IS THE UNION OF TWO COMPLEMENTARY CRITERIA, reconstructed here from
# the channel_mask sheet of the workbook (see README):
#   (1) dark-fraction — a subject x channel cell whose recording spends
#       >= DARK_FRACTION_CUT of its samples at the instrument dark floor
#       (frac_dark), plus the single rail-saturated cell a dark-floor test
#       cannot see (flat_not_dark == 1: P10 channel 12, pinned at 1.000000);
#   (2) detector exclusion — detector 8 (T8, channels 20/21) is dropped for
#       every subject: dark in 29 of 69 sessions and non-neurovascular when it
#       registers (HbO/HbR coupling +0.67 against -0.46 at detector 10 in the
#       same subjects).
# The sheet carries the resulting flags (mask_excluded, detector_excluded,
# use_in_aggregation); load_channel_mask() RECOMPUTES them from frac_dark /
# flat_not_dark / detector and stops if the recomputation disagrees with the
# deposited flags, so the mask logic is code, not a taken-on-trust column.
DARK_FRACTION_CUT <- 0.05
UNUSABLE_DETECTORS <- 8L

#' The union channel mask: one row per subject x channel cell, with the
#' recomputed flags verified against the deposited ones.
load_channel_mask <- function() {
  m <- read_sheet("channel_mask")
  recomputed <- m %>%
    transmute(subject_id, channel_index,
              mask_excluded_r = as.integer(frac_dark >= DARK_FRACTION_CUT |
                                             flat_not_dark == 1),
              detector_excluded_r = as.integer(detector %in% UNUSABLE_DETECTORS))
  chk <- m %>%
    select(subject_id, channel_index, mask_excluded, detector_excluded,
           active_pair, use_in_aggregation) %>%
    mutate(across(c(mask_excluded, detector_excluded, active_pair,
                    use_in_aggregation), as.integer)) %>%
    left_join(recomputed, by = c("subject_id", "channel_index"))
  # The deposited flags come back from the CSV round-trip as double; compare
  # as integers so the check tests values, not storage type.
  stopifnot(
    identical(as.integer(chk$mask_excluded), chk$mask_excluded_r),
    identical(as.integer(chk$detector_excluded), chk$detector_excluded_r),
    identical(as.integer(chk$use_in_aggregation),
              as.integer(chk$active_pair == 1 & chk$mask_excluded == 0 &
                           chk$detector_excluded == 0))
  )
  m
}

#' Excluded subject x channel cells under the union mask.
load_masked_out_cells <- function() {
  load_channel_mask() %>%
    filter(active_pair == 1, use_in_aggregation == 0) %>%
    distinct(subject_id, channel_index)
}

#' The superseded session-CV mask of 2026-07-30 (raw-intensity CV < 0.002 per
#' subject x channel, from the raw recordings — shipped as an intermediate).
#' Retained because the Fig 2c mirror-pair counts and the A23 audit were
#' declared against it.
load_dead_channels <- function() {
  load_intermediate("qa_detector_dropout_0730.csv") %>%
    filter(.data$dead_raw == 1) %>%
    select(subject_id, channel_index, roi, detector, raw_intensity_cv)
}

# Fixed dead-channel vector from the earliest pass (channels 18/19, detector 9,
# dark in all 69 subjects). Kept for the adjudication modules that declared it.
FLAT_CHANNELS <- c(18L, 19L)

#' Rebuild subject x condition ROI means from the condition x channel layer,
#' keeping only each subject's own usable channels (union mask applied).
#'
#' THE DEFAULT for any ROI quantity reported under the mask. The aggregation
#' (mean over the cell's usable channels) reproduces the delivered ROI columns
#' exactly when nothing is excluded — check_roi_rebuild() gates this.
#'
#' Returns the delivered columns (HbO_<roi>) alongside the masked ones for all
#' three windows (HbOm_<roi> difference, HbOm_resp_<roi>, HbOm_base_<roi>),
#' plus n_live_<roi>.
load_p24_roi_masked <- function(chromophore = "HbO") {
  dead <- load_masked_out_cells()
  map <- load_channels() %>% filter(active_pair == 1) %>%
    distinct(subject_id, channel_index, roi)
  agg <- load_condition_metrics(active_only = TRUE) %>%
    filter(.data$chromophore == !!chromophore) %>%
    left_join(dead %>% mutate(.dead = TRUE), by = c("subject_id", "channel_index")) %>%
    filter(is.na(.data$.dead)) %>%
    left_join(map, by = c("subject_id", "channel_index")) %>%
    filter(!is.na(roi)) %>%
    group_by(subject_id, marker, roi) %>%
    summarise(
      diff_m = mean(mean_response_minus_baseline_m5_30),
      resp_m = mean(mean_response_m5_30),
      base_m = mean(mean_baseline_m10_0),
      n_live = n(), .groups = "drop")
  wide <- function(v, prefix) {
    agg %>% mutate(roi = paste0(prefix, roi)) %>%
      select(subject_id, marker, roi, all_of(v)) %>%
      pivot_wider(names_from = roi, values_from = all_of(v))
  }
  load_p24() %>%
    left_join(wide("diff_m", paste0(chromophore, "m_")), by = c("subject_id", "marker")) %>%
    left_join(wide("resp_m", paste0(chromophore, "m_resp_")), by = c("subject_id", "marker")) %>%
    left_join(wide("base_m", paste0(chromophore, "m_base_")), by = c("subject_id", "marker")) %>%
    left_join(agg %>% mutate(roi = paste0("n_live_", roi)) %>%
                select(subject_id, marker, roi, n_live) %>%
                pivot_wider(names_from = roi, values_from = n_live),
              by = c("subject_id", "marker"))
}

#' Rebuild ROI means excluding a FIXED channel set (superseded; kept for the
#' channel-10 adjudication modules that declared it).
load_p24_roi_corrected <- function(exclude_channels = FLAT_CHANNELS,
                                   chromophore = "HbO") {
  map <- load_channels() %>% filter(active_pair == 1) %>%
    distinct(subject_id, channel_index, roi)
  agg <- load_condition_metrics(active_only = TRUE) %>%
    filter(.data$chromophore == !!chromophore) %>%
    filter(!channel_index %in% exclude_channels) %>%
    left_join(map, by = c("subject_id", "channel_index")) %>%
    filter(!is.na(roi)) %>%
    group_by(subject_id, marker, roi) %>%
    summarise(value = mean(mean_response_minus_baseline_m5_30), .groups = "drop") %>%
    mutate(roi = paste0(chromophore, "c_", roi)) %>%
    pivot_wider(names_from = roi, values_from = value)
  load_p24() %>% left_join(agg, by = c("subject_id", "marker"))
}

#' Verify that the condition-layer ROI rebuild reproduces the delivered ROI
#' columns when nothing is excluded. Returns the worst absolute discrepancy
#' per ROI (uM); run_all.R asserts it stays below 1e-9.
check_roi_rebuild <- function(chromophore = "HbO") {
  map <- load_channels() %>% filter(active_pair == 1) %>%
    distinct(subject_id, channel_index, roi)
  agg <- load_condition_metrics(active_only = TRUE) %>%
    filter(.data$chromophore == !!chromophore) %>%
    left_join(map, by = c("subject_id", "channel_index")) %>%
    filter(!is.na(roi)) %>%
    group_by(subject_id, marker, roi) %>%
    summarise(value = mean(mean_response_minus_baseline_m5_30), .groups = "drop")
  d <- load_p24() %>%
    select(subject_id, marker, all_of(paste0(chromophore, "_", ROIS))) %>%
    pivot_longer(all_of(paste0(chromophore, "_", ROIS)),
                 names_to = "roi", values_to = "delivered",
                 names_prefix = paste0(chromophore, "_")) %>%
    left_join(agg, by = c("subject_id", "marker", "roi"))
  map_dfr(ROIS, function(rr) {
    dd <- d %>% filter(roi == rr)
    keep <- !is.na(dd$delivered) & !is.na(dd$value)
    tibble(roi = rr, n_compared = sum(keep), n_delivered_na = sum(is.na(dd$delivered)),
           max_abs_diff = max(abs(dd$delivered[keep] - dd$value[keep])))
  })
}

# ------------------------------------------------------------------- modelling

# Model backbone: y ~ <fixed> + (1 | subject_id).
p24_lmm <- function(data, fixed, y, keep_clip = FALSE, reml = TRUE) {
  re <- if (keep_clip) "(1 | subject_id) + (1 | clip_id)" else "(1 | subject_id)"
  f <- as.formula(paste(y, "~", fixed, "+", re))
  lmerTest::lmer(f, data = data, REML = reml,
                 control = lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4)))
}

#' z-score, NA-safe, returning a plain vector.
z <- function(x) as.vector(scale(x))

#' Benjamini-Hochberg q values. Family size is always stated by the caller.
bh <- function(p) p.adjust(p, method = "BH")

#' Cluster bootstrap by subject: resample subjects with replacement, refit,
#' return percentile CI for the named fixed-effect terms. BOOT_N = 1,000.
boot_lmm <- function(data, fixed, y, terms, n = BOOT_N, seed = 24) {
  set.seed(seed)
  subs <- unique(data$subject_id)
  out <- matrix(NA_real_, nrow = n, ncol = length(terms), dimnames = list(NULL, terms))
  for (b in seq_len(n)) {
    pick <- sample(subs, length(subs), replace = TRUE)
    dd <- map_dfr(seq_along(pick), function(i) {
      data %>% filter(subject_id == pick[i]) %>% mutate(subject_id = paste0(pick[i], "_b", i))
    })
    fit <- try(suppressMessages(suppressWarnings(p24_lmm(dd, fixed, y))), silent = TRUE)
    if (inherits(fit, "try-error")) next
    co <- lme4::fixef(fit)
    for (tm in terms) if (tm %in% names(co)) out[b, tm] <- co[[tm]]
  }
  tibble(term = terms,
         boot_n_ok = apply(out, 2, function(v) sum(!is.na(v))),
         ci_lo = apply(out, 2, quantile, 0.025, na.rm = TRUE),
         ci_hi = apply(out, 2, quantile, 0.975, na.rm = TRUE))
}

#' Write a CSV into the repo's output/ tree.
write_outcome <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(x, path)
  cat("  wrote", path, "\n")
}

#' Tee console output to a report txt while still printing it.
open_log <- function(path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- file(path, open = "wt")
  sink(con, split = TRUE)
  con
}
close_log <- function(con) { sink(); close(con) }

#' Compact fixed-effect table from an lmerTest fit.
fixef_table <- function(fit, keep = NULL) {
  s <- as.data.frame(summary(fit)$coefficients)
  s <- tibble(term = rownames(s), estimate = s[[1]], se = s[[2]],
              df = s[[3]], t = s[[4]], p = s[[5]])
  if (!is.null(keep)) s <- s %>% filter(term %in% keep)
  s
}
