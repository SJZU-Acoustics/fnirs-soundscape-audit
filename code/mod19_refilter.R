# mod19_refilter.R — Does the 0.01 Hz high-pass create the campaign's null?
#
# Port of exploration/analysis_19_refilter for the release repo. The re-filter
# pipeline (code/python/refilter_pipeline.py, a Python re-implementation of the
# Homer3 chain) re-derives the concentration series from the raw SNIRF at four
# high-pass settings (0.01 / 0.005 / 0.002 Hz / none) under two motion-
# correction arms (none / wavelet), holding everything else fixed. The raw
# recordings are NOT part of the Mendeley deposit, and the pipeline's per-event
# channel-level output is too large to ship, so this module reads three compact
# intermediates derived from it (see README):
#
#   refilter_condition_metrics.csv  subject x marker x ROI x arm x hpf, HbO
#     window means. The three ROIs are aggregated the way the frozen builder
#     aggregates (per-event mean over the ROI's active channels, then the mean
#     over the cell's presentations); `head` is the mean over all of the
#     cell's channel x presentation rows, the aggregation the whole-head
#     summary uses.
#   refilter_event_metrics_roi.csv  subject x presentation x ROI x arm x hpf
#     per-event ROI means with channel counts, unmasked and under the declared
#     dropout mask (session-CV dead cells + detector 8; mod21's gate-C
#     decomposition recovers its pooled mean as the count-weighted mean of
#     these event means).
#   refilter_channel10_metrics.csv  subject x marker x arm x hpf at channel 10.
#
# UNITS. Haemoglobin columns in these intermediates are already in
# micromol/L (uM) — converted once, at build time. Do NOT rescale them.
#
# THREE GATES. Gate A (pruning) and gate B (waveform agreement) are builder-
# side and are verified here from their shipped CSVs. Gate C (functional) is
# this module's business: does the hpf = 0.01 branch reproduce the campaign's
# results — the response-window null, the baseline carrying the comfort slope,
# the cycle-position contrast, channel 10?
#
# Displays fed: SI Table S3 (five gate-C rows), SI Table S2 row A19
# (channel 10 across the four corners), Table 1 row 5 and Table 3 row 2
# (waveform median r = 0.84 no-MC / 0.66 wavelet, gate B, HbO).

source("code/load_data.R")

OUT <- file.path("output", "analysis_19_refilter")
con <- open_log(file.path(OUT, "refilter_analysis_report.txt"))
rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")

AXES <- c("comfort_wi", "richness_wi")

# ---------------------------------------------------------------- 0. gates A and B
rule("0. GATES A AND B — builder-side reproduction checks (shipped CSVs)")

gate_a <- load_intermediate("gate_a_pruning.csv")
cat("Gate A (pair pruning vs the delivered pairActive): ",
    nrow(gate_a), "/", nrow(gate_a), " subjects reproduced, ",
    sum(gate_a$n_disagreements), " channel disagreements, ",
    sum(gate_a$n_active_mine == gate_a$n_active_delivered),
    " matching active counts\n", sep = "")
stopifnot(nrow(gate_a) == 69L, sum(gate_a$n_disagreements) == 0L)

gate_b <- load_intermediate("gate_b_reproduction.csv") %>%
  filter(chromophore == "HbO")
wave_r <- gate_b %>% group_by(motion_correction) %>%
  summarise(median_r = median(r), n = n(), .groups = "drop")
cat("Gate B (waveform agreement with the delivered series at hpf = 0.01, HbO):\n")
print(as.data.frame(wave_r), digits = 3, row.names = FALSE)
cat("  -> Table 1 row 5 / Table 3 row 2: median r = 0.84 (no MC) / 0.66 (wavelet)\n")

# ---------------------------------------------------------------- 1. assemble
rule("1. ASSEMBLE — re-filtered condition cells at every filter setting")

cond <- load_intermediate("refilter_condition_metrics.csv")
design <- .analysis_design("full")
cell_roi <- cond %>% filter(roi != "head") %>%
  inner_join(design, by = c("subject_id", "marker"))

cat("settings:", n_distinct(paste(cond$motion_correction, cond$hpf)),
    " cells per setting (3 ROIs):",
    nrow(cell_roi) / n_distinct(paste(cond$motion_correction, cond$hpf)) / 3, "\n")

# ------------------------------------------------------- 2. declared families
rule("2. THE DECLARED PRIMARY FAMILY AT EVERY SETTING (2 axes x 3 ROIs = 6, BH)")

WINDOWS <- c(difference = "response_minus_baseline_m5_30",
             response   = "response_m5_30",
             baseline   = "baseline_m10_0")

fam <- expand_grid(mc = unique(cell_roi$motion_correction),
                   hpf = sort(unique(cell_roi$hpf)),
                   window = names(WINDOWS)) %>%
  pmap_dfr(function(mc, hpf, window) {
    d0 <- cell_roi %>% filter(motion_correction == mc, hpf == !!hpf)
    map_dfr(ROIS, function(rr) {
      d <- d0 %>% filter(roi == rr)
      m <- p24_lmm(d, paste(AXES, collapse = " + "), WINDOWS[[window]])
      fixef_table(m, keep = AXES) %>%
        transmute(mc, hpf, window, roi = rr, axis = term,
                  estimate, se, p)
    })
  }) %>%
  group_by(mc, hpf, window) %>% mutate(q_BH = bh(p)) %>% ungroup()

for (w in names(WINDOWS)) {
  cat("\n---", w, "window — comfort axis, all settings ---\n")
  print(as.data.frame(fam %>% filter(window == w, axis == "comfort_wi") %>%
    select(mc, hpf, roi, estimate, se, p, q_BH) %>% arrange(mc, hpf, roi)),
    digits = 3, row.names = FALSE)
}
cat("\nMinimum q in each declared family of 6:\n")
print(as.data.frame(fam %>% group_by(mc, hpf, window) %>%
  summarise(min_q = min(q_BH), n_survivors = sum(q_BH < 0.05), .groups = "drop") %>%
  arrange(window, mc, hpf)), digits = 3, row.names = FALSE)
write_outcome(fam, file.path(OUT, "declared_family_by_setting.csv"))

# ------------------------------------------- 3. GATE C: functional equivalence
rule("3. GATE C — does hpf = 0.01 reproduce the campaign's frozen results?")

frozen <- load_p24()
# frozen reference values, recomputed here from the workbook so the comparison
# is like-for-like (same model, same sample, same axes)
frozen_ref <- map_dfr(ROIS, function(rr) {
  m <- p24_lmm(frozen, paste(AXES, collapse = " + "), paste0("HbO_", rr))
  fixef_table(m, keep = AXES) %>% transmute(roi = rr, axis = term,
                                            frozen_estimate = estimate, frozen_p = p)
})
cmp <- fam %>% filter(window == "difference", hpf == 0.01) %>%
  select(mc, roi, axis, refilter_estimate = estimate, refilter_p = p) %>%
  inner_join(frozen_ref, by = c("roi", "axis"))
cat("Frozen difference score vs the re-filtered hpf = 0.01 branch:\n")
print(as.data.frame(cmp), digits = 3, row.names = FALSE)
cat("\nSign agreement:", sum(sign(cmp$refilter_estimate) == sign(cmp$frozen_estimate)),
    "/", nrow(cmp), "\n")
cat("Correlation of estimates across the 12 cells:",
    round(cor(cmp$refilter_estimate, cmp$frozen_estimate), 3), "\n")
write_outcome(cmp, file.path(OUT, "gate_c_frozen_vs_refilter.csv"))

# SI Table S3 gate-C rows, PFC comfort: difference -0.052 (p = 0.034) frozen vs
# -0.057 (p = 0.042) rebuild; baseline +0.031 (p = 0.083) frozen (the campaign's
# frozen baseline window) vs +0.037 (p = 0.064) rebuild. The frozen difference
# row is recomputed above; the frozen baseline anchor is quoted from the frozen
# campaign (the analysis table carries difference scores only).
rebuild_base <- fam %>% filter(window == "baseline", hpf == 0.01, mc == "none",
                               roi == "PFC_Frontal", axis == "comfort_wi")
cat("\nSI Table S3 PFC comfort rows — rebuild side at hpf = 0.01 (no MC):",
    "difference", round(cmp$refilter_estimate[cmp$mc == "none" &
                        cmp$roi == "PFC_Frontal" & cmp$axis == "comfort_wi"], 3),
    "(p =", signif(cmp$refilter_p[cmp$mc == "none" & cmp$roi == "PFC_Frontal" &
                                    cmp$axis == "comfort_wi"], 2), "),",
    "baseline", round(rebuild_base$estimate, 3),
    "(p =", signif(rebuild_base$p, 2), ")\n")

# --------------------------- 4. the diagnosis: where does the slope live?
rule("4. THE DIAGNOSIS — baseline share of the difference-score comfort slope")

# On the frozen data 68% of the prefrontal comfort slope is pre-stimulus. If the
# 0.01 Hz high-pass creates that, relaxing it should move the slope out of the
# baseline and into the response window.
share <- fam %>% filter(axis == "comfort_wi") %>%
  select(mc, hpf, window, roi, estimate) %>%
  pivot_wider(names_from = window, values_from = estimate) %>%
  mutate(baseline_share = baseline / (baseline - response),
         check_difference = difference)
cat("Comfort slope by window (uM per appraisal unit), and the baseline's share\n")
cat("of the difference score. difference = response - baseline by construction.\n\n")
print(as.data.frame(share %>% arrange(mc, roi, hpf)), digits = 3, row.names = FALSE)
write_outcome(share, file.path(OUT, "baseline_share_by_setting.csv"))
cat("Baseline share, PFC comfort: 0.68 frozen (quoted from the frozen\n",
    "campaign) vs",
    round(share$baseline_share[share$mc == "none" & abs(share$hpf - 0.01) < 1e-12 &
                               share$roi == "PFC_Frontal"], 3),
    "re-filtered at hpf = 0.01 (no MC)\n")

# ------------------------------------------- 5. the cycle-position contrast
rule("5. CYCLE POSITION — run-openers vs followers (the frozen p = 3.4e-6)")

# Events 1, 5 and 9 open a run and have their baseline in genuine rest; the
# other nine follow a clip by 60 s. If the high-pass creates the offset, the
# contrast should shrink as the corner is lowered.
blk <- load_intermediate("refilter_event_metrics_roi.csv") %>%
  filter(roi == "PFC_Frontal") %>%
  rename(diff = response_minus_baseline_m5_30, base = baseline_m10_0) %>%
  mutate(run_opener = as.integer(event_index %in% c(1, 5, 9)))

cyc <- blk %>% group_by(motion_correction, hpf) %>%
  group_modify(~ {
    m <- lmerTest::lmer(diff ~ run_opener + (1 | subject_id), data = .x)
    ft <- fixef_table(m, keep = "run_opener")
    mb <- lmerTest::lmer(base ~ run_opener + (1 | subject_id), data = .x)
    fb <- fixef_table(mb, keep = "run_opener")
    tibble(diff_contrast = ft$estimate, diff_p = ft$p,
           base_contrast = fb$estimate, base_p = fb$p,
           mean_diff_opener = mean(.x$diff[.x$run_opener == 1]),
           mean_diff_follower = mean(.x$diff[.x$run_opener == 0]))
  }) %>% ungroup()
print(as.data.frame(cyc), digits = 3, row.names = FALSE)
write_outcome(cyc, file.path(OUT, "cycle_position_by_setting.csv"))
cat("\nSI Table S3 cycle-position row: +0.065 (p = 3.4e-6) frozen vs",
    round(cyc$diff_contrast[cyc$motion_correction == "none" & abs(cyc$hpf - 0.01) < 1e-12], 3),
    "(p =", signif(cyc$diff_p[cyc$motion_correction == "none" & abs(cyc$hpf - 0.01) < 1e-12], 2),
    ") rebuild at hpf = 0.01\n")

# ------------------------------------------------------ 6. channel 10
rule("6. CHANNEL 10 — the campaign's one local candidate, at every setting")

ch10 <- load_intermediate("refilter_channel10_metrics.csv") %>%
  inner_join(design, by = c("subject_id", "marker")) %>%
  group_by(motion_correction, hpf) %>%
  group_modify(~ {
    mr <- p24_lmm(.x, paste(AXES, collapse = " + "), "response_m5_30")
    md <- p24_lmm(.x, paste(AXES, collapse = " + "), "response_minus_baseline_m5_30")
    fr <- fixef_table(mr, keep = "comfort_wi"); fd <- fixef_table(md, keep = "comfort_wi")
    tibble(n_cells = nrow(.x),
           response_est = fr$estimate, response_p = fr$p,
           difference_est = fd$estimate, difference_p = fd$p)
  }) %>% ungroup()
print(as.data.frame(ch10), digits = 3, row.names = FALSE)
write_outcome(ch10, file.path(OUT, "channel10_by_setting.csv"))
cn <- ch10 %>% filter(motion_correction == "none") %>% arrange(desc(hpf))
cat("\nSI Table S2 row A19 (no-MC arm, response window, corners 0.01/0.005/0.002/none):",
    paste(round(cn$response_est, 3), collapse = " / "), "uM; p =",
    paste(signif(cn$response_p, 2), collapse = " / "),
    "\n  and SI Table S3 row 5: -0.051 (p = 7.8e-4) frozen vs",
    round(cn$response_est[abs(cn$hpf - 0.01) < 1e-12], 3),
    "(p =", signif(cn$response_p[abs(cn$hpf - 0.01) < 1e-12], 2), ") rebuild\n")

# ---------------------------------------------- 7. whole-head, for completeness
rule("7. WHOLE-HEAD COMFORT SLOPE BY SETTING (declared target of this analysis)")

head_tab <- cond %>% filter(roi == "head") %>%
  inner_join(design, by = c("subject_id", "marker")) %>%
  group_by(motion_correction, hpf) %>%
  group_modify(~ {
    map_dfr(c("response_m5_30", "baseline_m10_0", "response_minus_baseline_m5_30"),
            function(y) {
              m <- p24_lmm(.x, paste(AXES, collapse = " + "), y)
              fixef_table(m, keep = "comfort_wi") %>%
                transmute(window = y, estimate, se, p)
            })
  }) %>% ungroup()
print(as.data.frame(head_tab), digits = 3, row.names = FALSE)
write_outcome(head_tab, file.path(OUT, "whole_head_by_setting.csv"))

rule("DONE")
close_log(con)
