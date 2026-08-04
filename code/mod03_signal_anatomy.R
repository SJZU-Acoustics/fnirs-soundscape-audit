# =============================================================================
# mod03_signal_anatomy.R — the cycle-position anatomy of the negative offset.
#
# Main-text Results (Section "No response-window appraisal association
# survived the declared analyses"): run-opening clips differ from followers by
# -0.023 against -0.088 uM (p = 3.4e-6). Recomputed here from the workbook's
# block_level sheet (per-presentation HbO head difference score) and checked
# against the Fig 1a anchor shipped in data/intermediates/.
# =============================================================================

source("code/load_data.R")
OUT <- "output/analysis_03_signal_anatomy"

b <- load_p24_block() %>%
  mutate(cycle_position = if_else(event_index %in% c(1, 5, 9),
                                  "opens a run", "follows a 90-s cycle"),
         opens = as.integer(cycle_position == "opens a run"))

cyc_cmp <- b %>% group_by(cycle_position) %>%
  summarise(n = n(), difference = mean(HbO_head),
            se_difference = sd(HbO_head) / sqrt(n()), .groups = "drop")
print(as.data.frame(cyc_cmp), digits = 4)
write_outcome(cyc_cmp, file.path(OUT, "cycle_position_comparison.csv"))

m_pos <- p24_lmm(b, "opens", "HbO_head")
cat("\nHbO_head (difference score) ~ opens_a_run + (1|subject):\n")
pos_tab <- fixef_table(m_pos)
print(as.data.frame(pos_tab), digits = 3)
write_outcome(pos_tab, file.path(OUT, "cycle_position_lmm.csv"))

# Cross-check against the Fig 1a anchor (raw-recording-era recompute).
anchor <- load_intermediate("fig1a_difference_score_check.csv")
chk <- cyc_cmp %>%
  left_join(anchor %>% select(cycle_position, a03_mean), by = "cycle_position") %>%
  mutate(abs_diff = abs(difference - a03_mean))
cat("\nFig 1a difference-score anchor check:\n")
print(as.data.frame(chk), digits = 6)
stopifnot(all(chk$abs_diff < 1e-9))

# --- Window decomposition: how much of the difference-score slope is baseline?
# Main text: "the baseline carried 68% of the comfort slope of the difference
# score". On the raw uM scale the decomposition is exact
# (response-minus-baseline = response - baseline), so each slope of rmb splits
# into the response-window slope minus the pre-stimulus baseline slope.
bm <- load_block_metrics(active_only = TRUE)
d  <- load_p24()
head_blocks <- bm %>%
  filter(chromophore == "HbO") %>%
  group_by(subject_id, event_index, marker) %>%
  summarise(rmb  = mean(response_minus_baseline_m5_30),
            raw  = mean(response_m5_30),
            base = mean(baseline_m10_0), .groups = "drop") %>%
  group_by(subject_id, marker) %>%
  summarise(across(c(rmb, raw, base), mean), .groups = "drop") %>%
  left_join(d %>% select(subject_id, marker, comfort_wi, richness_wi),
            by = c("subject_id", "marker"))

win_raw <- map_dfr(c("rmb", "raw", "base"), function(y) {
  m <- p24_lmm(head_blocks, "comfort_wi + richness_wi", y)
  fixef_table(m, keep = c("comfort_wi", "richness_wi")) %>% mutate(dv = y, .before = 1)
})
write_outcome(win_raw, file.path(OUT, "window_sensitivity_raw_uM.csv"))

est <- win_raw %>% select(dv, term, estimate) %>%
  pivot_wider(names_from = dv, values_from = estimate)
cat("\nShares of each difference-score slope (baseline enters with a minus sign):\n")
print(as.data.frame(est %>% mutate(
  share_response = raw / rmb,
  share_baseline = -base / rmb)), digits = 4)
cat("mod03 done.\n")
