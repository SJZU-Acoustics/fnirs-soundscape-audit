# mod21_correction_sweep.R — completeness sweep: does anything ELSE in the
# campaign depend on the compromised channels?
#
# Port of exploration/analysis_21_correction_sweep for the release repo. Four
# groups, as declared:
#   1. head-mean / global-regressor results — predicted to be affine-invariant
#      (coefficient rescales, t and p unchanged). Tested, not assumed.
#   2. families whose BH set spans the temporal ROIs (A14) — prefrontal
#      q-values can move even though prefrontal data do not.
#   3. nominal-only temporal seeds quoted in the synthesis (A10).
#   4. channel-level families — dead channels are structurally VOID tests.
#      Reported as a declared sensitivity that MAY NEVER CARRY A CLAIM; the
#      campaign's MC policy forbids shrinking a family to improve a q.
# plus a fifth found on review: the gate-C agreement statistic of the re-filter
# analysis (mod19) depends on the right-temporal ROI definition.
#
# MASK. A20 and A21 were DECLARED against the session-CV dropout determination
# of 2026-07-30 (raw-intensity CV < 0.002 per subject x channel, shipped as
# data/intermediates/qa_detector_dropout_0730.csv) plus detector 8, which is
# why that file ships. That mask — not the union mask of load_channel_mask() —
# is used here, so the declared numbers reproduce exactly.
#
# A15 INPUTS. The channel-family sensitivity re-runs the BH correction on A15's
# adjusted channel family (21 channels x 2 axes). Those per-channel fits derive
# from the deposited condition layer, so they are recomputed here from the
# workbook rather than shipped — the CR2 arm requires the clubSandwich package.
#
# Displays fed: SI Table S2 row A21 (channel 10 in the 34-test live-channel
# family: estimate -0.035 uM, p = 0.003, q = 0.0925); Table 3 row 3 (gate-C
# agreement r = 0.337 across the pooled 12 cells, 0.235 on the valid arm
# alone, 0.643 on the corrected ROI); the A10 right-temporal sex seed
# re-run (+0.039 uM, p = 0.023 as published -> +0.027 uM, p = 0.034 masked).

source("code/load_data.R")
suppressPackageStartupMessages(library(clubSandwich))

OUT <- file.path("output", "analysis_21_correction_sweep")
con <- open_log(file.path(OUT, "correction_sweep_report.txt"))
rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")
AXES <- c("comfort_wi", "richness_wi")

cat("mod21 — completeness sweep after the detector-dropout correction.\n")
cat("All concentrations in uM.\n")

dead      <- load_dead_channels()
mask_keys <- dead %>% distinct(subject_id, channel_index) %>% mutate(.dead = TRUE)
det8      <- load_channels() %>% filter(detector == 8L) %>% distinct(subject_id, channel_index)
drop_all  <- bind_rows(mask_keys %>% select(subject_id, channel_index), det8) %>%
  distinct(subject_id, channel_index) %>% mutate(.drop = TRUE)

bm <- load_block_metrics(active_only = TRUE) %>% filter(chromophore == "HbO")
map <- load_channels() %>% filter(active_pair == 1) %>%
  distinct(subject_id, channel_index, roi)

# ============================================ 1. is the head mean really affine?
rule("1. THE HEAD MEAN — is the affine prediction true?")
cat("A dark channel contributes a near-constant to a mean, so head_frozen should\n")
cat("be an affine transform of head_masked and every slope test on it should keep\n")
cat("its t and p while the coefficient rescales. Predicted before running.\n\n")

head_cells <- function(exclude) {
  b <- bm
  if (nrow(exclude)) b <- b %>% left_join(exclude, by = c("subject_id", "channel_index")) %>%
    filter(is.na(.drop) | !.drop) %>% select(-any_of(".drop"))
  b %>% group_by(subject_id, event_index, marker) %>%
    summarise(difference = mean(response_minus_baseline_m5_30),
              response = mean(response_m5_30), baseline = mean(baseline_m10_0),
              .groups = "drop") %>%
    group_by(subject_id, marker) %>%
    summarise(across(c(difference, response, baseline), mean), .groups = "drop")
}
# THREE head means, because the prediction is only about the dark channels:
#   f = frozen (everything)
#   m = mask only  -> removes near-constants ONLY. This is the affine test.
#   d = defensible -> also removes detector 8 where it was LIVE, which is a real
#       signal change, not a constant, so the prediction does not apply to it.
h_frozen <- head_cells(tibble(subject_id = character(), channel_index = integer(),
                              .drop = logical()))
h_mask   <- head_cells(mask_keys %>% transmute(subject_id, channel_index, .drop = TRUE))
h_def    <- head_cells(drop_all)
hd <- h_frozen %>% rename(diff_f = difference, resp_f = response, base_f = baseline) %>%
  left_join(h_mask %>% rename(diff_m = difference, resp_m = response, base_m = baseline),
            by = c("subject_id", "marker")) %>%
  left_join(h_def %>% rename(diff_d = difference, resp_d = response, base_d = baseline),
            by = c("subject_id", "marker")) %>%
  left_join(load_p24() %>% select(subject_id, marker, all_of(AXES)),
            by = c("subject_id", "marker"))

lin <- map_dfr(c("diff", "resp", "base"), function(w) {
  map_dfr(c(m = "mask only (affine test)", d = "also dropping detector 8"), function(lbl) {
    a <- hd[[paste0(w, "_f")]]; b <- hd[[paste0(w, "_", names(which(c(m = "mask only (affine test)", d = "also dropping detector 8") == lbl)))]]
    keep <- !is.na(a) & !is.na(b)
    tibble(window = w, versus = lbl, n = sum(keep), r = cor(a[keep], b[keep]),
           slope = coef(lm(a[keep] ~ b[keep]))[2],
           resid_sd_uM = sd(residuals(lm(a[keep] ~ b[keep]))))
  })
})
cat("head_frozen regressed on each alternative:\n")
print(as.data.frame(lin), digits = 4)
write_outcome(lin, file.path(OUT, "head_mean_affine_test.csv"))

cat("\nand the consequence for inference — the same comfort/richness models on both:\n")
VNAME <- c(f = "frozen head", m = "mask-only head", d = "defensible head (no det 8)")
hm <- map_dfr(c("diff", "resp", "base"), function(w) {
  map_dfr(c("f", "m", "d"), function(v) {
    y <- paste0(w, "_", v)
    fixef_table(p24_lmm(hd, paste(AXES, collapse = " + "), y), keep = AXES) %>%
      mutate(window = w, version = unname(VNAME[v]), .before = 1)
  })
}) %>% select(window, version, term, estimate, se, t, p)
print(as.data.frame(hm), digits = 4)
write_outcome(hm, file.path(OUT, "head_mean_models.csv"))

worst_t <- hm %>% select(window, version, term, t) %>%
  pivot_wider(names_from = version, values_from = t) %>%
  mutate(dt_mask = abs(`frozen head` - `mask-only head`),
         dt_def  = abs(`frozen head` - `defensible head (no det 8)`))
cat("\nchange in t against the frozen head:\n")
print(as.data.frame(worst_t %>% select(window, term, dt_mask, dt_def)), digits = 3)
cat("\nmask only (the affine prediction): max |dt| =", signif(max(worst_t$dt_mask), 3), "\n")
cat(if (max(worst_t$dt_mask) < 0.05)
  "  PREDICTION HELD — removing the dark channels rescales the head mean and\n  leaves every t and p where they were.\n" else
  "  PREDICTION FAILED — removing the dark channels alone moves inference.\n")
cat("dropping detector 8 as well: max |dt| =", signif(max(worst_t$dt_def), 3),
    "— this is a REAL signal change, not a\n  constant, so the prediction never applied to it; it is reported as the size\n  of the correction, not as a failure.\n")

# ================================ 2. A14's family of 6 spans the temporal ROIs
rule("2. A14's FAMILY OF 6 — its BH set spans the temporal ROIs")
cat("A14 fitted comfort slopes on sequence-detrended series, 2 windows x 3 ROIs,\n")
cat("BH across 6 (published min q = 0.061, no survivors). Two temporal p-values\n")
cat("enter that correction, so the prefrontal q can move even though prefrontal\n")
cat("data cannot.\n\n")

roi_blocks <- bm %>%
  left_join(drop_all, by = c("subject_id", "channel_index")) %>%
  filter(is.na(.drop)) %>%
  left_join(map, by = c("subject_id", "channel_index")) %>%
  filter(!is.na(roi)) %>%
  group_by(subject_id, sequence, event_index, marker, roi) %>%
  summarise(baseline = mean(baseline_m10_0), response = mean(response_m5_30),
            .groups = "drop") %>%
  group_by(subject_id, sequence, roi) %>%                    # A14's detrending
  mutate(across(c(baseline, response), ~ .x - mean(.x))) %>%
  ungroup() %>%
  left_join(load_p24() %>% select(subject_id, marker, all_of(AXES)),
            by = c("subject_id", "marker"))

a14 <- map_dfr(ROIS, function(rr) {
  map_dfr(c("baseline", "response"), function(w) {
    dd <- roi_blocks %>% filter(roi == rr)
    fixef_table(p24_lmm(dd, paste(AXES, collapse = " + "), w), keep = "comfort_wi") %>%
      mutate(roi = rr, window = w, .before = 1)
  })
}) %>% mutate(q_BH = bh(p), family_size = n())
print(as.data.frame(a14 %>% select(window, roi, estimate, se, p, q_BH)), digits = 3)
write_outcome(a14, file.path(OUT, "a14_family_masked.csv"))
cat("\nA14 as published: min q = 0.061, 0 survivors | under the mask: min q =",
    signif(min(a14$q_BH), 3), ", survivors:", sum(a14$q_BH < 0.05), "\n")

# ============================ 3. A10's right-temporal seed quoted in the synthesis
rule("3. A10's RIGHT-TEMPORAL SEED — sex difference in reactivity")
cat("A10 reported females more right-temporally reactive (+0.039 uM, p = 0.023),\n")
cat("nominal only, and the synthesis carries it as a seed observation. It is a\n")
cat("right-temporal quantity, so it is re-run on the corrected ROI.\n\n")

resp_cells <- bm %>%
  left_join(drop_all, by = c("subject_id", "channel_index")) %>%
  filter(is.na(.drop)) %>%
  left_join(map, by = c("subject_id", "channel_index")) %>%
  filter(!is.na(roi)) %>%
  group_by(subject_id, marker, roi) %>%
  summarise(response = mean(response_m5_30), .groups = "drop") %>%
  left_join(load_p24() %>% select(subject_id, marker, sex, age), by = c("subject_id", "marker"))

a10 <- map_dfr(ROIS, function(rr) {
  dd <- resp_cells %>% filter(roi == rr, !is.na(sex))
  fixef_table(p24_lmm(dd, "sex", "response")) %>%
    mutate(roi = rr, n_subjects = n_distinct(dd$subject_id), .before = 1)
}) %>% filter(!grepl("Intercept", term))
print(as.data.frame(a10 %>% select(roi, n_subjects, term, estimate, se, p)), digits = 3)
write_outcome(a10, file.path(OUT, "a10_sex_reactivity_masked.csv"))
cat("\nRight-temporal seed under the mask: +",
    round(a10$estimate[a10$roi == "Right_Temporal"], 4), " uM, p =",
    signif(a10$p[a10$roi == "Right_Temporal"], 3),
    "(quoted in the main text as +0.039, p = 0.023 -> corrected here)\n")

# ================== 4. channel families: the structurally void tests (SENSITIVITY)
rule("4. CHANNEL FAMILIES — the structurally void tests (DECLARED SENSITIVITY)")
cat("Dead channels cannot carry a hypothesis, so tests on them are structurally\n")
cat("void and inflate the BH denominator. Removing them LOWERS every other\n")
cat("channel's q. THIS MAY NOT CARRY A CLAIM — the MC policy forbids shrinking a\n")
cat("family to improve a q, and the decision to keep 18/19 was correct at the\n")
cat("time and is not revised. The purpose is to bound how much of channel 10's q\n")
cat("is family composition.\n\n")

# The campaign quotes ch10 at q = 0.114, which is the CR2 family; the
# model-based family holds q = 0.074. Both are recomputed from the workbook
# (A15's leave-that-channel-out global adjustment, response window), and the
# CR2 row is the one that matters. CR2 is applied to the SAME random-intercept
# fit (A08's method): it recomputes the standard error under arbitrary
# within-subject covariance and leaves the point estimate untouched.
decomp2 <- function(d, v) {
  d %>% group_by(subject_id) %>%
    mutate(!!paste0(v, "_bw") := mean(.data[[v]], na.rm = TRUE)) %>% ungroup() %>%
    mutate(!!paste0(v, "_wi") := .data[[v]] - .data[[paste0(v, "_bw")]])
}
cd <- load_condition_analysis("HbO") %>%
  left_join(map, by = c("subject_id", "channel_index")) %>% filter(!is.na(roi))
YCOL_RESPONSE <- "mean_response_m5_30"
chans <- sort(unique(cd$channel_index))
ch_rows <- list(); cr2_rows <- list()
for (ci in chans) {
  cell <- cd %>% group_by(subject_id, marker) %>%
    summarise(y = mean(.data[[YCOL_RESPONSE]][channel_index == ci]),
              g = mean(.data[[YCOL_RESPONSE]][channel_index != ci]),
              has = any(channel_index == ci), .groups = "drop") %>%
    filter(has, is.finite(y), is.finite(g)) %>%
    inner_join(.analysis_design("full"), by = c("subject_id", "marker")) %>%
    decomp2("g")
  m_adj <- p24_lmm(cell, paste(c(AXES, "g_wi", "g_bw"), collapse = " + "), "y")
  fa <- fixef_table(m_adj, keep = AXES)
  for (ax in AXES) ch_rows[[length(ch_rows) + 1]] <- tibble(
    channel_index = ci, axis = ax,
    est_adjusted = fa$estimate[fa$term == ax], p_adjusted = fa$p[fa$term == ax])
  ct <- as.data.frame(clubSandwich::coef_test(
    m_adj, vcov = "CR2", cluster = cell$subject_id, test = "Satterthwaite"))
  ct$term <- rownames(ct)
  for (ax in AXES) {
    i <- which(ct$term == ax)
    cr2_rows[[length(cr2_rows) + 1]] <- tibble(
      channel_index = ci, axis = ax, est = ct$beta[i], p_cr2 = ct$p_Satt[i])
  }
}
a15    <- bind_rows(ch_rows)  %>% mutate(q_adjusted = bh(p_adjusted))
a15cr2 <- bind_rows(cr2_rows) %>% mutate(q_cr2 = bh(p_cr2))

void <- c(18L, 19L, 20L, 21L)
recompute <- function(d, keep_ch, label) {
  dd <- d %>% filter(channel_index %in% keep_ch)
  dd %>% mutate(q_new = bh(p_adjusted)) %>%
    filter(channel_index == 10, grepl("comfort", axis)) %>%
    transmute(family = label, family_size = nrow(dd), channel = 10,
              estimate = est_adjusted, p = p_adjusted,
              q_as_published = q_adjusted, q_recomputed = q_new)
}
recompute_cr2 <- function(d, keep_ch, label) {
  dd <- d %>% filter(channel_index %in% keep_ch)
  dd %>% mutate(q_new = bh(p_cr2)) %>%
    filter(channel_index == 10, grepl("comfort", axis)) %>%
    transmute(family = label, family_size = nrow(dd), channel = 10,
              estimate = est, p = p_cr2, q_as_published = q_cr2, q_recomputed = q_new)
}
all_ch <- sort(unique(a15$channel_index))
sens <- bind_rows(
  recompute_cr2(a15cr2, all_ch, "CR2, as published: 21 channels x 2 axes = 42"),
  recompute_cr2(a15cr2, setdiff(all_ch, c(18L, 19L)), "CR2, minus det-9 channels = 38"),
  recompute_cr2(a15cr2, setdiff(all_ch, void), "CR2, minus det-9 and det-8 = 34"),
  recompute(a15, all_ch, "model-based, as published = 42"),
  recompute(a15, setdiff(all_ch, void), "model-based, minus both detectors = 34"))
print(as.data.frame(sens), digits = 3)
write_outcome(sens, file.path(OUT, "channel_family_sensitivity.csv"))
cat("\nRead this as a bound, not a result: channel 10's q moves only because the\n")
cat("denominator shrinks; its estimate and p are identical in all three rows.\n")
cat("SI Table S2 row A21 (34 live channels): estimate",
    round(sens$estimate[sens$family_size == 34 & grepl("CR2", sens$family)], 3),
    "uM, p =", signif(sens$p[sens$family_size == 34 & grepl("CR2", sens$family)], 2),
    ", q =", signif(sens$q_recomputed[sens$family_size == 34 & grepl("CR2", sens$family)], 3), "\n")

cat("\nA16's ridge: dead channels are constant predictors within subject, so after\n")
cat("centring they carry zero variance and cannot receive weight. Verified:\n")
ch_var <- bm %>%
  group_by(subject_id, channel_index) %>%
  summarise(v = var(response_minus_baseline_m5_30), .groups = "drop") %>%
  left_join(mask_keys, by = c("subject_id", "channel_index")) %>%
  mutate(state = if_else(is.na(.dead), "live", "dark")) %>%
  group_by(state) %>%
  summarise(n_cells = n(), median_var_uM2 = median(v, na.rm = TRUE),
            max_var_uM2 = max(v, na.rm = TRUE), .groups = "drop")
print(as.data.frame(ch_var), digits = 3)
write_outcome(ch_var, file.path(OUT, "a16_dead_channel_variance.csv"))

cat("\nA18's reliability statistics are between-subject CORRELATIONS of time\n")
cat("courses. A constant channel in a mean is an affine transform, and\n")
cat("correlation is affine-invariant, so those results are unchanged by\n")
cat("construction — no re-run needed, and section 1 above is the empirical check\n")
cat("of the same argument.\n")

# ============ 5. mod19's gate-C reproducibility statistic depends on them too
rule("5. MOD19's GATE-C AGREEMENT — a fifth dependence, found on review")
cat("The gate-C agreement pools BOTH motion-correction arms, including the\n")
cat("wavelet arm that fails gate C, and its right-temporal cells were computed on\n")
cat("the ROI that averages detector 8. Decomposed here: pooled 12 cells vs the\n")
cat("valid arm alone, as published vs the corrected ROI definition.\n\n")

design <- load_p24() %>% select(subject_id, marker, all_of(AXES))

# Pooled (one-step) subject x condition ROI means of the difference score.
roi_est <- function(dat, use_mask) {
  d <- dat
  if (use_mask) d <- d %>% left_join(drop_all, by = c("subject_id", "channel_index")) %>%
    filter(is.na(.drop)) %>% select(-.drop)
  cells <- d %>% left_join(map, by = c("subject_id", "channel_index")) %>%
    filter(!is.na(roi)) %>%
    group_by(subject_id, marker, roi) %>%
    summarise(y = mean(response_minus_baseline_m5_30), .groups = "drop") %>%
    left_join(design, by = c("subject_id", "marker"))
  bind_rows(lapply(ROIS, function(rr) {
    tb <- fixef_table(p24_lmm(cells %>% filter(roi == rr),
                              paste(AXES, collapse = " + "), "y"), keep = AXES)
    tb$roi <- rr
    tb
  })) %>% select(roi, term, estimate, se, p)
}

# The re-filtered side at hpf = 0.01, from the event-level intermediate: the
# pooled mean over a cell's channel x presentation rows is the count-weighted
# mean of the per-event ROI means (masked variant uses the masked columns).
ev <- load_intermediate("refilter_event_metrics_roi.csv") %>%
  filter(abs(hpf - 0.01) < 1e-9)
rb_cells <- function(arm, msk) {
  d <- ev %>% filter(motion_correction == arm)
  if (msk) {
    d <- d %>% filter(!is.na(diff_xdead)) %>%
      group_by(subject_id, marker, roi) %>%
      summarise(y = weighted.mean(diff_xdead, n_live), .groups = "drop")
  } else {
    d <- d %>% group_by(subject_id, marker, roi) %>%
      summarise(y = weighted.mean(response_minus_baseline_m5_30, n_channels),
                .groups = "drop")
  }
  cells <- d %>% left_join(design, by = c("subject_id", "marker"))
  bind_rows(lapply(ROIS, function(rr) {
    tb <- fixef_table(p24_lmm(cells %>% filter(roi == rr),
                              paste(AXES, collapse = " + "), "y"), keep = AXES)
    tb$roi <- rr
    tb
  })) %>% select(roi, term, estimate, se, p)
}

gate_c <- map_dfr(c(FALSE, TRUE), function(msk) {
  frz <- roi_est(bm, msk)
  both <- map_dfr(c("none", "wavelet"), function(a) {
    rb_cells(a, msk) %>% mutate(mc = a)
  }) %>% inner_join(frz, by = c("roi", "term"), suffix = c("_ref", "_frz"))
  tibble(roi_definition = if (msk) "corrected" else "as published",
         r_both_arms_12 = cor(both$estimate_ref, both$estimate_frz),
         r_valid_arm_6 = with(filter(both, mc == "none"), cor(estimate_ref, estimate_frz)),
         r_wavelet_arm_6 = with(filter(both, mc == "wavelet"), cor(estimate_ref, estimate_frz)))
})
print(as.data.frame(gate_c), digits = 3)
write_outcome(gate_c, file.path(OUT, "a19_gate_c_decomposed.csv"))
cat("\nTable 3 row 3: published 0.337 (both arms) =",
    round(gate_c$r_valid_arm_6[1], 3), "valid arm alone; on the corrected ROI",
    "the valid arm is",
    round(gate_c$r_valid_arm_6[gate_c$roi_definition == "corrected"], 3), "\n")
cat("Two corrections run in opposite directions: pooling the rejected wavelet\n")
cat("arm flattered the published number, and a substantial part of the published\n")
cat("limit was the DETECTOR, not the rebuild. Neither removes the limitation:\n")
cat("0.643 over six cells is moderate and the rebuild remains a controlled\n")
cat("stand-in for Homer3 rather than a reproduction of it.\n")

rule("VERDICT")
cat("1. head mean: max |dt| =", signif(max(worst_t$dt_mask), 3), "(mask only),",
    signif(max(worst_t$dt_def), 3), "(also dropping det 8).",
    "\n   The head mean is a non-claim-carrying summary under MC policy rule 2, so\n",
    "  nothing rests on it; the numbers quoted for it must nonetheless be restated.\n")
cat("2. A14 family: min q 0.061 (published) ->", signif(min(a14$q_BH), 3), "\n")
cat("3. A10 sex seed, right temporal: p =",
    signif(a10$p[a10$roi == "Right_Temporal"], 3), "\n")
cat("4. ch10 q under the void-test sensitivity:",
    paste(signif(sens$q_recomputed, 3), collapse = " -> "), "(bound only)\n")
cat("5. gate C: published 0.337 (both arms) =",
    signif(gate_c$r_valid_arm_6[1], 3), "valid arm alone;",
    "\n   on the corrected ROI the valid arm is",
    signif(gate_c$r_valid_arm_6[gate_c$roi_definition == "corrected"], 3), "\n")

close_log(con)
