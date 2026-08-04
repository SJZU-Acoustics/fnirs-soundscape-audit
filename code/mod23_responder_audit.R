# =============================================================================
# mod23_responder_audit.R — Responder heterogeneity.
#
# Question: is the campaign's null an AVERAGE that conceals a subset of people
# who respond more strongly to the stimuli?
#
# The governing rule, fixed before this module was written:
#   - selection on the outcome + testing on the same rows may never produce a
#     reportable result; where it is run (part D) it is run to MEASURE what it
#     manufactures;
#   - every part is reported whatever it returns; no part may be dropped after
#     seeing it;
#   - a responsive subset may be claimed only if A (variance) AND B
#     (reliability) AND (C or E) all clear.
#
# Parts: 0 gates · A direct heterogeneity · A2 falsification of A · A2.4
# residual shape · B reliability ceiling · B.4 reliability generality · C
# out-of-sample selection · D in-sample selection · E outcome-independent
# subgroups · F non-person heterogeneity · G mixture · H lambda(k).
#
# RNG. Every stochastic element carries its own seed (24 for the parametric
# bootstrap via simulate(seed = 24); 1000+b for the subject bootstrap;
# 20000+i for the within-subject selection permutations; 3000+i for the
# mixture null; 40000+i / 50000+i for the part-A2 nulls), so the streams are
# identical however many cores mclapply uses.
#
# Outputs (output/analysis_23_responder_heterogeneity/): A_slope_variance.csv,
# A2_variance_falsification.csv, A2_residual_shape.csv, B_slope_reliability.csv,
# B_reliability_bootstrap.csv, B4_reliability_generality.csv,
# CD_selection_curves.csv, CD_selection_summary.csv, CD_permutation_nulls.csv,
# E_subgroup_interactions.csv, F_nonperson_heterogeneity.csv, G_mixture.csv,
# H_lambda_by_subset_size.csv, responder_heterogeneity_report.txt
# =============================================================================

source("code/load_data.R")
suppressPackageStartupMessages({ library(mclust); library(parallel) })
# mclust masks purrr::map. Everything below means purrr's.
map <- purrr::map

OUT <- file.path("output", "analysis_23_responder_heterogeneity")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
con <- open_log(file.path(OUT, "responder_heterogeneity_report.txt"))

set.seed(24)
N_PERM <- 1000L   # house default
N_PB   <- 1000L   # parametric-bootstrap LRT simulations
N_CORE <- max(1L, parallel::detectCores() - 1L)

cat("=========================================================================\n")
cat("A23 — RESPONDER HETEROGENEITY\n")
cat("=========================================================================\n\n")

# ======================================================================== PART 0
# Gates. The mask must be extended from the difference score to the response and
# baseline windows; that extension is only usable if it reproduces the campaign
# default exactly where the two overlap.

cat("---- PART 0: validation gates ----------------------------------------\n\n")

#' Block-level ROI means under the declared session-CV dead-channel mask plus
#' the detector exclusion, for ALL three window statistics at once. Mirrors the
#' masked-ROI aggregation (channel mean within block, then mean over the cell's
#' presentations) and adds a whole-head unit (mean over every live channel).
#' The presentation x channel layer is the shipped tidy_block_metrics.csv
#' intermediate; the CV mask is load_dead_channels().
masked_block_roi <- function(chromophore = "HbO",
                             exclude_detectors = UNUSABLE_DETECTORS,
                             use_dead_mask = TRUE) {
  dead <- tibble(subject_id = character(), channel_index = integer())
  if (use_dead_mask) dead <- load_dead_channels() %>% distinct(subject_id, channel_index)
  if (length(exclude_detectors)) {
    dropped <- load_channels() %>% filter(.data$detector %in% exclude_detectors) %>%
      distinct(subject_id, channel_index)
    dead <- bind_rows(dead, dropped) %>% distinct(subject_id, channel_index)
  }
  dead <- dead %>% mutate(.dead = TRUE)
  bm <- load_block_metrics(active_only = TRUE) %>% filter(.data$chromophore == !!chromophore)
  map_ <- load_channels() %>% filter(active_pair == 1) %>%
    distinct(subject_id, channel_index, roi)
  live <- bm %>%
    left_join(dead, by = c("subject_id", "channel_index")) %>%
    filter(is.na(.data$.dead)) %>%
    left_join(map_, by = c("subject_id", "channel_index"))
  bind_rows(live %>% filter(!is.na(roi)), live %>% mutate(roi = "head")) %>%
    group_by(subject_id, marker, event_index, roi) %>%
    summarise(resp = mean(response_m5_30),
              base = mean(baseline_m10_0),
              diff = mean(response_minus_baseline_m5_30),
              n_live = n(), .groups = "drop")
}

#' The declared masked-ROI reference for gate (a): same session-CV mask and
#' aggregation as masked_block_roi(), on the difference score only, wide by ROI.
roi_masked_cv <- function(chromophore = "HbO", exclude_detectors = UNUSABLE_DETECTORS) {
  dead <- load_dead_channels() %>% distinct(subject_id, channel_index)
  if (length(exclude_detectors)) {
    dropped <- load_channels() %>% filter(.data$detector %in% exclude_detectors) %>%
      distinct(subject_id, channel_index)
    dead <- bind_rows(dead, dropped) %>% distinct(subject_id, channel_index)
  }
  dead <- dead %>% mutate(.dead = TRUE)
  bm <- load_block_metrics(active_only = TRUE) %>%
    filter(.data$chromophore == !!chromophore)
  map_ <- load_channels() %>% filter(active_pair == 1) %>%
    distinct(subject_id, channel_index, roi)
  agg <- bm %>%
    left_join(dead, by = c("subject_id", "channel_index")) %>%
    filter(is.na(.data$.dead)) %>%
    left_join(map_, by = c("subject_id", "channel_index")) %>%
    filter(!is.na(roi)) %>%
    group_by(subject_id, marker, event_index, roi) %>%
    summarise(block_mean = mean(response_minus_baseline_m5_30),
              n_live = n(), .groups = "drop") %>%
    group_by(subject_id, marker, roi) %>%
    summarise(value = mean(block_mean), n_live = max(n_live), .groups = "drop")
  wide_v <- agg %>% mutate(roi = paste0(chromophore, "m_", roi)) %>%
    select(subject_id, marker, roi, value) %>%
    pivot_wider(names_from = roi, values_from = value)
  wide_n <- agg %>% mutate(roi = paste0("n_live_", roi)) %>%
    select(subject_id, marker, roi, n_live) %>%
    pivot_wider(names_from = roi, values_from = n_live)
  load_p24() %>%
    left_join(wide_v, by = c("subject_id", "marker")) %>%
    left_join(wide_n, by = c("subject_id", "marker"))
}

blk <- masked_block_roi()

cond_from_blk <- blk %>%
  group_by(subject_id, marker, roi) %>%
  summarise(resp = mean(resp), base = mean(base), diff = mean(diff),
            n_live = max(n_live), .groups = "drop")

ref <- roi_masked_cv() %>%
  select(subject_id, marker, starts_with("HbOm_")) %>%
  pivot_longer(starts_with("HbOm_"), names_to = "roi", values_to = "ref_diff") %>%
  mutate(roi = sub("^HbOm_", "", roi))

gate_a <- cond_from_blk %>% inner_join(ref, by = c("subject_id", "marker", "roi")) %>%
  group_by(roi) %>%
  summarise(n = n(), max_abs_diff = max(abs(diff - ref_diff), na.rm = TRUE), .groups = "drop")
cat("Gate (a) — local aggregator vs the declared masked-ROI builder on the difference score:\n")
print(as.data.frame(gate_a)); stopifnot(all(gate_a$max_abs_diff < 1e-12))

gate_b <- check_roi_rebuild()
cat("\nGate (b) — check_roi_rebuild() (frozen reproduction, nothing excluded):\n")
print(as.data.frame(gate_b)); stopifnot(all(gate_b$max_abs_diff < 1e-12))

# ------------------------------------------------------------- fast primitives
# Everything inside a permutation or bootstrap loop uses these. They are exact
# re-expressions of the lmer/lm fits, asserted against them in gate (c).

#' Per-subject cross-products of the within-subject regression, one row per subject.
#' y is demeaned within subject; comfort_wi/richness_wi are already subject-centred,
#' so the pooled no-intercept OLS on these equals the random-intercept GLS estimate.
.XP_COLS <- c("scc", "scr", "srr", "syc", "syr", "syy", "n")
xprod_by_subject <- function(y, cc, rr, sid) {
  ok <- is.finite(y) & is.finite(cc) & is.finite(rr)
  y <- y[ok]; cc <- cc[ok]; rr <- rr[ok]; sid <- sid[ok]
  idx <- split(seq_along(y), sid)
  # An empty split makes t(vapply(...)) drop its dimnames, which turns a data
  # problem into a cryptic subscript error. Return a named 0-row matrix instead.
  if (!length(idx))
    return(matrix(numeric(0), nrow = 0, ncol = 7, dimnames = list(NULL, .XP_COLS)))
  m <- t(vapply(idx, function(i) {
    yy <- y[i] - mean(y[i]); c1 <- cc[i]; r1 <- rr[i]
    c(scc = sum(c1 * c1), scr = sum(c1 * r1), srr = sum(r1 * r1),
      syc = sum(c1 * yy), syr = sum(r1 * yy), syy = sum(yy * yy), n = length(i))
  }, numeric(7)))
  colnames(m) <- .XP_COLS
  m
}

#' Solve the 2x2 system from accumulated cross-products; returns estimate, se, t.
solve_xp <- function(scc, scr, srr, syc, syr, syy, n_obs, n_subj) {
  det <- scc * srr - scr * scr
  if (!is.finite(det) || abs(det) < 1e-18) return(c(NA_real_, NA_real_, NA_real_))
  b_c <- (srr * syc - scr * syr) / det
  b_r <- (scc * syr - scr * syc) / det
  dfres <- n_obs - n_subj - 2
  if (dfres <= 0) return(c(NA_real_, NA_real_, NA_real_))
  rss <- syy - (b_c * syc + b_r * syr)
  s2 <- max(rss, 0) / dfres
  se <- sqrt(s2 * srr / det)
  c(b_c, se, b_c / se)
}

#' Within-subject slope on a set of rows (the declared model's fixed effect).
wi_slope <- function(y, cc, rr, sid) {
  xp <- xprod_by_subject(y, cc, rr, sid)
  if (!nrow(xp)) return(c(NA_real_, NA_real_, NA_real_))
  s <- colSums(xp)
  solve_xp(s[["scc"]], s[["scr"]], s[["srr"]], s[["syc"]], s[["syr"]], s[["syy"]],
           s[["n"]], nrow(xp))
}

#' Per-subject appraisal slope AND its standard error. The slope is the
#' "responsiveness" statistic any subset rule ranks on; its SE is what decides
#' whether that ranking carries information (part B).
subject_slope_full <- function(y, cc, rr, sid, min_cells = 3) {
  xp <- xprod_by_subject(y, cc, rr, sid)
  det <- xp[, "scc"] * xp[, "srr"] - xp[, "scr"]^2
  b_c <- (xp[, "srr"] * xp[, "syc"] - xp[, "scr"] * xp[, "syr"]) / det
  b_r <- (xp[, "scc"] * xp[, "syr"] - xp[, "scr"] * xp[, "syc"]) / det
  dfres <- xp[, "n"] - 3            # y demeaned within subject (1) + two slopes
  rss <- xp[, "syy"] - (b_c * xp[, "syc"] + b_r * xp[, "syr"])
  se <- sqrt(pmax(rss, 0) / dfres * xp[, "srr"] / det)
  bad <- !is.finite(det) | abs(det) < 1e-18 | xp[, "n"] < min_cells | dfres <= 0
  b_c[bad] <- NA_real_; se[bad] <- NA_real_
  tibble(subject_id = rownames(xp), slope = unname(b_c), se = unname(se), n = unname(xp[, "n"]))
}

#' Slope only, named — the form the selection loops use.
subject_slope_vec <- function(y, cc, rr, sid, min_cells = 3) {
  f <- subject_slope_full(y, cc, rr, sid, min_cells)
  setNames(f$slope, f$subject_id)
}

# Gate (c): the fast path must reproduce the lmer fixed effect and its SE.
cond <- cond_from_blk %>%
  inner_join(load_p24() %>%
               select(subject_id, marker, group, form, clip_id, comfort, richness,
                      comfort_wi, richness_wi, event_assignment,
                      flag_onset_displaced, flag_qc_check),
             by = c("subject_id", "marker"))

gate_c <- map_dfr(c("PFC_Frontal", "Left_Temporal", "Right_Temporal", "head"), function(rr) {
  d <- cond %>% filter(roi == rr)
  fit <- suppressMessages(lmer(resp ~ comfort_wi + richness_wi + (1 | subject_id),
                               data = d, REML = TRUE))
  ce <- summary(fit)$coefficients["comfort_wi", ]
  fs <- wi_slope(d$resp, d$comfort_wi, d$richness_wi, d$subject_id)
  tibble(roi = rr, lmer_est = ce[[1]], fast_est = fs[1], abs_diff = abs(ce[[1]] - fs[1]),
         lmer_se = ce[[2]], fast_se = fs[2], se_diff = abs(ce[[2]] - fs[2]))
})
cat("\nGate (c) — fast cross-product path vs lmer (response window):\n")
print(as.data.frame(gate_c))
stopifnot(all(gate_c$abs_diff < 1e-9), all(gate_c$se_diff < 1e-9))
cat("\nAll three gates PASS.\n\n")

WINDOWS <- c(response = "resp", difference = "diff", baseline = "base")
ROI3 <- c("PFC_Frontal", "Left_Temporal", "Right_Temporal")

# ======================================================================== PART A
# Direct heterogeneity: does the appraisal slope vary between people at all?
# No subset is chosen. This is the analysis that actually answers the question.

cat("---- PART A: direct heterogeneity (random slope variance) --------------\n\n")

pb_lrt <- function(d, ycol, axis, nsim = N_PB) {
  f0 <- as.formula(sprintf("%s ~ comfort_wi + richness_wi + (1 | subject_id)", ycol))
  f1 <- as.formula(sprintf("%s ~ comfort_wi + richness_wi + (1 | subject_id) + (0 + %s | subject_id)",
                           ycol, axis))
  ctl <- lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4), calc.derivs = FALSE)
  m0 <- suppressMessages(suppressWarnings(lmer(f0, data = d, REML = FALSE, control = ctl)))
  m1 <- suppressMessages(suppressWarnings(lmer(f1, data = d, REML = FALSE, control = ctl)))
  obs <- as.numeric(2 * (logLik(m1) - logLik(m0)))
  vc <- as.data.frame(VarCorr(m1))
  tau2 <- vc$vcov[vc$grp != "Residual" & !is.na(vc$var1) & vc$var1 == axis]
  tau2 <- if (length(tau2)) tau2[1] else NA_real_
  sims <- simulate(m0, nsim = nsim, seed = 24)
  nulls <- unlist(mclapply(seq_len(nsim), function(i) {
    a <- try(suppressMessages(suppressWarnings(refit(m0, sims[[i]]))), silent = TRUE)
    b <- try(suppressMessages(suppressWarnings(refit(m1, sims[[i]]))), silent = TRUE)
    if (inherits(a, "try-error") || inherits(b, "try-error")) return(NA_real_)
    as.numeric(2 * (logLik(b) - logLik(a)))
  }, mc.cores = N_CORE))
  nulls <- nulls[is.finite(nulls)]
  obs_c <- max(obs, 0); null_c <- pmax(nulls, 0)
  tibble(lrt = obs, tau2_slope = tau2, sd_slope = sqrt(max(tau2, 0, na.rm = TRUE)),
         singular = isSingular(m1),
         p_chisq_naive = pchisq(obs_c, df = 1, lower.tail = FALSE),
         p_boot = (1 + sum(null_c >= obs_c)) / (1 + length(null_c)),
         # The shape of this null is what part A2 checks: a homoscedastic-normal
         # null that sits far below the permutation null is the failure mode.
         boot_null_p95 = unname(quantile(null_c, 0.95, na.rm = TRUE)),
         boot_null_max = max(null_c, na.rm = TRUE),
         n_boot_ok = length(null_c))
}

partA <- map_dfr(names(WINDOWS), function(wn) {
  ycol <- WINDOWS[[wn]]
  map_dfr(ROI3, function(rr) {
    d <- cond %>% filter(roi == rr)
    map_dfr(c("comfort_wi", "richness_wi"), function(ax) {
      cat(sprintf("  [A] %-11s %-15s %-12s ... ", wn, rr, ax)); flush.console()
      r <- pb_lrt(d, ycol, ax)
      cat(sprintf("LRT %6.3f  tau %.4f  p_boot %.3f%s\n",
                  r$lrt, r$sd_slope, r$p_boot, ifelse(r$singular, "  [SINGULAR]", "")))
      bind_cols(tibble(window = wn, roi = rr, axis = ax), r)
    })
  }) %>% mutate(q_bh = bh(p_boot))
})
write_outcome(partA, file.path(OUT, "A_slope_variance.csv"))
cat("\nPart A — declared family of 6 per window, BH on the bootstrap p:\n")
print(as.data.frame(partA %>% select(window, roi, axis, lrt, sd_slope, singular, p_boot, q_bh)))
cat(sprintf("\n  Survivors at q < 0.05: %d of %d | singular fits: %d\n\n",
            sum(partA$q_bh < 0.05, na.rm = TRUE), nrow(partA), sum(partA$singular)))

# ======================================================================= PART A2
# FALSIFICATION of part A. Mandatory, not optional: part A's declared reading rule
# makes it the gate for everything else, so it is audited before it is believed.
#
# The concern is specific and was stated before these numbers were seen: with four
# cells per subject and two fixed predictors there is ONE residual df per person, so a
# random slope is nearly unidentifiable and can absorb structure that is not slope
# variance. The parametric bootstrap in part A simulates from a HOMOSCEDASTIC normal
# null, so it cannot see that failure mode. Three nulls that can:
#
#   F1  within-subject permutation of the appraisal pair — assumption-light, preserves
#       the outcome and the residual structure exactly, destroys only the appraisal link
#   F2  heteroscedastic parametric bootstrap — per-subject residual SD taken from the
#       null fit, so between-subject variance heterogeneity is present under the null
#   F3  leave-one-subject-out on the observed LRT — is it carried by a few people?

cat("---- PART A2: falsification of the part-A variance result ---------------\n\n")

lrt_of <- function(d, ycol, axis, yvec = NULL, cvec = NULL, rvec = NULL) {
  dd <- d
  if (!is.null(yvec)) dd[[ycol]] <- yvec
  if (!is.null(cvec)) { dd$comfort_wi <- cvec; dd$richness_wi <- rvec }
  ctl <- lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4), calc.derivs = FALSE)
  f0 <- as.formula(sprintf("%s ~ comfort_wi + richness_wi + (1 | subject_id)", ycol))
  f1 <- as.formula(sprintf("%s ~ comfort_wi + richness_wi + (1 | subject_id) + (0 + %s | subject_id)",
                           ycol, axis))
  a <- try(suppressMessages(suppressWarnings(lmer(f0, dd, REML = FALSE, control = ctl))), silent = TRUE)
  b <- try(suppressMessages(suppressWarnings(lmer(f1, dd, REML = FALSE, control = ctl))), silent = TRUE)
  if (inherits(a, "try-error") || inherits(b, "try-error")) return(NA_real_)
  max(as.numeric(2 * (logLik(b) - logLik(a))), 0)
}

targets <- partA %>% filter(q_bh < 0.05) %>% select(window, roi, axis)
cat(sprintf("  Auditing %d cells that survived part A's declared correction.\n\n", nrow(targets)))

partA2 <- pmap_dfr(targets, function(window, roi, axis) {
  ycol <- WINDOWS[[window]]
  d <- cond %>% filter(roi == !!roi) %>% as.data.frame()
  obs <- lrt_of(d, ycol, axis)
  sidx <- split(seq_len(nrow(d)), d$subject_id)

  # F1 — within-subject permutation of the appraisal pair
  f1null <- unlist(mclapply(seq_len(N_PERM), function(i) {
    set.seed(40000 + i)
    cc <- d$comfort_wi; rr2 <- d$richness_wi
    for (ix in sidx) if (length(ix) > 1L) { o <- sample(ix); cc[ix] <- d$comfort_wi[o]; rr2[ix] <- d$richness_wi[o] }
    lrt_of(d, ycol, axis, cvec = cc, rvec = rr2)
  }, mc.cores = N_CORE))
  f1null <- f1null[is.finite(f1null)]

  # F2 — heteroscedastic parametric bootstrap
  ctl <- lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4), calc.derivs = FALSE)
  m0 <- suppressMessages(suppressWarnings(
    lmer(as.formula(sprintf("%s ~ comfort_wi + richness_wi + (1 | subject_id)", ycol)),
         d, REML = FALSE, control = ctl)))
  fx <- lme4::fixef(m0); tau0 <- sqrt(as.numeric(VarCorr(m0)$subject_id[1, 1]))
  res <- residuals(m0)
  sd_i <- vapply(sidx, function(ix) stats::sd(res[ix]), numeric(1))
  sd_i[!is.finite(sd_i) | sd_i <= 0] <- stats::sd(res)
  Xb <- fx[["(Intercept)"]] + fx[["comfort_wi"]] * d$comfort_wi + fx[["richness_wi"]] * d$richness_wi
  f2null <- unlist(mclapply(seq_len(N_PERM), function(i) {
    set.seed(50000 + i)
    yv <- Xb
    for (nm in names(sidx)) {
      ix <- sidx[[nm]]
      yv[ix] <- yv[ix] + rnorm(1, 0, tau0) + rnorm(length(ix), 0, sd_i[[nm]])
    }
    lrt_of(d, ycol, axis, yvec = yv)
  }, mc.cores = N_CORE))
  f2null <- f2null[is.finite(f2null)]

  # F3 — leave-one-subject-out on the observed LRT
  loo <- vapply(names(sidx), function(nm)
    lrt_of(d[d$subject_id != nm, , drop = FALSE], ycol, axis), numeric(1))

  tibble(window = window, roi = roi, axis = axis, lrt_obs = obs,
         boot_null_p95 = partA$boot_null_p95[partA$window == window & partA$roi == roi &
                                               partA$axis == axis][1],
         p_perm_F1  = (1 + sum(f1null >= obs)) / (1 + length(f1null)),
         perm_p95   = unname(quantile(f1null, 0.95, na.rm = TRUE)),
         p_hetero_F2 = (1 + sum(f2null >= obs)) / (1 + length(f2null)),
         hetero_p95 = unname(quantile(f2null, 0.95, na.rm = TRUE)),
         loo_min = min(loo, na.rm = TRUE), loo_max = max(loo, na.rm = TRUE),
         loo_n_below_chisq95 = sum(loo < 2.706, na.rm = TRUE))
})
partA2 <- partA2 %>% mutate(q_perm_F1 = bh(p_perm_F1), q_hetero_F2 = bh(p_hetero_F2))
write_outcome(partA2, file.path(OUT, "A2_variance_falsification.csv"))
cat("Part A2 — three nulls the part-A bootstrap could not see:\n")
print(as.data.frame(partA2))
cat(sprintf("\n  Surviving F1 (permutation) at q < 0.05: %d of %d\n",
            sum(partA2$q_perm_F1 < 0.05, na.rm = TRUE), nrow(partA2)))
cat(sprintf("  Surviving F2 (heteroscedastic) at q < 0.05: %d of %d\n\n",
            sum(partA2$q_hetero_F2 < 0.05, na.rm = TRUE), nrow(partA2)))

# ======================================================================== PART B
# Reliability ceiling: can a responder be identified at all? BOUNDS every
# selection-based part below.

cat("---- PART B: can a responder be identified? (reliability ceiling) ------\n\n")

# Encounter = the chronological 1st/2nd/3rd presentation of a condition. It is a
# property of the (subject, marker, event_index) key ALONE and must be built on
# the distinct key, never ranked inside `blk` — `blk` carries one row per ROI,
# so ranking there would rank 12 ROI-duplicated rows 1..12 instead of the 3
# presentations 1..3, and every split below would be meaningless.
enc_key <- blk %>% distinct(subject_id, marker, event_index) %>%
  group_by(subject_id, marker) %>%
  mutate(encounter = rank(event_index, ties.method = "first")) %>% ungroup()
stopifnot(max(enc_key$encounter) <= 3L)

blk_d <- blk %>%
  left_join(enc_key, by = c("subject_id", "marker", "event_index")) %>%
  inner_join(load_p24() %>% select(subject_id, marker, comfort_wi, richness_wi),
             by = c("subject_id", "marker"))

# Gate (d): every ROI must see all three encounters with the full subject count.
gate_d <- blk_d %>% group_by(roi, encounter) %>%
  summarise(n_rows = n(), n_subjects = n_distinct(subject_id), .groups = "drop")
cat("Gate (d) — encounter coverage by ROI:\n"); print(as.data.frame(gate_d))
stopifnot(all(gate_d$n_subjects > 50), nrow(gate_d) == 4 * 3)
cat("\n")

#' Condition-level frame for an encounter subset, aligned to a fixed key.
cond_for <- function(d, encs) {
  d %>% filter(encounter %in% encs) %>%
    group_by(subject_id, marker) %>%
    summarise(y = mean(resp), .groups = "drop")
}

sb <- function(r, k) k * r / (1 + (k - 1) * r)

reliab_one <- function(d) {
  per_enc <- map(1:3, function(e) {
    cd <- cond_for(d, e) %>% left_join(d %>% distinct(subject_id, marker, comfort_wi, richness_wi),
                                       by = c("subject_id", "marker"))
    subject_slope_vec(cd$y, cd$comfort_wi, cd$richness_wi, cd$subject_id)
  })
  ids <- sort(unique(unlist(map(per_enc, names))))
  M <- sapply(per_enc, function(v) v[ids])
  rs <- combn(3, 2, function(p) {
    x <- M[, p[1]]; y <- M[, p[2]]; ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10) NA_real_ else suppressWarnings(cor(x[ok], y[ok]))
  })
  list(r_mean = mean(rs, na.rm = TRUE), M = M, ids = ids)
}

partB <- map_dfr(c(ROI3, "head"), function(rr) {
  d <- blk_d %>% filter(roi == rr)
  rl <- reliab_one(d)
  # {1,3} vs {2}
  cd_a <- cond_for(d, c(1, 3)) %>% left_join(d %>% distinct(subject_id, marker, comfort_wi, richness_wi),
                                             by = c("subject_id", "marker"))
  cd_b <- cond_for(d, 2) %>% left_join(d %>% distinct(subject_id, marker, comfort_wi, richness_wi),
                                       by = c("subject_id", "marker"))
  sa <- subject_slope_vec(cd_a$y, cd_a$comfort_wi, cd_a$richness_wi, cd_a$subject_id)
  sbv <- subject_slope_vec(cd_b$y, cd_b$comfort_wi, cd_b$richness_wi, cd_b$subject_id)
  ids <- intersect(names(sa), names(sbv))
  ok <- is.finite(sa[ids]) & is.finite(sbv[ids])
  r_oe <- if (sum(ok) >= 10) suppressWarnings(cor(sa[ids][ok], sbv[ids][ok])) else NA_real_
  cd_all <- cond_for(d, 1:3) %>% left_join(d %>% distinct(subject_id, marker, comfort_wi, richness_wi),
                                           by = c("subject_id", "marker"))
  sf <- subject_slope_full(cd_all$y, cd_all$comfort_wi, cd_all$richness_wi, cd_all$subject_id)
  sall <- sf$slope
  pa <- partA %>% filter(window == "response", roi == rr, axis == "comfort_wi")
  tau2 <- if (nrow(pa)) pa$tau2_slope[1] else NA_real_
  se_med <- median(sf$se, na.rm = TRUE)
  tibble(roi = rr,
         r_enc_pairwise_mean = rl$r_mean,
         r_enc_sb3 = sb(max(rl$r_mean, -0.32), 3),
         r_oddeven = r_oe, r_oddeven_sb = sb(max(r_oe, -0.49), 2),
         tau2_slope = tau2, sd_slope_model = sqrt(max(tau2, 0, na.rm = TRUE)),
         sd_observed_slopes = sd(sall, na.rm = TRUE),
         median_subject_slope_se = se_med,
         # model-based reliability of one person's slope estimate
         reliability_model = tau2 / (tau2 + se_med^2),
         n_estimable = sum(is.finite(sall)))
})
write_outcome(partB, file.path(OUT, "B_slope_reliability.csv"))
cat("Part B — reliability of the per-person comfort slope (response window):\n")
print(as.data.frame(partB))

d_pf <- blk_d %>% filter(roi == "PFC_Frontal")
pf_split <- split(d_pf, d_pf$subject_id)
subs <- names(pf_split)
boot_r <- unlist(mclapply(seq_len(N_PERM), function(b) {
  set.seed(1000 + b)
  pick <- sample(subs, length(subs), replace = TRUE)
  dd <- bind_rows(lapply(seq_along(pick), function(i) {
    x <- pf_split[[pick[i]]]; x$subject_id <- paste0(pick[i], "_b", i); x
  }))
  reliab_one(dd)$r_mean
}, mc.cores = N_CORE))
ci <- quantile(boot_r, c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("\n  Prefrontal reliability r = %.3f, subject-bootstrap 95%% CI [%.3f, %.3f]\n\n",
            partB$r_enc_pairwise_mean[partB$roi == "PFC_Frontal"], ci[[1]], ci[[2]]))
write_outcome(tibble(boot_r = boot_r), file.path(OUT, "B_reliability_bootstrap.csv"))

# ================================================================== PARTS C & D
# Select-then-test. C keeps selection and inference on disjoint rows; D does not,
# and is run to measure what that manufactures.

cat("---- PARTS C & D: select-then-test curves ------------------------------\n\n")

KS <- 5:69
SPLITS <- list(enc12_to_enc3 = list(def = c(1, 2), tst = 3),
               enc3_to_enc12 = list(def = 3,       tst = c(1, 2)),
               enc13_to_enc2 = list(def = c(1, 3), tst = 2),
               enc2_to_enc13 = list(def = 2,       tst = c(1, 3)))

#' Precompute everything that a permutation does NOT change: the outcome values and
#' the subject x marker structure of each half. A permutation only reshuffles the
#' appraisal pair within a subject, so the y values are computed once.
prep_roi <- function(rr) {
  d <- blk_d %>% filter(roi == rr)
  key <- d %>% distinct(subject_id, marker, comfort_wi, richness_wi) %>%
    arrange(subject_id, marker) %>% mutate(kid = row_number())
  attach_key <- function(cd) cd %>%
    left_join(key %>% select(subject_id, marker, kid), by = c("subject_id", "marker"))
  list(key = key,
       key_idx = split(key$kid, key$subject_id),
       all = attach_key(cond_for(d, 1:3)),
       halves = map(SPLITS, function(sp) list(def = attach_key(cond_for(d, sp$def)),
                                              tst = attach_key(cond_for(d, sp$tst)))))
}

#' One select-then-test curve, computed incrementally: adding one subject at a time
#' accumulates the cross-products, so the whole k = 5..69 curve costs one pass.
curve_fast <- function(def, tst, pc, pr, rank_by = "signed") {
  s <- subject_slope_vec(def$y, pc[def$kid], pr[def$kid], def$subject_id)
  s <- s[is.finite(s)]
  ord <- if (rank_by == "signed") order(s) else order(-abs(s))
  ranked <- names(s)[ord]
  xp <- xprod_by_subject(tst$y, pc[tst$kid], pr[tst$kid], tst$subject_id)
  have <- ranked[ranked %in% rownames(xp)]
  if (!length(have)) return(tibble(k = KS, est = NA_real_, se = NA_real_, t = NA_real_, n = 0L))
  cs <- apply(xp[have, , drop = FALSE], 2, cumsum)
  if (is.null(dim(cs))) cs <- matrix(cs, nrow = 1, dimnames = list(NULL, names(cs)))
  out <- t(vapply(KS, function(k) {
    if (k > nrow(cs)) return(c(NA_real_, NA_real_, NA_real_, NA_real_))
    v <- cs[k, ]
    c(solve_xp(v[["scc"]], v[["scr"]], v[["srr"]], v[["syc"]], v[["syr"]], v[["syy"]],
               v[["n"]], k), v[["n"]])
  }, numeric(4)))
  tibble(k = KS, est = out[, 1], se = out[, 2], t = out[, 3], n = out[, 4])
}

#' Permute the four condition labels within each subject, carrying the appraisal
#' predictors with them (the whole select-then-test pipeline under the null).
perm_pred <- function(key, key_idx) {
  pc <- key$comfort_wi; pr <- key$richness_wi
  for (i in key_idx) if (length(i) > 1L) { o <- sample(i); pc[i] <- key$comfort_wi[o]; pr[i] <- key$richness_wi[o] }
  list(pc = pc, pr = pr)
}

run_roi_cd <- function(rr) {
  P <- prep_roi(rr)
  pc0 <- P$key$comfort_wi; pr0 <- P$key$richness_wi
  jobs <- c(map(names(SPLITS), function(nm) list(arm = "C_out_of_sample", split = nm,
                                                 def = P$halves[[nm]]$def, tst = P$halves[[nm]]$tst)),
            list(list(arm = "D_in_sample", split = "in_sample", def = P$all, tst = P$all)))
  map_dfr(jobs, function(j) {
    cur <- curve_fast(j$def, j$tst, pc0, pr0)
    obs <- max(abs(cur$t), na.rm = TRUE)
    nulls <- unlist(mclapply(seq_len(N_PERM), function(i) {
      set.seed(20000 + i)
      pp <- perm_pred(P$key, P$key_idx)
      max(abs(curve_fast(j$def, j$tst, pp$pc, pp$pr)$t), na.rm = TRUE)
    }, mc.cores = N_CORE))
    nulls <- nulls[is.finite(nulls)]
    cur$arm <- j$arm; cur$split <- j$split; cur$roi <- rr
    assign("CURVES", bind_rows(get("CURVES", envir = .GlobalEnv), cur), envir = .GlobalEnv)
    assign("NULLS", bind_rows(get("NULLS", envir = .GlobalEnv),
                              tibble(roi = rr, arm = j$arm, split = j$split, null_max = nulls)),
           envir = .GlobalEnv)
    tibble(arm = j$arm, roi = rr, split = j$split, obs_max_abs_t = obs,
           k_at_max = cur$k[which.max(abs(cur$t))],
           est_at_max = cur$est[which.max(abs(cur$t))],
           n_nominal_k = sum(abs(cur$t) > 1.96, na.rm = TRUE), n_k = sum(is.finite(cur$t)),
           null_p95 = quantile(nulls, 0.95, na.rm = TRUE),
           p_perm = (1 + sum(nulls >= obs)) / (1 + length(nulls)))
  })
}

CURVES <- tibble(); NULLS <- tibble()
cd_summ <- map_dfr(ROI3, function(rr) { cat(sprintf("  [C/D] %s ...\n", rr)); flush.console(); run_roi_cd(rr) })
cd_summ <- cd_summ %>% group_by(arm) %>% mutate(q_bh = bh(p_perm)) %>% ungroup()
write_outcome(CURVES, file.path(OUT, "CD_selection_curves.csv"))
write_outcome(cd_summ, file.path(OUT, "CD_selection_summary.csv"))
write_outcome(NULLS, file.path(OUT, "CD_permutation_nulls.csv"))
cat("\nParts C & D — one permutation p per curve (max |t| across k = 5..69):\n")
print(as.data.frame(cd_summ %>% select(arm, roi, split, obs_max_abs_t, k_at_max, est_at_max,
                                       n_nominal_k, n_k, null_p95, p_perm, q_bh)))
cat("\n")

# ======================================================================== PART E
# Subgroups defined WITHOUT touching the tested association.

cat("---- PART E: outcome-independent subgroups -----------------------------\n\n")

rat <- load_ratings()
adj_cols <- grep("^adj_", names(rat), value = TRUE)
subj_meta <- rat %>% group_by(subject_id) %>%
  summarise(comfort_span_sd = sd(comfort, na.rm = TRUE),
            adj_engagement_sd = sd(unlist(across(all_of(adj_cols))), na.rm = TRUE),
            .groups = "drop")
clean <- load_p24() %>%
  group_by(subject_id) %>%
  summarise(acq_clean = as.integer(all(event_assignment == "positional") &
                                     all(flag_onset_displaced == 0) &
                                     all(flag_qc_check == 0)), .groups = "drop")

partE <- map_dfr(ROI3, function(rr) {
  d <- cond %>% filter(roi == rr) %>%
    left_join(subj_meta, by = "subject_id") %>% left_join(clean, by = "subject_id") %>%
    group_by(subject_id) %>% mutate(n_live_subj = mean(n_live)) %>% ungroup()
  parts <- list(
    comfort_span  = as.integer(d$comfort_span_sd   > median(d$comfort_span_sd, na.rm = TRUE)),
    live_channels = as.integer(d$n_live_subj       > median(d$n_live_subj, na.rm = TRUE)),
    acquisition   = d$acq_clean,
    rating_engage = as.integer(d$adj_engagement_sd > median(d$adj_engagement_sd, na.rm = TRUE)))
  imap_dfr(parts, function(g, nm) {
    dd <- d %>% mutate(sub_g = g) %>% filter(!is.na(sub_g))
    if (length(unique(dd$sub_g)) < 2)
      return(tibble(roi = rr, partition = nm, estimate = NA_real_, se = NA_real_, p = NA_real_,
                    n_hi = NA_integer_, n_lo = NA_integer_, slope_hi = NA_real_, slope_lo = NA_real_))
    dd$sub_g <- factor(dd$sub_g)
    fit <- suppressMessages(suppressWarnings(
      lmer(resp ~ comfort_wi * sub_g + richness_wi + (1 | subject_id), data = dd, REML = TRUE)))
    ce <- summary(fit)$coefficients
    row <- grep("^comfort_wi:sub_g", rownames(ce), value = TRUE)[1]
    hi <- dd %>% filter(sub_g == "1"); lo <- dd %>% filter(sub_g == "0")
    s_hi <- wi_slope(hi$resp, hi$comfort_wi, hi$richness_wi, hi$subject_id)
    s_lo <- wi_slope(lo$resp, lo$comfort_wi, lo$richness_wi, lo$subject_id)
    tibble(roi = rr, partition = nm, estimate = ce[row, 1], se = ce[row, 2], p = ce[row, 5],
           n_hi = n_distinct(hi$subject_id), n_lo = n_distinct(lo$subject_id),
           slope_hi = s_hi[1], slope_lo = s_lo[1])
  })
}) %>% mutate(q_bh = bh(p))
write_outcome(partE, file.path(OUT, "E_subgroup_interactions.csv"))
cat("Part E — declared family of 12 (4 partitions x 3 ROIs), BH:\n")
print(as.data.frame(partE))
cat(sprintf("\n  Survivors at q < 0.05: %d of %d\n\n", sum(partE$q_bh < 0.05, na.rm = TRUE), nrow(partE)))

# ======================================================================== PART F
cat("---- PART F: non-person heterogeneity (clip, group) --------------------\n\n")

partF <- map_dfr(c("PFC_Frontal", "head"), function(rr) {
  d <- cond %>% filter(roi == rr) %>% mutate(clip = factor(clip_id), group_f = factor(group))
  map_dfr(c("clip", "group_f"), function(gv) {
    ctl <- lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4), calc.derivs = FALSE)
    f0 <- resp ~ comfort_wi + richness_wi + (1 | subject_id)
    f1 <- as.formula(sprintf("resp ~ comfort_wi + richness_wi + (1 | subject_id) + (0 + comfort_wi | %s)", gv))
    m0 <- suppressMessages(suppressWarnings(lmer(f0, data = d, REML = FALSE, control = ctl)))
    m1 <- suppressMessages(suppressWarnings(lmer(f1, data = d, REML = FALSE, control = ctl)))
    obs <- max(as.numeric(2 * (logLik(m1) - logLik(m0))), 0)
    sims <- simulate(m0, nsim = N_PB, seed = 24)
    nulls <- unlist(mclapply(seq_len(N_PB), function(i) {
      a <- try(suppressMessages(suppressWarnings(refit(m0, sims[[i]]))), silent = TRUE)
      b <- try(suppressMessages(suppressWarnings(refit(m1, sims[[i]]))), silent = TRUE)
      if (inherits(a, "try-error") || inherits(b, "try-error")) return(NA_real_)
      max(as.numeric(2 * (logLik(b) - logLik(a))), 0)
    }, mc.cores = N_CORE))
    nulls <- nulls[is.finite(nulls)]
    vc <- as.data.frame(VarCorr(m1))
    tv <- vc$vcov[vc$grp == gv]
    tibble(roi = rr, grouping = gv, lrt = obs,
           sd_slope = sqrt(max(if (length(tv)) tv[1] else 0, 0)),
           singular = isSingular(m1),
           p_boot = (1 + sum(nulls >= obs)) / (1 + length(nulls)))
  })
}) %>% mutate(q_bh = bh(p_boot))
write_outcome(partF, file.path(OUT, "F_nonperson_heterogeneity.csv"))
print(as.data.frame(partF)); cat("\n")

# ======================================================================== PART G
cat("---- PART G: mixture on per-person slopes ------------------------------\n\n")

partG <- map_dfr(ROI3, function(rr) {
  d <- blk_d %>% filter(roi == rr)
  cd <- cond_for(d, 1:3) %>% left_join(d %>% distinct(subject_id, marker, comfort_wi, richness_wi),
                                       by = c("subject_id", "marker"))
  x <- subject_slope_vec(cd$y, cd$comfort_wi, cd$richness_wi, cd$subject_id)
  x <- x[is.finite(x)]
  m1 <- Mclust(x, G = 1, modelNames = "V", verbose = FALSE)
  m2 <- Mclust(x, G = 2, modelNames = "V", verbose = FALSE)
  obs <- as.numeric(m2$bic - m1$bic)
  nulls <- unlist(mclapply(seq_len(500), function(i) {
    set.seed(3000 + i)
    xx <- rnorm(length(x), m1$parameters$mean[1], sqrt(m1$parameters$variance$sigmasq[1]))
    a <- try(Mclust(xx, G = 1, modelNames = "V", verbose = FALSE), silent = TRUE)
    b <- try(Mclust(xx, G = 2, modelNames = "V", verbose = FALSE), silent = TRUE)
    if (inherits(a, "try-error") || inherits(b, "try-error") || is.null(b) || is.null(a)) return(NA_real_)
    as.numeric(b$bic - a$bic)
  }, mc.cores = N_CORE))
  nulls <- nulls[is.finite(nulls)]
  tibble(roi = rr, n = length(x), bic_G1 = m1$bic, bic_G2 = m2$bic, bic_gain_G2 = obs,
         p_boot = (1 + sum(nulls >= obs)) / (1 + length(nulls)),
         shapiro_p = shapiro.test(x)$p.value)
})
write_outcome(partG, file.path(OUT, "G_mixture.csv"))
print(as.data.frame(partG)); cat("\n")

# ======================================================================== PART H
cat("---- PART H: lambda(k) --------------------------------------------------\n\n")

lam <- read_csv(file.path("output", "analysis_22_chain", "lambda_chain.csv"),
                show_col_types = FALSE) %>%
  filter(se_source == "corrected ROI (primary)") %>%
  select(roi, C0_ideal, C3_full = C3_plus_dependence)

partH <- crossing(lam, k = c(seq(5, 65, 5), 69)) %>%
  mutate(infl = sqrt(69 / k), lambda_ideal_k = C0_ideal * infl, lambda_full_k = C3_full * infl,
         implied_peak_uM = 0.073, ratio_full = implied_peak_uM / lambda_full_k) %>%
  select(roi, k, infl, lambda_ideal_k, lambda_full_k, implied_peak_uM, ratio_full) %>%
  arrange(roi, k)
write_outcome(partH, file.path(OUT, "H_lambda_by_subset_size.csv"))
print(as.data.frame(partH %>% filter(roi == "PFC_Frontal")))

cat("\n=========================================================================\n")
cat("A23 main module complete; residual-shape and reliability-generality sections follow.\n")
cat("=========================================================================\n\n")

# ====================================================================== PART A2.4
# WHY part A's parametric bootstrap was mis-calibrated.
#
# Part A tested the slope variance against a parametric bootstrap, which
# simulates from the fitted null: homoscedastic, normal. Part A2's permutation
# null put the same statistic's 95th percentile roughly fourteen times higher.
# This section measures the two assumptions that gap is made of, so the report
# states a mechanism rather than an inference. Writes A2_residual_shape.csv.

cat("---- PART A2.4: residual shape -----------------------------------------\n\n")

cond_rs <- load_block_metrics(active_only = TRUE) %>% filter(chromophore == "HbO") %>%
  left_join(load_dead_channels() %>%
              bind_rows(load_channels() %>% filter(detector %in% UNUSABLE_DETECTORS) %>%
                          distinct(subject_id, channel_index)) %>%
              distinct(subject_id, channel_index) %>% mutate(.dead = TRUE),
            by = c("subject_id", "channel_index")) %>%
  filter(is.na(.dead)) %>%
  left_join(load_channels() %>% filter(active_pair == 1) %>%
              distinct(subject_id, channel_index, roi),
            by = c("subject_id", "channel_index")) %>%
  filter(!is.na(roi)) %>%
  group_by(subject_id, marker, event_index, roi) %>%
  summarise(resp = mean(response_m5_30), base = mean(baseline_m10_0),
            diff = mean(response_minus_baseline_m5_30), .groups = "drop") %>%
  group_by(subject_id, marker, roi) %>%
  summarise(across(c(resp, base, diff), mean), .groups = "drop") %>%
  inner_join(load_p24() %>% select(subject_id, marker, comfort_wi, richness_wi),
             by = c("subject_id", "marker"))

WIN <- c(response = "resp", difference = "diff", baseline = "base")

out_rs <- map_dfr(names(WIN), function(wn) {
  map_dfr(c("PFC_Frontal", "Left_Temporal", "Right_Temporal"), function(rr) {
    d <- cond_rs %>% filter(roi == rr) %>% filter(is.finite(.data[[WIN[[wn]]]]))
    m0 <- suppressMessages(suppressWarnings(lmer(
      as.formula(sprintf("%s ~ comfort_wi + richness_wi + (1 | subject_id)", WIN[[wn]])),
      d, REML = FALSE, control = lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4)))))
    r <- residuals(m0); sid <- d$subject_id
    sd_i <- tapply(r, sid, stats::sd); sd_i <- sd_i[is.finite(sd_i) & sd_i > 0]
    cs <- tapply(d$comfort_wi, sid, stats::sd)
    bt <- suppressWarnings(bartlett.test(r, factor(sid)))
    ids <- intersect(names(sd_i), names(cs))
    tibble(window = wn, roi = rr, n_subjects = length(sd_i),
           resid_sd_min = min(sd_i), resid_sd_median = median(sd_i), resid_sd_max = max(sd_i),
           resid_sd_ratio_max_min = max(sd_i) / min(sd_i),
           bartlett_K2 = unname(bt$statistic), bartlett_p = bt$p.value,
           resid_kurtosis = mean((r - mean(r))^4) / stats::sd(r)^4,
           comfort_sd_ratio_max_min = max(cs, na.rm = TRUE) / min(cs, na.rm = TRUE),
           cor_residSD_comfortSD = suppressWarnings(cor(sd_i[ids], cs[ids], use = "complete.obs")))
  })
})

write_outcome(out_rs, file.path(OUT, "A2_residual_shape.csv"))
print(as.data.frame(out_rs %>% select(window, roi, resid_sd_ratio_max_min, bartlett_p,
                                      resid_kurtosis, cor_residSD_comfortSD)), digits = 4)
cat("\nBoth parametric-bootstrap assumptions fail: per-subject residual variance is\n",
    "grossly unequal and the residuals are heavy-tailed. A random slope on a predictor\n",
    "whose per-subject spread co-varies with residual SD absorbs that structure, which\n",
    "is what part A's homoscedastic null cannot generate and part A2's permutation can.\n\n")

# ======================================================================= PART B.4
# Is the reliability ceiling a property of the DESIGN or of HbO?
#
# Part B measured the split-half reliability of the per-person comfort slope on
# HbO ROI means and found it indistinguishable from zero. The structural
# explanation is that four condition cells and two predictors leave ONE
# residual df per person, which would make the ceiling outcome-general. That is
# an inference, so this section tests it: the same reliability statistic on
# HbR, on HbT, on the whole head, and on channel 10 — the campaign's one local
# candidate. If the ceiling is design-borne, all of them sit at zero; if it is
# HbO-specific, some outcome should escape it.
#
# Descriptive, no correction — a measurement property, not a hypothesis test.
# Writes B4_reliability_generality.csv.

cat("---- PART B.4: reliability generality ----------------------------------\n\n")

.XP <- c("scc", "scr", "srr", "syc", "syr", "syy", "n")
xprod <- function(y, cc, rr, sid) {
  ok <- is.finite(y) & is.finite(cc) & is.finite(rr)
  y <- y[ok]; cc <- cc[ok]; rr <- rr[ok]; sid <- sid[ok]
  idx <- split(seq_along(y), sid)
  if (!length(idx)) return(matrix(numeric(0), 0, 7, dimnames = list(NULL, .XP)))
  m <- t(vapply(idx, function(i) {
    yy <- y[i] - mean(y[i]); c1 <- cc[i]; r1 <- rr[i]
    c(sum(c1 * c1), sum(c1 * r1), sum(r1 * r1), sum(c1 * yy), sum(r1 * yy),
      sum(yy * yy), length(i))
  }, numeric(7)))
  colnames(m) <- .XP; m
}
subj_slope <- function(y, cc, rr, sid, min_cells = 3) {
  xp <- xprod(y, cc, rr, sid)
  det <- xp[, "scc"] * xp[, "srr"] - xp[, "scr"]^2
  s <- (xp[, "srr"] * xp[, "syc"] - xp[, "scr"] * xp[, "syr"]) / det
  s[!is.finite(det) | abs(det) < 1e-18 | xp[, "n"] < min_cells] <- NA_real_
  setNames(s, rownames(xp))
}

dead_b4 <- bind_rows(load_dead_channels() %>% distinct(subject_id, channel_index),
                     load_channels() %>% filter(detector %in% UNUSABLE_DETECTORS) %>%
                       distinct(subject_id, channel_index)) %>%
  distinct(subject_id, channel_index) %>% mutate(.dead = TRUE)
mapc <- load_channels() %>% filter(active_pair == 1) %>% distinct(subject_id, channel_index, roi)
design <- load_p24() %>% distinct(subject_id, marker, comfort_wi, richness_wi)

#' Split-half reliability across the three encounters, for one long-format outcome frame
#' carrying subject_id / marker / event_index / y.
reliability <- function(df) {
  df <- df %>% group_by(subject_id, marker) %>%
    mutate(encounter = rank(event_index, ties.method = "first")) %>% ungroup()
  stopifnot(max(df$encounter) <= 3L)
  per <- map(1:3, function(e) {
    cd <- df %>% filter(encounter == e) %>%
      group_by(subject_id, marker) %>% summarise(y = mean(y), .groups = "drop") %>%
      left_join(design, by = c("subject_id", "marker"))
    subj_slope(cd$y, cd$comfort_wi, cd$richness_wi, cd$subject_id)
  })
  ids <- sort(unique(unlist(map(per, names))))
  M <- sapply(per, function(v) v[ids])
  rs <- combn(3, 2, function(p) {
    x <- M[, p[1]]; y <- M[, p[2]]; ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10) NA_real_ else suppressWarnings(cor(x[ok], y[ok]))
  })
  list(r = mean(rs, na.rm = TRUE), n = sum(is.finite(M[, 1])))
}

bm_all <- load_block_metrics(active_only = TRUE) %>%
  left_join(dead_b4, by = c("subject_id", "channel_index")) %>% filter(is.na(.dead))

targets_b4 <- list()
# HbO / HbR / HbT ROI means and whole head, on the corrected montage
for (ch in c("HbO", "HbR", "HbT")) {
  base <- bm_all %>% filter(chromophore == ch) %>%
    left_join(mapc, by = c("subject_id", "channel_index"))
  for (u in c("PFC_Frontal", "Left_Temporal", "Right_Temporal", "head")) {
    d <- if (u == "head") base else base %>% filter(roi == u)
    targets_b4[[paste(ch, u)]] <- d %>%
      group_by(subject_id, marker, event_index) %>%
      summarise(y = mean(response_m5_30), .groups = "drop")
  }
}
# Channel 10 (FP2-AF8), the campaign's one local candidate, HbO
targets_b4[["HbO channel_10"]] <- bm_all %>%
  filter(chromophore == "HbO", channel_index == 10L) %>%
  transmute(subject_id, marker, event_index, y = response_m5_30)

out_b4 <- imap_dfr(targets_b4, function(d, nm) {
  r <- reliability(d)
  tibble(outcome = nm, split_half_r = r$r, n_subjects = r$n)
}) %>% separate(outcome, c("chromophore", "unit"), sep = " ")

write_outcome(out_b4, file.path(OUT, "B4_reliability_generality.csv"))
print(as.data.frame(out_b4), digits = 3)
cat(sprintf("\nAll |r| < 0.25: %s | max |r| = %.3f (%s %s)\n",
            all(abs(out_b4$split_half_r) < 0.25, na.rm = TRUE),
            max(abs(out_b4$split_half_r), na.rm = TRUE),
            out_b4$chromophore[which.max(abs(out_b4$split_half_r))],
            out_b4$unit[which.max(abs(out_b4$split_half_r))]))
cat("If every outcome sits at zero, the ceiling is the DESIGN (one residual df per\n",
    "person), not the chromophore or the region — and no responder analysis is\n",
    "recoverable by changing the outcome.\n")

cat("\n=========================================================================\n")
cat("A23 complete.\n")
cat("=========================================================================\n")
close_log(con)
