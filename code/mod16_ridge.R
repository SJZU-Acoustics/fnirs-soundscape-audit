# Module 16 — Multivariate distributed readout of appraisal from the channel pattern
#
# Port of exploration/analysis_16_multivariate_pattern (batch-5 declaration):
#
# Every spatial analysis in A01-A15 is mass-univariate: one channel at a time,
# BH across the family. A08 closed the single-channel route under robust
# correction and A11 found the effect distributed. A weak DISTRIBUTED effect is
# exactly what per-channel BH is worst at detecting and what a cross-validated
# multivariate readout is best at — it asks one question, of the whole pattern,
# and answers it with one permutation p-value, so there is no family to correct
# and no opportunity to gerrymander one.
#
# Declared family (fixed before the original analysis ran):
#   PRIMARY   4 tests = 2 windows (response [5,30); frozen difference score)
#             x 2 pattern definitions (within-subject-centred channel values;
#             the same with each cell's across-channel mean removed = pure
#             spatial pattern, no global level), BH.
#   SECONDARY richness, the same 4, BH separately.
#
# Channel set: the 13 prefrontal channels (1-15 less the two distance-excluded
# 6 and 13). 64 of 69 subjects have all 13 active, against 50 of 69 for the full
# 21. The choice is made on COMPLETENESS and on prefrontal being the a-priori
# unit of A04/A06/A15 — not on any outcome. The full-21 complete-case set
# (n = 50) is reported alongside as a declared robustness variant.
#
# Estimator: ridge, nested leave-one-subject-out. lambda is chosen by an inner
# 5-fold-over-subjects CV *within the training subjects only* — never using the
# held-out subject, and re-selected inside every permutation, so the null runs
# the identical pipeline. Everything is linear in y, so all 1,000 permutations
# are carried as columns of one matrix and the whole nested CV is a handful of
# small matrix products.
#
# Run from the repository root (Rscript code/mod16_ridge.R). Reads only the
# Mendeley workbook via load_data.R; writes
# output/analysis_16_multivariate_pattern/. Runtime is long (the nested CV is
# re-run inside each of 1,000 permutations for every declared test).

source("code/load_data.R")

OUT <- file.path("output", "analysis_16_multivariate_pattern")
con <- open_log(file.path(OUT, "multivariate_pattern_report.txt"))
rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")

set.seed(24)
N_PERM  <- BOOT_N                      # 1,000, house default
LAMBDAS <- 10^seq(-3, 3, length.out = 13)
PFC13   <- setdiff(1:15, c(6, 13))
ALL21   <- setdiff(1:23, c(6, 13))
WINDOWS <- c(response   = "mean_response_m5_30",
             difference = "mean_response_minus_baseline_m5_30")

cd <- load_condition_analysis("HbO")

# ------------------------------------------------------------- matrix builder
# Returns X (cells x channels) and y, both within-subject centred, restricted to
# subjects with every requested channel active. `spatial = TRUE` first removes
# each CELL's across-channel mean, leaving pure spatial pattern with no global
# level; the within-subject centring then follows.
center_within <- function(M, sid) {
  for (s in unique(sid)) {
    i <- which(sid == s)
    M[i, ] <- sweep(M[i, , drop = FALSE], 2, colMeans(M[i, , drop = FALSE]), "-")
  }
  M
}
build_xy <- function(channels, ycol, target, spatial) {
  w <- cd %>% filter(channel_index %in% channels) %>%
    select(subject_id, marker, channel_index, value = all_of(ycol)) %>%
    pivot_wider(names_from = channel_index, values_from = value, names_prefix = "ch")
  chcols <- paste0("ch", channels)
  keep <- w %>% filter(if_all(all_of(chcols), ~ !is.na(.x))) %>%
    count(subject_id) %>% filter(n == 4) %>% pull(subject_id)
  w <- w %>% filter(subject_id %in% keep) %>% arrange(subject_id, marker) %>%
    inner_join(.analysis_design("full") %>% select(subject_id, marker, all_of(target)),
               by = c("subject_id", "marker"))
  X <- as.matrix(w[, chcols])
  if (spatial) X <- X - rowMeans(X)
  X <- center_within(X, w$subject_id)
  list(X = X, y = as.numeric(w[[target]]), subject = w$subject_id,
       n_sub = n_distinct(w$subject_id), n_cell = nrow(w))
}

# ------------------------------------------------------- nested LOSO ridge + permutation
run_mvpa <- function(dat, n_perm = N_PERM, seed = 24) {
  X <- dat$X; y <- dat$y; sid <- dat$subject
  n <- nrow(X); p <- ncol(X); subs <- unique(sid)
  nS <- length(subs); nL <- length(LAMBDAS)

  # Y: column 1 observed, columns 2..n_perm+1 within-subject permutations.
  # Permuting the 4 condition labels within a subject is the exact null: it
  # destroys the pattern-appraisal link while preserving every subject's own
  # pattern set, their 4 ratings, and the within-subject centring.
  set.seed(seed)
  Y <- matrix(0, n, n_perm + 1); Y[, 1] <- y
  idx_by_sub <- lapply(subs, function(s) which(sid == s))
  for (j in seq_len(n_perm)) {
    yy <- y
    for (i in idx_by_sub) yy[i] <- y[i][sample.int(length(i))]
    Y[, j + 1] <- yy
  }

  # inner folds: subjects split into 5 groups, fixed across permutations
  inner_group <- setNames(sample(rep_len(1:5, nS)), subs)

  PRED  <- array(0, dim = c(n, n_perm + 1, nL))       # outer predictions per lambda
  IMSE  <- array(0, dim = c(nS, n_perm + 1, nL))      # inner CV error per outer fold

  for (f in seq_len(nS)) {
    te <- which(sid == subs[f]); tr <- setdiff(seq_len(n), te)
    Xtr <- X[tr, , drop = FALSE]; Xte <- X[te, , drop = FALSE]
    Ytr <- Y[tr, , drop = FALSE]
    G <- crossprod(Xtr); XtY <- crossprod(Xtr, Ytr)

    for (l in seq_len(nL)) {
      B <- solve(G + LAMBDAS[l] * diag(p), XtY)        # p x (n_perm+1)
      PRED[te, , l] <- Xte %*% B
    }

    # inner CV, training subjects only
    tr_subs <- subs[subs != subs[f]]
    for (g in 1:5) {
      ite_sub <- tr_subs[inner_group[tr_subs] == g]
      if (!length(ite_sub)) next
      ite <- which(sid %in% ite_sub); itr <- setdiff(tr, ite)
      Xitr <- X[itr, , drop = FALSE]; Xite <- X[ite, , drop = FALSE]
      Gi <- crossprod(Xitr); XtYi <- crossprod(Xitr, Y[itr, , drop = FALSE])
      for (l in seq_len(nL)) {
        Bi <- solve(Gi + LAMBDAS[l] * diag(p), XtYi)
        R <- Y[ite, , drop = FALSE] - Xite %*% Bi
        IMSE[f, , l] <- IMSE[f, , l] + colSums(R^2)
      }
    }
  }

  # per outer fold and per permutation, take the inner-CV-optimal lambda
  pred <- matrix(0, n, n_perm + 1)
  best_lambda <- matrix(0, nS, n_perm + 1)
  for (f in seq_len(nS)) {
    te <- which(sid == subs[f])
    bl <- max.col(-IMSE[f, , , drop = TRUE], ties.method = "first")   # (n_perm+1)
    best_lambda[f, ] <- LAMBDAS[bl]
    for (l in unique(bl)) {
      j <- which(bl == l)
      pred[te, j] <- PRED[te, j, l]
    }
  }

  stat <- vapply(seq_len(n_perm + 1), function(j) {
    s <- suppressWarnings(cor(pred[, j], Y[, j])); if (is.na(s)) 0 else s }, numeric(1))
  obs <- stat[1]; null <- stat[-1]
  list(r_obs = obs,
       p_perm = (1 + sum(null >= obs)) / (n_perm + 1),      # one-sided: prediction
       null_mean = mean(null), null_sd = sd(null),
       null_q95 = unname(quantile(null, 0.95)),
       median_lambda = median(best_lambda[, 1]),
       n_sub = dat$n_sub, n_cell = dat$n_cell, p = p)
}

# ------------------------------------------------------------------ 1. primary
rule("1. PRIMARY FAMILY — comfort, 13 prefrontal channels, 4 tests, BH")

grid <- expand_grid(window = names(WINDOWS), pattern = c("channels", "spatial"))
prim <- pmap_dfr(grid, function(window, pattern) {
  dat <- build_xy(PFC13, WINDOWS[[window]], "comfort_wi", spatial = (pattern == "spatial"))
  r <- run_mvpa(dat)
  tibble(target = "comfort_wi", window = window, pattern = pattern, !!!r)
}) %>% mutate(q_BH = bh(p_perm))
print(as.data.frame(prim %>% select(window, pattern, n_sub, n_cell, p, r_obs,
                                    null_mean, null_q95, p_perm, q_BH)),
      digits = 3, row.names = FALSE)
write_outcome(prim, file.path(OUT, "primary_comfort_pfc13.csv"))

# ---------------------------------------------------------------- 2. secondary
rule("2. SECONDARY FAMILY — richness, same 4 tests, BH separately")

sec <- pmap_dfr(grid, function(window, pattern) {
  dat <- build_xy(PFC13, WINDOWS[[window]], "richness_wi", spatial = (pattern == "spatial"))
  r <- run_mvpa(dat)
  tibble(target = "richness_wi", window = window, pattern = pattern, !!!r)
}) %>% mutate(q_BH = bh(p_perm))
print(as.data.frame(sec %>% select(window, pattern, n_sub, n_cell, r_obs,
                                   null_q95, p_perm, q_BH)), digits = 3, row.names = FALSE)
write_outcome(sec, file.path(OUT, "secondary_richness_pfc13.csv"))

# ------------------------------------------------- 3. robustness: all 21, n=50
rule("3. ROBUSTNESS — all 21 channels, complete cases (n = 50 subjects)")

rob <- pmap_dfr(grid, function(window, pattern) {
  dat <- build_xy(ALL21, WINDOWS[[window]], "comfort_wi", spatial = (pattern == "spatial"))
  r <- run_mvpa(dat)
  tibble(target = "comfort_wi", window = window, pattern = pattern, !!!r)
})
print(as.data.frame(rob %>% select(window, pattern, n_sub, n_cell, p, r_obs,
                                   null_q95, p_perm)), digits = 3, row.names = FALSE)
write_outcome(rob, file.path(OUT, "robustness_all21_comfort.csv"))

# ----------------------------------------------------- 4. univariate benchmark
rule("4. BENCHMARK — what the best SINGLE channel achieves under the same CV")

# The point of comparison that matters: does the pattern beat its own best
# channel? Same nested LOSO, same permutation null, one predictor at a time.
bench <- map_dfr(c("response", "difference"), function(wn) {
  dat <- build_xy(PFC13, WINDOWS[[wn]], "comfort_wi", spatial = FALSE)
  map_dfr(seq_len(ncol(dat$X)), function(k) {
    d1 <- list(X = dat$X[, k, drop = FALSE], y = dat$y, subject = dat$subject,
               n_sub = dat$n_sub, n_cell = dat$n_cell)
    r <- run_mvpa(d1, n_perm = 200)
    tibble(window = wn, channel_index = PFC13[k], r_obs = r$r_obs, p_perm = r$p_perm)
  })
})
cat("\nBest single channels per window (nested LOSO, 200 perms):\n")
print(as.data.frame(bench %>% group_by(window) %>% slice_max(r_obs, n = 3) %>% ungroup()),
      digits = 3, row.names = FALSE)
cat("\nPattern vs best single channel:\n")
print(as.data.frame(
  prim %>% filter(pattern == "channels") %>% select(window, r_pattern = r_obs) %>%
    left_join(bench %>% group_by(window) %>% summarise(r_best_single = max(r_obs),
              best_channel = channel_index[which.max(r_obs)], .groups = "drop"),
              by = "window")), digits = 3, row.names = FALSE)
write_outcome(bench, file.path(OUT, "single_channel_benchmark.csv"))

rule("DONE")
close_log(con)
