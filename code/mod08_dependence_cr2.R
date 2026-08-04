# mod08_dependence_cr2.R — Dependence-structure and estimand audit.
#
# Release port of the stage-2 exploration pass A08
# (exploration/analysis_08_dependence_structure_audit).
#
# Why this extra pass exists: module 07's leave-one-subject-out and cluster
# bootstrap showed that channel 10 was not a single-person artefact, but a
# direct cluster-sandwich calculation materially widened its SE. The module 06
# model has four condition rows per subject and a random intercept only; its
# model-based SE assumes the remaining within-subject covariance is correctly
# specified.
#
# This module therefore repeats module 06's entire, unchanged 42-test family
# per window with clubSandwich CR2 SEs clustered by subject and Satterthwaite
# df. This is a post-discovery audit. No family is narrowed, and no
# replacement winner is selected.
#
# Run from the repository root:  Rscript code/mod08_dependence_cr2.R

source("code/load_data.R")
suppressPackageStartupMessages(library(clubSandwich))

OUT <- file.path("output", "analysis_08_dependence_structure_audit")
con <- open_log(file.path(OUT, "dependence_structure_audit_report.txt"))

WINDOWS <- c(
  difference = "mean_response_minus_baseline_m5_30",
  response   = "mean_response_m5_30",
  baseline   = "mean_baseline_m10_0"
)
AXES <- c("comfort_wi", "richness_wi")
TARGET_CHANNEL <- 10L
TARGET_DV <- "mean_response_m5_30"

rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n",
                        strrep("=", 74), "\n", sep = "")

d <- load_condition_analysis("HbO")
channels <- sort(unique(d$channel_index))
stopifnot(length(channels) == 21L)

cat("Module 08 — dependence-structure and estimand audit\n")
cat("Post-discovery audit: same 42 tests/window as module 06, CR2 clustered",
    "by subject.\n")

# ------------------------------------------------------ 1. robust full screen
rule("1. REPEAT THE FULL CHANNEL SCREEN WITH CR2 SUBJECT-CLUSTER INFERENCE")

robust_scan <- map_dfr(names(WINDOWS), function(win) {
  y <- unname(WINDOWS[[win]])
  map_dfr(channels, function(cc) {
    dd <- d %>% filter(channel_index == cc)
    fit <- p24_lmm(dd, paste(AXES, collapse = " + "), y)
    ct <- as.data.frame(
      clubSandwich::coef_test(
        fit,
        vcov = "CR2",
        cluster = dd$subject_id,
        test = "Satterthwaite"
      )
    )
    ct$term <- rownames(ct)
    as_tibble(ct) %>%
      filter(term %in% AXES) %>%
      transmute(
        window = win,
        channel_index = cc,
        n_cells = nrow(dd),
        n_subjects = n_distinct(dd$subject_id),
        term,
        estimate = beta,
        se_CR2 = SE,
        df_CR2 = df_Satt,
        t_CR2 = tstat,
        p_CR2 = p_Satt
      )
  })
}) %>%
  group_by(window) %>%
  mutate(
    q_CR2_BH = bh(p_CR2),
    family_size = n(),
    rank_CR2 = min_rank(p_CR2)
  ) %>%
  ungroup()
stopifnot(all(robust_scan %>% count(window) %>% pull(n) == 42L))

a06_file <- file.path("output", "analysis_06_channelwise_spatial_scan",
                      "channel_scan_all.csv")
if (!file.exists(a06_file)) {
  stop("Module 06 output not found at ", a06_file, "\n",
       "Run code/mod06_channel_scan.R first.", call. = FALSE)
}
a06 <- readr::read_csv(a06_file, show_col_types = FALSE) %>%
  select(
    window, channel_index, term,
    se_model = se, df_model = df, t_model = t,
    p_model = p, q_model_BH = q_BH
  )

comparison <- robust_scan %>%
  left_join(a06, by = c("window", "channel_index", "term")) %>%
  mutate(
    se_inflation_CR2_over_model = se_CR2 / se_model,
    p_ratio_CR2_over_model = p_CR2 / p_model
  ) %>%
  arrange(factor(window, levels = names(WINDOWS)), p_CR2)

family_summary <- comparison %>%
  group_by(window) %>%
  summarise(
    family_size = first(family_size),
    model_based_min_q = min(q_model_BH),
    model_based_survivors = sum(q_model_BH < 0.05),
    CR2_min_p = min(p_CR2),
    CR2_min_q = min(q_CR2_BH),
    CR2_survivors = sum(q_CR2_BH < 0.05),
    CR2_best_channel = channel_index[which.min(p_CR2)],
    CR2_best_axis = term[which.min(p_CR2)],
    .groups = "drop"
  )

print(as.data.frame(family_summary), digits = 4)
cat("\nTop CR2 tests in each family:\n")
print(as.data.frame(
  comparison %>%
    group_by(window) %>%
    slice_min(p_CR2, n = 5, with_ties = FALSE) %>%
    ungroup() %>%
    select(window, channel_index, term, estimate, se_model, se_CR2,
           p_model, q_model_BH, p_CR2, q_CR2_BH)
), digits = 4)

write_outcome(comparison, file.path(OUT, "channel_scan_CR2_comparison.csv"))
write_outcome(family_summary, file.path(OUT, "family_summary_CR2.csv"))

# ---------------------------------------------------- 2. target estimand panel
rule("2. CHANNEL-10 ESTIMATOR TRIANGULATION")

cc <- d %>% filter(channel_index == TARGET_CHANNEL)
ri <- p24_lmm(cc, paste(AXES, collapse = " + "), TARGET_DV)
ri_conv <- fixef_table(ri, keep = "comfort_wi")
ri_cr2_raw <- as.data.frame(
  clubSandwich::coef_test(
    ri,
    vcov = "CR2",
    cluster = cc$subject_id,
    test = "Satterthwaite"
  )
)
ri_cr2_raw$term <- rownames(ri_cr2_raw)
ri_cr2 <- as_tibble(ri_cr2_raw) %>% filter(term == "comfort_wi")

# Random comfort slope: a different covariance model, estimable with four
# conditions per subject. It is a targeted sensitivity model, not a new screen.
rs <- lmerTest::lmer(
  mean_response_m5_30 ~ comfort_wi + richness_wi +
    (1 + comfort_wi | subject_id),
  data = cc,
  REML = TRUE,
  control = lmerControl(
    check.conv.singular = .makeCC("ignore", tol = 1e-4)
  )
)
rs_tab <- fixef_table(rs, keep = "comfort_wi")

# Subject fixed-effect coefficient with a CR1 cluster sandwich, implemented on
# within-subject demeaned rows. This reproduces the pooled coefficient without
# relying on the random-intercept covariance model.
within <- cc %>%
  group_by(subject_id) %>%
  mutate(
    y_dm = .data[[TARGET_DV]] - mean(.data[[TARGET_DV]]),
    comfort_dm = comfort_wi - mean(comfort_wi),
    richness_dm = richness_wi - mean(richness_wi)
  ) %>%
  ungroup()
X <- as.matrix(within %>% select(comfort_dm, richness_dm))
y <- within$y_dm
beta <- solve(crossprod(X), crossprod(X, y))
resid <- as.vector(y - X %*% beta)
cluster_ids <- unique(within$subject_id)
meat <- matrix(0, nrow = ncol(X), ncol = ncol(X))
for (sid in cluster_ids) {
  ii <- within$subject_id == sid
  score <- crossprod(X[ii, , drop = FALSE], resid[ii])
  meat <- meat + tcrossprod(score)
}
N <- nrow(X)
K <- ncol(X)
G <- length(cluster_ids)
bread <- solve(crossprod(X))
V_CR1 <- (G / (G - 1)) * ((N - 1) / (N - K)) *
  bread %*% meat %*% bread
se_CR1 <- sqrt(diag(V_CR1))
t_CR1 <- beta / se_CR1
p_CR1 <- 2 * pt(abs(t_CR1), df = G - 1, lower.tail = FALSE)

# Average of separately estimated person-specific slopes. With four conditions
# and two predictors each slope has one residual df, so this intentionally
# gives every participant equal weight and is a different, very noisy estimand.
subject_slopes <- cc %>%
  group_by(subject_id) %>%
  group_modify(~ {
    fit <- lm(mean_response_m5_30 ~ comfort_wi + richness_wi, data = .x)
    tibble(
      comfort_slope = unname(coef(fit)["comfort_wi"]),
      comfort_sd = sd(.x$comfort_wi),
      design_kappa = kappa(model.matrix(fit))
    )
  }) %>%
  ungroup()
subject_slope_test <- t.test(subject_slopes$comfort_slope)

a07_file <- file.path("output", "analysis_07_focal_signal_adjudication",
                      "bootstrap_ci.csv")
if (!file.exists(a07_file)) {
  stop("Module 07 output not found at ", a07_file, "\n",
       "Run code/mod07_channel10_adjudication.R first.", call. = FALSE)
}
boot <- readr::read_csv(a07_file, show_col_types = FALSE) %>%
  filter(term == "comfort_wi")

method_panel <- bind_rows(
  tibble(
    method = "random-intercept LMM, model SE",
    estimand = "pooled within-subject slope",
    estimate = ri_conv$estimate,
    se = ri_conv$se,
    df = ri_conv$df,
    p = ri_conv$p,
    ci_lo = ri_conv$estimate - qt(0.975, ri_conv$df) * ri_conv$se,
    ci_hi = ri_conv$estimate + qt(0.975, ri_conv$df) * ri_conv$se
  ),
  tibble(
    method = "random-intercept LMM, CR2 SE",
    estimand = "pooled within-subject slope",
    estimate = ri_cr2$beta,
    se = ri_cr2$SE,
    df = ri_cr2$df_Satt,
    p = ri_cr2$p_Satt,
    ci_lo = ri_cr2$beta - qt(0.975, ri_cr2$df_Satt) * ri_cr2$SE,
    ci_hi = ri_cr2$beta + qt(0.975, ri_cr2$df_Satt) * ri_cr2$SE
  ),
  tibble(
    method = "random-comfort-slope LMM",
    estimand = "population mean of varying slopes",
    estimate = rs_tab$estimate,
    se = rs_tab$se,
    df = rs_tab$df,
    p = rs_tab$p,
    ci_lo = rs_tab$estimate - qt(0.975, rs_tab$df) * rs_tab$se,
    ci_hi = rs_tab$estimate + qt(0.975, rs_tab$df) * rs_tab$se
  ),
  tibble(
    method = "subject fixed effects, CR1 SE",
    estimand = "pooled within-subject slope",
    estimate = unname(beta[1]),
    se = unname(se_CR1[1]),
    df = G - 1,
    p = unname(p_CR1[1]),
    ci_lo = unname(beta[1] - qt(0.975, G - 1) * se_CR1[1]),
    ci_hi = unname(beta[1] + qt(0.975, G - 1) * se_CR1[1])
  ),
  tibble(
    method = "unweighted mean of 68 subject slopes",
    estimand = "mean person-specific slope",
    estimate = mean(subject_slopes$comfort_slope),
    se = sd(subject_slopes$comfort_slope) / sqrt(nrow(subject_slopes)),
    df = unname(subject_slope_test$parameter),
    p = subject_slope_test$p.value,
    ci_lo = subject_slope_test$conf.int[1],
    ci_hi = subject_slope_test$conf.int[2]
  ),
  tibble(
    method = "subject-cluster bootstrap, 1,000",
    estimand = "pooled within-subject slope",
    estimate = ri_conv$estimate,
    se = NA_real_,
    df = NA_real_,
    p = NA_real_,
    ci_lo = boot$ci_lo,
    ci_hi = boot$ci_hi
  )
)

print(as.data.frame(method_panel), digits = 4)
cat("\nRandom-slope model singular:", lme4::isSingular(rs), "\n")
cat("Subject-specific slopes negative:",
    sum(subject_slopes$comfort_slope < 0), "of", nrow(subject_slopes),
    "; median", round(median(subject_slopes$comfort_slope), 4), "\n")

write_outcome(method_panel, file.path(OUT, "channel10_method_comparison.csv"))
write_outcome(subject_slopes, file.path(OUT, "channel10_subject_slopes.csv"))

# ------------------------------------------ 3. effect of the SE assumption
rule("3. MATERIALITY OF THE DEPENDENCE ASSUMPTION")

target_compare <- comparison %>%
  filter(
    window == "response",
    channel_index == TARGET_CHANNEL,
    term == "comfort_wi"
  ) %>%
  transmute(
    estimate,
    model_se = se_model,
    CR2_se = se_CR2,
    SE_increase_pct = 100 * (se_CR2 / se_model - 1),
    model_p = p_model,
    model_q = q_model_BH,
    CR2_p = p_CR2,
    CR2_q = q_CR2_BH
  )
print(as.data.frame(target_compare), digits = 4)
write_outcome(target_compare, file.path(OUT, "channel10_inference_shift.csv"))

close_log(con)
cat("Module 08 done.\n")
