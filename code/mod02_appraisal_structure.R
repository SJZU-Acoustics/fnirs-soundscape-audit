# A02 — Structure of the subjective side.
#
# Do the eight adjectives really span a two-dimensional space with Comfort and
# Richness as its axes? How much of each axis is a property of the clip (all
# listeners agree) versus of the listener (idiosyncratic appraisal)? And does any
# listener fail to use the scale, which would make them uninformative for the
# within-subject slopes A04 rests on?

source("code/load_data.R")

OUT  <- file.path("output", "analysis_02_appraisal_structure")
con  <- open_log(file.path(OUT, "appraisal_structure_report.txt"))

d <- load_p24()
r <- load_ratings()
subj <- load_subjects()

rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")
cat("A02 — structure of the subjective appraisal space\n")

adj_cols <- paste0("adj_", ADJECTIVES$zh)
stopifnot(all(adj_cols %in% names(r)))
A <- r %>% select(all_of(adj_cols)) %>% rename_with(~ ADJECTIVES$en, all_of(adj_cols))

# ------------------------------------------------- 1. dimensionality of the 8
rule("1. IS THE ADJECTIVE SPACE TWO-DIMENSIONAL?")
cat("n =", nrow(A), "clip ratings; response scale 1-5.\n\n")
cat("item means and SDs:\n")
print(as.data.frame(A %>% pivot_longer(everything()) %>% group_by(name) %>%
  summarise(mean = mean(value), sd = sd(value), min = min(value), max = max(value),
            .groups = "drop")), digits = 3)

pc <- prcomp(A, scale. = TRUE)
ev <- pc$sdev^2
cat("\neigenvalues / variance explained:\n")
print(as.data.frame(tibble(component = seq_along(ev), eigenvalue = ev,
       pct = 100 * ev / sum(ev), cum_pct = cumsum(100 * ev / sum(ev)))), digits = 3)

cat("\nloadings on the first two components, with each item's correlation with the\n")
cat("study's Comfort and Richness axes:\n")
load2 <- tibble(item = colnames(A), PC1 = pc$rotation[, 1], PC2 = pc$rotation[, 2],
                r_comfort = as.vector(cor(A, r$comfort)),
                r_richness = as.vector(cor(A, r$richness)))
print(as.data.frame(load2), digits = 3)
write_outcome(load2, file.path(OUT, "adjective_loadings.csv"))

# Is the PC plane the same plane as the Comfort/Richness plane? Regress each
# axis on the two PC scores: if the axes are a rotation of the first two PCs,
# R2 is essentially 1.
sc <- as_tibble(pc$x[, 1:2]); names(sc) <- c("pc1", "pc2")
rot <- tibble(
  axis = c("comfort", "richness"),
  R2_on_first_two_PCs = c(summary(lm(r$comfort ~ sc$pc1 + sc$pc2))$r.squared,
                          summary(lm(r$richness ~ sc$pc1 + sc$pc2))$r.squared))
cat("\nare Comfort/Richness a rotation of the first two components?\n")
print(as.data.frame(rot), digits = 5)
write_outcome(rot, file.path(OUT, "axis_pc_rotation.csv"))
cat("r(comfort, richness) =", round(cor(r$comfort, r$richness), 4), "\n")

# --------------------------------------- 2. consensus: clip versus listener
rule("2. CONSENSUS — HOW MUCH OF AN AXIS IS THE CLIP, HOW MUCH THE LISTENER?")
cat("Variance components from y ~ 1 + (1|subject_id) + (1|clip_id).\n")
cat("Each clip was heard by 16-18 listeners; each listener rated 4 clips.\n\n")

vc <- map_dfr(c("comfort", "richness"), function(v) {
  m <- suppressMessages(lmer(as.formula(paste(v, "~ 1 + (1|subject_id) + (1|clip_id)")),
                             data = r, REML = TRUE))
  vs <- as.data.frame(VarCorr(m))
  s2_sub  <- vs$vcov[vs$grp == "subject_id"]
  s2_clip <- vs$vcov[vs$grp == "clip_id"]
  s2_res  <- vs$vcov[vs$grp == "Residual"]
  tot <- s2_sub + s2_clip + s2_res
  tibble(axis = v, var_subject = s2_sub, var_clip = s2_clip, var_residual = s2_res,
         pct_subject = 100 * s2_sub / tot, pct_clip = 100 * s2_clip / tot,
         pct_residual = 100 * s2_res / tot)
})
print(as.data.frame(vc), digits = 3)
write_outcome(vc, file.path(OUT, "axis_variance_components.csv"))

cat("\nreliability of a clip's mean position (ICC(2,k), k = mean raters per clip):\n")
k <- nrow(r) / n_distinct(r$clip_id)
icc2k <- vc %>% mutate(k = k,
  ICC_1 = var_clip / (var_clip + var_subject + var_residual),
  ICC_2k = (var_clip) / (var_clip + (var_subject + var_residual) / k)) %>%
  select(axis, k, ICC_1, ICC_2k)
print(as.data.frame(icc2k), digits = 3)
write_outcome(icc2k, file.path(OUT, "clip_position_reliability.csv"))
cat("\nICC_1 is how much a single listener's rating tells you about the clip;\n")
cat("ICC_2k is how reliable the clip's MEAN position is with k raters — the\n")
cat("quantity that matters, because the stimulus component of the within-subject\n")
cat("variation (comfort_stim / richness_stim) is built from those clip means.\n")

# ------------------------------------------------ 3. circumplex geometry check
rule("3. CIRCUMPLEX GEOMETRY OF THE FOUR NOMINAL QUADRANTS")
geo <- d %>% group_by(quadrant) %>%
  summarise(comfort = mean(comfort), richness = mean(richness), .groups = "drop") %>%
  mutate(angle_deg = (atan2(richness, comfort) * 180 / pi) %% 360,
         radius = sqrt(comfort^2 + richness^2))
print(as.data.frame(geo), digits = 3)
cat("\nintended angles: Q1 45, Q2 135, Q3 225, Q4 315 degrees.\n")
cat("observed gaps between consecutive quadrants (degrees):",
    paste(round(diff(c(geo$angle_deg, geo$angle_deg[1] + 360)), 1), collapse = ", "), "\n")
write_outcome(geo, file.path(OUT, "quadrant_geometry.csv"))

cat("\nsame, per clip (16 cells) — the three known off-quadrant clips should show:\n")
geo_clip <- d %>% group_by(group, quadrant, clip_id, flag_off_quadrant_stimulus) %>%
  summarise(comfort = mean(comfort), richness = mean(richness), .groups = "drop") %>%
  mutate(angle_deg = (atan2(richness, comfort) * 180 / pi) %% 360,
         radius = sqrt(comfort^2 + richness^2),
         intended_angle = c(Q1 = 45, Q2 = 135, Q3 = 225, Q4 = 315)[as.character(quadrant)],
         angle_error = ((angle_deg - intended_angle + 180) %% 360) - 180)
print(as.data.frame(geo_clip %>% arrange(desc(abs(angle_error)))), digits = 3)
write_outcome(geo_clip, file.path(OUT, "clip_geometry.csv"))

# ------------------------------------------- 4. does every listener use the scale?
rule("4. DOES EVERY LISTENER DISCRIMINATE? (uninformative subjects for A04)")
disc <- r %>% group_by(subject_id) %>%
  summarise(sd_comfort = sd(comfort), sd_richness = sd(richness),
            range_comfort = diff(range(comfort)), range_richness = diff(range(richness)),
            n_distinct_comfort = n_distinct(round(comfort, 6)), .groups = "drop")
cat("per-subject spread across their four clips:\n")
print(as.data.frame(disc %>% summarise(across(-subject_id,
      list(min = min, median = median, max = max)))), digits = 3)
flat <- disc %>% filter(sd_comfort < 0.10 | sd_richness < 0.10)
cat("\nsubjects with near-flat use of an axis (SD < 0.10 across their four clips):",
    nrow(flat), "\n")
if (nrow(flat) > 0) print(as.data.frame(flat), digits = 3)
write_outcome(disc, file.path(OUT, "listener_discrimination.csv"))

# ------------------------------------------ 5. do person covariates shape ratings?
rule("5. DO PERSON COVARIATES SHAPE THE RATINGS? (family = 8, BH)")
pm <- r %>% group_by(subject_id) %>%
  summarise(comfort_mean = mean(comfort), richness_mean = mean(richness),
            comfort_sd = sd(comfort), richness_sd = sd(richness), .groups = "drop") %>%
  left_join(subj %>% select(subject_id, sex_code, age, who5_total_0_25, nss_total_keyed),
            by = "subject_id")
tests <- expand_grid(outcome = c("comfort_mean", "richness_mean"),
                     predictor = c("age", "sex_code", "who5_total_0_25", "nss_total_keyed"))
res <- pmap_dfr(tests, function(outcome, predictor) {
  dd <- pm %>% select(y = all_of(outcome), x = all_of(predictor)) %>% drop_na()
  ct <- cor.test(dd$y, dd$x)
  tibble(outcome = outcome, predictor = predictor, n = nrow(dd),
         r = unname(ct$estimate), p = ct$p.value)
}) %>% mutate(q_BH = bh(p)) %>% arrange(p)
print(as.data.frame(res), digits = 3)
write_outcome(res, file.path(OUT, "covariate_rating_scan.csv"))
cat("\nfamily size = 8 (2 rating summaries x 4 person variables), BH corrected.\n")

close_log(con)
cat("A02 done.\n")
