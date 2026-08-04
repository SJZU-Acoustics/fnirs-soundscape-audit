# mod10_person_side.R — The person side: who responds, and who responds with appraisal?
#
# Person covariates were previously scanned against the RATINGS (nothing
# survived); this module scans them against the HAEMODYNAMICS. Two questions:
#   (i) between-person: does a trait (age, sex, WHO-5, noise sensitivity)
#       predict a person's mean haemodynamic reactivity?
#   (ii) moderation: does a trait moderate the within-subject appraisal->HbO
#       slope (a cross-level interaction)?
#
# Declared families: response window [5,30) is primary (the candidate-bearing
# window, free of the baseline contamination); the frozen difference score is
# descriptive alongside, no correction.
#   (i) 4 traits x 3 ROIs = 12 tests, BH. Age/sex n = 69; WHO-5/NSS n = 59.
#   (ii) 4 traits x 2 axes x 3 ROIs = 24 interaction tests, BH.
#
# AGE CAVEAT. The deposited workbook anonymises age to five-year bands and the
# loader's `age` column is the band lower bound; the frozen pipeline used exact
# integer ages. Age-dependent rows therefore differ from the frozen run, and
# because BH q values are computed on the whole family p-vector, every q in the
# two families shifts slightly. The module uses the loader's age as-is; the
# report log prints both the recomputed (band-age) and the frozen (exact-age)
# family minima.
#
# Outputs (output/analysis_10_person_side/):
#   between_person_reactivity.csv,
#   between_person_reactivity_difference_descriptive.csv,
#   trait_moderation_response.csv,
#   trait_moderation_difference_descriptive.csv
# Display anchors: main-text Table 2 and SI Table S1 quote the family minima
# (min q 0.273 reactivity / 0.277 moderation, frozen exact-age values).

if (!exists("load_p24")) source("code/load_data.R")

OUT <- file.path("output", "analysis_10_person_side")
con <- open_log(file.path(OUT, "person_side_report.txt"))

rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n", strrep("=", 74), "\n", sep = "")

# subject x condition x ROI means for one chromophore, joined to the canonical
# analysis design — built here from the public loaders.
condition_roi_long <- function(chromophore) {
  map <- load_channels() %>%
    filter(active_pair == 1) %>%
    distinct(subject_id, channel_index, roi)
  load_condition_metrics(active_only = TRUE) %>%
    filter(.data$chromophore == !!chromophore) %>%
    left_join(map, by = c("subject_id", "channel_index")) %>%
    filter(!is.na(roi)) %>%
    group_by(subject_id, marker, roi) %>%
    summarise(
      across(c(mean_baseline_m10_0, mean_stim_m0_30, mean_response_m5_30,
               mean_response_minus_baseline_m5_30), mean),
      n_channels = n(),
      .groups = "drop"
    ) %>%
    inner_join(.analysis_design("full"), by = c("subject_id", "marker"))
}

roi_long <- condition_roi_long("HbO") %>%
  mutate(sex = factor(sex_code, levels = c(1, 2), labels = c("male", "female")))
d_diff <- load_p24() %>%
  mutate(sex = factor(sex_code, levels = c(1, 2), labels = c("male", "female")))

TRAITS <- c("age", "sex", "who5", "nss")

# ================================================================== (i) between
rule("(i) BETWEEN-PERSON REACTIVITY ~ TRAITS — response window, family of 12")

# subject-level mean reactivity per ROI (equal weight over the four cells)
subj_react <- roi_long %>%
  group_by(subject_id, roi) %>%
  summarise(reactivity = mean(mean_response_m5_30, na.rm = TRUE), .groups = "drop") %>%
  left_join(load_subjects() %>%
              transmute(subject_id, sex = factor(sex_code, c(1, 2), c("male", "female")),
                        age, who5_total_0_25, nss_total_keyed, has_covariates),
            by = "subject_id")

between_rows <- map_dfr(ROIS, function(rr) {
  dd <- subj_react %>% filter(roi == rr, !is.na(reactivity))
  m_demo <- lm(reactivity ~ z(age) + sex, data = dd)
  m_quest <- lm(reactivity ~ z(who5_total_0_25) + z(nss_total_keyed),
                data = dd %>% filter(has_covariates == 1))
  bind_rows(
    broom::tidy(m_demo) %>%
      filter(term %in% c("z(age)", "sexfemale")) %>%
      mutate(trait = c("age", "sex"), n = nrow(dd)),
    broom::tidy(m_quest) %>%
      filter(term %in% c("z(who5_total_0_25)", "z(nss_total_keyed)")) %>%
      mutate(trait = c("who5", "nss"), n = sum(dd$has_covariates == 1))
  ) %>% mutate(roi = rr, .before = 1)
}) %>% mutate(q_BH = bh(p.value), family_size = n())

print(as.data.frame(between_rows %>%
      select(roi, trait, estimate, std.error, statistic, p.value, q_BH, n)), digits = 3)
write_outcome(between_rows, file.path(OUT, "between_person_reactivity.csv"))

# descriptive difference-score version (no correction)
between_diff <- map_dfr(ROIS, function(rr) {
  y <- paste0("HbO_", rr)
  dd <- d_diff %>%
    group_by(subject_id) %>%
    summarise(reactivity = mean(.data[[y]], na.rm = TRUE),
              sex = first(sex), age = first(age),
              who5 = first(who5_total_0_25), nss = first(nss_total_keyed),
              has_covariates = first(has_covariates), .groups = "drop") %>%
    filter(!is.na(reactivity))
  m1 <- lm(reactivity ~ z(age) + sex, data = dd)
  m2 <- lm(reactivity ~ z(who5) + z(nss), data = dd %>% filter(has_covariates == 1))
  bind_rows(
    broom::tidy(m1) %>% filter(term %in% c("z(age)", "sexfemale")) %>%
      mutate(trait = c("age", "sex")),
    broom::tidy(m2) %>% filter(term %in% c("z(who5)", "z(nss)")) %>%
      mutate(trait = c("who5", "nss"))
  ) %>% mutate(roi = rr, .before = 1)
}) %>% mutate(window = "difference (descriptive, uncorrected)")
write_outcome(between_diff, file.path(OUT, "between_person_reactivity_difference_descriptive.csv"))

# ================================================================== (ii) moderation
rule("(ii) TRAIT x APPRAISAL MODERATION — response window, family of 24")

trait_col <- list(age = "age", sex = "sex",
                  who5 = "who5_total_0_25", nss = "nss_total_keyed")

mod_rows <- map_dfr(TRAITS, function(tr) {
  col <- trait_col[[tr]]
  map_dfr(ROIS, function(rr) {
    dd <- roi_long %>% filter(roi == rr, !is.na(mean_response_m5_30))
    if (tr %in% c("who5", "nss")) dd <- dd %>% filter(has_covariates == 1)
    dd <- dd %>% mutate(tv = if (tr == "sex") as.numeric(sex == "female") * 1.0
                            else z(.data[[col]]))
    f <- "comfort_wi * tv + richness_wi * tv + comfort_bw + richness_bw"
    m <- p24_lmm(dd, f, "mean_response_m5_30")
    fixef_table(m, keep = c("comfort_wi:tv", "richness_wi:tv", "tv:richness_wi")) %>%
      mutate(term = recode(term, "tv:richness_wi" = "richness_wi:tv"),
             trait = tr, roi = rr, n_rows = nrow(dd),
             n_subjects = n_distinct(dd$subject_id), .before = 1)
  })
}) %>% mutate(q_BH = bh(p), family_size = n())

print(as.data.frame(mod_rows %>%
      select(trait, roi, term, estimate, se, t, p, q_BH, n_subjects)), digits = 3)
write_outcome(mod_rows, file.path(OUT, "trait_moderation_response.csv"))

# descriptive difference-score version (no correction)
mod_diff <- map_dfr(TRAITS, function(tr) {
  col <- trait_col[[tr]]
  map_dfr(ROIS, function(rr) {
    y <- paste0("HbO_", rr)
    dd <- d_diff %>% filter(!is.na(.data[[y]]))
    if (tr %in% c("who5", "nss")) dd <- dd %>% filter(has_covariates == 1)
    dd <- dd %>% mutate(tv = if (tr == "sex") as.numeric(sex == "female") * 1.0
                            else z(.data[[col]]))
    f <- "comfort_wi * tv + richness_wi * tv + comfort_bw + richness_bw"
    m <- p24_lmm(dd, f, y)
    fixef_table(m, keep = c("comfort_wi:tv", "richness_wi:tv", "tv:richness_wi")) %>%
      mutate(term = recode(term, "tv:richness_wi" = "richness_wi:tv"),
             trait = tr, roi = rr, n_subjects = n_distinct(dd$subject_id), .before = 1)
  })
}) %>% mutate(window = "difference (descriptive, uncorrected)")
write_outcome(mod_diff, file.path(OUT, "trait_moderation_difference_descriptive.csv"))

cat("\nfamily sizes: between", unique(between_rows$family_size),
    "; moderation", unique(mod_rows$family_size), "\n")
cat("surviving q < 0.05: between", sum(between_rows$q_BH < 0.05),
    "; moderation", sum(mod_rows$q_BH < 0.05), "\n")
cat("family minima (five-year-band age): between", min(between_rows$q_BH),
    "; moderation", min(mod_rows$q_BH), "\n")
cat("frozen exact-age minima quoted in Table 2 / Table S1: 0.273 / 0.277\n")

close_log(con)
cat("\nmod10 done.\n")
