# =============================================================================
# Recompute of A18's in-memory permutation nulls (position contrast), so that
# Fig 5a can draw the null envelope for the position-minus-neither contrast.
#
# Method: exact replica of isc_one() in code/mod18_clip_reliability.R
# Same seed (24), same unit order (PFC, LT, RT, head) -> identical RNG stream,
# so the recomputed permutation p values must equal the module's stored ones in
#   output/analysis_18_clip_reliability_and_adjectives/neural_clip_reliability.csv
# That equality is the verification gate of this script; run mod18 first.
# Output: output/data_lock/fig5a_position_contrast_null.csv
#         output/data_lock/fig5a_contrast_observed.csv
# =============================================================================

source("code/load_data.R")
source("code/style.R")

set.seed(24)
N_PERM <- BOOT_N  # 1000
BINS_POST <- c("[0,5)", "[5,10)", "[10,15)", "[15,20)", "[20,25)", "[25,30)",
               "[30,35)", "[35,40)", "[40,45)", "[45,50)", "[50,55)", "[55,60)")
UNITS <- c("PFC_Frontal", "Left_Temporal (descriptive)",
           "Right_Temporal (descriptive)", "head (descriptive)")
UNIT_LABEL <- c("PFC_Frontal" = "Prefrontal", "Left_Temporal (descriptive)" = "Left temporal",
                "Right_Temporal (descriptive)" = "Right temporal",
                "head (descriptive)" = "Whole head")

bins <- load_intermediate("condition_bin_values_hbo.csv")

subj_group <- load_subjects() %>% select(subject_id, group)
wave_id <- load_p24() %>%
  select(subject_id, marker, clip_id, flag_shared_q4_waveform) %>%
  mutate(waveform = if_else(flag_shared_q4_waveform == 1, "shared_Q4_wdt16", clip_id))

isc_one <- function(unit) {
  w <- bins %>% filter(target == unit, bin %in% BINS_POST) %>%
    select(subject_id, marker, bin, referenced) %>%
    mutate(bin = factor(bin, levels = BINS_POST)) %>%
    arrange(subject_id, marker, bin) %>%
    pivot_wider(names_from = bin, values_from = referenced) %>%
    inner_join(subj_group, by = "subject_id") %>%
    inner_join(wave_id %>% select(subject_id, marker, waveform),
               by = c("subject_id", "marker"))
  M <- as.matrix(w[, BINS_POST])
  ok <- stats::complete.cases(M) & apply(M, 1, function(v) sd(v) > 0)
  w <- w[ok, ]; M <- M[ok, , drop = FALSE]

  Z <- t(scale(t(M)))
  C <- (Z %*% t(Z)) / (ncol(Z) - 1)

  n <- nrow(w)
  same_sub <- outer(w$subject_id, w$subject_id, "==")
  pairs    <- !same_sub & upper.tri(C)

  stat_for <- function(lab_wave, lab_mark) {
    same_wave <- outer(lab_wave, lab_wave, "==")
    same_mark <- outer(lab_mark, lab_mark, "==")
    a  <- mean(C[pairs & same_wave])
    pc <- mean(C[pairs & same_mark & !same_wave])
    b  <- mean(C[pairs & !same_mark])
    c(audio_vs_position = a - pc, audio_vs_diff = a - b, position_vs_diff = pc - b,
      r_same_wave = a, r_pos_ctrl = pc, r_diff = b)
  }
  obs <- stat_for(w$waveform, w$marker)

  idx_by_sub <- split(seq_len(n), w$subject_id)
  null <- matrix(NA_real_, N_PERM, 3,
                 dimnames = list(NULL, c("audio_vs_position", "audio_vs_diff",
                                         "position_vs_diff")))
  for (j in seq_len(N_PERM)) {
    lw <- w$waveform; lm_ <- w$marker
    for (i in idx_by_sub) {
      o <- sample.int(length(i)); lw[i] <- w$waveform[i][o]; lm_[i] <- w$marker[i][o]
    }
    null[j, ] <- stat_for(lw, lm_)[1:3]
  }
  list(unit = unname(UNIT_LABEL[unit]), obs = obs, null = null,
       p_position_vs_diff = (1 + sum(null[, "position_vs_diff"] >=
                                       obs[["position_vs_diff"]])) / (N_PERM + 1),
       p_audio_vs_position = (1 + sum(null[, "audio_vs_position"] >=
                                        obs[["audio_vs_position"]])) / (N_PERM + 1))
}

res <- lapply(UNITS, isc_one)  # same order as the module run -> same RNG stream

# ---- verification gate: recomputed p must equal stored p --------------------
stored <- readr::read_csv(
  file.path("output", "analysis_18_clip_reliability_and_adjectives",
            "neural_clip_reliability.csv"), show_col_types = FALSE)

chk <- tibble(
  unit = vapply(res, `[[`, "", "unit"),
  p_position_recomputed = vapply(res, `[[`, 0, "p_position_vs_diff"),
  p_audio_recomputed    = vapply(res, `[[`, 0, "p_audio_vs_position")
) %>% left_join(stored %>% select(unit, p_position_vs_diff, p_audio_vs_position),
                by = "unit")
print(as.data.frame(chk), digits = 4)
stopifnot(isTRUE(all.equal(chk$p_position_recomputed, chk$p_position_vs_diff)),
          isTRUE(all.equal(chk$p_audio_recomputed, chk$p_audio_vs_position)))
cat("VERIFICATION PASSED: recomputed permutation p values match the stored table.\n")

# ---- dump the plotted null (position contrast) + observed values ------------
null_df <- bind_rows(lapply(res, function(r) {
  tibble(unit = r$unit, perm_id = seq_len(N_PERM),
         null_position_vs_diff = r$null[, "position_vs_diff"])
}))
obs_df <- tibble(
  unit = vapply(res, `[[`, "", "unit"),
  observed_position_vs_diff = vapply(res, function(r) r$obs[["position_vs_diff"]], 0),
  observed_audio_vs_position = vapply(res, function(r) r$obs[["audio_vs_position"]], 0),
  p_position_vs_diff = vapply(res, `[[`, 0, "p_position_vs_diff"),
  p_audio_vs_position = vapply(res, `[[`, 0, "p_audio_vs_position"))

write_csv(null_df, file.path(LOCKDIR, "fig5a_position_contrast_null.csv"))
write_csv(obs_df,  file.path(LOCKDIR, "fig5a_contrast_observed.csv"))
cat("wrote", file.path(LOCKDIR, "fig5a_position_contrast_null.csv"), "\n")
