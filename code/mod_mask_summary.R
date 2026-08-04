# Mask summary — the headline dead-channel mask facts, recomputed from the
# workbook's channel_mask sheet via load_channel_mask() (which recomputes and
# verifies the deposited flags). Backs the mask row of the main results table
# and the Methods mask paragraph.
#
# THE MASK IS THE UNION OF TWO CRITERIA (see load_data.R):
#   (1) dark-fraction — a subject x channel cell whose recording spends
#       >= DARK_FRACTION_CUT of its samples at the instrument dark floor, plus
#       the single rail-saturated cell a dark-floor test cannot see;
#   (2) detector exclusion — detector 8 (T8, channels 20/21) dropped for every
#       subject: dark in 29 of 69 sessions and non-neurovascular when it
#       registers.

source("code/load_data.R")

OUT <- file.path("output")

m  <- load_channel_mask()
cv <- load_dead_channels() %>% distinct(subject_id, channel_index)

# The analysis montage: channels active in at least one subject (channels 6
# and 13 are inactive for all 69, so the montage is 21 channels).
montage_channels <- m %>% group_by(channel_index) %>%
  summarise(any_active = any(active_pair == 1), .groups = "drop") %>%
  filter(any_active) %>% pull(channel_index)

n_sheet    <- nrow(m)
n_montage  <- sum(m$channel_index %in% montage_channels)          # 21 x 69
n_active   <- sum(m$active_pair == 1)
n_dark     <- sum(m$mask_excluded == 1)
n_dark_active   <- sum(m$mask_excluded == 1 & m$active_pair == 1)
n_dark_inactive <- sum(m$mask_excluded == 1 & m$active_pair == 0)
n_union_active  <- sum(m$active_pair == 1 & m$use_in_aggregation == 0)

# Per-ROI usable-channel counts (union mask, active cells only).
usable <- m %>% filter(use_in_aggregation == 1) %>%
  count(subject_id, roi, name = "n_usable") %>%
  complete(subject_id, roi, fill = list(n_usable = 0L))
roi_stat <- function(rr) {
  v <- usable %>% filter(roi == rr) %>% pull(n_usable)
  c(full = sum(v == max(v)), any = sum(v >= 1), empty = sum(v == 0))
}
pfc <- roi_stat("PFC_Frontal"); lt <- roi_stat("Left_Temporal"); rt <- roi_stat("Right_Temporal")

# Against the superseded session-CV mask.
dark_cells  <- m %>% filter(mask_excluded == 1) %>% distinct(subject_id, channel_index)
union_cells <- m %>% filter(active_pair == 1, use_in_aggregation == 0) %>%
  distinct(subject_id, channel_index)
n_newly_dark  <- nrow(anti_join(dark_cells, cv, by = c("subject_id", "channel_index")))
n_lost_dark   <- nrow(anti_join(cv, dark_cells, by = c("subject_id", "channel_index")))
n_newly_union <- nrow(anti_join(union_cells, cv, by = c("subject_id", "channel_index")))
n_lost_union  <- nrow(anti_join(cv, union_cells, by = c("subject_id", "channel_index")))

rail <- m %>% filter(flat_not_dark == 1)

summary_tab <- tribble(
  ~metric, ~value,
  "cells_in_channel_mask_sheet",                    n_sheet,
  "montage_channels",                               length(montage_channels),
  "montage_cells_channels_x_subjects",              n_montage,
  "delivered_active_cells",                         n_active,
  "excluded_cells_dark_fraction_criterion",         n_dark,
  "  of_which_retained_by_active_pair",             n_dark_active,
  "  of_which_flagged_inactive_by_active_pair",     n_dark_inactive,
  "excluded_cells_union_mask_active",               n_union_active,
  "pfc_subjects_all_13_channels_usable",            pfc[["full"]],
  "left_temporal_subjects_any_channel_usable",      lt[["any"]],
  "left_temporal_subjects_no_channel_usable",       lt[["empty"]],
  "right_temporal_subjects_any_channel_usable",     rt[["any"]],
  "right_temporal_subjects_no_channel_usable",      rt[["empty"]],
  "newly_excluded_vs_cv_mask_dark_fraction",        n_newly_dark,
  "cv_cells_lost_dark_fraction_strict_superset",    n_lost_dark,
  "newly_excluded_vs_cv_mask_union_active",         n_newly_union,
  "cv_cells_not_in_union_active_set",               n_lost_union)

cat("mask summary (recomputed from the channel_mask sheet):\n\n")
print(as.data.frame(summary_tab))
cat("\nthe rail-saturated cell a dark-floor test cannot see:",
    rail$subject_id, "channel", rail$channel_index,
    "(detector", rail$detector, ",", rail$roi, ")\n")

write_outcome(summary_tab, file.path(OUT, "mask_summary.csv"))
