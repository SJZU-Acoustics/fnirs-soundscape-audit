# data/

Place the Mendeley Data workbook here:

- `Soundscape_appraisal_fNIRS_prefrontal_temporal_data.xlsx` — DOI [10.17632/cr2tpv999j](https://doi.org/10.17632/cr2tpv999j) (CC BY 4.0). Sheets: `analysis_table`, `block_level`, `subjects`, `subjective_ratings`, `stimuli`, `events`, `channels`, `channel_mask`, `channel_condition`, plus the variable and value dictionaries.

The workbook is not redistributed with this repository (it is the data deposit itself); download it from Mendeley Data.

## intermediates/

Small derived tables shipped with the repository. Each exists because its source — the raw SNIRF recordings or the delivered continuous concentration series — is **not** part of the Mendeley deposit. The code that produced them from those sources is in `code/python/` (raw-recording QA and re-filter, requiring the raw recordings) or was run against the delivered series during the analysis campaign; the modules in `code/` recompute every displayed statistic *from* these tables, so no display item is taken on trust.

| File | Provenance | Feeds |
|---|---|---|
| `qa_detector_dropout_0730.csv` | Session-CV channel screen on the raw recordings (`code/python/qa_detector_dropout.py`); the superseded mask the corrected-ROI rechecks and the Fig 2c mirror-pair counts were declared against | `mod20`, `mod23`, `fig2` |
| `tidy_block_metrics.csv` | Participant × presentation × channel × chromophore window statistics, averaged from the delivered concentration series (the deposit carries only the participant × condition layer) | `mod23` |
| `condition_bin_values_hbo.csv` | Participant × condition 5-s event-referenced HbO bins, epoched from the delivered concentration series | `mod11`, `mod18`, `fig5` |
| `clip_acoustic_features.csv` | Per-clip acoustic features computed from the stimulus WAV files (not redistributed) | `mod13` |
| `refilter_condition_metrics.csv` | Condition-level ROI/head HbO window means of the re-filtered series (four high-pass settings × two motion-correction arms), from `code/python/refilter_pipeline.py` on the raw recordings; µmol/L | `mod19`, `mod21` |
| `refilter_event_metrics_roi.csv` | Per-presentation ROI window means of the same re-filtered series (the cycle-position contrast's layer) | `mod19`, `mod21` |
| `refilter_channel10_metrics.csv` | Channel-10 condition-level values of the same re-filtered series | `mod19` |
| `gate_a_pruning.csv`, `gate_b_reproduction.csv` | Re-filter gates: channel pruning agreement and per-participant waveform correlations between the rebuilt and delivered series | `mod19` |
| `a17_glm_subject_betas.csv` | Per-participant canonical-HRF GLM modulator betas, fitted on the delivered concentration series | `mod17` |
| `fig1a_*.csv`, `fig1c_refilter.csv` | Fig 1 panels a (cycle anatomy, recomputed from the delivered series on the frozen onsets) and c (re-filter contrasts) | `fig1` |
| `fig2a_*.csv`, `fig2b_*.csv` | Fig 2 panels a (raw intensity traces of the three failure kinds) and b (per-session dropout calendar, rest intervals, release onsets), from `code/python/qa_detector_dropout.py` / `qa_dropout_onset_0803.py` on the raw recordings | `fig2` |

Fig 1b (the filter amplitude response) needs no intermediate: it is recomputed from the filter definition in `code/mod_fig1b_filter_gain.R`. Fig 2c–d are recomputed from the workbook's `channels` and `channel_mask` sheets.
