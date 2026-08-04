R code for reproducing the statistical analyses, figures and tables for the manuscript "When design timing, preprocessing and channel quality undermine environmental-appraisal fNIRS".

## Requirements

- R 4.5+ (developed and verified on R 4.5.3)

- CRAN packages: `tidyverse`, `readxl`, `lme4`, `lmerTest`, `broom`, `broom.mixed`, `clubSandwich`, `emmeans`, `mclust`, `signal`, `ggplot2`, `ggtext`, `patchwork`

- Install:

  ```r
  install.packages(c("tidyverse", "readxl", "lme4", "lmerTest", "broom",
                     "broom.mixed", "clubSandwich", "emmeans", "mclust",
                     "signal", "ggplot2", "ggtext", "patchwork"))
  ```

- No non-standard hardware is required. The full pipeline (below) is compute-heavy: the cluster bootstraps, permutation audits and cross-validated fits mean a complete run takes on the order of one to a few hours on a normal desktop.

## Data

The analysis reads a single input: the Mendeley Data workbook
`Soundscape_appraisal_fNIRS_prefrontal_temporal_data.xlsx`.

1. Download it from Mendeley Data, **DOI [10.17632/cr2tpv999j](https://doi.org/10.17632/cr2tpv999j)** (CC BY 4.0).
2. Place the `.xlsx` file in the `data/` folder (see `data/README.md`).

The workbook holds the analysis-ready tables: one row per participant and condition, the per-presentation record, the eight adjective ratings, the stimulus inventory, the acquisition event record, the channel-to-region mapping, the per-session channel-quality mask, and per-channel condition concentrations for all three chromophores.

The raw SNIRF recordings and the delivered continuous concentration series are **not** part of the deposit. The small derived tables that the raw-dependent display items need are shipped in `data/intermediates/` (each documented in `data/README.md`), and the Python scripts that produced them from the raw recordings ship in `code/python/` for inspection — they are documentation-grade and are not part of the standard reproduction path (`run_all.R` calls them only if a `data/raw_snirf/` folder is supplied locally).

## File structure

- `run_all.R` — master script: checks the ROI-rebuild gate, runs every analysis module, renders every manuscript figure, and finishes with a verification pass comparing the regenerated display values against the manuscript's quoted numbers.
- `code/load_data.R` — single data entry: reads the workbook sheets (`col_types = "text"`, types re-inferred), converts haemoglobin from mol/L to µmol/L once, reconstructs the union channel mask from the `channel_mask` sheet, and provides the shared loaders and modelling helpers every module uses.
- `code/style.R` — shared publication plot styling.
- `code/mod*.R` — analysis modules, one per declared analysis family of the paper (appraisal structure, confirmatory ROI families, geometry variants, channel scan and adjudication, cluster-robust audit, chromophore triangulation, trait moderation, cycle time course, repetition structure, acoustic drivers, detrended coupling, global adjustment, multivariate pattern, HRF GLM, clip reliability, re-filter gates, corrected-ROI rechecks, correction sweep, the λ chain, the responder-heterogeneity audit, and the union-mask families).
- `code/fig*.R`, `code/sifig1_inference_audit.R` — the manuscript display items (Figures 1–6 and Supplementary Figure S1).
- `code/make_tables.R`, `code/verify_displays.R` — regenerate the numbers behind Tables 1–3 and Supplementary Tables S1–S3, and check every quoted display value.
- `code/python/` — `qa_detector_dropout.py`, `qa_dropout_onset_0803.py`, `refilter_pipeline.py`: the raw-recording QA and re-filter pipelines (documentation-grade; require the raw SNIRF recordings, which are not deposited — contact the authors).
- `data/intermediates/` — small derived tables whose source is not deposited (see `data/README.md`).

## Usage

From the repository root:

```bash
Rscript run_all.R
```

Outputs are written to:

- `output/figures/` — Figures 1–6 and Supplementary Figure S1 (PNG, 600 dpi)
- `output/data_lock/` — the plotted values behind every figure panel (CSV)
- `output/tables/` — the numbers behind Tables 1–3 and Supplementary Tables S1–S3
- `output/<analysis>/` — the full regenerated result tables of each analysis module

The console log reports every family minimum, gate check and anchor comparison. To save it alongside the outputs:

```bash
Rscript run_all.R 2>&1 | tee output/run_log.txt
```

## Notes

- Region-of-interest means must be aggregated through the per-session channel mask: the delivered quality flags pass dead channels. `load_data.R` reconstructs the union mask (absolute dark-floor criterion ≥ 5% of samples, a flatness-off-floor criterion catching one rail-saturated channel, and the exclusion of detector 8 for every participant) from the `channel_mask` sheet and verifies the reconstruction against the deposited `use_in_aggregation` flag on every run; `run_all.R` additionally gates that the masked aggregation reproduces the delivered ROI columns exactly when nothing is excluded.
- Deposited values are stored at Excel's 15-significant-digit precision. Against the working pipeline, regenerated point estimates agree to ≲1e-7 relative; inferential columns (standard errors, degrees of freedom, p and q values) wiggle up to ~1e-7 relative through the mixed-model fits, and the one filter-numerics-sensitive file (`recovery_per_subject.csv`) to 6.3e-7. Differences at that scale are round-trip noise, not discrepancies; every quoted display value is reproduced at its quoted precision (see below).
- `output/` is produced at run time and is safe to delete.

## Verification

`code/verify_displays.R` replays every number quoted in the main text, Tables 1–3, Supplementary Tables S1–S3 and the figure anchors against the regenerated outputs — 196 checks, all passing at the quoted precision — and writes `output/verification/display_number_check.csv`. Quantities whose source is the raw recording QA (not deposited) are recorded as documented constants rather than recomputed, and the caveats below are recorded as flags, never silently passed.

## Reconciliation with the manuscript (2026-08-04)

At the initial release, four deposit-vs-manuscript deviations were on record here. The current manuscript draft agrees with this release on all four; each keeps a footnote or parenthetical documenting the original value:

- **Trait-reactivity family minimum (Table 2 / Table S1)**: age is anonymised to five-year bands in the deposit (the working pipeline used exact integer ages), and the family's minimum-setting test is co-modelled with z(age), so the deposit regenerates the family minimum as q = 0.576 while the exact-age analysis table gives 0.273. The manuscript now quotes 0.576 and footnotes the exact-age 0.273. The companion trait-moderation family (0.277) reproduces exactly, and both values are far-from-significance non-detections — no conclusion changes. Age-conditioned descriptive rows in the covariate scan shift for the same reason.
- **Frozen difference-slope rounding (Table S3 / Section 3.4)**: the frozen prefrontal difference-score comfort slope regenerates as −0.05147 µM; the manuscript quoted −0.052 (carried from a four-decimal lock) and now quotes the deposit-regenerable −0.051.
- **Mask counting (Methods / Table 1 / Section 3.2)**: the union mask excludes 273 of the 1,393 delivered active subject × channel cells; 289 cells meet the criteria over the full montage sheet, 44 of them already pruned by the delivered screen, and the delivered screen had retained 245 of the 289. The manuscript previously said "289 of 1,449" and "every dark channel had been retained"; it now states the 273/1,393 counting with the 289/44 parenthetical and the 245/289 retained count. `output/mask_summary.csv` carries all of these counts.
- **SI survivor label**: the single no-high-pass difference-score survivor (q = 0.035) regenerates as *prefrontal* richness; the supplementary text said right-temporal and now says prefrontal. The q value was correct throughout.

## Remaining caveats

- **Fig 2b binomial**: the quoted P = 2.9e-4 reproduces as P(X ≥ 10) over 14 onsets at the quoted 24.5% rest share (2.86e-4); the audit's per-subject averaging variant (2.7e-4) is the value the DISPLAY_LOCK records — both are documented in the verification ledger.
- The dark-channel SNR ≈ 6,700 and the dRange factor of 13 are documented values from the raw-data QA (raw recordings not deposited); the shipped excerpt trace reproduces the same order (mean/SD ≈ 5.9 × 10³).
- The HbR cycle-bin companion file of the cycle-time-course module is not portable (HbO-only bins are shipped); it is exploratory and appears in no display item.

## License

Code in this repository is released under the MIT License (see `LICENSE`). The input data are archived separately under CC BY 4.0 at Mendeley Data (DOI [10.17632/cr2tpv999j](https://doi.org/10.17632/cr2tpv999j)).
