"""
Re-filter pipeline — P24 stage-1 rebuild from raw SNIRF (authorised 2026-07-30).

WHY THIS EXISTS
---------------
The collector's preprocessing high-passes at 0.01 Hz (cutoff period 100 s) while
the presentation cycle is 90 s (0.0111 Hz). The design's fundamental therefore
sits on the filter's shoulder; conventional block-design guidance wants the
cutoff period at least twice the task cycle (here >= 180 s, i.e. <= 0.0056 Hz).
A03/A11/A15/A17/A18 all converge on that interaction, and A17 showed no
estimator downstream can escape it because the artefact is at the design's own
frequency. Only re-deriving the concentration series with a different high-pass
can test it, and that cannot be done from the delivered (already filtered)
series because the sub-0.01 Hz content is gone.

WHAT IT DOES
------------
Reimplements the collector's declared Homer3 pipeline in Python and re-runs it
from the raw SNIRF at several high-pass settings:

    PruneChannels -> Intensity2OD -> MotionCorrectWavelet -> BandpassFilt -> OD2Conc

Everything except the high-pass corner is held fixed, and the OD stage
(including wavelet motion correction) is computed ONCE per subject and shared
across the filter settings. The filter contrast is therefore exact and
internally valid **regardless of whether this reimplementation matches Homer3
bit-for-bit** — the only thing that differs between branches is hpf.

TWO GATES
---------
A. Pruning gate — does this pipeline reproduce the delivered `pairActive`?
B. Reproduction gate — at hpf = 0.01 (the collector's own setting), how closely
   does this pipeline reproduce the delivered concentration series `dcData`?
   This is what calibrates how far conclusions transfer to the frozen tables.
   It has never been tested: stage 1 proved our *epoching* is the collector's,
   never that their stated preprocessing reproduces their delivered series.

WHAT IT DOES NOT DO
-------------------
It does not touch, overwrite or re-freeze anything in data/processed/. Outputs
go to data/processed_refilter/ and are explicitly an alternative, not a
replacement. The frozen tables remain the campaign's reference.

Reference material, all from the 2026-07-28 delivery (declared rebuilds):
  05 batch_preprocess_fnirs_homer.m   — the parameter set and pipeline order
  01 hmrR_OD2Conc.m                   — the modified Beer-Lambert step
  02 GetExtinctions.m                 — the Prahl extinction table (case 0)

Run:  ~/Documents/FullVE/bin/python data/code/refilter_pipeline.py
"""

from __future__ import annotations

import re
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pywt
from scipy.signal import butter, filtfilt

# NOTE (release copy): documentation-grade. This pipeline ran over the raw
# SNIRF recordings and the delivered MATLAB exports, which are NOT part of the
# Mendeley deposit (contact the authors). The paths below reproduce the
# analysis campaign's layout and are kept for the record; its outputs are the
# shipped intermediates refilter_condition_metrics.csv, gate_a_pruning.csv and
# gate_b_reproduction.csv (see data/README.md).
PROJECT = Path(__file__).resolve().parents[2]
RAW_DIR = PROJECT / "data/raw_snirf/from-author-2026-07-18/01 原始未处理snirf数据"
MAT_DIR = PROJECT / "data/raw_snirf/from-author-2026-07-18/03 预处理后的 MATLAB 数据"
REF_DIR = PROJECT / "data/raw_snirf/from-author-2026-07-28"
PROC = PROJECT / "data/raw_snirf/processed"
OUT = PROJECT / "data/raw_snirf/processed_refilter"
OUT.mkdir(parents=True, exist_ok=True)

# ---- the collector's declared parameters (05 batch_preprocess_fnirs_homer.m) ----
D_RANGE = (1e-3, 1e9)
SNR_THRESH = 2.0
SD_RANGE = (0.0, 45.0)
PPF = np.array([6.0, 6.0])
WAVELET_IQR = 1.5
LPF = 0.10
HPF_REFERENCE = 0.01

# the experiment: high-pass settings to compare. lpf is held at 0.10 throughout.
#   0.01   the collector's setting (cutoff period 100 s vs a 90 s cycle)
#   0.005  cutoff period 200 s  — clears the >=2x task-cycle rule of thumb
#   0.002  cutoff period 500 s  — comfortably clear
#   0.0    no high-pass at all; drift left in, handled by the window difference
HPF_SETTINGS = [0.01, 0.005, 0.002, 0.0]

# Motion correction is the one step that cannot be reproduced from the delivered
# material: Homer3's hmrR_MotionCorrectWavelet padding and level choice are not
# recoverable (the original scripts were not retained), and empirically neither
# this implementation nor any wavelet/threshold variant tried reproduces the
# delivered series better than omitting the step (median r 0.82 with, 0.87
# without). Rather than pick one and hope, the filter experiment is run under
# BOTH, so the unverifiable step cancels out of the comparison: if the two arms
# agree about the filter, the conclusion does not depend on it.
MC_SETTINGS = ["none", "wavelet"]

WINDOWS = {
    "baseline_m10_0": (-10.0, 0.0),
    "stim_m0_30": (0.0, 30.0),
    "response_m5_30": (5.0, 30.0),
}


# --------------------------------------------------------------- extinctions
def load_extinctions(wavelengths: np.ndarray) -> np.ndarray:
    """Prahl table, case 0 of GetExtinctions.m, x2.303, then /10 for /cm -> /mm.

    Returns e (n_wavelength x 2) with columns [HbO, HbR], matching
    `e = GetExtinctions(Lambda); e = e(:,1:2) / 10;` in 01 hmrR_OD2Conc.m.
    """
    lines = (REF_DIR / "02 GetExtinctions.m").read_text(
        encoding="utf-8", errors="replace"
    ).split("\n")
    start = next(i for i, l in enumerate(lines) if "vLambdaHbOHb = [" in l)
    rows = []
    for line in lines[start + 1:]:
        m = re.match(r"\s*([\d.]+)\s+([\d.]+)\s+([\d.]+);", line)
        if m:
            rows.append([float(x) for x in m.groups()])
        elif "];" in line and rows:
            break
    t = np.array(rows)
    e = np.column_stack([
        np.interp(wavelengths, t[:, 0], t[:, 1]) * 2.303,
        np.interp(wavelengths, t[:, 0], t[:, 2]) * 2.303,
    ])
    return e / 10.0


# ------------------------------------------------------------- SNIRF reading
def read_snirf(path: Path):
    with h5py.File(path, "r") as f:
        d = f["nirs/data1"]
        y = d["dataTimeSeries"][()]                       # (nT, 46)
        t = np.asarray(d["time"][()]).ravel()
        n_ml = sum(1 for k in d if k.startswith("measurementList"))
        ml = np.array([
            [int(d[f"measurementList{i}"][k][()])
             for k in ("sourceIndex", "detectorIndex", "wavelengthIndex")]
            for i in range(1, n_ml + 1)
        ])
        p = f["nirs/probe"]
        wl = np.asarray(p["wavelengths"][()]).ravel()
        sp = p["sourcePos3D"][()]
        dp = p["detectorPos3D"][()]
    return y, t, ml, wl, sp, dp


# ------------------------------------------------------------- pipeline steps
def prune_channels(y, ml, rho_pair):
    """hmrR_PruneChannels: dRange on the mean, SNR = mean/std, S-D distance."""
    dm = np.abs(y).mean(axis=0)
    ds = y.std(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = np.where(ds > 0, dm / ds, 0.0)
    ok_meas = (dm >= D_RANGE[0]) & (dm <= D_RANGE[1]) & (snr >= SNR_THRESH)
    n_pair = ml[ml[:, 2] == 1].shape[0]
    ok_pair = np.ones(n_pair, dtype=bool)
    for k in range(n_pair):
        cols = np.where((ml[:, 0] == ml[k, 0]) & (ml[:, 1] == ml[k, 1]))[0]
        ok_pair[k] = ok_meas[cols].all()
    ok_pair &= (rho_pair >= SD_RANGE[0]) & (rho_pair <= SD_RANGE[1])
    return ok_pair


def intensity_to_od(y):
    """hmrR_Intensity2OD: dod = -log(|d| / mean(|d|))."""
    dm = np.abs(y).mean(axis=0)
    dm = np.where(dm == 0, 1.0, dm)
    return -np.log(np.abs(y) / dm)


def motion_correct_wavelet(dod, iqr_factor=WAVELET_IQR, wavelet="db2"):
    """Molavi & Dumont (2012) wavelet motion correction, as Homer3 applies it.

    Per channel: discrete wavelet decomposition (periodised Daubechies-2, the
    4-coefficient Daubechies filter Homer's WaveLab call uses); at each detail
    level, coefficients outside [Q1 - k*IQR, Q3 + k*IQR] are treated as motion
    and set to zero; inverse transform.

    Homer3's exact padding and level choice are not recoverable from the
    delivered material (the original scripts were not retained). This is
    therefore a faithful implementation of the published method rather than a
    bit-exact port — which is why the filter contrast is run as a within-pipeline
    comparison sharing this same output.
    """
    out = np.empty_like(dod)
    n = dod.shape[0]
    level = pywt.dwt_max_level(n, pywt.Wavelet(wavelet).dec_len)
    for c in range(dod.shape[1]):
        coeffs = pywt.wavedec(dod[:, c], wavelet, mode="periodization", level=level)
        new = [coeffs[0]]
        for det in coeffs[1:]:
            q1, q3 = np.percentile(det, [25, 75])
            iqr = q3 - q1
            lo, hi = q1 - iqr_factor * iqr, q3 + iqr_factor * iqr
            new.append(np.where((det < lo) | (det > hi), 0.0, det))
        rec = pywt.waverec(new, wavelet, mode="periodization")
        out[:, c] = rec[:n]
    return out


def bandpass(y, fs, hpf, lpf):
    """hmrBandpassFilt: 3rd-order Butterworth low-pass then high-pass, filtfilt.

    Homer applies the two stages separately rather than designing one band-pass;
    hpf = 0 skips the high-pass entirely.
    """
    out = y
    if lpf and lpf > 0:
        b, a = butter(3, lpf * 2.0 / fs, btype="low")
        out = filtfilt(b, a, out, axis=0)
    if hpf and hpf > 0:
        b, a = butter(3, hpf * 2.0 / fs, btype="high")
        out = filtfilt(b, a, out, axis=0)
    return out


def pair_columns(ml):
    """For each wl1 pair, the (wl1_col, wl2_col) indices into the data matrix.

    Source-detector distance uses the 3D probe positions: 2D positions give the
    same waveform but a median SD ratio of 1.65 against the delivered series
    versus 1.20 for 3D, so 3D is what the collector used.
    """
    idx1 = np.where(ml[:, 2] == 1)[0]
    cols = []
    for k in idx1:
        cols.append([
            np.where((ml[:, 0] == ml[k, 0]) & (ml[:, 1] == ml[k, 1]) &
                     (ml[:, 2] == w))[0][0]
            for w in (1, 2)
        ])
    return np.array(cols)


def od_to_conc(dod, cols, rho_pair, e):
    """01 hmrR_OD2Conc.m with ppf = [6, 6] (so the rho*ppf branch is taken)."""
    einv = np.linalg.inv(e.T @ e) @ e.T          # 2 x n_wavelength
    n_pair = cols.shape[0]
    nT = dod.shape[0]
    hbo = np.empty((nT, n_pair))
    hbr = np.empty((nT, n_pair))
    for k in range(n_pair):
        scaled = dod[:, cols[k]] / (rho_pair[k] * PPF)   # (nT, 2)
        conc = scaled @ einv.T                           # (nT, 2) -> [HbO, HbR]
        hbo[:, k], hbr[:, k] = conc[:, 0], conc[:, 1]
    return hbo, hbr


# ------------------------------------------------------------------ epoching
def window_means(series, t, onsets):
    """Mean of `series` in each declared window around each onset.

    series: (nT, n_pair). Returns dict window -> (n_event, n_pair).
    """
    out = {}
    for name, (lo, hi) in WINDOWS.items():
        vals = np.full((len(onsets), series.shape[1]), np.nan)
        for i, o in enumerate(onsets):
            m = (t >= o + lo) & (t < o + hi)
            if m.any():
                vals[i] = series[m].mean(axis=0)
        out[name] = vals
    return out


# ---------------------------------------------------------------------- main
def main():
    t_start = time.time()
    crosswalk = pd.read_csv(PROJECT / "data/private/subject_crosswalk.csv")
    events = pd.read_csv(PROC / "tidy_events.csv")
    events = events[events["onset_sec"].notna()]
    frozen_ch = pd.read_csv(PROC / "tidy_channels.csv")

    prune_rows, gate_rows, metric_rows = [], [], []

    for i, r in crosswalk.iterrows():
        pid, sid = r["subject_id"], r["fnirs_id"]
        y, t, ml, wl, sp, dp = read_snirf(RAW_DIR / f"{sid}.snirf")
        fs = 1.0 / np.median(np.diff(t))
        e = load_extinctions(wl)

        idx1 = np.where(ml[:, 2] == 1)[0]
        rho = np.array([
            np.linalg.norm(sp[ml[k, 0] - 1] - dp[ml[k, 1] - 1]) for k in idx1
        ])

        # ---- gate A: pruning
        ok_pair = prune_channels(y, ml, rho)
        with h5py.File(MAT_DIR / f"{sid}.mat", "r") as f:
            delivered_active = np.asarray(f["pairActive"][()]).ravel().astype(int)
            dc_delivered = np.asarray(f["dcData"][()])          # (69, nT) or (nT, 69)
        if dc_delivered.shape[0] != 69:
            dc_delivered = dc_delivered.T
        prune_rows.append({
            "subject_id": pid,
            "n_active_mine": int(ok_pair.sum()),
            "n_active_delivered": int(delivered_active.sum()),
            "n_disagreements": int((ok_pair.astype(int) != delivered_active).sum()),
        })

        # ---- OD once; the two motion-correction arms once. Both are shared
        # across every filter setting, so within an arm the ONLY thing that
        # differs between branches is the high-pass corner.
        cols = pair_columns(ml)
        dod_raw = intensity_to_od(y)
        dod_by_mc = {"none": dod_raw, "wavelet": motion_correct_wavelet(dod_raw)}

        ev = events[events["subject_id"] == pid].sort_values("event_index")
        onsets = ev["onset_sec"].to_numpy()
        markers = ev["marker"].to_numpy()
        eidxs = ev["event_index"].to_numpy()
        d_hbo = dc_delivered[0::3, :].T                        # (nT, 23)
        d_hbr = dc_delivered[1::3, :].T

        for mc in MC_SETTINGS:
            for hpf in HPF_SETTINGS:
                dod_f = bandpass(dod_by_mc[mc], fs, hpf, LPF)
                hbo, hbr = od_to_conc(dod_f, cols, rho, e)

                # ---- gate B: at the collector's own setting, reproduce dcData
                if hpf == HPF_REFERENCE:
                    for k in range(hbo.shape[1]):
                        if not ok_pair[k]:
                            continue
                        for nm, mine, theirs in (("HbO", hbo[:, k], d_hbo[:, k]),
                                                 ("HbR", hbr[:, k], d_hbr[:, k])):
                            sd = theirs.std()
                            gate_rows.append({
                                "subject_id": pid, "motion_correction": mc,
                                "channel_index": k + 1, "chromophore": nm,
                                "r": float(np.corrcoef(mine, theirs)[0, 1]),
                                "rmse_over_sd": float(
                                    np.sqrt(((mine - theirs) ** 2).mean()) / sd
                                ) if sd > 0 else np.nan,
                                "sd_mine": float(mine.std()), "sd_theirs": float(sd),
                            })

                # ---- window statistics on the frozen onsets
                for nm, series in (("HbO", hbo), ("HbR", hbr)):
                    wm = window_means(series, t, onsets)
                    base, stim, resp = (wm["baseline_m10_0"], wm["stim_m0_30"],
                                        wm["response_m5_30"])
                    for j in range(len(eidxs)):
                        for k in range(series.shape[1]):
                            if not ok_pair[k]:
                                continue
                            metric_rows.append({
                                "subject_id": pid, "event_index": int(eidxs[j]),
                                "marker": markers[j], "channel_index": k + 1,
                                "chromophore": nm, "motion_correction": mc,
                                "hpf": hpf,
                                "baseline_m10_0": base[j, k],
                                "stim_m0_30": stim[j, k],
                                "response_m5_30": resp[j, k],
                                "response_minus_baseline_m5_30":
                                    resp[j, k] - base[j, k],
                            })

        if (i + 1) % 10 == 0:
            print(f"  {i + 1}/{len(crosswalk)} subjects "
                  f"({time.time() - t_start:.0f}s)", flush=True)

    pd.DataFrame(prune_rows).to_csv(OUT / "gate_a_pruning.csv", index=False)
    gate = pd.DataFrame(gate_rows)
    gate.to_csv(OUT / "gate_b_reproduction.csv", index=False)
    met = pd.DataFrame(metric_rows)
    met.to_csv(OUT / "refilter_block_metrics.csv", index=False)

    print("\n=== GATE A — pruning ===")
    pa = pd.DataFrame(prune_rows)
    print(f"  subjects with identical pairActive: "
          f"{(pa.n_disagreements == 0).sum()} / {len(pa)}")
    print(f"  total channel disagreements: {pa.n_disagreements.sum()}")

    print("\n=== GATE B — reproduction of the delivered series at hpf = 0.01 ===")
    for (mc, nm), g in gate.groupby(["motion_correction", "chromophore"]):
        print(f"  mc={mc:8s} {nm}: median r = {g.r.median():.4f}   "
              f"10th pct r = {g.r.quantile(0.10):.4f}   "
              f"median RMSE/SD = {g.rmse_over_sd.median():.4f}   "
              f"median SD ratio = {(g.sd_mine / g.sd_theirs).median():.4f}")
    print("  NOTE: this is a waveform gate. It says how well this Python stand-in")
    print("  reproduces Homer3, not whether the filter experiment is valid — that")
    print("  rests on the functional gate in refilter_analysis.R.")

    print(f"\nblock metrics written: {len(met):,} rows, "
          f"{met.hpf.nunique()} filter settings")
    print(f"total {time.time() - t_start:.0f}s")


if __name__ == "__main__":
    main()
