"""QA — when within a session does the dark-channel light disappear? (2026-08-03)

Trigger. The round-07 reply (delivery 2026-08-04) explains both detector failures
as a loose fibre connector at those optode positions, and adds a testable claim:
*"每次实验开始前都是确保全部探头有信号了才开始的，可能是接头维持不了多久"* — every
optode was verified to carry signal before each session began, so the connector
must have released shortly afterwards. Item 3 of the same reply states that
Cortiview showed per-channel signal quality and that it was checked every session.

Those two statements are falsifiable against material already held. If a
connector held at the operator's pre-session check and released during the
recording, the affected channels must carry real light for some initial stretch
of the raw intensity series and only then fall to the dark floor. If the series
sits at the dark floor from its first sample, the loss preceded the recording's
first sample.

Criterion. `qa_detector_dropout.py` used the session-level coefficient of
variation, which is a *session*-scale statistic: over a 10-s window a live
channel is often smooth enough to look flat by that test, so it cannot time an
onset. This script uses the absolute dark floor instead. DATA_AUDIT_08 section 4
established that a dark cell reads a constant ~0.01312 at both wavelengths; on
the raw samples the band [0.0130, 0.0133] holds 99.99% of samples from
session-dark columns and 4.8% of samples from live ones (the latter being the
live columns' own dropout stretches, which this script then measures). A sample
is dark when **both wavelengths** sit in the band.

Four passes:
  1. Onset — per subject x detector, the first and last non-dark sample, and
     whether any dark session opens with live light.
  2. Partial sessions — for detectors 8 and 9, how much of each session is dark
     and where the dark stretches sit. This is where an intra-session release
     would be visible if it happens at all.
  3. Detector-8 day clustering — the reply says nothing happened between 1 and
     3 October. Permutes session -> date to test per-session independence.
  4. Specificity — the same statistics on the four live temporal channels.

Reads only. Writes `data/processed/qa_dropout_onset_0803.csv` and
`qa_dropout_onset_0803_report.txt`. Nothing frozen is touched.

Usage:  ~/Documents/FullVE/bin/python data/code/qa_dropout_onset_0803.py
"""

from __future__ import annotations

import collections
import csv
import random
import statistics
from pathlib import Path

import h5py
import numpy as np

# NOTE (release copy): documentation-grade. This QA ran over the raw SNIRF
# recordings, which are NOT part of the Mendeley deposit (contact the authors).
# The paths below reproduce the analysis campaign's layout and are kept for the
# record; its per-session determinations are deposited in the workbook's
# channel_mask sheet. The subject crosswalk is private and not shipped.
ROOT = Path(__file__).resolve().parents[2]
RAW_SNIRF = ROOT / "data" / "raw_snirf" / "from-author-2026-07-18" / "01 原始未处理snirf数据"
PROCESSED = ROOT / "data" / "raw_snirf" / "processed"
CROSSWALK = ROOT / "data" / "raw_snirf" / "private" / "subject_crosswalk.csv"

DARK_LO, DARK_HI = 0.0130, 0.0133   # the instrument dark floor, both wavelengths
WL_AGREE = 1e-4                      # a dark sample reads the same at both wavelengths
FLAT_CV = 2e-3                       # a wavelength this flat carries no pulsation
FS = 10.0                            # Hz
DARK_SESSION_FRACTION = 0.99         # a session is "dark" for a detector at this
PARTIAL_FRACTION = 0.05              # a cell is materially contaminated at this
MIN_RUN_S = 5.0                      # ignore runs shorter than this when timing
FAULTY = (9, 8)
MASK_SOURCE = "qa_detector_dropout_0730.csv"   # the determination now in force
SEED = 20260803
N_PERM = 20000


def load_channel_map():
    """(source, detector) -> channel, channel -> (detector, roi), detector -> channels."""
    per_subject = collections.defaultdict(dict)
    meta = {}
    with open(PROCESSED / "tidy_channels.csv", encoding="utf8") as fh:
        for row in csv.DictReader(fh):
            ch = int(row["channel_index"])
            per_subject[row["subject_id"]][(int(row["source"]), int(row["detector"]))] = ch
            meta[ch] = (int(row["detector"]), row["roi"])
    maps = list(per_subject.values())
    assert all(m == maps[0] for m in maps), "montage differs between subjects"
    by_detector = collections.defaultdict(list)
    for ch, (det, _) in sorted(meta.items()):
        by_detector[det].append(ch)
    return maps[0], meta, {d: tuple(c) for d, c in by_detector.items()}


def load_mask():
    """The determination in force: subject x channel -> flagged dark by pass A."""
    flagged = set()
    path = PROCESSED / MASK_SOURCE
    with open(path, encoding="utf8") as fh:
        for row in csv.DictReader(fh):
            if row.get("dead_raw", row.get("flagged", "0")) in ("1", "True", "true"):
                flagged.add((row["subject_id"], int(row["channel_index"])))
    return flagged


def load_crosswalk():
    with open(CROSSWALK, encoding="utf8") as fh:
        return {row["fnirs_id"]: row["subject_id"] for row in csv.DictReader(fh)}


def runs(mask: np.ndarray):
    """Yield (start, stop, value) for maximal constant runs of a boolean mask."""
    if mask.size == 0:
        return
    edges = np.flatnonzero(np.diff(mask.astype(np.int8))) + 1
    bounds = np.concatenate(([0], edges, [mask.size]))
    for a, b in zip(bounds[:-1], bounds[1:]):
        yield int(a), int(b), bool(mask[a])


def scan():
    sd_to_channel, meta, by_detector = load_channel_map()
    channel_to_sd = {v: k for k, v in sd_to_channel.items()}
    crosswalk = load_crosswalk()
    rows = []
    sessions = {}
    grids = {}

    for path in sorted(RAW_SNIRF.glob("*.snirf")):
        sid = crosswalk.get(path.stem, path.stem)
        with h5py.File(path, "r") as h5:
            tags = h5["/nirs/metaDataTags"]

            def tag(name: str) -> str:
                try:
                    return np.array(tags[name]).tobytes().decode("utf8", "ignore").strip("\x00")
                except Exception:
                    return "?"

            sessions[sid] = tag("MeasurementDate")
            # condition markers, read by the stim group's *name* — never its
            # index (builder rule; a positional read produced the withdrawn
            # 07-28 finding). Names "3".."6" are the M3-M6 conditions.
            cond = []
            for group in h5["nirs"].keys():
                if not group.startswith("stim"):
                    continue
                name = np.array(h5[f"nirs/{group}/name"]).tobytes()\
                    .decode("utf8", "ignore").strip("\x00")
                if name in ("3", "4", "5", "6"):
                    block = np.atleast_2d(np.array(h5[f"nirs/{group}/data"]))
                    cond.extend(block[:, 0].tolist())
            grids[sid] = sorted(cond)
            data = h5["/nirs/data1"]
            series = np.array(data["dataTimeSeries"])
            columns = collections.defaultdict(list)
            for key in data.keys():
                if not key.startswith("measurementList"):
                    continue
                entry = data[key]
                idx = int(key.replace("measurementList", "")) - 1
                pair = (int(np.array(entry["sourceIndex"]).ravel()[0]),
                        int(np.array(entry["detectorIndex"]).ravel()[0]))
                columns[pair].append(idx)

        n = series.shape[0]
        for ch in sorted(meta):
            detector, roi = meta[ch]
            cols = sorted(columns[channel_to_sd[ch]])
            wl = [series[:, c] for c in cols]
            in_band = np.logical_and.reduce([(x >= DARK_LO) & (x <= DARK_HI) for x in wl])
            # a dark sample also reads the *same* value at both wavelengths;
            # tissue attenuates 760 and 850 nm differently, so this separates a
            # dark floor from a dim but real channel
            agree = np.abs(wl[0] - wl[1]) < WL_AGREE
            dark = in_band & agree
            # A channel can also be flat somewhere OTHER than the dark floor —
            # a wavelength pinned at the converter's upper rail carries no
            # cardiac pulsation either, and no dark-floor test can see it. One
            # such cell exists (P10 channel 12, wavelength 2 at exactly 1.0),
            # and it is the one cell the superseded session-CV mask caught that
            # this criterion does not. Recorded per cell so the mask can be the
            # union of the two rather than a replacement of one by the other.
            flat_not_dark = int(any(
                (x.std() / abs(x.mean())) < FLAT_CV and
                not (DARK_LO <= x.mean() <= DARK_HI) for x in wl))
            live_idx = np.flatnonzero(~dark)
            frac = float(dark.mean())
            lead = int(np.argmax(dark)) if dark.any() else n
            long_dark = [(a, b) for a, b, v in runs(dark)
                         if v and (b - a) >= MIN_RUN_S * FS]
            rows.append({
                "subject_id": sid,
                "date": sessions[sid],
                "detector": detector,
                "channel": ch,
                "roi": roi,
                "n_samples": n,
                "frac_dark": f"{frac:.4f}",
                "frac_in_band_only": f"{float(in_band.mean()) - frac:.4f}",
                "flat_not_dark": flat_not_dark,
                "session_dark": int(frac >= DARK_SESSION_FRACTION),
                "partial_dark": int(PARTIAL_FRACTION <= frac < DARK_SESSION_FRACTION),
                "leading_live_seconds": f"{lead / FS:.1f}",
                "last_live_sample_s": f"{live_idx[-1] / FS:.1f}" if live_idx.size else "",
                "n_dark_runs_ge5s": len(long_dark),
                "first_dark_run_start_s": f"{long_dark[0][0] / FS:.1f}" if long_dark else "",
            })
    return rows, sessions, meta, by_detector, grids


def day_clustering_test(rows, sessions, detector, by_detector):
    per_subject = {}
    for r in rows:
        if r["detector"] == detector and r["channel"] == by_detector[detector][0]:
            per_subject[r["subject_id"]] = int(r["session_dark"])
    subjects = sorted(per_subject)
    dark = [per_subject[s] for s in subjects]
    dates = [sessions[s] for s in subjects]
    p = sum(dark) / len(dark)
    by_date = collections.defaultdict(lambda: [0, 0])
    for d, x in zip(dates, dark):
        by_date[d][0] += x
        by_date[d][1] += 1

    def chi(assignment):
        acc = collections.defaultdict(lambda: [0, 0])
        for d, x in zip(dates, assignment):
            acc[d][0] += x
            acc[d][1] += 1
        return sum((k - m * p) ** 2 / (m * p * (1 - p)) for k, m in acc.values() if m)

    observed = chi(dark)
    rng = random.Random(SEED)
    pool = list(dark)
    hits = sum(1 for _ in range(N_PERM) if (rng.shuffle(pool) or chi(pool)) >= observed)
    return by_date, p, observed, (hits + 1) / (N_PERM + 1)


def main():
    rows, sessions, meta, by_detector, grids = scan()
    out_csv = PROCESSED / "qa_dropout_onset_0803.csv"
    with open(out_csv, "w", newline="", encoding="utf8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    flagged = load_mask()

    lines = []
    add = lines.append
    add("QA — dark-channel onset and its extent across the montage (2026-08-03)")
    add("=" * 74)
    add(f"dark sample: both wavelengths in [{DARK_LO}, {DARK_HI}] AND agreeing to "
        f"<{WL_AGREE}; session dark at >={DARK_SESSION_FRACTION:.0%} of samples, "
        f"partial at >={PARTIAL_FRACTION:.0%}")
    add("")

    add("1. The montage, by detector")
    add("-" * 74)
    add("  det  channels          roi              wholly dark   partial   clean")
    for det in sorted(by_detector):
        chans = by_detector[det]
        ch = chans[0]
        sub = [r for r in rows if r["channel"] == ch]
        whole = sum(int(r["session_dark"]) for r in sub)
        part = sum(int(r["partial_dark"]) for r in sub)
        add(f"  {det:3d}  {str(chans):16s}  {meta[ch][1]:14s} "
            f"{whole:8d}   {part:8d}  {len(sub) - whole - part:6d}")
    add("")

    add("2. Onset — does any wholly dark session open with live light?")
    add("-" * 74)
    for det in FAULTY:
        for ch in by_detector[det]:
            sub = [r for r in rows if r["channel"] == ch and int(r["session_dark"])]
            lead = [r for r in sub if float(r["leading_live_seconds"]) >= MIN_RUN_S]
            add(f"  detector {det}, channel {ch}: {len(sub)} wholly dark sessions; "
                f"{len(lead)} open with >={MIN_RUN_S:.0f} s of live light")
    add("")

    add("3. Partial sessions — is a release ever seen inside a recording?")
    add("-" * 74)
    for det in sorted(by_detector):
        ch = by_detector[det][0]
        sub = [r for r in rows if r["channel"] == ch]
        partial = [r for r in sub if int(r["partial_dark"])]
        if not partial:
            continue
        mid = [r for r in partial if float(r["leading_live_seconds"]) >= 30]
        add(f"  detector {det} ({meta[ch][1]}): {len(partial)} partial sessions, "
            f"{len(mid)} of them live for >=30 s before the first dark sample")
        for r in sorted(partial, key=lambda r: -float(r["frac_dark"])):
            add(f"      {r['subject_id']} {r['date']}: dark {float(r['frac_dark']):.1%}, "
                f"live until {r['leading_live_seconds']} s, "
                f"{r['n_dark_runs_ge5s']} dark runs >=5 s")
    add("")

    add("4. Detector 8 — by-date pattern against per-session independence")
    add("-" * 74)
    by_date, p, observed, pval = day_clustering_test(rows, sessions, 8, by_detector)
    add(f"  overall wholly-dark rate p = {p:.3f}")
    add("    date         wholly dark   partial   mean dark fraction")
    for d in sorted(by_date):
        k, m = by_date[d]
        day_rows = [r for r in rows
                    if r["channel"] == by_detector[8][0] and r["date"] == d]
        part = sum(int(r["partial_dark"]) for r in day_rows)
        mean_frac = statistics.mean(float(r["frac_dark"]) for r in day_rows)
        add(f"    {d}   {k:2d} / {m:2d}        {part:3d}       {mean_frac:.3f}")
    add(f"  statistic {observed:.2f}; permutation p = {pval:.5f} "
        f"({N_PERM} permutations, session->date shuffled)")
    add("")

    add("5. What the mask in force does with this")
    add("-" * 74)
    add(f"  mask source: {MASK_SOURCE} (session-level CV criterion)")
    missed = [r for r in rows
              if int(r["partial_dark"]) or (int(r["session_dark"])
                                            and (r["subject_id"], r["channel"]) not in flagged)]
    missed = [r for r in missed if (r["subject_id"], r["channel"]) not in flagged]
    by_roi = collections.Counter(r["roi"] for r in missed)
    add(f"  subject x channel cells with >={PARTIAL_FRACTION:.0%} dark time that the "
        f"mask retains: {len(missed)}")
    for roi, k in by_roi.most_common():
        add(f"      {roi:16s} {k}")
    add("  the worst of them:")
    for r in sorted(missed, key=lambda r: -float(r["frac_dark"]))[:20]:
        add(f"      {r['subject_id']} ch{r['channel']:2d} (det {r['detector']:2d}, "
            f"{r['roi']:14s}) dark {float(r['frac_dark']):.1%}")
    add("")

    add("6. Where in the session do the dropouts start?")
    add("-" * 74)
    add("  Each session runs 12 condition events: three sequences of four, ~90 s")
    add("  apart within a sequence, with a rest interval between sequences. Tests")
    add("  whether dropout onsets fall in the rest intervals more often than their")
    add("  share of session time.")
    seen, in_gap, out_gap, share = set(), [], [], []
    for r in sorted(rows, key=lambda r: (r["subject_id"], r["detector"])):
        t = r["first_dark_run_start_s"]
        if not t or float(t) <= 0 or float(r["frac_dark"]) < PARTIAL_FRACTION:
            continue
        key = (r["subject_id"], r["detector"])
        if key in seen:
            continue
        seen.add(key)
        t = float(t)
        grid = grids[r["subject_id"]]
        if len(grid) != 12:
            continue
        gaps = [(grid[3], grid[4]), (grid[7], grid[8])]
        duration = r["n_samples"] and int(r["n_samples"]) / FS
        share.append(sum(b - a for a, b in gaps) / duration)
        (in_gap if any(a <= t <= b for a, b in gaps) else out_gap).append(
            (r["subject_id"], r["detector"], t))
    n = len(in_gap) + len(out_gap)
    p_share = statistics.mean(share) if share else 0.0
    add(f"  rest intervals occupy {p_share:.1%} of the recording")
    add(f"  dropout onsets inside a rest interval: {len(in_gap)} / {n} "
        f"(expected {n * p_share:.1f})")
    # exact binomial tail
    from math import comb
    tail = sum(comb(n, k) * p_share ** k * (1 - p_share) ** (n - k)
               for k in range(len(in_gap), n + 1))
    add(f"  binomial P(X >= {len(in_gap)}) = {tail:.2e}")
    add("  inside:  " + ", ".join(f"{s} det{d} {t:.0f}s" for s, d, t in in_gap))
    add("  outside: " + ", ".join(f"{s} det{d} {t:.0f}s" for s, d, t in out_gap))
    add("")
    add(f"per-cell output: {out_csv.relative_to(ROOT)}")

    text = "\n".join(lines)
    (PROCESSED / "qa_dropout_onset_0803_report.txt").write_text(text + "\n", encoding="utf8")
    print(text)


if __name__ == "__main__":
    main()
