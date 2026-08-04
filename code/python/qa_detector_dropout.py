"""QA — detector-level dropout in the temporal montage (2026-07-30).

Trigger. A03 section 3b established that channels 18 and 19 are dead in all 69
subjects and read that as two dead channels. This script asks the prior question:
*which optode* they share, and whether the same failure reaches any other channel.

Two independent passes, deliberately kept separate:

  Pass A — raw SNIRF (69 files, delivery 2026-07-18). Per subject x channel, the
  coefficient of variation of the raw intensity time series, min over the two
  wavelengths. Light that has traversed tissue always carries a cardiac pulsation,
  so a channel with no relative fluctuation received no tissue light. This pass is
  the physical criterion and owns the diagnosis.

  Pass B — the frozen block layer (`tidy_block_metrics.csv`, digest-verified by
  the analysis helpers). Per subject x channel, the SD of the HbO response metric
  across that subject's blocks, relative to the campaign-wide median. This pass is
  what stage 2 can compute without leaving the freeze, and it is the criterion
  `helpers.R` implements.

The passes are compared, not merged: pass B is usable only because pass A shows
it recovers the same subject x channel set.

Reads only. Writes `data/processed/qa_detector_dropout_0730.csv` (one row per
subject x channel) and `qa_detector_dropout_0730_report.txt`. Nothing frozen is
touched; the outputs are QA artefacts and are not part of the freeze.

Usage:  ~/Documents/FullVE/bin/python data/code/qa_detector_dropout.py
"""

from __future__ import annotations

import collections
import csv
import statistics
from pathlib import Path

import h5py
import numpy as np

# NOTE (release copy): documentation-grade. This QA ran over the raw SNIRF
# recordings, which are NOT part of the Mendeley deposit (contact the authors).
# The paths below reproduce the analysis campaign's layout and are kept for the
# record; its output is the shipped intermediate qa_detector_dropout_0730.csv
# (see data/README.md). The subject crosswalk is private and not shipped.
ROOT = Path(__file__).resolve().parents[2]
RAW_SNIRF = ROOT / "data" / "raw_snirf" / "from-author-2026-07-18" / "01 原始未处理snirf数据"
PROCESSED = ROOT / "data" / "raw_snirf" / "processed"
CROSSWALK = ROOT / "data" / "raw_snirf" / "private" / "subject_crosswalk.csv"

# Pass A: a channel whose raw intensity varies by less than this fraction of its
# own mean carried no cardiac pulsation. Live channels sit at 0.017-0.077 and the
# dead cluster at ~1e-4, so the threshold falls in a two-order-of-magnitude gap;
# its sensitivity is reported rather than assumed.
CV_FLAT = 0.002
CV_GRID = (0.0005, 0.001, 0.002, 0.005, 0.010)

# Pass B: an equivalent floor on the frozen layer, as a fraction of the median
# per-subject-channel response SD.
FROZEN_FLAT_FRACTION = 0.01

TEMPORAL_ROIS = {"Left_Temporal": (16, 17, 18, 19), "Right_Temporal": (20, 21, 22, 23)}


def load_channel_map() -> tuple[dict, dict, dict]:
    """(source, detector) -> channel_index, plus per-subject roi and active flags.

    The montage is identical for all 69 subjects; that is asserted, not assumed.
    """
    per_subject: dict[str, dict] = collections.defaultdict(dict)
    active: dict[str, dict] = collections.defaultdict(dict)
    meta: dict[int, tuple] = {}
    with open(PROCESSED / "tidy_channels.csv", encoding="utf8") as fh:
        for row in csv.DictReader(fh):
            ch = int(row["channel_index"])
            key = (int(row["source"]), int(row["detector"]))
            per_subject[row["subject_id"]][key] = ch
            active[row["subject_id"]][ch] = int(row["active_pair"])
            meta[ch] = (key[0], key[1], row["channel_name"], row["roi"],
                        float(row["sd_distance_mm"]))
    maps = list(per_subject.values())
    assert all(m == maps[0] for m in maps), "montage differs between subjects"
    return maps[0], active, meta


def load_crosswalk() -> dict[str, str]:
    """Original fNIRS filename stem -> anonymous subject_id."""
    out = {}
    with open(CROSSWALK, encoding="utf8") as fh:
        for row in csv.DictReader(fh):
            out[row["fnirs_id"]] = row["subject_id"]
    return out


def scan_raw(sd_to_channel: dict) -> tuple[dict, dict, dict]:
    """Pass A. Returns {subject: {channel: cv}}, {subject: {channel: (mean per
    wavelength)}} and {subject: (date, time)}."""
    cv_by_subject: dict[str, dict] = {}
    level_by_subject: dict[str, dict] = {}
    session: dict[str, tuple] = {}
    crosswalk = load_crosswalk()
    for path in sorted(RAW_SNIRF.glob("*.snirf")):
        stem = path.stem
        sid = crosswalk.get(stem, stem)
        with h5py.File(path, "r") as h5:
            tags = h5["/nirs/metaDataTags"]

            def tag(name: str) -> str:
                try:
                    return np.array(tags[name]).tobytes().decode("utf8", "ignore").strip("\x00")
                except Exception:
                    return "?"

            session[sid] = (tag("MeasurementDate"), tag("MeasurementTime"))
            data = h5["/nirs/data1"]
            series = np.array(data["dataTimeSeries"])
            columns: dict[tuple, list] = collections.defaultdict(list)
            for key in data.keys():
                if not key.startswith("measurementList"):
                    continue
                entry = data[key]
                index = int(key.replace("measurementList", "")) - 1
                pair = (int(np.array(entry["sourceIndex"]).ravel()[0]),
                        int(np.array(entry["detectorIndex"]).ravel()[0]))
                columns[pair].append(index)
        cvs, levels = {}, {}
        for pair, cols in columns.items():
            ch = sd_to_channel.get(pair)
            if ch is None:  # pair not in the analysis montage
                continue
            cvs[ch] = min(
                float(series[:, c].std() / abs(series[:, c].mean())) for c in cols
            )
            # mean raw intensity at each wavelength — a dead channel reads a
            # constant dark floor, identical at both, which is the diagnosis
            levels[ch] = tuple(float(series[:, c].mean()) for c in sorted(cols))
        cv_by_subject[sid] = cvs
        level_by_subject[sid] = levels
    return cv_by_subject, level_by_subject, session


def scan_frozen() -> tuple[dict, float]:
    """Pass B. Returns {subject_id: {channel: response SD}} and the global median."""
    acc: dict[tuple, list] = collections.defaultdict(list)
    with open(PROCESSED / "tidy_block_metrics.csv", encoding="utf8") as fh:
        for row in csv.DictReader(fh):
            if row["chromophore"] != "HbO":
                continue
            try:
                value = float(row["response_minus_baseline_m5_30"])
            except (TypeError, ValueError):
                continue
            acc[(row["subject_id"], int(row["channel_index"]))].append(value)
    sd_by_subject: dict[str, dict] = collections.defaultdict(dict)
    all_sds = []
    for (sid, ch), values in acc.items():
        if len(values) < 3:
            continue
        sd = statistics.pstdev(values)
        sd_by_subject[sid][ch] = sd
        all_sds.append(sd)
    return sd_by_subject, statistics.median(all_sds)


def main() -> None:
    sd_to_channel, active, meta = load_channel_map()
    cv_raw, level_raw, session = scan_raw(sd_to_channel)
    sd_frozen, median_sd = scan_frozen()
    subjects = sorted(cv_raw)
    floor = FROZEN_FLAT_FRACTION * median_sd

    rows = []
    for sid in subjects:
        for ch in sorted(meta):
            cv = cv_raw[sid].get(ch)
            sd = sd_frozen.get(sid, {}).get(ch)
            level = level_raw[sid].get(ch, ())
            src, det, name, roi, dist = meta[ch]
            rows.append({
                "subject_id": sid,
                "session_date": session[sid][0],
                "session_time": session[sid][1],
                "channel_index": ch,
                "channel_name": name,
                "roi": roi,
                "source": src,
                "detector": det,
                "sd_distance_mm": f"{dist:.4f}",
                "active_pair": active[sid][ch],
                "raw_intensity_cv": "" if cv is None else f"{cv:.8f}",
                "raw_mean_wl1": f"{level[0]:.6f}" if len(level) > 0 else "",
                "raw_mean_wl2": f"{level[1]:.6f}" if len(level) > 1 else "",
                "frozen_response_sd_uM": "" if sd is None else f"{sd * 1e6:.8f}",
                "dead_raw": "" if cv is None else int(cv < CV_FLAT),
                "dead_frozen": "" if sd is None else int(sd < floor),
            })
    with open(PROCESSED / "qa_detector_dropout_0730.csv", "w", newline="", encoding="utf8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    out = []

    def say(line: str = "") -> None:
        out.append(line)

    say("QA — detector-level dropout in the temporal montage")
    say("=" * 72)
    n_active = sum(1 for ch in meta if any(active[s][ch] for s in subjects))
    say(f"subjects scanned: {len(subjects)}   montage pairs: {len(meta)} "
        f"({n_active} ever active; 6 and 13 are distance-excluded)")
    say(f"pass A criterion: raw-intensity CV < {CV_FLAT}")
    say(f"pass B criterion: frozen response SD < {FROZEN_FLAT_FRACTION:.0%} of the "
        f"campaign median ({median_sd * 1e6:.4f} uM -> floor {floor * 1e6:.6f} uM)")
    say()

    say("1. Which optodes do the dead channels share?")
    say("-" * 72)
    by_detector = collections.defaultdict(list)
    for ch, (src, det, name, roi, dist) in meta.items():
        by_detector[det].append(ch)
    say(f"{'detector':>9}{'channels':>12}{'roi':>16}{'sessions dead (both ch)':>26}")
    for det in sorted(by_detector):
        chs = sorted(by_detector[det])
        n_dead = sum(
            1 for sid in subjects
            if all(cv_raw[sid].get(c, 9) < CV_FLAT for c in chs)
        )
        roi = meta[chs[0]][3]
        say(f"{det:>9}{str(chs):>12}{roi:>16}{n_dead:>26}")
    say()

    say("2. Per-channel dead counts, and the sensitivity of the threshold")
    say("-" * 72)
    header = f"{'ch':>4}{'name':>14}{'roi':>16}" + "".join(f"{t:>9}" for t in CV_GRID)
    say(header + f"{'frozen':>9}")
    for ch in sorted(meta):
        line = f"{ch:>4}{meta[ch][2]:>14}{meta[ch][3]:>16}"
        for thr in CV_GRID:
            line += f"{sum(1 for s in subjects if cv_raw[s].get(ch, 9) < thr):>9}"
        line += f"{sum(1 for s in subjects if sd_frozen.get(s, {}).get(ch, 9e9) < floor):>9}"
        say(line)
    say()
    say("The `frozen` column is pass B; the five preceding columns are pass A at")
    say("five thresholds spanning a factor of twenty.")
    say()

    say("3. Agreement between the two passes")
    say("-" * 72)
    set_a = {(s, c) for s in subjects for c, v in cv_raw[s].items() if v < CV_FLAT}
    set_b = {(s, c) for s in subjects for c, v in sd_frozen.get(s, {}).items() if v < floor}
    inter = set_a & set_b
    say(f"pass A flags {len(set_a)} subject x channel cells; pass B flags {len(set_b)}")
    say(f"intersection {len(inter)}   A only {len(set_a - set_b)}   B only {len(set_b - set_a)}")
    denom = len(set_a | set_b)
    say(f"Jaccard agreement {len(inter) / denom:.3f}" if denom else "no cells flagged")
    say()

    say("4. The dark floor — what a dead channel reads")
    say("-" * 72)
    live = [v for s in subjects for v in cv_raw[s].values() if v >= CV_FLAT]
    dead = [v for s in subjects for v in cv_raw[s].values() if v < CV_FLAT]
    say(f"live cells  n = {len(live):>5}   CV min {min(live):.5f}  median {statistics.median(live):.5f}")
    say(f"dead cells  n = {len(dead):>5}   CV max {max(dead):.5f}  median {statistics.median(dead):.5f}")
    say()
    lvl_dead = [level_raw[s][c] for s in subjects for c, v in cv_raw[s].items()
                if v < CV_FLAT and len(level_raw[s].get(c, ())) == 2]
    lvl_live = [level_raw[s][c] for s in subjects for c, v in cv_raw[s].items()
                if v >= CV_FLAT and len(level_raw[s].get(c, ())) == 2]
    say("mean raw intensity, by wavelength (arbitrary instrument units):")
    for name, pool in (("dead", lvl_dead), ("live", lvl_live)):
        wl1 = [p[0] for p in pool]
        wl2 = [p[1] for p in pool]
        say(f"  {name:>4}: wl1 median {statistics.median(wl1):.5f} "
            f"(IQR {statistics.quantiles(wl1)[0]:.5f}-{statistics.quantiles(wl1)[2]:.5f})   "
            f"wl2 median {statistics.median(wl2):.5f} "
            f"(IQR {statistics.quantiles(wl2)[0]:.5f}-{statistics.quantiles(wl2)[2]:.5f})")
    say("A dead cell reads the same constant at BOTH wavelengths. Tissue")
    say("attenuates 760 and 850 nm differently, so equality at the two is the")
    say("signature of light that never entered tissue.")
    say(f"dRange lower bound is 0.001, so the floor passes it by a factor of "
        f"{statistics.median([p[0] for p in lvl_dead]) / 0.001:.0f}.")
    say()

    say("5. Does `active_pair` catch any of it?")
    say("-" * 72)
    say(f"{'ch':>4}{'dead (pass A)':>16}{'of which active_pair == 1':>28}")
    for ch in sorted(meta):
        flagged = [s for s in subjects if cv_raw[s].get(ch, 9) < CV_FLAT]
        if not flagged:
            continue
        kept = [s for s in flagged if active[s][ch] == 1]
        say(f"{ch:>4}{len(flagged):>16}{len(kept):>28}")
    say()

    say("6. Live channels per temporal ROI, per subject (active and not dead)")
    say("-" * 72)
    for roi, chs in TEMPORAL_ROIS.items():
        dist = collections.Counter()
        for sid in subjects:
            live_n = sum(
                1 for c in chs
                if active[sid][c] == 1 and cv_raw[sid].get(c, 9) >= CV_FLAT
            )
            dist[live_n] += 1
        say(f"{roi}: " + ", ".join(
            f"{k} live -> {v} subjects" for k, v in sorted(dist.items())))
    say()

    say("7. Session pattern — permanent fault or intermittent connection?")
    say("-" * 72)
    say(f"{'date':<12}{'n':>4}" + "".join(
        f"{'det' + str(d):>10}" for d in sorted(by_detector)))
    by_date = collections.defaultdict(list)
    for sid in subjects:
        by_date[session[sid][0]].append(sid)
    for date in sorted(by_date):
        line = f"{date:<12}{len(by_date[date]):>4}"
        for det in sorted(by_detector):
            chs = by_detector[det]
            n = sum(1 for s in by_date[date]
                    if all(cv_raw[s].get(c, 9) < CV_FLAT for c in chs))
            line += f"{n:>10}"
        say(line)
    say()
    say("Counts are sessions in which every channel of that detector was dead.")

    text = "\n".join(out) + "\n"
    (PROCESSED / "qa_detector_dropout_0730_report.txt").write_text(text, encoding="utf8")
    print(text)


if __name__ == "__main__":
    main()
