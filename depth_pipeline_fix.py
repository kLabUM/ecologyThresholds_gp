"""
Drop-in replacement for the per-section depth loading in dataSorting.ipynb (cells 2a/2b/2c).

Fixes:
  1. Each section reads ONLY its own sensor files. The notebook currently leaves
     `depth_data_mm_25` bound to the MIDDLE 2025 file and concatenates it into the
     Upper and Lower series.
  2. Uses the correct stage datum PER GAUGE, instead of hard-coding ARB012's
     7313.2 for all three sites.
  3. Drops MaxBotix error codes (_value <= 500 and _value >= 9000) before aggregating.
  4. Applies the outlier trim to *stage*, after the datum is removed, so a bimodal
     raw series can no longer silently disable it (which is what happens to Upper
     today: 1489 points in, 1489 out, nothing trimmed).

WHY A CONSTANT DATUM AND NOT THE `offset` COLUMN
------------------------------------------------
The `offset` column changes mid-record at two sites, but in both cases `device_id`
is unchanged across the change -- these are datum revisions entered retrospectively,
not physical remountings. Applying them per row injects a step that isn't real:

    Middle, 2025-06-09 (device 074A both sides):  per-row +1.27 ft  vs  constant -0.22 ft
    Lower,  2026-06-08 (device 067A both sides):  per-row +2.62 ft  vs  constant +0.12 ft

A correct datum removes steps; these create them. So: one constant per gauge.

Lower's genuine sensor swap is later (067A -> 147A, between 2026-07-02 and 07-06),
with offset 6705.6 on both sides. That transition shows +0.48 ft regardless of
treatment and may be real.

>>> CONFIRM ARB112's datum against your deployment notes before publishing. The
>>> choice between 5943.6 and 6705.6 moves Lower's absolute stage by 2.5 ft. It
>>> does NOT change the fitted GP (a constant shift cancels under StandardScaler),
>>> only the reported values and any absolute-stage threshold.

Usage in the notebook, replacing cells 2a/2b/2c:

    from depth_pipeline_fix import SECTION_GAUGES, load_depth_series, attach_upstream_depth

    for section, frame in [('Middle', complete_data_Middle),
                           ('Upper',  complete_data_Upper),
                           ('Lower',  complete_data_Lower)]:
        stage = load_depth_series(section)
        attach_upstream_depth(frame, stage, label=section)
"""

import numpy as np
import pandas as pd

MM_TO_FT = 0.00328084
RESAMPLE_MINUTES = 300
OUTLIER_SD = 2.0

# MaxBotix out-of-range / error codes
ERROR_LOW = 500
ERROR_HIGH = 9000

# One gauge per section. `datum_mm` is the stage reference for THAT gauge.
SECTION_GAUGES = {
    "Middle": dict(
        site_id="ARB012",
        datum_mm=7313.2,
        paths=["data/25_depthSensorData_Middle.csv",
               "data/26_depthSensorData_Middle.csv"],
    ),
    "Upper": dict(
        site_id="ARB017",
        datum_mm=1587.5,
        paths=["data/26_depthSensorData_Upper.csv"],
    ),
    "Lower": dict(
        site_id="ARB112",
        datum_mm=6705.6,          # <-- CONFIRM against deployment notes (see header)
        paths=["data/26_depthSensorData_Lower.csv"],
    ),
}


def load_depth_series(section, gauges=None, resample_minutes=RESAMPLE_MINUTES,
                      outlier_sd=OUTLIER_SD, verbose=True):
    """Return a stage series in FEET for one section, indexed by time."""
    cfg = (gauges or SECTION_GAUGES)[section]
    datum, paths = cfg["datum_mm"], cfg["paths"]

    frames = []
    for path in paths:
        raw = pd.read_csv(path, skiprows=3)
        missing = {"_time", "_value"} - set(raw.columns)
        if missing:
            raise ValueError(f"{path}: missing required column(s) {sorted(missing)}")
        frames.append(raw)

    depth = pd.concat(frames, ignore_index=True)
    n_raw = len(depth)

    # Guard against the exact bug this function replaces: one section, one gauge.
    if "site_id" in depth.columns:
        sites = sorted(depth["site_id"].dropna().unique())
        if len(sites) > 1:
            raise ValueError(
                f"{section}: files span {len(sites)} sites {sites} -- these are not "
                f"all from the same gauge. Check the file list."
            )
        if sites and sites[0] != cfg["site_id"]:
            raise ValueError(
                f"{section}: expected gauge {cfg['site_id']} but files contain "
                f"{sites[0]}. The datum {datum} mm belongs to {cfg['site_id']}."
            )

    # Drop sensor error codes BEFORE aggregating.
    bad = (depth["_value"] <= ERROR_LOW) | (depth["_value"] >= ERROR_HIGH)
    depth = depth[~bad]

    depth["_time"] = pd.to_datetime(depth["_time"])
    depth = depth.set_index("_time").sort_index()

    # Stage above the gauge datum, then aggregate, then trim on stage.
    stage_mm = (datum - depth["_value"]).resample(f"{resample_minutes}min").median().dropna()
    n_pre = len(stage_mm)
    mean, sd = stage_mm.mean(), stage_mm.std()
    stage_mm = stage_mm[np.abs(stage_mm - mean) <= outlier_sd * sd]

    stage_ft = (stage_mm * MM_TO_FT).sort_index()

    if verbose:
        print(f"[{section}] gauge {cfg['site_id']}, datum {datum} mm")
        print(f"[{section}]   {n_raw} raw rows -> {int(bad.sum())} error codes dropped "
              f"-> {n_pre} resampled -> {len(stage_ft)} after {outlier_sd}sd trim")
        print(f"[{section}]   stage {stage_ft.min():.3f} to {stage_ft.max():.3f} ft "
              f"({stage_ft.index.min().date()} to {stage_ft.index.max().date()})")

    return stage_ft


def attach_upstream_depth(frame, stage_ft, label="", date_col="Collection Date",
                          out_col="upstream_depth", verbose=True):
    """Add a daily-mean `upstream_depth` column (feet) to `frame`, in place.

    Observations on days with no sensor record get NaN -- which is correct, and is
    what should now happen to Lower's 50 observations from 2025, which are currently
    matched against the Middle gauge ~19 river miles upstream.
    """
    dates = pd.to_datetime(frame[date_col], format="mixed", dayfirst=False).dt.date
    daily = stage_ft.groupby(stage_ft.index.date).mean()
    frame[out_col] = dates.map(daily).astype(float)

    if verbose:
        matched = int(frame[out_col].notna().sum())
        print(f"[{label}] matched {matched}/{len(frame)} observations "
              f"({len(frame) - matched} with no same-day sensor record)")
        if matched:
            vals = frame[out_col].dropna()
            print(f"[{label}]   upstream_depth {vals.min() * 0.3048:.3f} to "
                  f"{vals.max() * 0.3048:.3f} m across {vals.nunique()} distinct survey days")
            if vals.nunique() < 30:
                print(f"[{label}]   NOTE: only {vals.nunique()} distinct values -- this "
                      f"feature is effectively a survey-day index. Use GroupKFold on date.")
    return frame


def audit_gauges(gauges=None):
    """Print deployment segments per gauge so datum changes are visible up front."""
    for section, cfg in (gauges or SECTION_GAUGES).items():
        depth = pd.concat([pd.read_csv(p, skiprows=3) for p in cfg["paths"]],
                          ignore_index=True)
        depth["_time"] = pd.to_datetime(depth["_time"])
        cols = [c for c in ["site_id", "device_id", "offset"] if c in depth.columns]
        seg = depth.groupby(cols)["_time"].agg(["min", "max", "count"])
        print(f"--- {section} (using datum {cfg['datum_mm']} mm) ---")
        print(seg.to_string())
        if "offset" in depth.columns and depth["offset"].nunique() > 1:
            per_dev = depth.groupby("device_id")["offset"].nunique()
            revised = per_dev[per_dev > 1]
            if len(revised):
                print(f"  ! offset changes with device_id UNCHANGED for "
                      f"{list(revised.index)} -- datum revision, not a remounting. "
                      f"Do not apply the offset column per row.")
        print()


if __name__ == "__main__":
    audit_gauges()
    for section in SECTION_GAUGES:
        load_depth_series(section)
        print()
