# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "influxdb-client>=1.40",
#   "matplotlib>=3.8",
#   "pandas>=2.0",
# ]
# ///
"""
Water Surface Elevation Ridgeline — Digital Water Lab (dark mode).

One ridge per site, fixed top→bottom order, water-surface elevation (feet),
per-site normalized and anchored at the streambed, smoothed with a 6-point
moving median. Main-stem sites use a colorblind-safe ColorBrewer YlGnBu ramp
(reed-green → river-blue, upstream→downstream); tributaries are muted slate.
Mirrors the Cowork artifact.

Run:
  uv run --script ridgeline_plot.py
  (or:  pip install influxdb-client matplotlib pandas; python ridgeline_plot.py)

Produces water_surface_elevation_ridgeline.png. Set INFLUX_TOKEN as an env var
(or paste into TOKEN_DEFAULT).
"""

from __future__ import annotations

import os
from datetime import datetime, timedelta, timezone

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from influxdb_client import InfluxDBClient
from matplotlib.patches import Patch

# ────────────────────────────────────────────────────────────────────────────
# CONFIG
# ────────────────────────────────────────────────────────────────────────────

from env import INFLUX_URL, INFLUX_ORG, TOKEN_DEFAULT, INFLUX_TOKEN

# Fixed display order, top → bottom (manually specified; no auto-sort).
SITES: list[str] = [
    "ARB034", "ARB026", "ARB077", "ARB013", "ARB017", "ARB061", "ARB107",
    "ARB082", "ARB083", "ARB109", "ARB012", "ARB027", "ARB028", "ARB118",
    "ARB025", "ARB062", "ARB018", "ARB103", "ARB063", "ARB031", "ARB120",
    "ARB059", "ARB112",
]

MAIN_STEM_SITES: set[str] = {
    "ARB034", "ARB026", "ARB013", "ARB012", "ARB027", "ARB025", "ARB031",
    "ARB059", "ARB109", "ARB077", "ARB112", "ARB118", "ARB120",
}

# Time range.
STOP: datetime = datetime.now(timezone.utc)
START: datetime = STOP - timedelta(days=30)

# Display: elevation always in feet; smoothing fixed to a 6-point moving median;
# tight (overlapping) layout. (These match the artifact, which no longer exposes
# unit / smoothing / overlap toggles.)
FT_PER_MM = 0.00328084
SMOOTHING_WINDOW = 6

# Axes-fraction x where the elevation-range "–" sits; the title centers here.
DASH_AX = 1.05

# "all" — every site in SITES. "main" — only MAIN_STEM_SITES.
SHOW = "all"

# Dark-mode palette.
BG = "#000000"
FG_LIGHT = "#e6edf3"
FG_MUTED = "#8b949e"
GRID = "#30363d"

OUT_PATH = "water_surface_elevation_ridgeline.png"

# ────────────────────────────────────────────────────────────────────────────
# COLOR SYSTEM (matches the artifact)
# ────────────────────────────────────────────────────────────────────────────

def _hex_to_rgb(h): h = h.lstrip("#"); return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
def _rgb_to_hex(rgb): return "#%02x%02x%02x" % tuple(max(0, min(255, int(round(c)))) for c in rgb)
def _lerp(a, b, t): return tuple(a[i] + (b[i] - a[i]) * t for i in range(3))
def _lighten(rgb, f): return _lerp(rgb, (255, 255, 255), f)

def _ramp(stops, t):
    if t <= 0: return stops[0]
    if t >= 1: return stops[-1]
    seg = t * (len(stops) - 1); i = int(seg)
    return _lerp(stops[i], stops[i + 1], seg - i)

_GRAD_STOPS = [_hex_to_rgb(h) for h in ("#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8")]
_TRIB_RGB = _hex_to_rgb("#9fb3c8")

def build_site_colors():
    """fill (data color), stroke/label (lightened so it reads on black)."""
    colors = {}
    main_order = [s for s in SITES if s in MAIN_STEM_SITES]
    m = len(main_order)
    for i, s in enumerate(main_order):
        rgb = _ramp(_GRAD_STOPS, i / (m - 1) if m > 1 else 0)
        colors[s] = {"fill": _rgb_to_hex(rgb), "fill_alpha": 0.85,
                     "stroke": _rgb_to_hex(_lighten(rgb, 0.30))}
    for s in (x for x in SITES if x not in MAIN_STEM_SITES):
        colors[s] = {"fill": _rgb_to_hex(_TRIB_RGB), "fill_alpha": 0.55,
                     "stroke": _rgb_to_hex(_lighten(_TRIB_RGB, 0.30))}
    return colors

SITE_COLORS = build_site_colors()

# ────────────────────────────────────────────────────────────────────────────
# QUERIES
# ────────────────────────────────────────────────────────────────────────────

def pick_aggregate(days: float) -> str:
    if days <= 1:   return "5m"
    if days <= 7:   return "15m"
    if days <= 30:  return "1h"
    if days <= 90:  return "3h"
    if days <= 365: return "6h"
    return "1d"

def build_depth_query(sites, start, stop) -> str:
    sites_literal = "[" + ", ".join(f'"{s}"' for s in sites) + "]"
    window = pick_aggregate((stop - start).total_seconds() / 86400)
    return f"""
from(bucket: "depth")
  |> range(start: {start.isoformat().replace('+00:00', 'Z')}, stop: {stop.isoformat().replace('+00:00', 'Z')})
  |> filter(fn: (r) => r._measurement == "maxbotix_depth")
  |> filter(fn: (r) => contains(value: r.site_id, set: {sites_literal}))
  |> aggregateWindow(every: {window}, fn: mean, createEmpty: false)
  |> keep(columns: ["_time", "_value", "site_id"])
""".strip()

def build_name_query(sites) -> str:
    sites_literal = "[" + ", ".join(f'"{s}"' for s in sites) + "]"
    return f"""
from(bucket: "meta")
  |> range(start: -180d)
  |> filter(fn: (r) => r._field == "gps")
  |> filter(fn: (r) => contains(value: r.site_id, set: {sites_literal}))
  |> last()
  |> keep(columns: ["site_id", "site_name"])
""".strip()

def build_site_config_query(sites) -> str:
    sites_literal = "[" + ", ".join(f'"{s}"' for s in sites) + "]"
    return f"""
from(bucket: "sites")
  |> range(start: -365d)
  |> filter(fn: (r) => r._field == "offset" or r._field == "elevation")
  |> filter(fn: (r) => r._value > 0.0)
  |> filter(fn: (r) => contains(value: r._measurement, set: {sites_literal}))
  |> group(columns: ["_measurement", "_field"])
  |> last()
  |> keep(columns: ["_measurement", "_field", "_value"])
""".strip()

def fetch_data():
    if not INFLUX_TOKEN:
        raise RuntimeError("INFLUX_TOKEN not set. Export it as an env var or paste into TOKEN_DEFAULT.")
    with InfluxDBClient(url=INFLUX_URL, token=INFLUX_TOKEN, org=INFLUX_ORG) as client:
        qapi = client.query_api()
        depth_df = qapi.query_data_frame(build_depth_query(SITES, START, STOP))
        name_df = qapi.query_data_frame(build_name_query(SITES))
        config_df = qapi.query_data_frame(build_site_config_query(SITES))

    if isinstance(depth_df, list):
        depth_df = pd.concat(depth_df, ignore_index=True) if depth_df else pd.DataFrame()
    if isinstance(name_df, list):
        name_df = pd.concat(name_df, ignore_index=True) if name_df else pd.DataFrame()
    if isinstance(config_df, list):
        config_df = pd.concat(config_df, ignore_index=True) if config_df else pd.DataFrame()

    if depth_df.empty:
        raise RuntimeError("No depth records returned for the selected sites/range.")

    depth_df = depth_df[["_time", "_value", "site_id"]].copy()
    depth_df["_time"] = pd.to_datetime(depth_df["_time"], utc=True)
    depth_df["_value"] = pd.to_numeric(depth_df["_value"], errors="coerce")
    depth_df = depth_df.dropna(subset=["_value"])

    names = {}
    if not name_df.empty and {"site_id", "site_name"}.issubset(name_df.columns):
        for _, r in name_df.iterrows():
            sid, sname = r.get("site_id"), r.get("site_name")
            if isinstance(sid, str) and isinstance(sname, str):
                names[sid] = sname

    offsets, elevations = {}, {}
    if not config_df.empty and {"_measurement", "_field", "_value"}.issubset(config_df.columns):
        for _, r in config_df.iterrows():
            sid, field = r.get("_measurement"), r.get("_field")
            v = pd.to_numeric(r.get("_value"), errors="coerce")
            if not isinstance(sid, str) or not pd.notna(v) or v <= 0:
                continue
            if field == "offset":
                offsets[sid] = float(v)
            elif field == "elevation":
                elevations[sid] = float(v) * 1000.0
    return depth_df, names, offsets, elevations

# ────────────────────────────────────────────────────────────────────────────
# RENDER
# ────────────────────────────────────────────────────────────────────────────

def render(depth_df, names, offsets, elevations) -> None:
    def streambed_mm(s):
        o, e = offsets.get(s, float("nan")), elevations.get(s, float("nan"))
        return (e - o) if (np.isfinite(o) and np.isfinite(e)) else float("nan")

    visible_sites = [s for s in SITES if s in MAIN_STEM_SITES] if SHOW == "main" else list(SITES)

    # Per-site elevation series (mm), 6-point moving median.
    series = {}
    for s in visible_sites:
        sdf = depth_df[depth_df["site_id"] == s].sort_values("_time").reset_index(drop=True)
        elev_mm = elevations.get(s, float("nan"))
        if not sdf.empty and np.isfinite(elev_mm):
            sdf["_display"] = elev_mm - sdf["_value"]
            sdf = sdf.dropna(subset=["_display"])
            if not sdf.empty:
                sdf["_smoothed"] = sdf["_display"].rolling(SMOOTHING_WINDOW, min_periods=1).median()
                sdf = sdf.dropna(subset=["_smoothed"])
        else:
            sdf = sdf.iloc[0:0]
        series[s] = sdf

    row_spacing, ridge_height = 0.75, 1.4   # tight (overlapping)
    n = len(visible_sites)

    fig, ax = plt.subplots(figsize=(13, max(7, 0.5 * n + 2)))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)

    y_positions = {s: i * row_spacing for i, s in enumerate(visible_sites)}

    for site in visible_sites:
        sdf = series[site]
        y_base = y_positions[site]
        if sdf.empty or len(sdf) < 2:
            continue

        col = SITE_COLORS[site]
        data_min = float(sdf["_smoothed"].min())
        sbed = streambed_mm(site)
        floor = min(sbed, data_min) if np.isfinite(sbed) else data_min
        v_min = floor
        v_max = max(float(sdf["_smoothed"].max()), v_min + 1)
        v_range = v_max - v_min

        y_top = y_base - (sdf["_smoothed"].clip(lower=v_min) - v_min) / v_range * ridge_height
        x = sdf["_time"].to_numpy()
        # Translucent fills so overlapping ridges blend; black top line keeps
        # each ridge's crest legible.
        ax.fill_between(x, y_top, np.full_like(y_top, y_base, dtype=float),
                        color=col["fill"], alpha=col["fill_alpha"], linewidth=0)
        ax.plot(x, y_top, color="#000000", linewidth=1.2)

        # Range in three aligned pieces so every "–" sits in one column (DASH_AX),
        # vertically centred on the ridge band, in feet.
        smax = float(sdf["_smoothed"].max())
        yc = y_base - ridge_height * 0.5
        ax.annotate(f"{floor*FT_PER_MM:.2f}", xy=(DASH_AX - 0.008, yc),
                    xycoords=("axes fraction", "data"), ha="right", va="center",
                    fontsize=10, color=col["stroke"])
        ax.annotate("–", xy=(DASH_AX, yc), xycoords=("axes fraction", "data"),
                    ha="center", va="center", fontsize=10, color=col["stroke"])
        ax.annotate(f"{smax*FT_PER_MM:.2f} ft", xy=(DASH_AX + 0.008, yc),
                    xycoords=("axes fraction", "data"), ha="left", va="center",
                    fontsize=10, color=col["stroke"])

    # Y axis: descriptive name (primary) + site ID, colored to match each ridge,
    # vertically centred on each ridge band.
    tick_labels = [f"{names.get(s, s)}\n{s}" for s in visible_sites]
    ax.set_yticks([y_positions[s] - ridge_height * 0.5 for s in visible_sites])
    ax.set_yticklabels(tick_labels, fontsize=10)
    for lbl, s in zip(ax.get_yticklabels(), visible_sites):
        lbl.set_color(SITE_COLORS[s]["stroke"])
    ax.invert_yaxis()
    ax.tick_params(axis="y", length=0, pad=6)

    locator = mdates.AutoDateLocator(minticks=4, maxticks=8)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
    ax.tick_params(axis="x", labelsize=12, colors=FG_MUTED)

    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.spines["bottom"].set_linewidth(0.6)

    # Banner (top-left) + the "Water Surface Elevation" title over the right column.
    fig.suptitle("Digital Water Lab — Huron River Network",
                 x=0.04, y=0.985, ha="left", fontsize=16, fontweight="bold", color=FG_LIGHT)
    fig.text(
        0.04, 0.955,
        f"per-site normalized · anchored at streambed (sites.elevation − sites.offset) · "
        f"moving median, 6 pts · {START.date()} → {STOP.date()}",
        fontsize=10, color=FG_MUTED, ha="left",
    )
    # Title centered over the dash column of the elevation ranges below.
    ax.annotate(
        "Water Surface\nElevation",
        xy=(DASH_AX, 1.0), xycoords="axes fraction", xytext=(0, 14), textcoords="offset points",
        ha="center", va="bottom", fontsize=12, fontweight="bold", color=FG_LIGHT,
    )

    leg = ax.legend(
        handles=[
            Patch(facecolor="#7fcdbb", edgecolor="#c7e9b4", label="Main stem (upstream→downstream)"),
            Patch(facecolor="#9fb3c8", edgecolor="#cfd8e3", label="Tributary"),
        ],
        loc="lower right", frameon=False, fontsize=10, bbox_to_anchor=(1.0, -0.06),
    )
    for txt in leg.get_texts():
        txt.set_color(FG_MUTED)

    plt.subplots_adjust(left=0.17, right=0.82, top=0.92, bottom=0.06)
    plt.savefig(OUT_PATH, dpi=150, facecolor=fig.get_facecolor())
    print(f"Saved → {os.path.abspath(OUT_PATH)}")

# ────────────────────────────────────────────────────────────────────────────
# MAIN
# ────────────────────────────────────────────────────────────────────────────

def main() -> None:
    depth_df, names, offsets, elevations = fetch_data()
    print(
        f"Fetched {len(depth_df):,} aggregated points across "
        f"{depth_df['site_id'].nunique()}/{len(SITES)} sites "
        f"({START.date()} → {STOP.date()}); "
        f"{len(offsets)} offsets, {len(elevations)} elevations."
    )
    render(depth_df, names, offsets, elevations)

if __name__ == "__main__":
    main()
