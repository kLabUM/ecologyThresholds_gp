# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "influxdb-client>=1.40",
#   "matplotlib>=3.8",
#   "pandas>=2.0",
# ]
# ///
"""
Main Stem Water Surface Elevation — overlay line plot (dark mode).

All 13 main-stem sites overlaid on a single set of axes (shared time, full
elevation range in feet). Each site is a single line, colored along the same
upstream → downstream ColorBrewer YlGnBu ramp as the ridge plot. Mirrors the
Cowork overlay artifact; the ridge plot files are untouched.

Run:
  uv run --script overlay_plot.py
  (or:  pip install influxdb-client matplotlib pandas; python overlay_plot.py)

Produces main_stem_elevation_overlay.png. Set INFLUX_TOKEN as an env var
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

# ────────────────────────────────────────────────────────────────────────────
# CONFIG
# ────────────────────────────────────────────────────────────────────────────

INFLUX_URL = os.environ.get("INFLUX_URL", "https://us-east-1-1.aws.cloud2.influxdata.com")
INFLUX_ORG = os.environ.get("INFLUX_ORG", "Digital Water Lab")
TOKEN_DEFAULT = ""
INFLUX_TOKEN = os.environ.get("INFLUX_TOKEN", TOKEN_DEFAULT)

# Main-stem sites in upstream → downstream order (same as in the ridge plot).
MAIN_SITES: list[str] = [
    "ARB034", "ARB026", "ARB077", "ARB013", "ARB109", "ARB012", "ARB027",
    "ARB118", "ARB025", "ARB031", "ARB120", "ARB059", "ARB112",
]

STOP: datetime = datetime.now(timezone.utc)
START: datetime = STOP - timedelta(days=30)

FT_PER_MM = 0.00328084
SMOOTHING_WINDOW = 6

# Dark-mode palette.
BG = "#000000"
FG_LIGHT = "#e6edf3"
FG_MUTED = "#8b949e"
GRID = "#21262d"
AXIS = "#30363d"

OUT_PATH = "main_stem_elevation_overlay.png"

# ────────────────────────────────────────────────────────────────────────────
# COLOR GRADIENT (matches the ridge plot)
# ────────────────────────────────────────────────────────────────────────────

def _hex_to_rgb(h): h = h.lstrip("#"); return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
def _rgb_to_hex(rgb): return "#%02x%02x%02x" % tuple(max(0, min(255, int(round(c)))) for c in rgb)
def _lerp(a, b, t): return tuple(a[i] + (b[i] - a[i]) * t for i in range(3))

def _ramp(stops, t):
    if t <= 0: return stops[0]
    if t >= 1: return stops[-1]
    seg = t * (len(stops) - 1); i = int(seg)
    return _lerp(stops[i], stops[i + 1], seg - i)

_GRAD_STOPS = [_hex_to_rgb(h) for h in ("#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8")]

SITE_COLORS = {
    s: _rgb_to_hex(_ramp(_GRAD_STOPS, i / (len(MAIN_SITES) - 1) if len(MAIN_SITES) > 1 else 0))
    for i, s in enumerate(MAIN_SITES)
}

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

def build_elevation_query(sites) -> str:
    sites_literal = "[" + ", ".join(f'"{s}"' for s in sites) + "]"
    return f"""
from(bucket: "sites")
  |> range(start: -365d)
  |> filter(fn: (r) => r._field == "elevation")
  |> filter(fn: (r) => r._value > 0.0)
  |> filter(fn: (r) => contains(value: r._measurement, set: {sites_literal}))
  |> group(columns: ["_measurement"])
  |> last()
  |> keep(columns: ["_measurement", "_value"])
""".strip()

def fetch_data():
    if not INFLUX_TOKEN:
        raise RuntimeError("INFLUX_TOKEN not set. Export it or paste into TOKEN_DEFAULT.")
    with InfluxDBClient(url=INFLUX_URL, token=INFLUX_TOKEN, org=INFLUX_ORG) as client:
        qapi = client.query_api()
        depth_df = qapi.query_data_frame(build_depth_query(MAIN_SITES, START, STOP))
        name_df = qapi.query_data_frame(build_name_query(MAIN_SITES))
        elev_df = qapi.query_data_frame(build_elevation_query(MAIN_SITES))
    if isinstance(depth_df, list):
        depth_df = pd.concat(depth_df, ignore_index=True) if depth_df else pd.DataFrame()
    if isinstance(name_df, list):
        name_df = pd.concat(name_df, ignore_index=True) if name_df else pd.DataFrame()
    if isinstance(elev_df, list):
        elev_df = pd.concat(elev_df, ignore_index=True) if elev_df else pd.DataFrame()
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

    elevations = {}
    if not elev_df.empty and {"_measurement", "_value"}.issubset(elev_df.columns):
        for _, r in elev_df.iterrows():
            sid = r.get("_measurement")
            v = pd.to_numeric(r.get("_value"), errors="coerce")
            if isinstance(sid, str) and pd.notna(v) and v > 0:
                elevations[sid] = float(v) * 1000.0
    return depth_df, names, elevations

# ────────────────────────────────────────────────────────────────────────────
# RENDER
# ────────────────────────────────────────────────────────────────────────────

def render(depth_df, names, elevations) -> None:
    # Build per-site smoothed elevation series in FEET.
    series = {}
    for s in MAIN_SITES:
        sdf = depth_df[depth_df["site_id"] == s].sort_values("_time").reset_index(drop=True)
        elev_mm = elevations.get(s, float("nan"))
        if sdf.empty or not np.isfinite(elev_mm):
            series[s] = pd.DataFrame(columns=["_time", "_ft"])
            continue
        elev_ft = (elev_mm - sdf["_value"]) * FT_PER_MM
        elev_ft = elev_ft.rolling(SMOOTHING_WINDOW, min_periods=1).median()
        series[s] = pd.DataFrame({"_time": sdf["_time"], "_ft": elev_ft}).dropna()

    visible = [s for s in MAIN_SITES if len(series[s]) > 1]

    fig, ax = plt.subplots(figsize=(13, 7.5))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)

    for s in visible:
        sdf = series[s]
        ax.plot(sdf["_time"], sdf["_ft"], color=SITE_COLORS[s], linewidth=1.6, label=names.get(s, s))

    # End-of-line labels placed just past each line in data coords. Collision
    # avoidance is done after the axes have settled so we can convert a pixel
    # gap into data units. ax.text with clip_on=False renders into the margin.
    ends = sorted(
        ({"s": s, "t": series[s]["_time"].iloc[-1], "v": float(series[s]["_ft"].iloc[-1])} for s in visible),
        key=lambda e: -e["v"],   # top first (highest elevation)
    )
    fig.canvas.draw()
    inv = ax.transData.inverted()
    one_px_v = abs(inv.transform((0, 1))[1] - inv.transform((0, 0))[1])
    min_gap_v = 14 * one_px_v
    last_x = max(e["t"] for e in ends)
    offset_t = pd.Timedelta(minutes=max(30, int((STOP - START).total_seconds() / 60 / 80)))
    placed_v = float("inf")
    for e in ends:
        v = e["v"]
        if placed_v - v < min_gap_v:
            v = placed_v - min_gap_v
        placed_v = v
        ax.text(
            last_x + offset_t, v, names.get(e["s"], e["s"]),
            ha="left", va="center", fontsize=11, fontweight="bold",
            color=SITE_COLORS[e["s"]], clip_on=False,
        )

    # Axes styling.
    ax.set_ylabel("Elevation (ft)", color=FG_LIGHT, fontsize=12, fontweight="bold")
    ax.tick_params(axis="y", labelsize=11, colors=FG_MUTED)
    ax.tick_params(axis="x", labelsize=11, colors=FG_MUTED)
    locator = mdates.AutoDateLocator(minticks=4, maxticks=8)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
    ax.grid(True, axis="y", color=GRID, linewidth=0.5, linestyle=(0, (2, 3)))
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(AXIS); ax.spines[side].set_linewidth(0.6)

    # Title.
    fig.suptitle("Main Stem Water Surface Elevation — Digital Water Lab",
                 x=0.04, y=0.97, ha="left", fontsize=16, fontweight="bold", color=FG_LIGHT)
    fig.text(
        0.04, 0.935,
        f"{len(visible)}/{len(MAIN_SITES)} main-stem sites · upstream → downstream gradient · "
        f"6-pt moving median · {START.date()} → {STOP.date()}",
        fontsize=10, color=FG_MUTED, ha="left",
    )

    plt.subplots_adjust(left=0.07, right=0.82, top=0.90, bottom=0.08)
    plt.savefig(OUT_PATH, dpi=150, facecolor=fig.get_facecolor())
    print(f"Saved → {os.path.abspath(OUT_PATH)}")

# ────────────────────────────────────────────────────────────────────────────
# MAIN
# ────────────────────────────────────────────────────────────────────────────

def main() -> None:
    depth_df, names, elevations = fetch_data()
    print(
        f"Fetched {len(depth_df):,} aggregated points across "
        f"{depth_df['site_id'].nunique()}/{len(MAIN_SITES)} main-stem sites "
        f"({START.date()} → {STOP.date()})."
    )
    render(depth_df, names, elevations)

if __name__ == "__main__":
    main()
