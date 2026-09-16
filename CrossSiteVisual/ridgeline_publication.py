# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "influxdb-client>=1.40",
#   "matplotlib>=3.8",
#   "pandas>=2.0",
# ]
# ///
"""
Publication-ready water-depth figure for the Digital Water Lab MAIN-STEM
network (13 sites), styled as small multiples with per-row y-axes so the scale
is visible right next to each curve. Compact vertical layout with no title or
subtitle (caption is added downstream, in the paper).

Also overlays daily precipitation (mm) per site as translucent bars hanging
from the top of each panel, using the shared-scale rainfall pulled from
Open-Meteo by fetch_rainfall.py.

Design conventions (chosen for print-quality reproduction):
  • Light mode: white background, black axes, dark ink.
  • SI units: water depth in meters (raw is mm; converted at render).
  • One stacked panel per site (13 total), shared x-axis. Each panel's y-scale
    runs 0 → ceil(site max depth) for the plotted window, with integer ticks
    at the two extremes.
  • Rainfall bars share a single mm scale across all panels, so bar heights
    are directly comparable between sites. The right axis matches the left
    (integer end ticks, flush labels) and is titled Precipitation (mm).
  • S1-S13 labels inside each panel (mid-left), ordered N → S
    (upstream → downstream).
  • X-axis dates formatted as YYYY-MM-DD (ISO 8601).
  • Outputs both PNG (300 DPI) and vector PDF.

Run:
  uv run --script ridgeline_publication.py
  (or: pip install influxdb-client matplotlib pandas; python ridgeline_publication.py)

Set INFLUX_TOKEN as an env var (or paste into TOKEN_DEFAULT).
Requires rainfall_by_site.csv next to this script (produced by fetch_rainfall.py).
"""

from __future__ import annotations

import csv
import math
import os
import warnings
from datetime import datetime, timedelta, timezone
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from influxdb_client import InfluxDBClient
from influxdb_client.client.warnings import MissingPivotFunction
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory

warnings.simplefilter("ignore", MissingPivotFunction)

# ────────────────────────────────────────────────────────────────────────────
# CONFIG
# ────────────────────────────────────────────────────────────────────────────

from env import INFLUX_URL, INFLUX_ORG, TOKEN_DEFAULT, INFLUX_TOKEN

# All 23 sites in N → S order; MAIN_STEM filters to the 13 shown in the figure.
SITES: list[str] = [
    "ARB034", "ARB026", "ARB077", "ARB013", "ARB017", "ARB061", "ARB107",
    "ARB082", "ARB083", "ARB109", "ARB012", "ARB027", "ARB028", "ARB118",
    "ARB025", "ARB062", "ARB018", "ARB103", "ARB063", "ARB031", "ARB120",
    "ARB059", "ARB112",
]
MAIN_STEM: set[str] = {
    "ARB034", "ARB026", "ARB013", "ARB012", "ARB027", "ARB025", "ARB031",
    "ARB059", "ARB109", "ARB077", "ARB112", "ARB118", "ARB120",
}
SITES_TO_PLOT: list[str] = [s for s in SITES if s in MAIN_STEM]

START: datetime = datetime(2026, 4, 1, tzinfo=timezone.utc)
STOP: datetime = datetime(2026, 7, 1, tzinfo=timezone.utc)

MM_TO_M = 0.001
SMOOTHING_WINDOW = 6           # points, on the aggregated series

# Publication palette (light mode).
BG = "#ffffff"
FG = "#111111"
GRID = "#d0d0d0"
FILL_COLOR   = "#2a3f5f"       # dark navy — water fill
STROKE_COLOR = "#0f1e35"
FILL_ALPHA   = 0.75

RAIN_COLOR = "#4a90d9"         # muted blue — rainfall bars
RAIN_ALPHA = 0.40              # translucent so overhang doesn't hide water
RAIN_CSV = Path("rainfall_by_site.csv")
INFLUX_TIMEOUT_MS = 120_000    # longer ranges need more than the 10 s default

# Yellow highlight bands spanning every panel (top of S1 → bottom of S13).
HIGHLIGHT_PERIODS: list[tuple[datetime, datetime]] = [
    (datetime(2026, 4, 15, tzinfo=timezone.utc),
     datetime(2026, 5, 15, tzinfo=timezone.utc)),
    (datetime(2026, 6, 9, tzinfo=timezone.utc),
     datetime(2026, 6, 13, tzinfo=timezone.utc)),
]
HIGHLIGHT_FACE = "#ffe14a"
HIGHLIGHT_ALPHA = 0.22

# Orange baseline extensions (1-based panel numbers S1…S13) drawn to the
# right of the stack, labeled "Dam".
DAM_PANELS = {1, 2, 3, 5, 8, 9, 10, 11, 12}
DAM_COLOR = "#e67e22"
DAM_LW = 1.3
# Axes-fraction: left spine, across the panel, through the precip labels.
DAM_X0 = 0.00
DAM_X1 = 1.28
DAM_LABEL_X = 1.20

OUT_PNG = "ridgeline_publication.png"
OUT_PDF = "ridgeline_publication.pdf"
OUT_KEY = "ridgeline_publication_key.csv"

# ────────────────────────────────────────────────────────────────────────────
# QUERIES  (same shape as ridgeline_plot.py)
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
  |> filter(fn: (r) => r._field == "offset")
  |> filter(fn: (r) => r._value > 0.0)
  |> filter(fn: (r) => contains(value: r._measurement, set: {sites_literal}))
  |> group(columns: ["_measurement"])
  |> last()
  |> keep(columns: ["_measurement", "_value"])
""".strip()

def fetch_data():
    if not INFLUX_TOKEN:
        raise RuntimeError("INFLUX_TOKEN not set. Export it or paste into TOKEN_DEFAULT.")
    with InfluxDBClient(
        url=INFLUX_URL, token=INFLUX_TOKEN, org=INFLUX_ORG,
        timeout=INFLUX_TIMEOUT_MS,
    ) as client:
        qapi = client.query_api()
        depth_df = qapi.query_data_frame(build_depth_query(SITES_TO_PLOT, START, STOP))
        name_df = qapi.query_data_frame(build_name_query(SITES_TO_PLOT))
        cfg_df = qapi.query_data_frame(build_site_config_query(SITES_TO_PLOT))
    if isinstance(depth_df, list):
        depth_df = pd.concat(depth_df, ignore_index=True) if depth_df else pd.DataFrame()
    if isinstance(name_df, list):
        name_df = pd.concat(name_df, ignore_index=True) if name_df else pd.DataFrame()
    if isinstance(cfg_df, list):
        cfg_df = pd.concat(cfg_df, ignore_index=True) if cfg_df else pd.DataFrame()
    if depth_df.empty:
        raise RuntimeError("No depth records returned for the selected sites/range.")

    depth_df = depth_df[["_time", "_value", "site_id"]].copy()
    depth_df["_time"] = pd.to_datetime(depth_df["_time"], utc=True)
    depth_df["_value"] = pd.to_numeric(depth_df["_value"], errors="coerce")
    depth_df = depth_df.dropna(subset=["_value"])

    names: dict[str, str] = {}
    if not name_df.empty and {"site_id", "site_name"}.issubset(name_df.columns):
        for _, r in name_df.iterrows():
            sid, sname = r.get("site_id"), r.get("site_name")
            if isinstance(sid, str) and isinstance(sname, str):
                names[sid] = sname

    offsets: dict[str, float] = {}
    if not cfg_df.empty and {"_measurement", "_value"}.issubset(cfg_df.columns):
        for _, r in cfg_df.iterrows():
            sid = r.get("_measurement")
            v = pd.to_numeric(r.get("_value"), errors="coerce")
            if isinstance(sid, str) and pd.notna(v) and v > 0:
                offsets[sid] = float(v)                # mm
    return depth_df, names, offsets

def load_rainfall() -> pd.DataFrame:
    """Load per-site daily precipitation (mm) produced by fetch_rainfall.py."""
    if not RAIN_CSV.exists():
        print(f"NOTE: {RAIN_CSV} not found — skipping rainfall overlay. "
              f"Run fetch_rainfall.py first.")
        return pd.DataFrame(columns=["site_id", "date", "precip_mm"])
    r = pd.read_csv(RAIN_CSV)
    r["date"] = pd.to_datetime(r["date"], utc=True)
    r["precip_mm"] = pd.to_numeric(r["precip_mm"], errors="coerce").fillna(0.0)
    # Trim to the plotting window.
    r = r[(r["date"] >= START) & (r["date"] < STOP + timedelta(days=1))]
    return r

# ────────────────────────────────────────────────────────────────────────────
# RENDER (small multiples)
# ────────────────────────────────────────────────────────────────────────────

def render(depth_df: pd.DataFrame, names: dict, offsets: dict, rain_df: pd.DataFrame) -> None:
    n = len(SITES_TO_PLOT)

    # Per-site depth (m) with a 6-pt rolling median; clipped at 0.
    # y-ceiling is the integer just above (or equal to) that site's max.
    series: dict[str, pd.DataFrame] = {}
    y_ceilings: dict[str, int] = {}
    for s in SITES_TO_PLOT:
        sdf = depth_df[depth_df["site_id"] == s].sort_values("_time").reset_index(drop=True)
        off_mm = offsets.get(s, float("nan"))
        if sdf.empty or not np.isfinite(off_mm):
            series[s] = pd.DataFrame(columns=["_time", "_m"])
            y_ceilings[s] = 1
            continue
        depth_m = (off_mm - sdf["_value"]) * MM_TO_M
        depth_m = depth_m.clip(lower=0)
        depth_m = depth_m.rolling(SMOOTHING_WINDOW, min_periods=1).median()
        series[s] = pd.DataFrame({"_time": sdf["_time"], "_m": depth_m}).dropna()
        site_max = float(series[s]["_m"].max()) if not series[s].empty else 0.0
        y_ceilings[s] = max(1, int(np.ceil(site_max)))

    # Shared rainfall scale across all sites. Round up to a tidy step so the
    # bars in every panel are directly comparable.
    if not rain_df.empty and rain_df["precip_mm"].max() > 0:
        raw_max = float(rain_df["precip_mm"].max())
        # Round up to the next multiple of 5 mm (min 10 for a sensible axis).
        rain_max = max(10.0, math.ceil(raw_max / 5.0) * 5.0)
    else:
        rain_max = 0.0

    # Stacked panels, shared x. Extra width on the right for Dam stubs.
    fig_w = 8.0
    fig_h = max(8.4, 0.56 * n + 1.15)
    fig, axes = plt.subplots(n, 1, sharex=True, figsize=(fig_w, fig_h))
    fig.patch.set_facecolor(BG)
    rain_axes: list = []

    for i, s in enumerate(SITES_TO_PLOT):
        ax = axes[i]
        ax.set_facecolor(BG)

        sdf = series[s]
        if not sdf.empty and len(sdf) > 1:
            ax.fill_between(sdf["_time"], sdf["_m"], 0,
                            color=FILL_COLOR, alpha=FILL_ALPHA, linewidth=0)
            ax.plot(sdf["_time"], sdf["_m"], color=STROKE_COLOR, linewidth=0.7)

        # Per-site y-scale: 0 → ceil(max depth in this window), integer ticks.
        # Extreme labels sit flush with their ticks (not vertically centered):
        # top of the max against the top tick, bottom of "0" against the bottom.
        y_ceiling = y_ceilings[s]
        y_ticks = [0, y_ceiling]
        ax.set_ylim(0, y_ceiling)
        ax.set_yticks(y_ticks)
        y_labels = ax.set_yticklabels([f"{t}" for t in y_ticks])
        y_labels[0].set_verticalalignment("bottom")
        y_labels[-1].set_verticalalignment("top")
        ax.tick_params(axis="y", labelsize=8, colors=FG, length=2)

        # Rainfall bars hanging from the top (shared mm scale).
        # Inverted ylim so 0 sits at the top tick and rain_max at the bottom.
        if rain_max > 0:
            rax = ax.twinx()
            rax.set_facecolor("none")
            rax.set_ylim(rain_max, 0)
            rdf = rain_df[rain_df["site_id"] == s]
            if not rdf.empty:
                # bar centers at noon UTC so each bar sits over its day.
                bar_x = rdf["date"] + pd.Timedelta(hours=12)
                rax.bar(
                    bar_x, rdf["precip_mm"],
                    width=0.8, align="center",
                    color=RAIN_COLOR, alpha=RAIN_ALPHA,
                    linewidth=0, zorder=3,
                )
            rain_ticks = [0, int(rain_max)]
            rax.set_yticks(rain_ticks)
            r_labels = rax.set_yticklabels([f"{t}" for t in rain_ticks])
            # Inverted axis: "0" is the top-most number, rain_max the bottom-most.
            r_labels[0].set_verticalalignment("top")
            r_labels[-1].set_verticalalignment("bottom")
            rax.tick_params(axis="y", labelsize=8, colors=FG, length=2)
            for side in ("top", "right", "bottom", "left"):
                rax.spines[side].set_visible(False)
            rain_axes.append(rax)
            # Match the water-axis x-range on the twin so bars don't drift.
            rax.set_xlim(START, STOP)
            rax.margins(x=0)

        # Site label just inside the plot, right of the y-axis and vertically
        # "centered" on the panel so the curve stays flush to the left spine.
        ax.text(0.006, 0.75, f"S{i + 1}", transform=ax.transAxes,
                ha="left", va="center", fontsize=10, fontweight="bold",
                color=FG, family="serif")

        # Clean spines: left + right visible for every panel; bottom on the last.
        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_color(FG); ax.spines["left"].set_linewidth(0.5)
        ax.spines["right"].set_visible(True)
        ax.spines["right"].set_color(FG); ax.spines["right"].set_linewidth(0.5)
        if i < n - 1:
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(axis="x", bottom=False, labelbottom=False)
        else:
            ax.spines["bottom"].set_color(FG); ax.spines["bottom"].set_linewidth(0.5)
            locator = mdates.AutoDateLocator(minticks=4, maxticks=6)
            ax.xaxis.set_major_locator(locator)
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
            ax.tick_params(axis="x", labelsize=8, colors=FG)

        # Shared x-range with no side margins so the hydrograph meets the
        # left spine (and the right edge) instead of floating inward.
        ax.set_xlim(START, STOP)
        ax.margins(x=0)

    # Figure-level y-axis labels (rotated, tight against the tick numbers).
    fig.text(0.005, 0.55, "Water Level (m)", rotation=90, va="center", ha="center",
             fontsize=11, color=FG, family="serif")
    fig.text(0.745, 0.55, "Precipitation (mm)", rotation=90, va="center", ha="center",
             fontsize=11, color=FG, family="serif")

    plt.subplots_adjust(left=0.05, right=0.70, top=0.96, bottom=0.11, hspace=0.30)
    # Draw once so tick Text objects exist, then pin extreme labels flush
    # with their ticks (draw/layout can otherwise restore the default va).
    fig.canvas.draw()

    # Yellow bands in date-x / figure-y so each cuts through every panel
    # (and the gaps) instead of being clipped to a single axes.
    periods = [(a, b) for a, b in HIGHLIGHT_PERIODS if b > a]
    if periods:
        for ax in np.atleast_1d(axes):
            ax.set_facecolor("none")
            ax.set_zorder(1)
        pos_top = axes[0].get_position()
        pos_bot = axes[-1].get_position()
        trans = blended_transform_factory(axes[0].transData, fig.transFigure)
        band_h = pos_top.y1 - pos_bot.y0
        for t0, t1 in periods:
            x0 = mdates.date2num(t0)
            x1 = mdates.date2num(t1)
            fig.add_artist(Rectangle(
                (x0, pos_bot.y0), x1 - x0, band_h,
                transform=trans,
                facecolor=HIGHLIGHT_FACE,
                alpha=HIGHLIGHT_ALPHA,
                edgecolor="none",
                linewidth=0,
                zorder=0,
                clip_on=False,
            ))

    # Bold orange baseline stubs to the right of selected panels, labeled Dam.
    for i, ax in enumerate(np.atleast_1d(axes)):
        if (i + 1) not in DAM_PANELS:
            continue
        trans = blended_transform_factory(ax.transAxes, ax.transData)
        ax.plot(
            [DAM_X0, DAM_X1], [0, 0],
            transform=trans,
            color=DAM_COLOR,
            linewidth=DAM_LW,
            solid_capstyle="butt",
            clip_on=False,
            zorder=5,
        )
        ax.annotate(
            "Dam",
            xy=(DAM_LABEL_X, 0.0),
            xycoords=trans,
            xytext=(0, 3),
            textcoords="offset points",
            ha="center", va="bottom",
            fontsize=8, fontweight="bold",
            color=DAM_COLOR, family="serif",
            annotation_clip=False,
            zorder=6,
        )

    # Column titles for the Dam strip: N (top) → S (bottom).
    pos0 = axes[0].get_position()
    posn = axes[-1].get_position()
    dam_x = pos0.x0 + DAM_LABEL_X * pos0.width
    fig.text(dam_x, pos0.y1 + 0.012, "Upstream",
             ha="center", va="bottom", fontsize=9, fontweight="bold",
             color=FG, family="serif")
    fig.text(dam_x, posn.y0 - 0.018, "Downstream",
             ha="center", va="top", fontsize=9, fontweight="bold",
             color=FG, family="serif")
    for ax in np.atleast_1d(axes):
        labels = ax.get_yticklabels()
        if len(labels) >= 2:
            labels[0].set_verticalalignment("bottom")
            labels[-1].set_verticalalignment("top")
    for rax in rain_axes:
        labels = rax.get_yticklabels()
        if len(labels) >= 2:
            # Inverted rainfall axis: 0 at the top tick, max at the bottom tick.
            labels[0].set_verticalalignment("top")
            labels[-1].set_verticalalignment("bottom")
    plt.savefig(OUT_PNG, dpi=300, facecolor=BG, bbox_inches="tight")
    plt.savefig(OUT_PDF, facecolor=BG, bbox_inches="tight")
    print(f"Saved → {os.path.abspath(OUT_PNG)}")
    print(f"Saved → {os.path.abspath(OUT_PDF)}")

# ────────────────────────────────────────────────────────────────────────────
# MAIN
# ────────────────────────────────────────────────────────────────────────────

def main() -> None:
    depth_df, names, offsets = fetch_data()
    rain_df = load_rainfall()
    print(
        f"Fetched {len(depth_df):,} aggregated points across "
        f"{depth_df['site_id'].nunique()}/{len(SITES_TO_PLOT)} main-stem sites "
        f"({START.date()} → {STOP.date()})."
    )
    if not rain_df.empty:
        print(
            f"Rainfall: {len(rain_df):,} daily records across "
            f"{rain_df['site_id'].nunique()} sites "
            f"(max = {rain_df['precip_mm'].max():.1f} mm/day)."
        )
    render(depth_df, names, offsets, rain_df)

if __name__ == "__main__":
    main()
