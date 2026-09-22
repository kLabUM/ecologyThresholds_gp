"""The two continuous inputs, whole deployment seasons, raw and smoothed, with the survey days marked.
Output: figures/00k_full_series.png"""
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"; F = HERE / "figures"
g = pd.read_csv(D / "series_10min.csv", parse_dates=["time_utc"]).set_index("time_utc"); g.index = pd.to_datetime(g.index, utc=True)
tr = pd.read_csv(D / "readings_transect.csv"); days = tr.groupby("date").year.first()
INK, AXIS, GRID, BLUE, ORANGE, C25, C26 = "#0b0b0b", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#2a78d6", "#1a9e74"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig, axes = plt.subplots(3, 1, figsize=(14, 10), gridspec_kw={"hspace": 0.38})
def survey_lines(ax, yr):
    for d, y in days.items():
        if y == yr: ax.axvline(pd.Timestamp(d, tz="UTC") + pd.Timedelta(hours=16), color=INK, lw=0.7, alpha=0.5, ls=(0, (3, 3)))
w = g[g.index.year == 2025]
ax = axes[0]; ax.plot(w.index, w.stage_raw_m, color=INK, lw=0.4, alpha=0.4, label="raw, 10-min"); ax.plot(w.index, w.stage_m, color=BLUE, lw=0.9, label="smoothed (what the model will see)"); survey_lines(ax, 2025)
ax.set_title("2025 stage at the Middle sensor (dashed lines: the nine survey days)", fontsize=9); ax.set_ylabel("stage (m)")
ax = axes[1]; ax.plot(w.index, w.adcp_raw_ms, color=INK, lw=0.4, alpha=0.4, label="raw, 10-min bins"); ax.plot(w.index, w.adcp_ms, color=ORANGE, lw=0.9, label="smoothed (what the model will see)"); survey_lines(ax, 2025)
ax.set_title("2025 velocity at the ADCP", fontsize=9); ax.set_ylabel("velocity (m/s)"); ax.set_xlim(axes[0].get_xlim())
w = g[g.index.year == 2026]
ax = axes[2]; ax.plot(w.index, w.stage_raw_m, color=INK, lw=0.4, alpha=0.4, label="raw, 10-min"); ax.plot(w.index, w.stage_m, color=BLUE, lw=0.9, label="smoothed (what the model will see)"); survey_lines(ax, 2026)
ax.set_title("2026 stage at the Middle sensor (dashed lines: the two survey days); no ADCP in 2026", fontsize=9); ax.set_ylabel("stage (m)")
for ax in axes: ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7, loc="upper right")
fig.savefig(F / "00k_full_series.png", bbox_inches="tight"); print("saved 00k")
