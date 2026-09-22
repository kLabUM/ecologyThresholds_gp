"""Step 0d. One 10-minute grid for the two continuous records, both smoothed the same way.
Stage: the corrected series from 05 (stage_final_m) with one more small centred median (2 h).
ADCP: raw readings put on the same 10-minute grid (median within each bin), then the same filters as the stage:
1 h median, 3 h median, 2 h median. Output: data_clean/series_10min.csv (raw and smoothed columns for both),
figures/00j_adcp_smooth.png"""
import json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import matplotlib.dates as mdates
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"; F = HERE / "figures"
sm = pd.read_csv(D / "stage_simple_10min.csv", parse_dates=["time_utc"]).set_index("time_utc"); sm.index = pd.to_datetime(sm.index, utc=True); sm = sm.sort_index()
grid = pd.date_range(sm.index.min().floor("10min"), sm.index.max().ceil("10min"), freq="10min", tz="UTC")
out = pd.DataFrame(index=grid); out.index.name = "time_utc"
out["stage_raw_m"] = sm.stage_raw_m.reindex(grid, method="nearest", tolerance=pd.Timedelta("5min"))
out["stage_corrected_m"] = sm.stage_final_m.reindex(grid, method="nearest", tolerance=pd.Timedelta("5min"))
out["stage_m"] = out.stage_corrected_m.rolling(13, center=True, min_periods=5).median()                   # one more small median, 2 h
RAW = HERE / "data_raw" if (HERE / "data_raw").exists() else ROOT / "data"
a = pd.read_csv(RAW / "25_ADCPData.csv", skiprows=3); a = a[a["_field"] == "value"]
ad = pd.Series(a["_value"].astype(float).values, index=pd.to_datetime(a["_time"], utc=True)).sort_index()
ad10 = ad.resample("10min").median()                                                                        # on the grid, median within each 10-min bin
out["adcp_raw_ms"] = ad10.reindex(grid)
v = out.adcp_raw_ms.rolling(7, center=True, min_periods=3).median()                                        # 1 h median (glitches)
v = v.rolling(19, center=True, min_periods=7).median()                                                     # 3 h median
v = v.rolling(13, center=True, min_periods=5).median()                                                     # 2 h median
out["adcp_ms"] = v.rolling(7, center=True, min_periods=3).mean()                                           # 1 h mean: the SonTek logs at 0.01 m/s steps, medians alone leave a staircase
out.loc[out.adcp_raw_ms.isna(), "adcp_ms"] = np.nan                                                        # no smoothed value where there was no reading
out.to_csv(D / "series_10min.csv")
# how much smoother
def roughness(x): return float((x - x.rolling(7, center=True).median()).std())
info = {"grid": [str(grid[0]), str(grid[-1]), len(grid)], "adcp_bins_with_data": int(out.adcp_raw_ms.notna().sum()),
        "stage_10min_roughness_cm": {"raw": roughness(out.stage_raw_m) * 100, "corrected(05)": roughness(out.stage_corrected_m) * 100, "final": roughness(out.stage_m) * 100},
        "adcp_10min_roughness_ms": {"raw": roughness(out.adcp_raw_ms), "final": roughness(out.adcp_ms)},
        "adcp_vs_stage_residual_sd_ms": {}}
for c in ["adcp_raw_ms", "adcp_ms"]:
    j = out[[c, "stage_m"]].dropna(); b = np.polyfit(j.stage_m, j[c], 2); info["adcp_vs_stage_residual_sd_ms"][c] = float((j[c] - np.polyval(b, j.stage_m)).std())
json.dump(info, open(D / "series_10min_info.json", "w"), indent=2); print(json.dumps(info, indent=2))
INK, INK2, AXIS, GRID, BLUE, ORANGE = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig = plt.figure(figsize=(13, 7.5)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.1], hspace=0.4, wspace=0.32)
w = out[out.adcp_raw_ms.notna()]
ax = fig.add_subplot(gs[0, :]); ax.plot(w.index, w.adcp_raw_ms, color=INK, lw=0.4, alpha=0.4, label="ADCP raw, 10-min bins"); ax.plot(w.index, w.adcp_ms, color=ORANGE, lw=0.9, label="ADCP smoothed (1 h, 3 h, 2 h medians, 1 h mean)")
ax.set_ylabel("velocity at the ADCP (m/s)"); ax.set_title("ADCP, whole record", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7, loc="upper right")
z = out[(out.index >= "2025-06-14") & (out.index < "2025-06-28")]
ax = fig.add_subplot(gs[1, 0]); ax.plot(z.index, z.adcp_raw_ms, color=INK, lw=0.5, alpha=0.5, label="ADCP raw"); ax.plot(z.index, z.adcp_ms, color=ORANGE, lw=1.1, label="ADCP smoothed")
ax2 = ax.twinx(); ax2.plot(z.index, z.stage_m, color=BLUE, lw=1.0, label="stage, final"); ax2.set_ylabel("stage (m)", color=BLUE); ax2.tick_params(axis="y", colors=BLUE)
ax.set_ylabel("velocity (m/s)"); ax.set_title("Two weeks in June 2025, both series on the same grid (stage in blue, right axis)", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.xaxis.set_major_locator(mdates.AutoDateLocator(maxticks=7)); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper left")
ax = fig.add_subplot(gs[1, 1]); j = out.dropna(subset=["adcp_raw_ms", "stage_m"]); ax.scatter(j.stage_m, j.adcp_raw_ms, s=3, color=INK, alpha=0.15, lw=0, label="raw 10-min"); ax.scatter(j.stage_m, j.adcp_ms, s=3, color=ORANGE, alpha=0.5, lw=0, label="smoothed")
ax.set_xlabel("stage, final (m)"); ax.set_ylabel("velocity at the ADCP (m/s)"); ax.set_title(f"Velocity against stage: scatter about a curve {info['adcp_vs_stage_residual_sd_ms']['adcp_raw_ms']*1000:.0f} mm/s raw, {info['adcp_vs_stage_residual_sd_ms']['adcp_ms']*1000:.0f} mm/s smoothed", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7, markerscale=4)
fig.savefig(F / "00j_adcp_smooth.png", bbox_inches="tight"); plt.close(fig); print("saved 00j")
