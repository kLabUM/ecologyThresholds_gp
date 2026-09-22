"""Step 0c, the simple way. Night = sun below the horizon at the site (42.38 N, 83.92 W). Keep the night readings as they
are, draw straight lines across the days. Output: data_clean/stage_night_10min.csv, figures/00g_stage_night.png"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"; F = HERE / "figures"
LAT, LON = 42.38, -83.92
st = pd.read_csv(D / "stage_10min.csv", parse_dates=["time_utc"]).set_index("time_utc").stage_m.dropna()
def sun_elevation(t_utc, lat=LAT, lon=LON):
    """solar elevation in degrees (NOAA approximation), t_utc a DatetimeIndex in UTC"""
    doy = t_utc.dayofyear.values; hour = (t_utc.hour + t_utc.minute / 60).values
    g = 2 * np.pi / 365 * (doy - 1 + (hour - 12) / 24)
    eqt = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g) - 0.014615 * np.cos(2 * g) - 0.040849 * np.sin(2 * g))
    decl = 0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g) + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g)
    tst = hour * 60 + eqt + 4 * lon; ha = np.deg2rad(tst / 4 - 180)
    return np.rad2deg(np.arcsin(np.sin(np.deg2rad(lat)) * np.sin(decl) + np.cos(np.deg2rad(lat)) * np.cos(decl) * np.cos(ha)))
elev = sun_elevation(st.index); night = elev < 0
out = pd.DataFrame({"stage_raw_m": st.values, "sun_elevation_deg": elev, "is_night": night}, index=st.index)
out["stage_night_m"] = out.stage_raw_m.where(night).interpolate(method="time")     # night readings kept, straight lines across the day
out.to_csv(D / "stage_night_10min.csv")
d = (out.stage_night_m - out.stage_raw_m) * 100
print(f"night = sun below the horizon: {night.mean()*100:.0f}% of readings; night length in June about {24*night[(out.index.month==6)].mean():.1f} h, in September {24*night[(out.index.month==9)].mean():.1f} h")
print(f"daytime readings: night-connected minus raw, mean {d[~night].mean():+.2f} cm, 5th to 95th percentile {d[~night].quantile(.05):+.2f} to {d[~night].quantile(.95):+.2f} cm")
INK, AXIS, GRID, BLUE, NIGHTC = "#0b0b0b", "#c3c2b7", "#e1e0d9", "#2a78d6", "#dcdcf5"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig, axes = plt.subplots(3, 1, figsize=(11, 8.5), gridspec_kw={"hspace": 0.4})
for ax, (t0, t1, title) in zip(axes, [("2025-07-27", "2025-08-06", "Ten dry days, 2025"), ("2025-06-14", "2025-06-28", "June 2025 storm"), ("2026-06-12", "2026-06-30", "June 2026, the two survey days")]):
    w = out[(out.index >= pd.Timestamp(t0, tz="UTC")) & (out.index < pd.Timestamp(t1, tz="UTC"))]
    ax.fill_between(w.index, w.stage_raw_m.min() - 0.05, w.stage_raw_m.max() + 0.05, where=w.is_night, color=NIGHTC, alpha=0.7, lw=0)
    ax.plot(w.index, w.stage_raw_m, color=INK, lw=0.5, alpha=0.6, label="raw, 10-min"); ax.plot(w.index, w.stage_night_m, color=BLUE, lw=1.1, label="night readings, connected")
    ax.set_ylim(w.stage_raw_m.min() - 0.03, w.stage_raw_m.max() + 0.03); ax.set_title(title + "  (shaded: sun below the horizon)", fontsize=9); ax.set_ylabel("stage (m)"); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7, loc="upper left")
fig.savefig(F / "00g_stage_night.png", bbox_inches="tight"); plt.close(fig); print("saved 00g")
