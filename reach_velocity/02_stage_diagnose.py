"""Step 0c (part 1). Diagnose the diurnal pattern in the MaxBotix stage record before cleaning it.
Writes figures/00e_stage_diurnal.png and prints the amplitude and timing by month."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"; F = HERE / "figures"
INK, INK2, AXIS, GRID, BLUE, ORANGE, NIGHT = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#dcdcf5"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
st = pd.read_csv(D / "stage_10min.csv", parse_dates=["time_utc"]).set_index("time_utc").stage_m.dropna()
loc = st.tz_convert("Etc/GMT+5")            # Eastern standard time all year (solar time, no DST jump)
df = pd.DataFrame({"stage": loc.values}, index=loc.index); df["hour"] = df.index.hour + df.index.minute / 60; df["day"] = df.index.date; df["month"] = df.index.to_period("M")
# anomaly against a 24 h centred rolling median (removes the hydrologic signal at the daily scale, keeps the within-day pattern)
df["base24"] = df.stage.rolling("24h", center=True).median(); df["anom"] = df.stage - df.base24
# dry days: little change day to day, so the within-day pattern is the sensor's, not the river's
drange = df.groupby("day").stage.agg(lambda s: s.max() - s.min()); nightmed = df[(df.hour >= 23) | (df.hour < 5)].groupby("day").stage.median()
d2d = nightmed.diff().abs(); dry = d2d[d2d < 0.02].index; df["dry"] = df.day.isin(dry)
prof = df[df.dry].groupby([df.month.astype(str), df.hour.round(0)]).anom.mean().unstack(0)
amp = (prof.max() - prof.min()) * 100; peak_hr = prof.idxmax(); trough_hr = prof.idxmin()
print("Diurnal pattern on dry days (anomaly from the 24 h rolling median), by month:")
for m in prof.columns: print(f"  {m}: peak-to-trough {amp[m]:.1f} cm, highest at {peak_hr[m]:02.0f}:00, lowest at {trough_hr[m]:02.0f}:00 EST; dry days {df[df.dry & (df.month.astype(str)==m)].day.nunique()}")
night = df[(df.hour >= 23) | (df.hour < 5)]; day_ = df[(df.hour >= 11) & (df.hour < 17)]
print(f"\nSpread within a night (23:00 to 05:00), dry days, median sd: {night[night.dry].groupby('day').stage.std().median()*100:.2f} cm; within an afternoon (11:00 to 17:00): {day_[day_.dry].groupby('day').stage.std().median()*100:.2f} cm")
print(f"Afternoon minus night, dry days: mean {(day_[day_.dry].groupby('day').stage.mean() - night[night.dry].groupby('day').stage.mean()).mean()*100:+.2f} cm (median {(day_[day_.dry].groupby('day').stage.mean() - night[night.dry].groupby('day').stage.mean()).median()*100:+.2f})")
# figure: (a) hour-of-day profile by month; (b) ten dry days in June 2025; (c) a storm window
fig = plt.figure(figsize=(11, 6.6)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.4, wspace=0.22)
ax = fig.add_subplot(gs[0, 0]); cmap = plt.get_cmap("Blues")
for i, m in enumerate(prof.columns): ax.plot(prof.index, prof[m] * 100, color=cmap(0.3 + 0.7 * i / max(1, len(prof.columns) - 1)), lw=1.2, label=m)
ax.axvspan(0, 5, color=NIGHT, alpha=0.6, lw=0); ax.axvspan(23, 24, color=NIGHT, alpha=0.6, lw=0); ax.axhline(0, color=INK, lw=0.5)
ax.set_xlim(0, 24); ax.set_xticks(range(0, 25, 3)); ax.set_xlabel("hour of day (Eastern standard time)"); ax.set_ylabel("stage minus 24 h rolling median (cm)"); ax.set_title("Within-day pattern on dry days, by month (shaded: 23:00 to 05:00)", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=6.5, ncol=2)
def window(ax, t0, t1, title):
    w = df[(df.index >= pd.Timestamp(t0, tz="Etc/GMT+5")) & (df.index < pd.Timestamp(t1, tz="Etc/GMT+5"))]
    for d in pd.date_range(t0, t1, freq="D", tz="Etc/GMT+5"): ax.axvspan(d - pd.Timedelta(hours=1), d + pd.Timedelta(hours=5), color=NIGHT, alpha=0.6, lw=0)
    ax.plot(w.index, w.stage, color=INK, lw=0.6); nm = w[(w.hour >= 23) | (w.hour < 5)].groupby("day").stage.median()
    ax.plot([pd.Timestamp(d, tz="Etc/GMT+5") + pd.Timedelta(hours=2) for d in nm.index], nm.values, "o", color=BLUE, markersize=4, label="night median (23:00 to 05:00)")
    ax.set_title(title, fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.set_ylabel("stage (m)"); ax.legend(fontsize=6.5, loc="upper right")
window(fig.add_subplot(gs[0, 1]), "2025-07-27", "2025-08-06", "Ten dry days, July to August 2025")
window(fig.add_subplot(gs[1, 0]), "2025-06-14", "2025-06-28", "A storm, June 2025")
window(fig.add_subplot(gs[1, 1]), "2026-06-12", "2026-06-30", "The two 2026 survey days (18 and 25 June)")
fig.suptitle("MaxBotix stage record: the diurnal sensor pattern", fontsize=9, y=0.96); fig.savefig(F / "00e_stage_diurnal.png", bbox_inches="tight"); plt.close(fig); print("saved 00e")
