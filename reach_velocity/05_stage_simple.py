"""Step 0c, simple diurnal correction, whole record.
1. median filter the raw 10-min series (1 h window) to take out short glitches
2. night line: median of the night readings (sun below the horizon) for each night, straight lines between nights
3. daily wobble = filtered minus night line, averaged by half hour of the day on dry days (night-to-night change under 2 cm),
   one curve per month; it is zero at night by construction
4. corrected = filtered minus the wobble for that month and time of day (a fixed few cm, so peaks stay)
5. final = 3 h centred rolling median of the corrected series
Output: data_clean/stage_simple_10min.csv, figures/00h_stage_full.png, 00i_stage_steps.png, stage_simple_checks.json"""
import json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import matplotlib.dates as mdates
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"; F = HERE / "figures"
TZ = "Etc/GMT+5"
st = pd.read_csv(D / "stage_10min.csv", parse_dates=["time_utc"]).set_index("time_utc").stage_m.dropna()
sun = pd.read_csv(D / "stage_night_10min.csv", parse_dates=["time_utc"]).set_index("time_utc").sun_elevation_deg
df = pd.DataFrame({"raw": st}); df["sun"] = sun.reindex(df.index); df = df.tz_convert(TZ)
df["filtered"] = df.raw.rolling(7, center=True, min_periods=3).median()                                  # 1. median filter, 1 h
night = df.sun < 0; night_id = ((~night).astype(int).diff().fillna(0) == -1).cumsum()                     # 2. one id per night
nl = df[night].groupby(night_id[night]).filtered.agg(["median"]); nl["t"] = df[night].groupby(night_id[night]).apply(lambda g: g.index[len(g) // 2])
df["night_line"] = np.interp(df.index.astype("int64"), nl.t.astype("int64").values, nl["median"].values)
df["hour"] = df.index.hour + df.index.minute / 60; df["hbin"] = df.hour // 0.5 * 0.5; df["month"] = df.index.strftime("%Y-%m"); df["night_date"] = (df.index + pd.Timedelta(hours=1)).normalize()
nm = df[night].groupby("night_date").filtered.median(); dry_nights = nm[nm.diff().abs() < 0.02].index; df["dry"] = df.night_date.isin(dry_nights)
df["wobble_obs"] = df.filtered - df.night_line                                                             # 3. daily wobble
tmpl_all = df[df.dry].groupby("hbin").wobble_obs.mean(); tmpl = {}
for m, g in df[df.dry].groupby("month"): tmpl[m] = g.groupby("hbin").wobble_obs.mean().reindex(tmpl_all.index).interpolate() if g.night_date.nunique() >= 6 else tmpl_all
df["wobble"] = [tmpl.get(m, tmpl_all).get(h, 0.0) for m, h in zip(df.month, df.hbin)]
df["corrected"] = df.filtered - df.wobble                                                                  # 4.
df["final"] = df.corrected.rolling("3h", center=True).median()                                            # 5.
out = df[["raw", "filtered", "night_line", "wobble", "corrected", "final", "dry"]].copy(); out.columns = ["stage_raw_m", "stage_medianfilter_m", "night_line_m", "daily_wobble_m", "stage_corrected_m", "stage_final_m", "dry_day"]
out.index.name = "time_est"; out["time_utc"] = out.index.tz_convert("UTC"); out.to_csv(D / "stage_simple_10min.csv")
pd.DataFrame(tmpl).assign(all_months=tmpl_all).rename_axis("hour_est").to_csv(D / "stage_simple_wobble_by_month.csv")
# checks
def hourprof(x): b = x.rolling("24h", center=True).median(); an = (x - b)[df.dry]; return an.groupby(df.hour[df.dry] // 1).mean() * 100
hp_raw, hp_fin = hourprof(df.raw), hourprof(df.final)
rise = df.raw.rolling("12h").max() - df.raw.rolling("12h").min(); ev = []
for t in rise.sort_values(ascending=False).index:
    if len(ev) >= 8: break
    if all(abs((t - e).total_seconds()) > 5 * 86400 for e in ev): ev.append(t)
pk = [{"time": str(t), "rise_m": float(rise[t]), "raw_max_m": float(df.raw[t - pd.Timedelta(hours=12): t + pd.Timedelta(hours=12)].max()), "final_max_m": float(df.final[t - pd.Timedelta(hours=12): t + pd.Timedelta(hours=12)].max())} for t in ev]
checks = {"dry_nights": int(len(dry_nights)), "dry_day_pattern_cm": {"raw": float(hp_raw.max() - hp_raw.min()), "final": float(hp_fin.max() - hp_fin.min())}, "wobble_cm": {"min": float(df.wobble.min() * 100), "max": float(df.wobble.max() * 100), "afternoon_mean_15_19h": float(df.wobble[(df.hour >= 15) & (df.hour < 19)].mean() * 100)},
          "peak_change_cm": {"mean": float(np.mean([(p["final_max_m"] - p["raw_max_m"]) * 100 for p in pk])), "worst": float(min((p["final_max_m"] - p["raw_max_m"]) * 100 for p in pk))}, "final_minus_night_line_at_night_cm_sd": float(((df.final - df.night_line)[night]).std() * 100), "peaks": pk}
json.dump(checks, open(D / "stage_simple_checks.json", "w"), indent=2); print(json.dumps({k: v for k, v in checks.items() if k != "peaks"}, indent=2)); print(pd.DataFrame(pk).round(3).to_string(index=False))
# figures
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, NIGHTC = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#dcdcf5"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig, axes = plt.subplots(2, 1, figsize=(13, 7), gridspec_kw={"hspace": 0.35})
for ax, yr in zip(axes, [2025, 2026]):
    w = df[df.index.year == yr]; ax.plot(w.index, w.raw, color=INK, lw=0.4, alpha=0.45, label="raw, 10-min"); ax.plot(w.index, w.final, color=BLUE, lw=0.9, label="final: median filter, daily wobble removed, 3 h median")
    ax.set_title(f"{yr}, whole record", fontsize=9); ax.set_ylabel("stage (m)"); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7, loc="upper right")
fig.savefig(F / "00h_stage_full.png", bbox_inches="tight"); plt.close(fig)
fig = plt.figure(figsize=(13, 8)); gs = fig.add_gridspec(3, 2, height_ratios=[1.4, 1.4, 1], hspace=0.45, wspace=0.18)
def shade(ax, w):
    ax.fill_between(w.index, -10, 10, where=w.sun < 0, color=NIGHTC, alpha=0.7, lw=0); ax.xaxis.set_major_locator(mdates.AutoDateLocator(maxticks=7)); ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
for j, (t0, t1, title) in enumerate([("2025-07-27", "2025-08-06", "ten dry days, 2025"), ("2026-06-12", "2026-06-30", "June 2026, storms and the two survey days")]):
    w = df[(df.index >= pd.Timestamp(t0, tz=TZ)) & (df.index < pd.Timestamp(t1, tz=TZ))]
    ax = fig.add_subplot(gs[0, j]); shade(ax, w); ax.plot(w.index, w.raw, color=INK, lw=0.5, alpha=0.5, label="raw"); ax.plot(w.index, w.filtered, color=INK2, lw=0.8, label="median filter (1 h)"); ax.plot(w.index, w.night_line, color=GREEN, lw=1.0, label="night line (nights connected)")
    ax.set_ylim(w.raw.min() - 0.03, w.raw.max() + 0.03); ax.set_title(f"Steps 1 and 2: {title}", fontsize=9); ax.set_ylabel("stage (m)"); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=6.5, loc="upper left")
    ax = fig.add_subplot(gs[1, j]); shade(ax, w); ax.plot(w.index, w.filtered, color=INK2, lw=0.6, alpha=0.7, label="median filter"); ax.plot(w.index, w.final, color=BLUE, lw=1.1, label="final (wobble removed, 3 h median)")
    ax.set_ylim(w.raw.min() - 0.03, w.raw.max() + 0.03); ax.set_title("Steps 4 and 5: the correction applied", fontsize=9); ax.set_ylabel("stage (m)"); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=6.5, loc="upper left")
ax = fig.add_subplot(gs[2, 0]); cmap = plt.get_cmap("Blues")
for i, (m, t) in enumerate(tmpl.items()): ax.plot(t.index, t.values * 100, color=cmap(0.3 + 0.7 * i / max(1, len(tmpl) - 1)), lw=1.1, label=m)
ax.axvspan(0, 5, color=NIGHTC, alpha=0.7, lw=0); ax.axvspan(21, 24, color=NIGHTC, alpha=0.7, lw=0); ax.axhline(0, color=INK, lw=0.5); ax.set_xlim(0, 24); ax.set_xticks(range(0, 25, 3)); ax.set_xlabel("hour of day (EST)"); ax.set_ylabel("wobble (cm)")
ax.set_title("Step 3: the daily wobble, filtered minus night line on dry days, by month", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=6, ncol=3)
ax = fig.add_subplot(gs[2, 1]); ax.plot(hp_raw.index, hp_raw.values, color=INK, lw=1.2, label="raw"); ax.plot(hp_fin.index, hp_fin.values, color=BLUE, lw=1.2, label="final"); ax.axhline(0, color=INK, lw=0.5); ax.axvspan(0, 5, color=NIGHTC, alpha=0.7, lw=0); ax.axvspan(21, 24, color=NIGHTC, alpha=0.7, lw=0)
ax.set_xlim(0, 24); ax.set_xticks(range(0, 25, 3)); ax.set_xlabel("hour of day (EST)"); ax.set_ylabel("stage minus 24 h median (cm)"); ax.set_title(f"Check: within-day pattern on dry days, {checks['dry_day_pattern_cm']['raw']:.1f} cm before, {checks['dry_day_pattern_cm']['final']:.1f} cm after; peaks change {checks['peak_change_cm']['mean']:+.1f} cm on average", fontsize=8.5); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7)
fig.savefig(F / "00i_stage_steps.png", bbox_inches="tight"); plt.close(fig); print("saved 00h, 00i")
