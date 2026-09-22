"""The two new MaxBotix records (Upper ARB017, Lower ARB112) diagnosed the way the Middle's was (02_stage_diagnose.py, entry 3):
the raw record, the within-day pattern on dry days by month (the sensor's daytime artefact), the night against the afternoon,
the spikes and implausible readings, and the offset tags. Nothing is cleaned here. Output: figures/16_other_sensors_diagnose.png"""
import numpy as np, pandas as pd
from pathlib import Path
RAW = Path(__file__).resolve().parent / "data_raw"
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import matplotlib.dates as mdates
from pathlib import Path
F = Path(__file__).resolve().parent / "figures"
INK, INK2, AXIS, GRID, BLUE, ORANGE, RED, NIGHT = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#c0392b", "#dcdcf5"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
SENS = {"Upper (ARB017)": ([str(RAW / "26_depthSensorData_Upper.csv")], 1587.5), "Lower (ARB112)": ([str(RAW / "25_depthSensorData_Lower.csv"), str(RAW / "26_depthSensorData_Lower.csv")], 6705.6)}
fig = plt.figure(figsize=(18, 11)); gs = fig.add_gridspec(2, 4, wspace=0.3, hspace=0.45, width_ratios=[2.2, 1, 1, 1])
for i, (name, (files, off)) in enumerate(SENS.items()):
    d = pd.concat([pd.read_csv(f, skiprows=3) for f in files]); d = d[d["_field"] == "value"].copy(); d["t"] = pd.to_datetime(d["_time"], utc=True); d["range_mm"] = d["_value"].astype(float); d["off_tag"] = d.offset.astype(float); d = d.sort_values("t").drop_duplicates("t")
    d["stage"] = (off - d.range_mm) / 1000; loc = d.set_index("t").stage.tz_convert("Etc/GMT+5"); df = pd.DataFrame({"stage": loc.values}, index=loc.index); df["hour"] = df.index.hour + df.index.minute / 60; df["day"] = df.index.date; df["month"] = df.index.to_period("M")
    plaus = df.stage.between(0.05, 3.0); med6 = df.stage.where(plaus).rolling("6h", center=True).median(); spike = plaus & ((df.stage - med6).abs() > 0.15)
    print(f"{name}: {len(df)} readings {df.index.min().date()} to {df.index.max().date()}; implausible (stage < 0.05 or > 3 m) {int((~plaus).sum())} ({(~plaus).mean()*100:.1f}%); spikes > 15 cm off a 6 h median {int(spike.sum())} ({spike.mean()*100:.1f}%); offset tags {sorted(d.off_tag.unique())}; median step {d.t.diff().dt.total_seconds().median()/60:.0f} min; gaps > 6 h: {int((d.t.diff().dt.total_seconds()/3600 > 6).sum())}")
    clean = df.stage.where(plaus & ~spike); anom = clean - clean.rolling("24h", center=True).median(); df["anom"] = anom
    nightmed = df[(df.hour >= 23) | (df.hour < 5)].groupby("day").stage.median(); d2d = nightmed.diff().abs(); dry = d2d[d2d < 0.02].index; df["dry"] = df.day.isin(dry)
    prof = df[df.dry].groupby([df.month.astype(str), df.hour.round(0)]).anom.mean().unstack(0); amp = (prof.max() - prof.min()) * 100
    night = df[((df.hour >= 23) | (df.hour < 5)) & df.dry].groupby("day").stage; aft = df[(df.hour >= 11) & (df.hour < 17) & df.dry].groupby("day").stage
    print(f"   dry days {len(dry)}; within-day pattern on dry days, peak to trough by month (cm): " + ", ".join(f"{m} {a:.1f}" for m, a in amp.items()) + f"; afternoon minus night on dry days {(aft.mean() - night.mean()).mean()*100:+.1f} cm; sd within a night {night.std().median()*100:.1f} cm")
    # (a) the raw record with the dropped readings
    ax = fig.add_subplot(gs[i, 0]); ax.plot(df.index, df.stage.where(plaus), color=INK, lw=0.4); ax.scatter(df.index[~plaus], np.full((~plaus).sum(), 0.0), s=3, color=RED, lw=0); ax.scatter(df.index[spike], df.stage[spike], s=3, color=ORANGE, lw=0)
    ax.set_ylabel("stage (m) = offset − range"); ax.set_title(f"{name}: the raw 10-min record (red at the axis: implausible, n {int((~plaus).sum())}; orange: spikes, n {int(spike.sum())})", fontsize=8.5); ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %y")); ax.grid(True); ax.set_axisbelow(True); ax.set_ylim(-0.05, min(3.0, np.nanquantile(df.stage.where(plaus), .999) + 0.3))
    # (b) the within-day pattern by month on dry days
    ax = fig.add_subplot(gs[i, 1]); cm = plt.get_cmap("viridis")
    for k, m in enumerate(prof.columns): ax.plot(prof.index, prof[m] * 100, color=cm(k / max(len(prof.columns) - 1, 1)), lw=1, label=f"{m} ({amp[m]:.1f} cm)")
    ax.axvspan(0, 5, color=NIGHT, alpha=0.5, lw=0); ax.axvspan(23, 24, color=NIGHT, alpha=0.5, lw=0); ax.set_xlim(0, 24); ax.set_xticks(range(0, 25, 6)); ax.set_xlabel("hour of day (EST)"); ax.set_ylabel("stage minus 24 h rolling median (cm)"); ax.set_title("Within-day pattern on dry days, by month\n(the Middle's sensor showed 2.2 cm, low in the afternoon)", fontsize=8); ax.legend(fontsize=5.5); ax.grid(True); ax.set_axisbelow(True)
    # (c) ten dry days zoomed
    ax = fig.add_subplot(gs[i, 2]); dd_ = sorted(dry)
    if len(dd_) >= 10:
        mid = dd_[len(dd_) // 2]; w = df[(df.index.date >= mid) & (df.index.date < mid + pd.Timedelta(days=10))]; ax.plot(w.index, w.stage.where(plaus.reindex(w.index)), color=INK, lw=0.6)
        for day in sorted(set(w.index.date)): ax.axvspan(pd.Timestamp(day, tz="Etc/GMT+5"), pd.Timestamp(day, tz="Etc/GMT+5") + pd.Timedelta(hours=5), color=NIGHT, alpha=0.5, lw=0)
        ax.set_title(f"Ten dry days from {mid} (shaded: 00:00 to 05:00)", fontsize=8); ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b")); ax.set_ylabel("stage (m)"); ax.grid(True); ax.set_axisbelow(True)
    # (d) the distribution of 10-min changes (spikes show as heavy tails)
    ax = fig.add_subplot(gs[i, 3]); dif = df.stage.where(plaus).diff() * 100; ax.hist(dif.clip(-60, 60).dropna(), bins=121, color=BLUE, lw=0); ax.set_yscale("log"); ax.set_xlabel("change between consecutive readings (cm)"); ax.set_ylabel("count (log)"); ax.set_title(f"Reading-to-reading changes: sd {dif.std():.1f} cm, beyond ±15 cm: {int((dif.abs() > 15).sum())}", fontsize=8); ax.grid(True); ax.set_axisbelow(True)
fig.suptitle("The two new MaxBotix records, diagnosed as the Middle's was: raw, within-day pattern, dry-day zoom, spikes", fontsize=10, y=0.95); fig.savefig(F / "16_other_sensors_diagnose.png", bbox_inches="tight"); plt.close(fig); print("saved 16")
