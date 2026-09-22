"""The surveyed cross-sections. The survey sheets carry an RTK 'Elevation' at every measurement point (ft, NAVD88 presumably)
and explicit 'bank height' points. Elevation is the bed: within a transect it falls as depth rises (r -0.9 to -1.0) and
elevation + depth (the water surface) is flat across the transect to a few cm. This script writes them out and compares
the water surface with the stage sensor along the reach, by survey day, to test whether the vertical datum is the same
across sessions. Outputs: data_clean/locations.csv (columns bed_rtk_m, wse_rtk_m), banks.csv, rtk_sessions.csv,
figures/00q_rtk_datum.png"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"; F = HERE / "figures"; RAW = HERE / "data_raw"; FT = 0.3048
raw = pd.concat([pd.read_csv(RAW / f, dtype=str).assign(year=int("20" + f[:2])) for f in ["25_ecohydroData.csv", "26_ecohydroData.csv"]], ignore_index=True)
for c in ["Elevation", "Northing", "Easting", "Total water depth"]: raw[c] = pd.to_numeric(raw[c], errors="coerce")
for c in ["Transect ID", "Measurement Point ID", "Type"]: raw[c] = raw[c].astype(str).str.strip()
raw["date"] = pd.to_datetime(raw["Collection Date"]).dt.strftime("%Y-%m-%d")
RELABEL = {("2025-07-14", "AT0710"): "AT07", ("2025-07-14", "AT0711"): "AT07", ("2025-07-14", "AT0714"): "AT07"}
for (d, pt), t in RELABEL.items(): raw.loc[(raw.date == d) & (raw["Measurement Point ID"] == pt), "Transect ID"] = t
mid = raw[(raw.Northing > 280000) & (raw.Northing < 323000)].copy()
loc = pd.read_csv(D / "locations.csv"); pts = mid[mid.Type == "Transect"].drop_duplicates(["date", "Transect ID", "Measurement Point ID"]).copy()
pts["location_id"] = pts.date + "|" + pts["Transect ID"] + "|" + pts["Measurement Point ID"]; pts = pts.set_index("location_id")
loc["bed_rtk_m"] = loc.location_id.map(pts.Elevation) * FT; loc["wse_rtk_m"] = loc.bed_rtk_m + loc.water_depth_m; loc.to_csv(D / "locations.csv", index=False)
banks = mid[mid["Point Notes"].fillna("").str.contains("bank height", case=False)][["date", "Transect ID", "Measurement Point ID", "Easting", "Northing", "Elevation", "Point Notes"]].copy(); banks["bank_rtk_m"] = banks.Elevation * FT; banks.to_csv(D / "banks.csv", index=False)
# sessions: one row per (date, transect)
ses = loc.groupby(["date", "transect"]).agg(year=("year", "first"), along_m=("along_m", "median"), stage_m=("stage_m", "median"), thalweg_m=("bed_rtk_m", "min"), mean_bed_m=("bed_rtk_m", "mean"), wse_m=("wse_rtk_m", "median"), wse_sd_m=("wse_rtk_m", "std"), n=("bed_rtk_m", "size")).reset_index()
ses["wse_minus_stage_m"] = ses.wse_m - ses.stage_m; ses = ses.sort_values("along_m"); ses.round(3).to_csv(D / "rtk_sessions.csv", index=False)
print(ses.round(2).to_string(index=False))
# ---- tie every session to the stage sensor: WSE - stage = offset(day) - slope * along, one slope for the reach
fit = ses[~((ses.date == "2025-07-14") & (ses.transect == "AT07"))].copy(); days = sorted(fit.date.unique()); A = np.zeros((len(fit), len(days) + 1))
for k, (_, r) in enumerate(fit.iterrows()): A[k, days.index(r.date)] = 1.0; A[k, -1] = -r.along_m
coef, *_ = np.linalg.lstsq(A, fit.wse_minus_stage_m.values, rcond=None); slope = float(coef[-1]); off = pd.Series(coef[:-1], index=days, name="offset_m"); resid = fit.wse_minus_stage_m.values - A @ coef
off.to_csv(D / "rtk_session_offsets.csv"); print(f"surface slope {slope*1e4:.1f} cm per 100 m; residual sd {resid.std()*100:.1f} cm; offsets by day:\n{off.round(3).to_string()}")
loc["bed_rtk_sensor_datum_m"] = loc.bed_rtk_m - loc.date.map(off) + slope * loc.along_m      # RTK bed on the stage sensor's datum, surface slope removed
loc["wse_rtk_sensor_datum_m"] = loc.wse_rtk_m - loc.date.map(off) + slope * loc.along_m
loc.loc[(loc.date == "2025-07-14") & (loc.transect == "AT07"), ["bed_rtk_m", "wse_rtk_m", "bed_rtk_sensor_datum_m", "wse_rtk_sensor_datum_m"]] = np.nan   # copied elevations, no RTK of their own
loc.to_csv(D / "locations.csv", index=False)
banks["date_iso"] = pd.to_datetime(banks.date).dt.strftime("%Y-%m-%d"); banks["along_m"] = banks.apply(lambda r: loc[(loc.date == r.date_iso) & (loc.transect == r["Transect ID"])].along_m.median(), axis=1)
banks["bank_sensor_datum_m"] = banks.bank_rtk_m - banks.date_iso.map(off) + slope * banks.along_m; banks.to_csv(D / "banks.csv", index=False)
ses["day_offset_m"] = ses.date.map(off); ses["wse_minus_stage_detrended_m"] = ses.wse_minus_stage_m - ses.day_offset_m; ses.round(3).to_csv(D / "rtk_sessions.csv", index=False)
INK, INK2, AXIS, GRID, C25, C26, ORANGE = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#1a9e74", "#eb6834"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig, axes = plt.subplots(3, 1, figsize=(11, 10.5), sharex=True, gridspec_kw={"hspace": 0.3})
ax = axes[0]
for yr, col in [(2025, C25), (2026, C26)]:
    s = ses[ses.year == yr]; ax.plot(s.along_m, s.wse_m, "o-", color=col, lw=0.8, ms=4, label=f"{yr} water surface at the visit"); ax.plot(s.along_m, s.thalweg_m, "v", color=col, ms=5, alpha=0.7, label=f"{yr} thalweg (lowest bed point)"); ax.plot(s.along_m, s.mean_bed_m, "_", color=col, ms=8, alpha=0.7)
    for _, r in s.iterrows(): ax.annotate(f"{r.transect}\n{r.date[5:]}", (r.along_m, r.wse_m), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=5.5, color=col)
bk = banks.copy(); bk["along_m"] = bk.apply(lambda r: loc[(loc.date == r.date) & (loc.transect == r["Transect ID"])].along_m.median(), axis=1); ax.plot(bk.along_m, bk.bank_rtk_m, "s", color=ORANGE, ms=4, label="bank height points")
ax.set_ylabel("elevation, RTK (m)"); ax.set_title("RTK elevations along the reach: water surface, thalweg, mean bed (dash), bank tops. The 2026 sessions sit 0.2 to 0.8 m above the 2025 ones at the same places", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7, ncol=3, loc="lower left")
ax = axes[1]
for yr, col in [(2025, C25), (2026, C26)]:
    s = ses[ses.year == yr]; ax.plot(s.along_m, s.wse_minus_stage_m, "o", color=col, ms=5, label=str(yr))
    for _, r in s.iterrows(): ax.annotate(r.date[5:], (r.along_m, r.wse_minus_stage_m), xytext=(0, 6), textcoords="offset points", ha="center", fontsize=5.5, color=col)
ax.set_xlabel("along the reach (m, zero at the ADCP, positive downstream)"); ax.set_ylabel("water surface minus sensor stage (m)"); ax.set_title("If every session shared one datum and the sensor tracked the reach, this would be one smooth falling line", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7)
ax = axes[2]; xx = np.linspace(ses.along_m.min() - 20, ses.along_m.max() + 20, 10); ax.plot(xx, -slope * xx, color=INK, lw=0.8, label=f"one water-surface slope for the reach: {slope*1e4:.1f} cm per 100 m")
for yr, col in [(2025, C25), (2026, C26)]:
    s_ = ses[(ses.year == yr) & ~((ses.date == "2025-07-14") & (ses.transect == "AT07"))]; ax.plot(s_.along_m, s_.wse_minus_stage_detrended_m, "o", color=col, ms=5, label=f"{yr}, each day's offset removed")
ax.set_xlabel("along the reach (m, zero at the ADCP, positive downstream)"); ax.set_ylabel("water surface minus stage, day offset removed (m)"); ax.set_title(f"With one offset per survey day (a total station in 2025, a Reach RS3 in 2026) every session falls on one line: residual {resid.std()*100:.1f} cm", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7)
fig.savefig(F / "00q_rtk_datum.png", bbox_inches="tight"); plt.close(fig); print("saved 00q")
