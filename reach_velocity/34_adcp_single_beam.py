"""The ADCP under the single-beam assumption (the PI, 22 September 2026): Vb is one slanted beam's along-beam velocity, so the along-stream
velocity it samples is Vb / sin(beam angle). With a 25 degree beam that is Vb / 0.423. The chain's prediction at the instrument's own cell
(bed 4 cm below the sensor datum, features from the 2 m grid) is then compared with the record through 2025: level, shape and coverage.
Outputs: adcp_single_beam.json, adcp_single_beam_series.csv, figures/09_adcp_single_beam.png"""
import json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import matplotlib.dates as mdates
from gpcore import *; from uncertainty import *
BEAM_DEG = 25.0; SINB = np.sin(np.radians(BEAM_DEG))
df = load_points(); df["bed_m"] = df.bed_rtk_sensor_datum_m.fillna(df.bed_elev_m); ch = Chain(df, banks=banks_for("Middle"))
grid = pd.read_csv(D / "grid.csv"); cell = grid.iloc[[int(np.argmin(np.hypot(grid.x_m, grid.y_m)))]].reset_index(drop=True); cell["bed_m"] = -0.043   # the surveyed bed 1.7 m away (adcp_recheck.json)
ad = pd.read_csv(D / "adcp.csv", parse_dates=["time_utc"]); ad["time_utc"] = pd.to_datetime(ad.time_utc, utc=True); a2 = ad.set_index("time_utc")[["velocity_ms", "stage_m"]].resample("2h").mean().dropna(); a2["u_beam_ms"] = a2.velocity_ms / SINB
S = sample_series(ch, a2.stage_m.values, cell, n=400, seed=3); v = S["mean_v"]   # one point: its depth-averaged velocity per draw and time
a2["v_pred"] = np.nanmedian(v, 0); a2["v_lo"] = np.nanquantile(v, .05, 0); a2["v_hi"] = np.nanquantile(v, .95, 0)
r = a2.u_beam_ms / a2.v_pred; inb = (a2.u_beam_ms >= a2.v_lo) & (a2.u_beam_ms <= a2.v_hi)
out = {"beam_angle_deg": BEAM_DEG, "factor": SINB, "n_2h": int(len(a2)), "period": [str(a2.index.min().date()), str(a2.index.max().date())], "stage_range_m": [float(a2.stage_m.min()), float(a2.stage_m.max())],
       "level_ratio_record_over_field_median": float(r.median()), "level_ratio_iqr": [float(r.quantile(.25)), float(r.quantile(.75))], "implied_beam_angle_if_level_matched_deg": float(np.degrees(np.arcsin(np.clip((a2.velocity_ms / a2.v_pred).median(), 0, 1)))),
       "MAE_m_s": float(np.mean(np.abs(a2.u_beam_ms - a2.v_pred))), "bias_m_s": float(np.mean(a2.u_beam_ms - a2.v_pred)), "coverage90": float(inb.mean()), "log_ratio_sd_over_season": float(np.log(r).std(ddof=1)),
       "MAE_above_1.0_m": float(np.mean(np.abs(a2.u_beam_ms - a2.v_pred)[a2.stage_m > 1.0])), "coverage90_above_1.0_m": float(inb[a2.stage_m > 1.0].mean()), "n_above_1.0_m": int((a2.stage_m > 1.0).sum()),
       "field_logA_at_cell": float(cell.logA.iloc[0]), "field_logA_sd_at_cell": float(cell.logA_sd.iloc[0])}
json.dump(out, open(D / "adcp_single_beam.json", "w"), indent=1); a2.to_csv(D / "adcp_single_beam_series.csv"); print(json.dumps(out, indent=1))
INK, AXIS, GRID, BLUE, ORANGE, MUTED = "#0b0b0b", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig, axes = plt.subplots(1, 2, figsize=(14, 4.2), gridspec_kw={"width_ratios": [2.2, 1]})
ax = axes[0]; ax.fill_between(a2.index, a2.v_lo, a2.v_hi, color=BLUE, alpha=0.2, lw=0, label="the chain at the instrument's cell, 90% band"); ax.plot(a2.index, a2.v_pred, color=BLUE, lw=0.8); ax.plot(a2.index, a2.u_beam_ms, color=INK, lw=0.7, label=f"ADCP record / sin({BEAM_DEG:.0f}°): the along-stream velocity if Vb is one beam")
ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y")); ax.set_ylabel("depth-averaged velocity at the ADCP (m/s)"); ax.set_title(f"Never used in the fit: level ratio {out['level_ratio_record_over_field_median']:.2f}, MAE {out['MAE_m_s']*1000:.0f} mm/s, 90% band covers {out['coverage90']*100:.0f}% of {len(a2)} two-hour means (stage 0.66 to 1.20 m)", fontsize=8.5); ax.legend(loc="upper right"); ax.grid(True); ax.set_axisbelow(True)
ax = axes[1]; ax.scatter(a2.stage_m, r, s=6, color=INK, alpha=0.5, lw=0); ax.axhline(1, color=BLUE, lw=0.8); ax.axvspan(0.74, 1.09, color=GRID, alpha=0.5, lw=0); ax.set_xlabel("stage (m)"); ax.set_ylabel("record / chain"); ax.set_ylim(0.4, 1.8); ax.set_title(f"The ratio against stage (grey: surveyed range); sd of its log {out['log_ratio_sd_over_season']:.2f}", fontsize=8.5); ax.grid(True); ax.set_axisbelow(True)
fig.suptitle(f"The ADCP as a single {BEAM_DEG:.0f}° beam, against the field and the law through 2025", fontsize=10, y=1.02); fig.savefig(F / "09_adcp_single_beam.png", bbox_inches="tight"); plt.close(fig); print("saved 09")
