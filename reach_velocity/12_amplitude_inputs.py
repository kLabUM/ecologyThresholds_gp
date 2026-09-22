"""Step 0f. The inputs for the spatial model, with the log-law amplitude as the target.
Per location: a = amplitude of v(z) = a ln(z / z0) (z0 = 11 mm for the reach, from 08), its sd from the readings' implied
values, the shear velocity u* = 0.41 a, the bed elevation relative to the stage datum = stage at the visit minus water depth.
ADCP: the instrument sits on the bed (the PI, 21 Sep 2026) with upward beams, so its velocity is taken as the depth-averaged
velocity of the column above it, v_adcp = a (ln(h / z0) - 1), and the amplitude is a = v_adcp / (ln(h / z0) - 1) with the local
depth h = stage minus the h0 fitted in step 1b (0.63 m, the stage at which the point runs dry). The earlier 1 ft point
placeholder is kept as a column for comparison. Still to confirm with the author: the SonTek model and what 'Vb' is. The scaling of the ADCP
velocity with local depth is tested against the two physical candidates: exponent 1/2 (shear velocity at fixed slope) and 2/3
(Manning). Outputs: locations.csv (columns added), adcp_amplitude_2h.csv, amplitude_inputs.json, figures/00n_amplitude_inputs.png"""
import json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"; F = HERE / "figures"; RAW = HERE / "data_raw"
KAPPA, RHO = 0.41, 1000.0; z0 = json.load(open(D / "field_qc.json"))["z0_m"]; ZMIN, VMIN = 0.03, 0.005
Z_ADCP = 0.3048   # the notebooks' old point assumption, kept only for comparison
ADCP_ON_BED = True  # the PI, 21 Sep 2026: the instrument sits on the bed; its velocity is treated as the column average
h0_point = json.load(open(D / "step1b_laws.json"))["manning_2_3"]["params"][1] if (D / "step1b_laws.json").exists() else float("nan")   # h0 from the point (Manning) reading of the ADCP, step 1b (superseded; recorded for the logbook)
from scipy.optimize import curve_fit
loc = pd.read_csv(D / "locations.csv"); tr = pd.read_csv(D / "readings_transect.csv")
# ---- per-location amplitude uncertainty from the readings' implied amplitudes
a_sd = {}
for lid, p in tr.groupby("location_id"):
    ok = p[(p.height_m >= ZMIN) & (p.velocity_ms > VMIN)]
    if len(ok) >= 2: a_sd[lid] = float((ok.velocity_ms / np.log(ok.height_m / z0)).std(ddof=1))
pooled = json.load(open(D / "field_qc.json"))["transfer_check"]["MAE_m_s"] / np.log(0.3048 / z0)
loc["a_sd"] = loc.location_id.map(a_sd).fillna(pooled); loc["u_star"] = KAPPA * loc.a; loc["tau_bed_Pa"] = RHO * loc.u_star ** 2
loc["bed_elev_m"] = loc.stage_m - loc.water_depth_m          # relative to the stage sensor's zero
loc.to_csv(D / "locations.csv", index=False)
# ---- ADCP
g = pd.read_csv(D / "series_10min.csv", parse_dates=["time_utc"]).set_index("time_utc"); g.index = pd.to_datetime(g.index, utc=True)
h2 = g.dropna(subset=["adcp_ms", "stage_m"])[["stage_m", "adcp_ms"]].resample("2h").mean().dropna()
# column model fitted jointly: v = c (stage - h0)^b (ln((stage - h0)/z0) - 1); h0 and b come from this fit, not from the point model
def column_b(x, c, h0, b): return c * np.clip(x - h0, 1e-6, None) ** b * (np.log(np.clip(x - h0, z0 * 2.72, None) / z0) - 1)
pc, _ = curve_fit(column_b, h2.stage_m.values, h2.adcp_ms.values, p0=[0.05, 0.5, 0.5], bounds=([0, -2, 0], [10, 0.66, 2]), maxfev=20000); h0_adcp, b_col = float(pc[1]), float(pc[2])
ad = pd.DataFrame({"stage_m": h2.stage_m, "velocity_ms": h2.adcp_ms}); ad["local_depth_m"] = ad.stage_m - h0_adcp
ad["a_point_1ft_old"] = ad.velocity_ms / np.log(Z_ADCP / z0)
den = np.log(np.clip(ad.local_depth_m, z0 * np.e * 1.01, None) / z0) - 1; ad["a"] = np.where(ad.local_depth_m > z0 * np.e * 1.5, ad.velocity_ms / den, np.nan)   # column average: v = a (ln(h/z0) - 1)
ad.index.name = "time_utc"; ad.to_csv(D / "adcp_amplitude_2h.csv")
m = ad[ad.local_depth_m > 0.02]; b_free, logc = np.polyfit(np.log(m.local_depth_m), np.log(m.velocity_ms), 1)
def mae_b(b): c = np.exp(np.mean(np.log(m.velocity_ms) - b * np.log(m.local_depth_m))); return float(np.mean(np.abs(m.velocity_ms - c * m.local_depth_m ** b)) * 1000), c
ma = ad.dropna(subset=["a"]); ba, _ = np.polyfit(np.log(ma.local_depth_m), np.log(ma.a), 1)
info = {"z0_m": z0, "kappa": KAPPA, "adcp_on_bed_column_average": ADCP_ON_BED, "z_adcp_point_old_m": Z_ADCP, "h0_adcp_m_column_fit": h0_adcp, "h0_adcp_m_point_fit_step1b": h0_point, "column_fit_exponent_b": b_col, "column_fit_MAE_mm_s": float(np.mean(np.abs(h2.adcp_ms.values - column_b(h2.stage_m.values, *pc))) * 1000), "adcp_amplitude_range": [float(ma.a.min()), float(ma.a.max())], "adcp_amplitude_vs_depth_exponent": float(ba), "adcp_depth_scaling": {"exponent_free": float(b_free), "MAE_mm_s_free": mae_b(b_free)[0], "MAE_mm_s_b_0.5": mae_b(0.5)[0], "MAE_mm_s_b_2/3": mae_b(2 / 3)[0]},
        "locations": {"n": int(len(loc)), "a_range": [float(loc.a.min()), float(loc.a.max())], "a_sd_median": float(loc.a_sd.median()), "u_star_range_m_s": [float(loc.u_star.min()), float(loc.u_star.max())], "bed_elev_range_m": [float(loc.bed_elev_m.min()), float(loc.bed_elev_m.max())], "corr_a_vs_depth": float(np.corrcoef(loc.a, loc.water_depth_m)[0, 1])},
        "assumptions": ["z0 = 11 mm for every location (fit flat between 7 and 16 mm)", "ADCP on the bed, upward beams: its velocity is the column average, a = v / (ln(h/z0) - 1); SonTek model and the meaning of Vb still to confirm", "ADCP local depth = stage - h0 with h0 fitted under the column model (0.55 m)", "bed elevation = stage at the visit minus surveyed water depth, on the stage sensor's datum"]}
json.dump(info, open(D / "amplitude_inputs.json", "w"), indent=2); print(json.dumps(info, indent=2))
# ---- figure
INK, INK2, AXIS, GRID, BLUE, ORANGE, C25, C26 = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#2a78d6", "#1a9e74"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig = plt.figure(figsize=(14, 9)); gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.3], hspace=0.35, wspace=0.3)
for j, (col, cmap, label, vmin, vmax) in enumerate([("a", "Blues", "log-law amplitude a (m/s)", 0, 0.2), ("bed_elev_m", "viridis_r", "bed elevation, stage datum (m)", None, None)]):
    ax = fig.add_subplot(gs[:, j]); sc = ax.scatter(loc.x_m, loc.y_m, c=loc[col], cmap=cmap, vmin=vmin, vmax=vmax, s=16, edgecolor=INK, lw=0.2); ax.scatter([0], [0], marker="*", s=150, color=ORANGE, edgecolor="white", lw=0.6, zorder=6)
    fig.colorbar(sc, ax=ax, shrink=0.5, pad=0.02).set_label(label); ax.set_aspect("equal"); ax.set_xlabel("east of the ADCP (m)"); ax.set_ylabel("north of the ADCP (m)"); ax.grid(True); ax.set_axisbelow(True); ax.set_title(["The target field: amplitude per location", "The second field: bed elevation"][j], fontsize=9)
ax = fig.add_subplot(gs[0, 2]); ax.scatter(m.local_depth_m, m.velocity_ms, s=4, color=INK, alpha=0.3, lw=0); hh = np.linspace(m.local_depth_m.min(), m.local_depth_m.max(), 50)
for b, col, lab in [(0.5, C26, "exponent 1/2 (shear velocity, fixed slope)"), (2 / 3, BLUE, "exponent 2/3 (Manning)"), (b_free, ORANGE, f"free fit, exponent {b_free:.2f}")]:
    mae, c = mae_b(b); ax.plot(hh, c * hh ** b, color=col, lw=1.3, label=f"{lab}: {mae:.1f} mm/s")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(f"local depth at the ADCP, stage minus {h0_adcp:.2f} m (m)"); ax.set_ylabel("ADCP velocity (m/s)"); ax.set_title("ADCP: how velocity scales with local depth", fontsize=9); ax.grid(True, which="both"); ax.set_axisbelow(True); ax.legend(fontsize=7)
ax = fig.add_subplot(gs[1, 2])
for yr, col in [(2025, C25), (2026, C26)]: p = loc[loc.year == yr]; ax.errorbar(p.water_depth_m, p.a, yerr=1.96 * p.a_sd, fmt="o", color=col, ms=3, elinewidth=0.5, alpha=0.7, label=f"{yr}")
ax.set_xlabel("water depth at the visit (m)"); ax.set_ylabel("amplitude a (m/s)"); ax.set_title(f"Amplitude against depth across locations (r = {info['locations']['corr_a_vs_depth']:.2f})", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7)
fig.suptitle("Inputs for the spatial model: amplitude and bed elevation per location, and the ADCP's depth scaling", fontsize=9, y=0.95); fig.savefig(F / "00n_amplitude_inputs.png", bbox_inches="tight"); plt.close(fig); print("saved 00n")
