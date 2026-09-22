"""What the model learns now, shown before any GP is fitted.
(a) Figure 00r: the measured amplitude at every location, the local depth at the visit, and the de-staged conveyance
    A = a / h^beta, which is the one field left to learn; A across each transect; and a test of the de-staging on pairs of
    neighbouring transects surveyed on different days at different stages: which beta makes their A profiles agree.
(b) Figure 00s: what one number per location generates, on one surveyed cross-section: the measured bed, the water surface
    at three stages, and the depth-averaged velocity across the section at those stages from A and the law, with the
    survey-day values; and the vertical profile at one point at the three stages.
Outputs: figures/00r_conveyance_field.png, 00s_what_it_generates.png, data_clean/destage_pair_check.csv"""
import json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"; F = HERE / "figures"
z0 = json.load(open(D / "field_qc.json"))["z0_m"]; BETA = 0.5
loc = pd.read_csv(D / "locations.csv"); loc = loc[loc.a > 0].copy()
loc["h"] = loc.water_depth_m; loc["A"] = loc.a / loc.h ** BETA; loc["logA"] = np.log(loc.A)
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, DARK, MUTED = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#0d366b", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
# ---------------------------------------------------------------- pair check: neighbouring transects, different days
pairs = [("0618-T01", "ACORNT01"), ("0618-T01", "AT03"), ("0618-T02", "AT03"), ("0618-T02", "AT0401"), ("0618-T03", "AT11"), ("0618-T03", "AT12"), ("AT03", "AT05"), ("AT06", "AT07"), ("AT06", "AT08")]
def profile(t, col):
    p = loc[loc.transect == t].sort_values("across_m"); return p.across_m.values, p[col].values, p.stage_m.median(), p.date.iloc[0]
rows = []
for t1, t2 in pairs:
    x1, _, s1, d1 = profile(t1, "a"); x2, _, s2, d2 = profile(t2, "a")
    if d1 == d2: continue
    lo, hi = max(x1.min(), x2.min()), min(x1.max(), x2.max()); xs = np.linspace(lo, hi, 15)
    for beta in [0.0, 0.5, 2 / 3, 1.0]:
        p1 = loc[loc.transect == t1].sort_values("across_m"); p2 = loc[loc.transect == t2].sort_values("across_m")
        A1 = np.interp(xs, p1.across_m, p1.a / p1.h ** beta); A2 = np.interp(xs, p2.across_m, p2.a / p2.h ** beta)
        rows.append({"pair": f"{t1} vs {t2}", "days": f"{d1} vs {d2}", "stage1": s1, "stage2": s2, "beta": beta, "log_ratio_median": float(np.median(np.log(np.clip(A1, 1e-3, None) / np.clip(A2, 1e-3, None)))), "abs_log_ratio_mean": float(np.mean(np.abs(np.log(np.clip(A1, 1e-3, None) / np.clip(A2, 1e-3, None)))))})
pc = pd.DataFrame(rows); pc.to_csv(D / "destage_pair_check.csv", index=False)
summ = pc.groupby("beta").agg(mean_abs_log_ratio=("abs_log_ratio_mean", "mean"), mean_log_ratio=("log_ratio_median", "mean"), n_pairs=("pair", "nunique")).round(3); print("pair check, how far apart neighbouring transects from different days are, by beta (log units; 0.1 = 10%):\n", summ.to_string())
# day effect on the field: does log a or log A still depend on the day's stage across locations?
for col in ["a", "A"]:
    b = np.polyfit(np.log(loc.stage_m), np.log(loc[col]), 1)[0]; print(f"slope of log {col} on log(stage at visit) across all locations: {b:.2f}")
# ---------------------------------------------------------------- figure 00r
fig = plt.figure(figsize=(15, 10)); gs = fig.add_gridspec(2, 4, height_ratios=[1.6, 1], hspace=0.35, wspace=0.35)
for j, (col, cmap, lab, vmax) in enumerate([("a", "Blues", "measured amplitude a (m/s)", 0.2), ("h", "viridis", "depth at the visit (m)", 1.0), ("A", "Blues", "conveyance A = a / sqrt(h)", 0.25)]):
    ax = fig.add_subplot(gs[0, j]); sc = ax.scatter(loc.x_m, loc.y_m, c=loc[col], cmap=cmap, vmin=0, vmax=vmax, s=14, edgecolor=INK, lw=0.15); ax.scatter([0], [0], marker="*", s=120, color=ORANGE, edgecolor="white", lw=0.5, zorder=5)
    fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02).set_label(lab); ax.set_aspect("equal"); ax.grid(True); ax.set_axisbelow(True); ax.set_title(["What was measured: amplitude on the day", "What is measured too: depth on the day", "What is learned: the conveyance field"][j], fontsize=9); ax.set_xlabel("east of the ADCP (m)")
    if j == 0: ax.set_ylabel("north of the ADCP (m)")
ax = fig.add_subplot(gs[0, 3])
for t, p in loc.sort_values("across_m").groupby("transect"): ax.plot(p.across_m, p.A, "-", lw=0.7, alpha=0.75, color=BLUE if p.year.iloc[0] == 2025 else GREEN)
ax.set_xlabel("across (m), + = left bank"); ax.set_ylabel("conveyance A"); ax.set_title("A across every transect (blue 2025, green 2026)", fontsize=9); ax.grid(True); ax.set_axisbelow(True)
ax = fig.add_subplot(gs[1, 0:2]); dd = loc.groupby("date").agg(stage=("stage_m", "median"), a=("a", "median"), A=("A", "median"), a25=("a", lambda x: x.quantile(.25)), a75=("a", lambda x: x.quantile(.75)), A25=("A", lambda x: x.quantile(.25)), A75=("A", lambda x: x.quantile(.75))).sort_values("stage")
ax.errorbar(dd.stage, dd.a, yerr=[dd.a - dd.a25, dd.a75 - dd.a], fmt="o", color=MUTED, ms=5, capsize=2, label="amplitude a, median and quartiles of the day's locations"); ax.errorbar(dd.stage + 0.004, dd.A, yerr=[dd.A - dd.A25, dd.A75 - dd.A], fmt="s", color=BLUE, ms=5, capsize=2, label="conveyance A = a / sqrt(h)")
ax.set_xlabel("stage on the survey day (m)"); ax.set_ylabel("m/s"); ax.set_title("By survey day: the amplitude carries the day's stage, the conveyance should not", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7)
ax = fig.add_subplot(gs[1, 2:4]); s_ = summ.reset_index(); ax.bar(range(len(s_)), s_.mean_abs_log_ratio * 100, color=[MUTED, BLUE, DARK, GREEN]); ax.set_xticks(range(len(s_))); ax.set_xticklabels([f"β = {b:g}" + ("\n(no de-staging)" if b == 0 else "\n(shear law)" if b == 0.5 else "\n(Manning)" if b > 0.6 and b < 0.7 else "") for b in s_.beta])
ax.set_ylabel("mean |log ratio| between neighbours (%)"); ax.set_title(f"Neighbouring transects surveyed on different days ({int(s_.n_pairs.iloc[0])} pairs): how far apart their fields are, by β", fontsize=9); ax.grid(True, axis="y"); ax.set_axisbelow(True)
fig.suptitle("What is being learned: not a number, a field. The measured amplitude carries the day's stage; dividing by the depth law leaves the conveyance, one value per place", fontsize=9, y=0.95)
fig.savefig(F / "00r_conveyance_field.png", bbox_inches="tight"); plt.close(fig); print("saved 00r")
# ---------------------------------------------------------------- figure 00s: what one number per location generates
T = "AT10"; p = loc[loc.transect == T].sort_values("across_m"); bed = p.bed_rtk_sensor_datum_m.fillna(p.bed_elev_m).values; x = p.across_m.values; s_visit = p.stage_m.median()
stages = [0.75, s_visit, 1.10]; cols = [GREEN, BLUE, DARK]
fig = plt.figure(figsize=(14, 5.2)); gs = fig.add_gridspec(1, 3, width_ratios=[1.4, 1.4, 1], wspace=0.3)
ax = fig.add_subplot(gs[0, 0]); ax.fill_between(x, bed.min() - 0.1, bed, color="#d9d3c7", lw=0); ax.plot(x, bed, color=INK, lw=1.2, label="measured bed (RTK, sensor datum)")
for s, c in zip(stages, cols): ax.axhline(s, color=c, lw=1.0, ls="--", label=f"water surface at stage {s:.2f} m" + (" (survey day)" if s == s_visit else ""))
ax.set_xlabel("across (m)"); ax.set_ylabel("elevation on the sensor datum (m)"); ax.set_title(f"Transect {T}: the measured cross-section and three stages", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7, loc="upper right"); ax.set_ylim(bed.min() - 0.1, 1.25)
ax = fig.add_subplot(gs[0, 1])
for s, c in zip(stages, cols):
    h = np.clip(s - bed, 0, None); a = p.A.values * h ** BETA; vmean = np.where(h > z0 * 2.72, a * (np.log(np.clip(h, z0 * 2.72, None) / z0) - 1), 0.0); ax.plot(x, vmean, "-o", color=c, ms=3, lw=1.2, label=f"stage {s:.2f} m")
ax.errorbar(x, p.v_mean, yerr=1.96 * p.v_mean_sd.fillna(0.03), fmt="s", color=INK, ms=4, elinewidth=0.6, capsize=1.5, label="measured on the survey day"); ax.axhline(0.2, color=INK2, lw=0.5, ls=":")
ax.set_xlabel("across (m)"); ax.set_ylabel("depth-averaged velocity (m/s)"); ax.set_title("What A generates: velocity across the section at each stage (law + measured bed)", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7)
ax = fig.add_subplot(gs[0, 2]); k = int(np.argmax(p.A.values)); Ak, bk = p.A.values[k], bed[k]
for s, c in zip(stages, cols):
    h = s - bk; a = Ak * h ** BETA; z = np.linspace(z0 * 1.05, h, 60); ax.plot(a * np.log(z / z0), z, color=c, lw=1.3, label=f"stage {s:.2f} m, depth {h:.2f} m")
q = pd.read_csv(D / "readings_transect.csv"); q = q[q.location_id == p.location_id.values[k]]; ax.scatter(q.velocity_ms, q.height_m, s=12, color=INK, zorder=4, label="readings on the survey day")
ax.set_xlabel("velocity (m/s)"); ax.set_ylabel("height above the bed (m)"); ax.set_title(f"One point ({p.point.values[k]}): the profile at each stage", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7)
fig.suptitle("From one learned number per place to everything the paper needs: the law and the measured bed do the rest", fontsize=9, y=1.0)
fig.savefig(F / "00s_what_it_generates.png", bbox_inches="tight"); plt.close(fig); print("saved 00s")
