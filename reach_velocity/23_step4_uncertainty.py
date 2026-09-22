"""Step 4. The uncertainty chain in use: calibration on held-out days, intervals on four surveyed sections at four stages,
and the attribution of the interval to each source. Outputs: uncertainty_calibration.json, step4_intervals.csv,
step4_attribution.csv, figures/04_step4_uncertainty.png"""
import json, sys, time, warnings
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from gpcore import *; from uncertainty import *
warnings.filterwarnings("ignore"); t0 = time.time(); log = lambda *a: (print(f"[{time.time()-t0:5.0f}s]", *a), sys.stdout.flush())
df = load_points(); df["bed_m"] = df.bed_rtk_sensor_datum_m.fillna(df.bed_elev_m); ch = Chain(df); json.dump({**ch.calib, "discrepancy_split": ch.disc_split}, open(D / "uncertainty_calibration.json", "w"), indent=2); log(f"calibration: {ch.calib}")
SECTIONS = ["AT10", "0625-T02", "AT0401", "AT13"]; STAGES = [0.75, 0.90, 1.10, 1.30]; N = 400; THR = 0.20
rows, att = [], []
for t in SECTIONS:
    pts = df[df.transect == t].sort_values("across_m").reset_index(drop=True)
    for st in STAGES:
        smp = ch.sample(st, pts, n=N, seed=1); mv, sh = section_stats(smp, THR)
        rows.append({"section": t, "stage_m": st, "mean_v_median": float(np.nanmedian(mv)), "mean_v_lo": float(np.nanquantile(mv, .05)), "mean_v_hi": float(np.nanquantile(mv, .95)), "share_median": float(np.nanmedian(sh)), "share_lo": float(np.nanquantile(sh, .05)), "share_hi": float(np.nanquantile(sh, .95)), "wet_points": int((smp["h"][0] > 0.02).sum())})
        for src in ALL_SOURCES:
            s1 = ch.sample(st, pts, n=N, sources=(src,), seed=2); mv1, sh1 = section_stats(s1, THR)
            att.append({"section": t, "stage_m": st, "source": src, "mean_v_iqr90": float(np.nanquantile(mv1, .95) - np.nanquantile(mv1, .05)), "share_iqr90": float(np.nanquantile(sh1, .95) - np.nanquantile(sh1, .05))})
        log(f"{t} at {st:.2f} m: mean v {rows[-1]['mean_v_median']:.3f} [{rows[-1]['mean_v_lo']:.3f}, {rows[-1]['mean_v_hi']:.3f}]; share < {THR}: {rows[-1]['share_median']:.2f} [{rows[-1]['share_lo']:.2f}, {rows[-1]['share_hi']:.2f}]")
iv = pd.DataFrame(rows); iv.to_csv(D / "step4_intervals.csv", index=False); at = pd.DataFrame(att); at.to_csv(D / "step4_attribution.csv", index=False)
# ---- figure
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, DARK, MUTED = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#0d366b", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig = plt.figure(figsize=(15, 9.5)); gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.32)
t = "AT10"; pts = df[df.transect == t].sort_values("across_m").reset_index(drop=True); cols = {0.75: GREEN, 0.90: BLUE, 1.10: DARK, 1.30: ORANGE}
ax = fig.add_subplot(gs[0, 0])
for st in STAGES:
    smp = ch.sample(st, pts, n=N, seed=3); vb = smp["vbar"]; ax.fill_between(pts.across_m, np.quantile(vb, .05, axis=0), np.quantile(vb, .95, axis=0), color=cols[st], alpha=0.15, lw=0); ax.plot(pts.across_m, np.median(vb, axis=0), color=cols[st], lw=1.3, label=f"stage {st:.2f} m" + (" (above the surveyed range)" if st > 1.1 else ""))
ax.errorbar(pts.across_m, pts.a * (np.log(pts.water_depth_m / Z0) - 1), fmt="s", color=INK, ms=3, label=f"measured on the survey day ({pts.stage_m.median():.2f} m)"); ax.axhline(THR, color=INK2, lw=0.5, ls=":")
ax.set_xlabel("across (m)"); ax.set_ylabel("depth-averaged velocity (m/s)"); ax.set_title(f"{t}: the chain's output with 90% bands, all sources on", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=6.5)
ax = fig.add_subplot(gs[0, 1]); w = 0.2
for k, st in enumerate(STAGES):
    a_ = at[(at.section == t) & (at.stage_m == st)].set_index("source"); ax.bar(np.arange(len(ALL_SOURCES)) + (k - 1.5) * w, [a_.mean_v_iqr90[s_] for s_ in ALL_SOURCES], w, color=cols[st], label=f"{st:.2f} m")
ax.set_xticks(range(len(ALL_SOURCES))); ax.set_xticklabels(["A field\n(GP band)", "discrepancy\n(days held out)", "bed\n(6 cm)", "stage\n(1.5 cm)", "law\n(h0, γ ridge)", "z0\n(7 to 16 mm)"], fontsize=7); ax.set_ylabel("90% interval width of the section mean velocity (m/s)"); ax.set_title(f"{t}: which source drives the band on the section mean velocity", fontsize=9); ax.grid(True, axis="y"); ax.set_axisbelow(True); ax.legend(fontsize=7, title="stage")
ax = fig.add_subplot(gs[0, 2])
for k, st in enumerate(STAGES):
    a_ = at[(at.section == t) & (at.stage_m == st)].set_index("source"); ax.bar(np.arange(len(ALL_SOURCES)) + (k - 1.5) * w, [a_.share_iqr90[s_] * 100 for s_ in ALL_SOURCES], w, color=cols[st])
ax.set_xticks(range(len(ALL_SOURCES))); ax.set_xticklabels(["A field", "discrepancy", "bed", "stage", "law", "z0"], fontsize=7); ax.set_ylabel(f"90% interval width of the share below {THR} m/s (points)"); ax.set_title(f"{t}: the same, for the ecological share", fontsize=9); ax.grid(True, axis="y"); ax.set_axisbelow(True)
ax = fig.add_subplot(gs[1, 0:2]); xs = np.arange(len(SECTIONS)); w2 = 0.2
for k, st in enumerate(STAGES):
    r = iv[iv.stage_m == st].set_index("section").loc[SECTIONS]; ax.errorbar(xs + (k - 1.5) * w2, r.share_median * 100, yerr=[(r.share_median - r.share_lo) * 100, (r.share_hi - r.share_median) * 100], fmt="o", color=cols[st], ms=5, capsize=3, label=f"stage {st:.2f} m")
ax.set_xticks(xs); ax.set_xticklabels(SECTIONS); ax.set_ylabel(f"share of the section below {THR} m/s (%), median and 90% interval"); ax.set_title("Four surveyed sections at four stages: the ecological share with its interval", fontsize=9); ax.grid(True, axis="y"); ax.set_axisbelow(True); ax.legend(fontsize=7)
ax = fig.add_subplot(gs[1, 2]); c = ch.calib; ax.text(0.02, 0.9, "Calibration on held-out days (log amplitude)", fontsize=9, transform=ax.transAxes, fontweight="bold"); ax.text(0.02, 0.72, f"coverage of the GP band alone: {c['held_out_coverage95_before']:.2f}\n95th percentile of |z|: {c['q95_abs_z']:.2f} (1.96 if perfect)\ninflation applied: {c['inflation']:.2f}\ndiscrepancy sd added: {c['disc_sd_log']:.2f} in log\non {c['n']} held-out points", fontsize=8.5, transform=ax.transAxes, va="top"); ax.axis("off")
fig.suptitle("Step 4: the uncertainty chain, calibrated and attributed", fontsize=10, y=0.96); fig.savefig(F / "04_step4_uncertainty.png", bbox_inches="tight"); plt.close(fig); log("saved 04")
print(iv.round(3).to_string(index=False)); print(at.pivot_table(index="source", columns="stage_m", values="mean_v_iqr90").round(3).to_string())
