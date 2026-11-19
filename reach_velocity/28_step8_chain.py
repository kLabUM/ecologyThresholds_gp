"""Step 8b. What transfers to the other two reaches, tested the same way as the Middle (and the Middle again under the same protocol).
(1) The spatial part: the conveyance field from log a with no day term, every transect held out in turn; scored in log a and in
velocity (the depth-averaged velocity at the survey-day depth) against the reach mean. (2) The day part: each day held out, the day's
mean residual is the level a stage law has to explain; those levels against the local stage from the RTK water surface (tied days) and
against the nearest gauge. (3) Bed matching between sections surveyed twice: the RTK datum check. No season series for the reaches
without a sensor: the RTK shows the gauge does not give their water level. Outputs: step8_summary.csv, step8_transect_out.csv,
step8_day_levels.csv, step8_bed_matching.csv, figures/08_step8_transfer.png
Bank features: set C (entry 30) for all three reaches, so this script runs after 37_ (the Upper and Lower sensor ties)."""
import json, sys, time, warnings
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from gpcore import *; from uncertainty import banks_for
warnings.filterwarnings("ignore"); t0 = time.time(); log = lambda *a: (print(f"[{time.time()-t0:5.0f}s]", *a), sys.stdout.flush())
THR = 0.20; Z0 = {"Middle": json.load(open(D / "field_qc.json"))["z0_m"]}
def target(df):
    d = df[df.a > 0].copy(); d["y"] = np.log(d.a); d["sd_y"] = np.clip(d.a_sd / d.a, 0.02, 1.0); return d.reset_index(drop=True)
reaches, daystage, gauge_name = {}, {}, {"Middle": "USGS Ann Arbor, 18 km downstream", "Upper": "USGS near Hamburg, 6 km downstream", "Lower": "USGS Ann Arbor, 42 km upstream, three dams between"}
mid = load_points(); reaches["Middle"] = target(mid); daystage["Middle"] = mid.groupby("date").stage_m.mean()
for r in ["Upper", "Lower"]:
    reaches[r] = target(banks_for(r).training(pd.read_csv(D / f"reach_{r}_points.csv")));   # bank features from set C at the survey-day sensor level (entry 30); needs 37_ first
    Z0[r] = json.load(open(D / f"reach_{r}_qc.json"))["z0_m"]; dd = pd.read_csv(D / f"reach_{r}_days.csv"); daystage[r] = dd.set_index("date").stage_local_m
summary, to_rows, dl_rows = [], [], []
for reach, d in reaches.items():
    z0 = Z0[reach]; days = sorted(d.date.unique()); vfac = np.log(np.clip(d.water_depth_m, np.e * z0, None) / z0) - 1
    # ---- (1) every transect held out; no day term
    preds = []
    for t in sorted(d.transect.unique()):
        te = d[d.transect == t]; tr = d[d.transect != t]; m = Model(n_restarts=1).fit(tr); mu, sd = m.predict(te); st = np.sqrt(sd ** 2 + te.sd_y.values ** 2)
        preds.append(pd.DataFrame({"reach": reach, "transect": t, "date": te.date.values, "location_id": te.location_id.values, "y": te.y.values, "mu": mu, "sd": st, "vfac": vfac[te.index].values, "h": te.water_depth_m.values}))
    P = pd.concat(preds).reset_index(drop=True); to_rows.append(P); rm = float(np.average(d.y, weights=1 / d.sd_y ** 2)); rsd = float(d.y.std())
    s_model = scores(P.y.values, P.mu.values, P.sd.values); s_reach = scores(P.y.values, np.full(len(P), rm), np.full(len(P), rsd))
    v_meas = np.exp(P.y) * P.vfac; v_pred = np.exp(P.mu) * P.vfac; v_reach = np.exp(rm) * P.vfac; wet = P.h > np.e * z0
    mae_v, mae_v_reach = float(np.mean(np.abs(v_meas - v_pred)[wet])), float(np.mean(np.abs(v_meas - v_reach)[wet])); below = v_meas < THR; hit = float(np.mean((v_pred < THR) == below)); hit_reach = float(np.mean((v_reach < THR) == below)); base = float(max(below.mean(), 1 - below.mean()))
    # ---- (2) every day held out; the day's mean residual under the spatial-only model
    for day in days:
        te = d[d.date == day]; tr = d[d.date != day]
        if len(tr) < 20: continue
        m = Model(n_restarts=1).fit(tr); mu, sd = m.predict(te); res = te.y.values - mu; w = 1 / (sd ** 2 + te.sd_y.values ** 2)
        dl_rows.append({"reach": reach, "date": day, "n_points": int(len(te)), "n_transects": int(te.transect.nunique()), "offset": float(np.average(res, weights=w)), "offset_se": float(np.sqrt(1 / w.sum())), "Q_cfs": float(te.Q_cfs.iloc[0]) if "Q_cfs" in te else np.nan,
                        "stage_local_m": float(daystage[reach].get(day, np.nan)), "tied": bool(te.tied.all()) if "tied" in te else True, "along_median_m": float(te.along_m.median()), "gauge": gauge_name[reach]})
    DL = pd.DataFrame([r_ for r_ in dl_rows if r_["reach"] == reach]); ok = DL.stage_local_m.notna() & DL.Q_cfs.notna()
    rs = spearmanr(DL.stage_local_m[ok], DL.offset[ok])[0] if ok.sum() >= 4 else np.nan; rq = spearmanr(np.log(DL.Q_cfs[ok]), DL.offset[ok])[0] if ok.sum() >= 4 else np.nan
    slope_st = np.polyfit(DL.stage_local_m[ok], DL.offset[ok], 1)[0] if ok.sum() >= 4 else np.nan   # d(log a)/d(stage): the Middle's law gives 1/(stage - h0) ~ 2 per m at 0.95 m
    mf = Model(n_restarts=2).fit(d); sc = mf.scales()
    summary.append({"reach": reach, "points": int(len(d)), "transects": int(d.transect.nunique()), "days": len(days), "z0_mm": z0 * 1000, "a_median": float(d.a.median()), "share_below_thr_measured": float(below.mean()),
                    "CRPS_model": s_model["CRPS"], "CRPS_reach_mean": s_reach["CRPS"], "MAE_log_model": s_model["MAE"], "MAE_log_reach_mean": s_reach["MAE"], "coverage95_model": s_model["coverage95"],
                    "MAE_v_model_m_s": mae_v, "MAE_v_reach_mean_m_s": mae_v_reach, "hit_rate_below_thr_model": hit, "hit_rate_reach_mean": hit_reach, "base_rate": base,
                    "sd_day_offsets": float(DL.offset.std()), "n_days_tied": int(ok.sum()), "spearman_offset_vs_local_stage": rs, "spearman_offset_vs_lnQ": rq, "dlog_a_per_m_of_stage": slope_st, "sd_logA_total": rsd, **{f"scale_{k}": v for k, v in sc.items()}})
    log(f"{reach}: transect-out CRPS {s_model['CRPS']:.3f} vs reach mean {s_reach['CRPS']:.3f}, coverage {s_model['coverage95']:.2f}; velocity MAE {mae_v:.3f} vs {mae_v_reach:.3f} m/s; hit rate below {THR}: {hit:.2f} vs reach mean {hit_reach:.2f} (base {base:.2f}); day offsets sd {DL.offset.std():.2f}, vs local stage rho {rs:.2f} (n {int(ok.sum())}, {slope_st:.1f} per m), vs ln Q rho {rq:.2f}; " + ", ".join(f"{k} {v:.2g}" for k, v in sc.items()))
S = pd.DataFrame(summary); S.to_csv(D / "step8_summary.csv", index=False); pd.concat(to_rows).to_csv(D / "step8_transect_out.csv", index=False); DL_all = pd.DataFrame(dl_rows); DL_all.to_csv(D / "step8_day_levels.csv", index=False)
# ---- (3) bed matching between sections surveyed twice (Upper, Lower)
bm = []
for reach in ["Upper", "Lower"]:
    loc = pd.read_csv(D / f"reach_{reach}_points.csv"); T = pd.read_csv(D / f"reach_{reach}_transects.csv").sort_values("along_m").reset_index(drop=True)
    for i in range(len(T)):
        for j in range(i + 1, len(T)):
            A, B = T.iloc[i], T.iloc[j]
            if A.date == B.date or abs(A.along_m - B.along_m) > 5: continue
            pa = loc[loc.transect == A.transect].sort_values("across_m"); pb = loc[loc.transect == B.transect].sort_values("across_m"); lo, hi = max(pa.across_m.min(), pb.across_m.min()), min(pa.across_m.max(), pb.across_m.max())
            if hi - lo < 3 or len(pb) < 3 or len(pa) < 3: continue
            q = pa[(pa.across_m >= lo) & (pa.across_m <= hi)]; dz = q.elev_rtk_m.values - np.interp(q.across_m, pb.across_m, pb.elev_rtk_m)
            bm.append({"reach": reach, "A": A.transect, "B": B.transect, "gap_m": abs(A.along_m - B.along_m), "n": len(q), "bed_A_minus_B_m": float(np.median(dz)), "mad_m": float(1.4826 * np.median(np.abs(dz - np.median(dz)))), "wse_A_minus_B_m": A.wse_rtk_m - B.wse_rtk_m, "Q_A": A.Q_cfs, "Q_B": B.Q_cfs})
BM = pd.DataFrame(bm); BM.to_csv(D / "step8_bed_matching.csv", index=False); print(BM.round(3).to_string(index=False))
# ---- figure
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, DARK, MUTED = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#0d366b", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig = plt.figure(figsize=(19, 12)); gs = fig.add_gridspec(3, 4, wspace=0.32, hspace=0.55, width_ratios=[1.5, 1, 1, 1]); P_all = pd.concat(to_rows)
for i, reach in enumerate(["Middle", "Upper", "Lower"]):
    d = reaches[reach]; days = sorted(d.date.unique()); row = S[S.reach == reach].iloc[0]; DL = DL_all[DL_all.reach == reach]; cmap = plt.get_cmap("tab10"); dcol = {day: cmap(k % 10) for k, day in enumerate(days)}
    ax = fig.add_subplot(gs[i, 0]); v = d.a * (np.log(np.clip(d.water_depth_m, np.e * Z0[reach], None) / Z0[reach]) - 1); sc_ = ax.scatter(d.along_m, d.across_m, c=v, cmap="Blues", s=9, vmin=0, vmax=0.6, edgecolor=INK, lw=0.15); fig.colorbar(sc_, ax=ax, shrink=0.8, pad=0.02).set_label("depth-averaged velocity on the survey day (m/s)")
    for t, g in d.groupby("transect"): ax.plot([g.along_m.median()] * 2, [g.across_m.min(), g.across_m.max()], color=dcol[g.date.iloc[0]], lw=0.6, alpha=0.6)
    ax.set_xlabel("along the reach (m, downstream)"); ax.set_ylabel("across (m)"); ax.set_title(f"{reach}, straightened: {len(d)} locations, {d.transect.nunique()} transects, {len(days)} days over {d.along_m.max()-d.along_m.min():.0f} m; z0 {Z0[reach]*1000:.0f} mm; {row.share_below_thr_measured*100:.0f}% of points below {THR} m/s", fontsize=8); ax.grid(True); ax.set_axisbelow(True)
    ax = fig.add_subplot(gs[i, 1]); P = P_all[P_all.reach == reach]; ax.errorbar(P.y, P.mu, yerr=1.96 * P.sd, fmt="o", ms=2.5, color=BLUE, ecolor=BLUE, elinewidth=0.3, alpha=0.55); lim = [P.y.min() - 0.3, P.y.max() + 0.3]; ax.plot(lim, lim, color=AXIS, lw=0.8); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("measured log a"); ax.set_ylabel("predicted with the transect held out (no day term)"); ax.set_title(f"Transect held out: CRPS {row.CRPS_model:.3f} against the reach mean {row.CRPS_reach_mean:.3f}; coverage {row.coverage95_model:.2f}\nvelocity MAE {row.MAE_v_model_m_s:.3f} against {row.MAE_v_reach_mean_m_s:.3f} m/s; below {THR} m/s called right {row.hit_rate_below_thr_model*100:.0f}% (base rate {row.base_rate*100:.0f}%)", fontsize=7.5); ax.grid(True); ax.set_axisbelow(True)
    ax = fig.add_subplot(gs[i, 2]); q = DL[DL.stage_local_m.notna()]
    for _, r_ in q.iterrows(): ax.errorbar(r_.stage_local_m, r_.offset, yerr=1.96 * r_.offset_se, fmt="o", color=dcol[r_.date], ms=5, elinewidth=0.6); ax.annotate(r_.date[5:], (r_.stage_local_m, r_.offset), fontsize=6, xytext=(4, 3), textcoords="offset points")
    if len(q) >= 4: xx = np.linspace(q.stage_local_m.min(), q.stage_local_m.max(), 10); b = np.polyfit(q.stage_local_m, q.offset, 1); ax.plot(xx, np.polyval(b, xx), color=INK, lw=0.7, ls="--")
    ax.set_xlabel("local water level on the day (m): the stage sensor" if reach == "Middle" else "local water level on the day from the RTK surface (m NAVD88; tied days)"); ax.set_ylabel("day held out: mean residual of log a")
    ax.set_title(f"The day level against the local water level: Spearman {row.spearman_offset_vs_local_stage:.2f}, {row.dlog_a_per_m_of_stage:.1f} per m (n {row.n_days_tied})\nday levels have sd {row.sd_day_offsets:.2f} in log a; the Middle's law gives about 2 per m", fontsize=7.5); ax.grid(True); ax.set_axisbelow(True)
    ax = fig.add_subplot(gs[i, 3]); q = DL[DL.Q_cfs.notna()]
    for _, r_ in q.iterrows(): ax.errorbar(r_.Q_cfs * 0.0283168, r_.offset, yerr=1.96 * r_.offset_se, fmt="o", mfc=dcol[r_.date] if r_.tied else "white", mec=dcol[r_.date], ecolor=dcol[r_.date], ms=5, elinewidth=0.6); ax.annotate(r_.date[5:], (r_.Q_cfs * 0.0283168, r_.offset), fontsize=6, xytext=(4, 3), textcoords="offset points")
    ax.set_xscale("log"); ax.set_xlabel(f"{gauge_name[reach]} (m³/s)"); ax.set_ylabel("day held out: mean residual of log a"); ax.set_title(f"The same day levels against the nearest gauge: Spearman {row.spearman_offset_vs_lnQ:.2f}\n(open: days whose RTK surface could not be tied)", fontsize=7.5); ax.grid(True, which="both"); ax.set_axisbelow(True)
fig.suptitle("Step 8: what transfers. The spatial field is testable in every reach; the day level needs the local water level, which the sensor gives and the gauge does not", fontsize=10, y=0.94); fig.savefig(F / "08_step8_transfer.png", bbox_inches="tight"); plt.close(fig); log("saved 08")
print(S.round(3).T.to_string())
