"""Step 6. Deliverable 1: the spatial rating curves. Every surveyed section's mean depth-averaged velocity and share below
the threshold, every two hours through both seasons from the stage record, with intervals; the rating curve of three
sections with the survey-day measurement on it; a continuity check (discharge through each section on its survey day,
measured and modelled, against the Ann Arbor gauge); and the validation summary gathered from steps 2 to 5.
Outputs: step6_section_series.csv, step6_section_draws.npz, step6_continuity.csv, step6_validation_table.csv, figures/06_step6_rating_curves.png"""
import json, sys, time, warnings
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import matplotlib.dates as mdates
from gpcore import *; from uncertainty import *
warnings.filterwarnings("ignore"); t0 = time.time(); log = lambda *a: (print(f"[{time.time()-t0:5.0f}s]", *a), sys.stdout.flush())
THR = 0.20; N = 100
df = load_points(); df["bed_m"] = df.bed_rtk_sensor_datum_m.fillna(df.bed_elev_m); ch = Chain(df)
sec_pts = {t: g.assign(_x=g.x_m.round(2), _y=g.y_m.round(2)).drop_duplicates(["_x", "_y"]).sort_values("across_m").reset_index(drop=True) for t, g in df.groupby("transect")}
order = sorted(sec_pts, key=lambda t: sec_pts[t].along_m.median()); along = {t: float(sec_pts[t].along_m.median()) for t in order}
g10 = pd.read_csv(D / "series_10min.csv", parse_dates=["time_utc"]).set_index("time_utc"); g10.index = pd.to_datetime(g10.index, utc=True); st2 = g10.stage_m.resample("2h").mean().dropna(); st2 = st2[(st2.index.year >= 2025)]
times, stages = st2.index, st2.values; log(f"{len(order)} sections, {len(times)} two-hourly stages")
secs_all = {t: sec_pts[t] for t in order}; R = sample_reach(ch, stages, secs_all, n=N, seed=7, threshold=THR); rows = []
for t in order:
    mv, sh = R["mean_v"][t], R["share"][t]
    rows.append(pd.DataFrame({"section": t, "time_utc": times, "stage_m": stages, "mean_v_med": np.nanmedian(mv, 0), "mean_v_lo": np.nanquantile(mv, .05, 0), "mean_v_hi": np.nanquantile(mv, .95, 0), "share_med": np.nanmedian(sh, 0), "share_lo": np.nanquantile(sh, .05, 0), "share_hi": np.nanquantile(sh, .95, 0), "wet_frac": R["wet_frac"][t].mean(0)}))
ser = pd.concat(rows, ignore_index=True); ser.to_csv(D / "step6_section_series.csv", index=False)
np.savez_compressed(D / "step6_section_draws.npz", times=times.astype(str).values, stages=stages, reach_share=R["reach_share"], length_weights=np.array([R["length_weights"][t] for t in order]), sections=np.array(order), **{f"mean_v__{t}": R["mean_v"][t] for t in order}, **{f"share__{t}": R["share"][t] for t in order})
log(f"joint reach draws: discrepancy common {R['disc']['common_sd_log']:.3f}, independent {R['disc']['independent_sd_log']:.3f} (log); reach share at the median stage of the record: {np.nanmedian(R['reach_share']):.3f}")
# ---- rating curves for three sections (from one joint draw so the sections' bands are consistent)
grid_st = np.linspace(0.65, 1.35, 36); RC = sample_reach(ch, grid_st, secs_all, n=300, seed=8, threshold=THR); rc = {t: (np.nanmedian(RC["mean_v"][t], 0), np.nanquantile(RC["mean_v"][t], .05, 0), np.nanquantile(RC["mean_v"][t], .95, 0)) for t in ["AT10", "0625-T02", "AT0401"]}
# ---- section-level check: each day held out, the chain refitted without it, its sections' measured mean velocity against the band
cov_rows = []
for day in sorted(df.date.unique()):
    tr_ = df[df.date != day]; ch_d = Chain(tr_, n_restarts=1); te_secs = {t: sec_pts[t] for t in order if sec_pts[t].date.iloc[0] == day}
    for t, p in te_secs.items():
        stage = float(p.stage_m.median()); s1 = sample_reach(ch_d, [stage], {t: p}, n=300, seed=21, threshold=THR); w = width_weights(p); h = p.water_depth_m.values; vb = p.a.values * (np.log(np.clip(h, np.e * Z0, None) / Z0) - 1)
        v_meas = float(np.sum(vb * w) / w.sum()); sh_meas = float(np.sum((vb < THR) * w) / w.sum()); mv = s1["mean_v"][t][:, 0]; sh = s1["share"][t][:, 0]
        cov_rows.append({"section": t, "date": day, "stage_m": stage, "mean_v_measured": v_meas, "mean_v_med": float(np.nanmedian(mv)), "mean_v_lo": float(np.nanquantile(mv, .05)), "mean_v_hi": float(np.nanquantile(mv, .95)), "share_measured": sh_meas, "share_med": float(np.nanmedian(sh)), "share_lo": float(np.nanquantile(sh, .05)), "share_hi": float(np.nanquantile(sh, .95))})
covs = pd.DataFrame(cov_rows); covs["mean_v_in_band"] = (covs.mean_v_measured >= covs.mean_v_lo) & (covs.mean_v_measured <= covs.mean_v_hi); covs["share_in_band"] = (covs.share_measured >= covs.share_lo) & (covs.share_measured <= covs.share_hi); covs.to_csv(D / "step6_section_coverage.csv", index=False)
log(f"section level, day held out: mean velocity inside the 90% band for {covs.mean_v_in_band.mean()*100:.0f}% of {len(covs)} sections (MAE {np.mean(np.abs(covs.mean_v_measured - covs.mean_v_med)):.3f} m/s, median relative error {np.median(np.abs(covs.mean_v_measured - covs.mean_v_med) / covs.mean_v_measured)*100:.0f}%); share below {THR} inside the band for {covs.share_in_band.mean()*100:.0f}% (MAE {np.mean(np.abs(covs.share_measured - covs.share_med))*100:.0f} points)")
# ---- continuity: discharge through each section on its survey day, measured from the profiles and modelled from the chain
qd = pd.read_csv(D / "usgs_daily.csv").set_index("date").discharge_cfs; cont = []
for t in order:
    p = sec_pts[t]; day = p.date.iloc[0]; stage = float(p.stage_m.median()); y = p.across_m.values; dy = np.diff(y); w = np.r_[dy[0] / 2 + min(p.dist_bank_m.iloc[0], 5), (dy[:-1] + dy[1:]) / 2, dy[-1] / 2 + min(p.dist_bank_m.iloc[-1], 5)]
    h = p.water_depth_m.values; vbar = p.a.values * (np.log(np.clip(h, np.e * Z0, None) / Z0) - 1); q_meas = float(np.sum(vbar * h * w))
    hm = np.clip(stage - p.bed_m.values, 0, None); sm = ch.sample(stage, p, n=300, seed=9); qd_ = (sm["vbar"] * hm[None, :] * w[None, :]).sum(1)   # modelled per draw: sum of v h w over the section's points
    q_mod, q_lo, q_hi = float(np.median(qd_)), float(np.quantile(qd_, .05)), float(np.quantile(qd_, .95))
    cont.append({"section": t, "date": day, "stage_m": stage, "Q_measured_m3s": q_meas, "Q_model_m3s": q_mod, "Q_model_lo": q_lo, "Q_model_hi": q_hi, "Q_usgs_m3s": float(qd.get(day, np.nan)) * 0.0283168, "width_m": float(w.sum())})
cont = pd.DataFrame(cont); cont.to_csv(D / "step6_continuity.csv", index=False); log("continuity:\n" + cont.round(2).to_string(index=False))
# the reach rating through the measured discharges: Q = c (stage - h0)^m with h0 from the primary law; its scatter is the place-free consistency of the 20 sections
h0r = ch.primary[0]; br_ = np.polyfit(np.log(cont.stage_m - h0r), np.log(cont.Q_measured_m3s), 1); rr_ = cont.Q_measured_m3s / np.exp(np.polyval(br_, np.log(cont.stage_m - h0r)))
json.dump({"h0": float(h0r), "exponent": float(br_[0]), "coef": float(np.exp(br_[1])), "scatter_pct": float(rr_.std(ddof=1) * 100), "n": int(len(cont))}, open(D / "step6_reach_rating.json", "w"), indent=2)
log(f"reach rating from the measured section discharges: Q = {np.exp(br_[1]):.1f} (stage - {h0r:.2f})^{br_[0]:.2f} m3/s, scatter {rr_.std(ddof=1)*100:.0f}% (n {len(cont)})")
# continuity at common stage: every section's modelled discharge at the same stages (sum v h w per draw), against the measured rating; the between-section spread at a common stage is the test
cs_rows = []; common_st = np.linspace(0.70, 1.30, 13)
for t in order:
    p = sec_pts[t]; y = p.across_m.values; dy = np.diff(y); w = np.r_[dy[0] / 2 + min(p.dist_bank_m.iloc[0], 5), (dy[:-1] + dy[1:]) / 2, dy[-1] / 2 + min(p.dist_bank_m.iloc[-1], 5)]
    for st_ in common_st:
        hm = np.clip(st_ - p.bed_m.values, 0, None); sm = ch.sample(st_, p, n=200, seed=10); qd_ = (sm["vbar"] * hm[None, :] * w[None, :]).sum(1); cs_rows.append({"section": t, "stage_m": st_, "Q_med": float(np.median(qd_)), "Q_lo": float(np.quantile(qd_, .05)), "Q_hi": float(np.quantile(qd_, .95))})
cs = pd.DataFrame(cs_rows); cs.to_csv(D / "step6_continuity_common_stage.csv", index=False); spread_cs = cs.groupby("stage_m").Q_med.agg(lambda x: np.std(np.log(x), ddof=1))
log("continuity at common stage, sd of log Q across the 20 sections: " + ", ".join(f"{st_:.2f} m: {v:.3f}" for st_, v in spread_cs.items()))
# readings on held-out days under the adopted law (the step 3 table quoted the Q law's numbers)
rd_all = pd.read_csv(D / "readings_transect.csv"); dl = with_target(df, ("stage", ch.primary[1], ch.primary[0])); RH = []
for day in sorted(dl.date.unique()):
    m = Model(n_restarts=1).fit(dl[dl.date != day]); te = dl[dl.date == day]; mu, sd = m.predict(te); st_ = np.sqrt(sd ** 2 + te.sd_y.values ** 2); aa = dict(zip(te.location_id, np.exp(mu) * np.clip(te.stage_m.values - ch.primary[0], 0.02, None) ** ch.primary[1])); ss = dict(zip(te.location_id, st_))
    r = rd_all[rd_all.location_id.isin(te.location_id) & (rd_all.height_m > Z0 * np.e)].copy(); r["v_pred"] = [aa[l] * np.log(h / Z0) for l, h in zip(r.location_id, r.height_m)]; r["v_lo"] = [aa[l] * np.exp(-1.96 * ss[l]) * np.log(h / Z0) for l, h in zip(r.location_id, r.height_m)]; r["v_hi"] = [aa[l] * np.exp(1.96 * ss[l]) * np.log(h / Z0) for l, h in zip(r.location_id, r.height_m)]; RH.append(r)
RH = pd.concat(RH); RH.to_csv(D / "step6_readings_heldout_adopted_law.csv", index=False); rh_mae = float(np.mean(np.abs(RH.velocity_ms - RH.v_pred))); rh_cov = float(np.mean((RH.velocity_ms >= RH.v_lo) & (RH.velocity_ms <= RH.v_hi))); rh_mad = float(np.mean(np.abs(RH.velocity_ms - RH.velocity_ms.median())))
below = RH.velocity_ms < THR; rh_hit = float(np.mean((RH.v_pred < THR) == below)); rh_base = float(max(below.mean(), 1 - below.mean())); log(f"readings on held-out days, adopted law: MAE {rh_mae:.3f} m/s (mean absolute deviation of the readings {rh_mad:.3f}), coverage {rh_cov:.2f}, below/above {THR} called right {rh_hit*100:.0f}% (base {rh_base*100:.0f}%), n {len(RH)}")
same_day = cont.groupby("date").filter(lambda g: len(g) > 1); spread = same_day.groupby("date").Q_measured_m3s.agg(lambda x: x.max() / x.min()); ratio = (cont.Q_measured_m3s / cont.Q_usgs_m3s)
log(f"same-day sections: measured discharge ratio max/min by day {spread.round(2).to_dict()}; measured/USGS median {ratio.median():.2f} (IQR {ratio.quantile(.25):.2f} to {ratio.quantile(.75):.2f}); model/USGS median {(cont.Q_model_m3s / cont.Q_usgs_m3s).median():.2f}")
# ---- validation table gathered
s2 = pd.read_csv(D / "step2_metrics.csv"); s3 = pd.read_csv(D / "step3_metrics.csv"); s3b = json.load(open(D / "step3b_adcp_shape.json")); s3s = json.load(open(D / "step3_summary.json")); cal = json.load(open(D / "uncertainty_calibration.json")); tot = pd.read_csv(D / "step5_reach_totals.csv")
def g2(model, scheme, col): return float(s2[(s2.model == model) & (s2.scheme == scheme)][col].iloc[0])
val = [{"test": "section held out: CRPS of log A, model / reach mean", "value": f"{g2('full features + notes', 'transect_out', 'CRPS_log'):.3f} / {g2('reach mean', 'transect_out', 'CRPS_log'):.3f}", "note": "20 sections, 257 points; 95% band covers " + f"{g2('full features + notes', 'transect_out', 'coverage95'):.2f}"},
       {"test": "day held out: CRPS of log a, stage law / no law / oracle", "value": f"{float(s3b['best_stage_law']['CRPS']):.3f} / {float(s3[(s3.law == 'none') & (s3.scheme == 'day_out')].CRPS.iloc[0]):.3f} / {float(s3[(s3.law == 'oracle') & (s3.scheme == 'day_out')].CRPS.iloc[0]) if 'oracle' in set(s3.law) else float('nan'):.3f}", "note": "11 days; six days hold out a single section, so this fold set overlaps the one above"},
       {"test": "readings on held-out days: MAE / their mean absolute deviation (m/s)", "value": f"{rh_mae:.3f} / {rh_mad:.3f}", "note": f"{len(RH):,} readings, adopted law; 95% band covers {rh_cov:.2f}; below/above {THR} m/s called right {rh_hit*100:.0f}% (base rate {rh_base*100:.0f}%)"},
       {"test": "sections on held-out days: mean velocity inside the 90% band", "value": f"{covs.mean_v_in_band.mean()*100:.0f}% of {len(covs)}", "note": f"MAE {np.mean(np.abs(covs.mean_v_measured - covs.mean_v_med)):.3f} m/s; share below {THR} inside its band {covs.share_in_band.mean()*100:.0f}%"},
       {"test": "ADCP through 2025: shape of the stage law, level fitted", "value": f"{[v for k, v in s3b['adcp_shape'].items() if 'reach' in k and 'point' in k][0]['log_ratio_sd_over_season']:.3f} (sd of log ratio)", "note": "point reading; the level is not predicted by the field; under a column reading the shape is depth^0.46 (instrument model open)"},
       {"test": "continuity: same-day sections' measured discharge, max/min", "value": ", ".join(f"{v:.2f}" for v in spread.values), "note": "a consistency check on the measurements, not on the model"},
       {"test": "continuity at common stage: sd of log Q across the 20 sections", "value": f"{spread_cs[0.90]:.3f} at 0.90 m, {spread_cs[1.30]:.3f} at 1.30 m", "note": "modelled discharge of every section at the same stage; the spread is the model's between-section inconsistency"},
       {"test": "continuity: measured section discharge / gauge, median", "value": f"{ratio.median():.2f}", "note": "USGS Ann Arbor 18 km downstream; two tributary-storm days at 0.4 to 0.5"},
       {"test": "reach share below 0.2 m/s at 0.90 m: map / sections", "value": f"{float(tot[(tot.stage_m == 0.9) & tot.source.str.startswith('map')].share_below_thr_median.iloc[0]):.3f} / {float(tot[(tot.stage_m == 0.9) & tot.source.str.startswith('sections')].share_below_thr_median.iloc[0]):.3f}", "note": "the map carries the margins to the bank line; the sections stop at their outermost points"}]
val = pd.DataFrame(val); val.to_csv(D / "step6_validation_table.csv", index=False)
# ---- figure
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, DARK, MUTED = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#0d366b", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig = plt.figure(figsize=(16, 12)); gs = fig.add_gridspec(3, 3, height_ratios=[1.1, 1.1, 1.2], hspace=0.5, wspace=0.28)
for r_, yr in enumerate([2025, 2026]):
    ax = fig.add_subplot(gs[r_, 0:2]); m = ser[ser.time_utc.dt.year == yr]; mat = m.pivot(index="section", columns="time_utc", values="mean_v_med").loc[order]
    im = ax.imshow(mat.values, aspect="auto", cmap="Blues", vmin=0, vmax=0.6, extent=[mdates.date2num(mat.columns[0]), mdates.date2num(mat.columns[-1]), len(order) - 0.5, -0.5]); ax.xaxis_date(); ax.set_yticks(range(len(order))); ax.set_yticklabels([f"{t} ({along[t]:+.0f} m)" for t in order], fontsize=6)
    ax2 = ax.twinx(); st_y = st2[st2.index.year == yr]; ax2.plot(st_y.index, st_y.values, color=ORANGE, lw=0.8); ax2.set_ylabel("stage (m)", color=ORANGE); ax2.set_ylim(0.5, 1.8)
    ax.set_title(f"{yr}: section mean velocity every two hours from the stage record (rows: the surveyed sections, upstream at the top)", fontsize=9); ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    if r_ == 1: fig.colorbar(im, ax=ax, shrink=0.8, pad=0.08).set_label("section mean depth-averaged velocity (m/s)")
ax = fig.add_subplot(gs[0:2, 2])
for t, col in zip(["AT10", "0625-T02", "AT0401"], [BLUE, GREEN, DARK]):
    md, lo, hi = rc[t]; ax.fill_between(grid_st, lo, hi, color=col, alpha=0.12, lw=0); ax.plot(grid_st, md, color=col, lw=1.4, label=f"{t}"); p = sec_pts[t]; vb = p.a.values * (np.log(np.clip(p.water_depth_m.values, np.e * Z0, None) / Z0) - 1); ax.plot(p.stage_m.median(), vb.mean(), "s", color=col, ms=7, markeredgecolor="white")
ax.axvspan(0.74, 1.09, color=GRID, alpha=0.5, lw=0); ax.text(0.745, 0.93, "surveyed stage range", fontsize=7, color=INK2); ax.set_xlabel("stage (m)"); ax.set_ylabel("section mean depth-averaged velocity (m/s)"); ax.set_title("Three sections' rating curves with 90% bands\n(squares: measured on the survey day)", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7, loc="upper left"); ax.set_ylim(0, 1.0)
ax = fig.add_subplot(gs[2, 0]); ax.plot([2, 30], [2, 30], color=AXIS, lw=0.8); ax.errorbar(cont.Q_usgs_m3s, cont.Q_model_m3s, yerr=[cont.Q_model_m3s - cont.Q_model_lo, cont.Q_model_hi - cont.Q_model_m3s], fmt="o", color=BLUE, ms=4, elinewidth=0.5, label="modelled at the survey stage"); ax.scatter(cont.Q_usgs_m3s, cont.Q_measured_m3s, s=18, color=INK, marker="s", label="measured from the profiles")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("USGS discharge at Ann Arbor on the day (m³/s), 18 km downstream"); ax.set_ylabel("discharge through the section (m³/s)"); ax.set_title(f"Continuity: measured / gauge median {ratio.median():.2f}; same-day sections agree to a factor {spread.max():.2f} at worst", fontsize=8.5); ax.grid(True, which="both"); ax.set_axisbelow(True); ax.legend(fontsize=7)
ax = fig.add_subplot(gs[2, 1:3]); ax.axis("off"); ax.set_title("Validation summary (from steps 2 to 6)", fontsize=9)
for i, r in val.iterrows(): ax.text(0.0, 0.95 - i * 0.11, f"{r.test}: {r.value}\n   {r.note}", fontsize=7, transform=ax.transAxes, va="top"); ax.text(0.0, 0.0, "", fontsize=7.5, transform=ax.transAxes, va="top", ha="left")
fig.suptitle("Step 6: the spatial rating curves, delivered from one stage sensor", fontsize=10, y=0.95); fig.savefig(F / "06_step6_rating_curves.png", bbox_inches="tight"); plt.close(fig); log("saved 06"); print(val.round(3).to_string(index=False))
