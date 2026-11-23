"""Step 3. Transfer across stage: which law carries A from one day to another? Leave one survey day out, fit the GP on the
other days' de-staged amplitudes, predict A at the held-out day's points, re-stage with that day's own depth (or
discharge), compare with the measured amplitudes and with the readings. Laws: none; local depth to a power beta; discharge
at Ann Arbor to a power m. Oracle: the spatial-only error that remains after each held-out day's own median offset is
removed, which bounds what any day law can achieve. Then the per-day offsets against stage, discharge and year; the ADCP's
season predicted from the field; the two revisited pairs.
Outputs: step3_metrics.csv, step3_day_offsets.csv, step3_predictions.csv, step3_adcp.csv, step3_summary.json, figures/03_step3_transfer.png"""
import json, sys, time, warnings
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import matplotlib.dates as mdates
from gpcore import *
warnings.filterwarnings("ignore"); t0 = time.time(); log = lambda *a: (print(f"[{time.time()-t0:5.0f}s]", *a), sys.stdout.flush())
df = load_points(); days = sorted(df.date.unique()); log(f"{len(df)} points, {len(days)} days")
LAWS = {"none": ("none",), "depth^0.25": ("depth", 0.25), "depth^0.5": ("depth", 0.5), "depth^0.67": ("depth", 2 / 3), "depth^1.0": ("depth", 1.0), "depth^1.5": ("depth", 1.5), "Q^0.25": ("Q", 0.25), "Q^0.5": ("Q", 0.5), "Q^0.75": ("Q", 0.75), "Q^1.0": ("Q", 1.0), "Q^1.5": ("Q", 1.5)}
rows, offs, preds = [], [], []
for lname, law in LAWS.items():
    d = with_target(df, law); MU, SD = np.full(len(d), np.nan), np.full(len(d), np.nan)
    for day in days:
        tr, te = np.where(d.date != day)[0], np.where(d.date == day)[0]; m = Model().fit(d.iloc[tr]); mu, sd = m.predict(d.iloc[te]); MU[te], SD[te] = mu, sd
        off = float(np.median(d.y.values[te] - mu)); offs.append({"law": lname, "date": day, "offset_log": off, "stage_m": float(d.stage_m.iloc[te].median()), "Q_cfs": float(d.Q_cfs.iloc[te].iloc[0]), "year": int(d.year.iloc[te].iloc[0]), "n": len(te)})
    sd_tot = np.sqrt(SD ** 2 + d.sd_y.values ** 2); y_a, mu_a = np.log(d.a.values), MU + np.log(d.g.values)     # back to the amplitude on the day
    sc = scores(y_a, mu_a, sd_tot); rows.append({"law": lname, "scheme": "day_out", **sc})
    o = pd.DataFrame([x for x in offs if x["law"] == lname]).set_index("date").offset_log; y_or = y_a - d.date.map(o).values; sc_or = scores(y_or, mu_a, sd_tot); rows.append({"law": lname, "scheme": "day_out_oracle_offsets", **sc_or})
    log(f"{lname:12s} day-out MAE log a {sc['MAE']:.3f} CRPS {sc['CRPS']:.3f} cov {sc['coverage95']:.2f} | after removing each day's offset (oracle): MAE {sc_or['MAE']:.3f}")
    for i in range(len(d)): preds.append((lname, d.location_id.iat[i], d.date.iat[i], d.transect.iat[i], float(y_a[i]), float(mu_a[i]), float(sd_tot[i])))
met = pd.DataFrame(rows); met.to_csv(D / "step3_metrics.csv", index=False); off_df = pd.DataFrame(offs); off_df.to_csv(D / "step3_day_offsets.csv", index=False)
pr = pd.DataFrame(preds, columns=["law", "location_id", "date", "transect", "log_a", "mu_log_a", "sd"]); pr.to_csv(D / "step3_predictions.csv", index=False)
best = met[met.scheme == "day_out"].sort_values("CRPS").law.iloc[0]; log(f"best law by day-out CRPS: {best}")
# how the day offsets relate to stage, Q and year, for 'none' and the best law
rel = {}
for lname in ["none", best]:
    o = off_df[off_df.law == lname]; rel[lname] = {"sd_of_day_offsets": float(o.offset_log.std()), "spearman_vs_stage": float(spearmanr(o.stage_m, o.offset_log)[0]), "spearman_vs_Q": float(spearmanr(o.Q_cfs, o.offset_log)[0]), "mean_offset_2025": float(o[o.year == 2025].offset_log.mean()), "mean_offset_2026": float(o[o.year == 2026].offset_log.mean())}
    log(f"day offsets, {lname}: sd {rel[lname]['sd_of_day_offsets']:.2f}; vs stage rho {rel[lname]['spearman_vs_stage']:+.2f}; vs Q rho {rel[lname]['spearman_vs_Q']:+.2f}; 2025 mean {rel[lname]['mean_offset_2025']:+.2f}, 2026 mean {rel[lname]['mean_offset_2026']:+.2f}")
# reading-level velocities on held-out days, best law and none
tr_ = pd.read_csv(D / "readings_transect.csv"); tr_ = tr_[(tr_.height_m >= 0.03) & (tr_.velocity_ms > 0.005)]
rl = {}
for lname in ["none", best]:
    p = pr[pr.law == lname].set_index("location_id"); q = tr_[tr_.location_id.isin(p.index)].copy(); q["v_pred"] = np.exp(p.loc[q.location_id, "mu_log_a"].values) * np.log(q.height_m / Z0); q["v_lo"] = np.exp(p.loc[q.location_id, "mu_log_a"].values - 1.96 * p.loc[q.location_id, "sd"].values) * np.log(q.height_m / Z0); q["v_hi"] = np.exp(p.loc[q.location_id, "mu_log_a"].values + 1.96 * p.loc[q.location_id, "sd"].values) * np.log(q.height_m / Z0)
    rl[lname] = {"MAE_m_s": float(np.mean(np.abs(q.v_pred - q.velocity_ms))), "bias_m_s": float(np.mean(q.v_pred - q.velocity_ms)), "coverage95": float(np.mean((q.velocity_ms >= q.v_lo) & (q.velocity_ms <= q.v_hi))), "n_readings": int(len(q)), "sd_of_readings": float(q.velocity_ms.std())}
    if lname == best: q.to_csv(D / "step3_readings_heldout.csv", index=False)
    log(f"readings on held-out days, {lname}: MAE {rl[lname]['MAE_m_s']:.3f} m/s, bias {rl[lname]['bias_m_s']:+.3f}, coverage {rl[lname]['coverage95']:.2f}, spread of readings {rl[lname]['sd_of_readings']:.3f}")
# ADCP: A at its position from the field fitted on all days (best law), then its season
law_b = LAWS[best]
near = df.iloc[(df.x_m ** 2 + df.y_m ** 2).argmin()]; feat_adcp = {f: float(near[f]) for f in FULL}; feat_adcp["along_m"] = 0.0
import banklines_3dep as B
SF = B.StageFeatures("Middle"); df = SF.training(df); d_all = with_target(df, law_b); m_all = Model().fit(d_all)   # bank set C (entries 30 and 31): fixed lines, so the bank features at the ADCP are the same at every stage
ad = pd.read_csv(D / "adcp_amplitude_2h.csv", parse_dates=["time_utc"]); h0 = json.load(open(D / "amplitude_inputs.json"))["h0_adcp_m_column_fit"]; ad["h"] = ad.stage_m - h0; ad = ad[ad.h > Z0 * 3].copy()
LV = SF.grid(float(ad.stage_m.min()), float(ad.stage_m.max())); MUs, SDs = [], []
for lv in LV:
    dist_, rel_, _ = SF.features([13244857.0], [319995.5], lv); f_ = dict(feat_adcp, dist_bank_m=float(dist_[0]), rel_across=float(rel_[0]))
    mu_, sd_ = m_all.predict(pd.DataFrame([{**f_, **{c: 0.0 for c in NOTE_COLS}}])); MUs.append(mu_[0]); SDs.append(sd_[0])
muA_t, sdA_t = np.interp(ad.stage_m, LV, MUs), np.interp(ad.stage_m, LV, SDs); A_adcp = float(np.exp(np.interp(0.90, LV, MUs))); sdA = float(np.interp(0.90, LV, SDs))   # summary numbers at the reference stage 0.90 m
feat_adcp["dist_bank_m"], feat_adcp["rel_across"] = [float(x[0]) for x in SF.features([13244857.0], [319995.5], 0.90)[:2]]
if law_b[0] == "depth": g_t = ad.h ** law_b[1]
elif law_b[0] == "Q":
    qd = pd.read_csv(D / "usgs_daily.csv"); qmap = qd.set_index("date").discharge_cfs; g_t = (ad.time_utc.dt.strftime("%Y-%m-%d").map(qmap).values / 350.0) ** law_b[1]
else: g_t = np.ones(len(ad))
ad["a_pred"] = np.exp(muA_t) * g_t; ad["v_pred"] = ad.a_pred * (np.log(ad.h / Z0) - 1); ad["v_lo"] = ad.v_pred * np.exp(-1.96 * sdA_t); ad["v_hi"] = ad.v_pred * np.exp(1.96 * sdA_t)
ad["a_obs"] = ad.velocity_ms / (np.log(ad.h / Z0) - 1); adcp_res = {"A_adcp_from_field": A_adcp, "sd_logA": sdA, "nearest_point": str(near.location_id), "dist_bank_m": feat_adcp["dist_bank_m"], "MAE_m_s": float(np.mean(np.abs(ad.v_pred - ad.velocity_ms))), "bias_m_s": float(np.mean(ad.v_pred - ad.velocity_ms)), "coverage95": float(np.mean((ad.velocity_ms >= ad.v_lo) & (ad.velocity_ms <= ad.v_hi))), "ratio_obs_over_pred_median": float(np.median(ad.velocity_ms / ad.v_pred))}
ad.to_csv(D / "step3_adcp.csv", index=False); log(f"ADCP from the field ({best}): A {A_adcp:.3f} (sd log {sdA:.2f}); season MAE {adcp_res['MAE_m_s']*1000:.0f} mm/s, bias {adcp_res['bias_m_s']*1000:+.0f}, coverage {adcp_res['coverage95']:.2f}, obs/pred median {adcp_res['ratio_obs_over_pred_median']:.2f}")
# revisited pairs under each law
pairs = {}
for (t1, t2) in [("AT03", "AT05")]:
    for lname, law in LAWS.items():
        d = with_target(df, law); p1 = d[d.transect == t1].sort_values("across_m"); p2 = d[d.transect == t2].sort_values("across_m"); xs = np.linspace(max(p1.across_m.min(), p2.across_m.min()), min(p1.across_m.max(), p2.across_m.max()), 15)
        pairs[f"{t1}_vs_{t2}_{lname}"] = float(np.median(np.log(np.interp(xs, p1.across_m, p1.A) / np.interp(xs, p2.across_m, p2.A))))
rv = df[df.point.isin(["AT0710", "AT0711"])]
for lname, law in LAWS.items():
    d = with_target(rv, law); r = []
    for pt, g in d.groupby("point"):
        g = g.sort_values("date"); r.append(np.log(g.A.iloc[0] / g.A.iloc[1]))
    pairs[f"AT0710_11_jul09_vs_jul14_{lname}"] = float(np.mean(r))
summ = {"best_law": best, "day_offsets_relations": rel, "readings_heldout": rl, "adcp": adcp_res, "pairs_log_ratio": pairs, "scales_best_all_days": m_all.scales()}; json.dump(summ, open(D / "step3_summary.json", "w"), indent=2)
# ---------------------------------------------------------------- figure
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, DARK, MUTED = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#0d366b", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig = plt.figure(figsize=(15, 10)); gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.3)
ax = fig.add_subplot(gs[0, 0]); m1 = met[met.scheme == "day_out"].set_index("law"); m2 = met[met.scheme == "day_out_oracle_offsets"].set_index("law")
dl_ = [k for k in LAWS if k.startswith("depth")]; ql = [k for k in LAWS if k.startswith("Q")]
ax.plot([0] + [LAWS[k][1] for k in dl_], [m1.CRPS["none"]] + [m1.CRPS[k] for k in dl_], "o-", color=BLUE, label="local depth law, exponent β"); ax.plot([0] + [LAWS[k][1] for k in ql], [m1.CRPS["none"]] + [m1.CRPS[k] for k in ql], "s-", color=ORANGE, label="discharge law, exponent m")
ax.axhline(m2.CRPS.min(), color=INK, lw=0.8, ls="--"); ax.text(0.02, m2.CRPS.min() + 0.005, "oracle: each day's own offset removed (spatial error only)", fontsize=7, color=INK2)
ax.set_xlabel("exponent"); ax.set_ylabel("day-out CRPS of log amplitude"); ax.set_title("Which law carries A across days?", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7)
for k, (lname, ttl) in enumerate([("none", "Day offsets with no law"), (best, f"Day offsets after the best law ({best})")]):
    ax = fig.add_subplot(gs[0, 1 + k]); o = off_df[off_df.law == lname]
    for yr, col in [(2025, BLUE), (2026, GREEN)]: oo = o[o.year == yr]; ax.scatter(oo.stage_m, oo.offset_log, s=oo.Q_cfs / 6, color=col, alpha=0.8, edgecolor="white", lw=0.5, label=f"{yr} (marker size: discharge)")
    for _, r in o.iterrows(): ax.annotate(r.date[5:], (r.stage_m, r.offset_log), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=6, color=INK2)
    ax.axhline(0, color=INK, lw=0.6); ax.set_xlabel("stage on the day (m)"); ax.set_ylabel("median log(measured / predicted) on the held-out day"); ax.set_title(f"{ttl}: sd {rel[lname]['sd_of_day_offsets']:.2f}, vs stage rho {rel[lname]['spearman_vs_stage']:+.2f}", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7)
ax = fig.add_subplot(gs[1, 0:2]); ax.plot(ad.time_utc, ad.velocity_ms, color=INK, lw=0.6, alpha=0.7, label="ADCP, smoothed record"); ax.fill_between(ad.time_utc, ad.v_lo, ad.v_hi, color=BLUE, alpha=0.15, lw=0); ax.plot(ad.time_utc, ad.v_pred, color=BLUE, lw=1.0, label=f"predicted from the field's A at the ADCP ({best}), 95% band")
ax.set_ylabel("column velocity at the ADCP (m/s)"); ax.set_title(f"The ADCP's season from the field: MAE {adcp_res['MAE_m_s']*1000:.0f} mm/s, bias {adcp_res['bias_m_s']*1000:+.0f} mm/s, coverage {adcp_res['coverage95']:.2f}, A at the ADCP {A_adcp:.3f}", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7); ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
ax = fig.add_subplot(gs[1, 2]); q = pd.read_csv(D / "step3_readings_heldout.csv"); ax.scatter(q.velocity_ms, q.v_pred, s=4, color=BLUE, alpha=0.35, lw=0); ax.plot([0, 0.9], [0, 0.9], color=AXIS, lw=0.8)
ax.set_xlabel("measured velocity, reading on a held-out day (m/s)"); ax.set_ylabel("predicted (m/s)"); ax.set_title(f"Readings on held-out days ({best}): MAE {rl[best]['MAE_m_s']:.3f} m/s, coverage {rl[best]['coverage95']:.2f}, n = {rl[best]['n_readings']}", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.set_xlim(0, 0.9); ax.set_ylim(0, 0.9)
fig.suptitle("Step 3: transfer across stage, one survey day held out at a time", fontsize=10, y=0.96); fig.savefig(F / "03_step3_transfer.png", bbox_inches="tight"); plt.close(fig); log("saved 03")
print(met.round(3).to_string(index=False)); print(off_df[off_df.law.isin(["none", best])].pivot(index="date", columns="law", values="offset_log").round(2).to_string())
