"""Step 9b. The chain on the Upper and Lower reaches with their sensors. The sensor sits outside each surveyed stretch in a damped spot, so
the sections' water level is tied to it by a linear 'sensor-to-section rating' fitted on the tied survey days (the RTK water surface,
slope-corrected to a reference position, against the sensor stage), with the fit's scatter as the stage error. The bed stays on the RTK
datum (NAVD88). The law, the chain and the validations are the Middle's, with the reach's own roughness scale and law set.
Outputs per reach: reach_{r}_chain_points.csv, reach_{r}_law_grid_local.csv, reach_{r}_section_series.csv, reach_{r}_reach_series.csv,
reach_{r}_continuity.csv, reach_{r}_chain_summary.json, figures/15_{r}_sensor_chain.png"""
import json, sys, time, warnings
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import matplotlib.dates as mdates
import uncertainty as U
from gpcore import *; from uncertainty import *
warnings.filterwarnings("ignore"); t0 = time.time(); log = lambda *a: (print(f"[{time.time()-t0:5.0f}s]", *a), sys.stdout.flush())
THR = 0.20; nb = json.load(open(D / "nest_bayes.json")); BV = nb["velocity_coef_per_m_s"]["posterior"]   # the Middle's nest-learned velocity effect, transferred
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, DARK, MUTED = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#0d366b", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
def run(reach):
    qc = json.load(open(D / f"reach_{reach}_qc.json")); sq = json.load(open(D / f"reach_{reach}_sensor_qc.json")); pts = pd.read_csv(D / f"reach_{reach}_points_sensor.csv"); T = pd.read_csv(D / f"reach_{reach}_tie.csv")
    g10 = pd.read_csv(D / f"reach_{reach}_stage_10min.csv", parse_dates=["time_utc"]).set_index("time_utc"); g10.index = pd.to_datetime(g10.index, utc=True); sm = g10.stage_m.dropna()
    z0 = qc["z0_m"]; U.Z0 = z0; U.Z0_RANGE = tuple(qc["z0_within_5pct_of_min_m"]); S = float(pts.slope_cm_per_100m.iloc[0]) / 1e4; along_ref = float(pts.along_m.median()); out = {"reach": reach, "z0_mm": z0 * 1000}
    # ---- the sensor-to-section rating: the sections' water level at the reference position (RTK, slope-corrected) against the sensor stage, tied days
    good = T[T.rtk_good & T.stage_m.notna()].copy(); good["L"] = good.wse_rtk_m + S * (good.along_m - along_ref); day = good.groupby("date").agg(s=("stage_m", "median"), L=("L", "mean"), n=("transect", "size"))
    b = np.polyfit(day.s, day.L, 1); day["L_fit"] = np.polyval(b, day.s); sdL = float((day.L - day.L_fit).std(ddof=2)); U.STAGE_SD = sdL
    out["sensor_to_section_rating"] = {"L = a + b s": [float(b[1]), float(b[0])], "scatter_sd_m": sdL, "n_days": int(len(day)), "days": {d: [float(r.s), float(r.L), float(r.L_fit)] for d, r in day.iterrows()}, "meaning": "L is the water level at the reference position on the RTK datum; s the sensor stage; the sections move b metres per metre of sensor stage"}
    log(f"{reach}: sections' level L = {b[1]:.2f} + {b[0]:.2f} s (n {len(day)} days, scatter {sdL*100:.0f} cm); z0 {z0*1000:.1f} mm; slope {S*1e4:.1f} cm/100 m")
    # ---- the points on the RTK datum: the level at each reading's time from the rating, the bed from the RTK where the fix is good
    pts["L_m"] = np.polyval(b, pts.stage_m); pts["bed_rtk_datum_m"] = np.where(pts.rtk_good, pts.elev_rtk_m, pts.L_m - S * (pts.along_m - along_ref) - pts.water_depth_m)
    pts["bed_m"] = pts.bed_rtk_datum_m + S * (pts.along_m - along_ref)   # relative to the reference position's level, so that depth = L - bed_m
    pts["stage_m"] = pts.L_m; ref = float(pts.L_m.median()); pts["h_ref_m"] = np.clip(ref - pts.bed_m, 0, None); gs_ = pts.groupby("transect"); pts["rel_depth"] = (pts.h_ref_m / gs_.h_ref_m.transform("max").replace(0, np.nan)).fillna(0); pts["section_mean_depth_m"] = gs_.h_ref_m.transform("mean")
    d = pts[pts.a > 0].reset_index(drop=True); d.to_csv(D / f"reach_{reach}_chain_points.csv", index=False)
    # ---- the law grid in L
    rows = []; Lmin = float(np.polyval(b, sm.min()))
    for h0 in np.round(np.arange(Lmin - 0.5, Lmin + 0.01, 0.1), 2):
        for gam in [0.5, 0.75, 1.0, 1.25]:
            dd = with_target(d, ("stage", gam, h0)); Y, MU, SD = [], [], []
            for dy in sorted(dd.date.unique()):
                tr, te = dd[dd.date != dy], dd[dd.date == dy]
                if len(tr) < 20: continue
                m = Model(n_restarts=1).fit(tr); mu, sd = m.predict(te); Y.append(te.y.values); MU.append(mu); SD.append(np.sqrt(sd ** 2 + te.sd_y.values ** 2))
            rows.append({"h0": float(h0), "gamma": gam, **scores(np.concatenate(Y), np.concatenate(MU), np.concatenate(SD))})
    lg = pd.DataFrame(rows).sort_values("CRPS"); lg.to_csv(D / f"reach_{reach}_law_grid_local.csv", index=False); best = lg.iloc[0]; ridge = [(float(r.h0), float(r.gamma)) for _, r in lg[lg.CRPS <= best.CRPS + 0.01].iterrows()]; primary = (float(best.h0), float(best.gamma))
    dn = with_target(d, ("none",)); Y, MU, SD = [], [], []
    for dy in sorted(dn.date.unique()):
        tr, te = dn[dn.date != dy], dn[dn.date == dy]
        if len(tr) < 20: continue
        m = Model(n_restarts=1).fit(tr); mu, sd = m.predict(te); Y.append(te.y.values); MU.append(mu); SD.append(np.sqrt(sd ** 2 + te.sd_y.values ** 2))
    crps_none = scores(np.concatenate(Y), np.concatenate(MU), np.concatenate(SD))["CRPS"]; out["law"] = {"primary": primary, "set": ridge, "CRPS_best": float(best.CRPS), "CRPS_none": float(crps_none), "L_range_survey": [float(d.L_m.min()), float(d.L_m.max())], "L_range_record": [float(np.polyval(b, sm.min())), float(np.polyval(b, sm.max()))]}
    log(f"   law in L: best {primary} CRPS {best.CRPS:.3f} vs none {crps_none:.3f}; set of {len(ridge)}")
    # ---- the chain: transect-out with the law, calibration, section-level day-out coverage
    ch = Chain(d, primary=primary, ridge=ridge, n_restarts=2, banks=banks_for(reach)); dl = with_target(d, ("stage", primary[1], primary[0])); Y, MU, SD = [], [], []
    for t in sorted(dl.transect.unique()):
        te = dl[dl.transect == t]; m = Model(n_restarts=1).fit(dl[dl.transect != t]); mu, sd = m.predict(te); Y.append(te.y.values); MU.append(mu); SD.append(np.sqrt(sd ** 2 + te.sd_y.values ** 2))
    Y, MU, SD = map(np.concatenate, (Y, MU, SD)); rm = float(np.average(dl.y, weights=1 / dl.sd_y ** 2)); s_to = scores(Y, MU, SD); s_rm = scores(Y, np.full(len(Y), rm), np.full(len(Y), float(dl.y.std())))
    out["transect_out"] = {"CRPS_model": s_to["CRPS"], "CRPS_reach_mean": s_rm["CRPS"], "coverage95": s_to["coverage95"], "n": int(len(Y))}; out["calibration"] = {**ch.calib, "discrepancy_split": ch.disc_split}
    secs = {t: g.assign(_x=g.x_m.round(2), _y=g.y_m.round(2)).drop_duplicates(["_x", "_y"]).sort_values("across_m").reset_index(drop=True) for t, g in d.groupby("transect")}; order = sorted(secs, key=lambda t: secs[t].along_m.median())
    cov_rows = []
    for dy in sorted(d.date.unique()):
        tr_ = d[d.date != dy]
        if len(tr_) < 20 or tr_.date.nunique() < 3: continue
        ch_d = Chain(tr_, primary=primary, ridge=ridge, n_restarts=1, banks=banks_for(reach))
        for t in [t for t in order if secs[t].date.iloc[0] == dy]:
            p = secs[t]; L = float(p.L_m.median()); s1 = sample_reach(ch_d, [L], {t: p}, n=300, seed=21, threshold=THR); w = width_weights(p); h = p.water_depth_m.values; vb = p.a.values * (np.log(np.clip(h, np.e * z0, None) / z0) - 1)
            cov_rows.append({"section": t, "date": dy, "L_m": L, "mean_v_measured": float(np.sum(vb * w) / w.sum()), "mean_v_med": float(np.nanmedian(s1["mean_v"][t])), "mean_v_lo": float(np.nanquantile(s1["mean_v"][t], .05)), "mean_v_hi": float(np.nanquantile(s1["mean_v"][t], .95))})
    cv = pd.DataFrame(cov_rows); cv["in_band"] = (cv.mean_v_measured >= cv.mean_v_lo) & (cv.mean_v_measured <= cv.mean_v_hi); cv.to_csv(D / f"reach_{reach}_section_coverage.csv", index=False)
    out["section_day_out"] = {"n": int(len(cv)), "in_band": float(cv.in_band.mean()), "MAE_m_s": float(np.mean(np.abs(cv.mean_v_measured - cv.mean_v_med))), "median_relative_error": float(np.median(np.abs(cv.mean_v_measured - cv.mean_v_med) / cv.mean_v_measured))}
    log(f"   transect-out CRPS {s_to['CRPS']:.3f} vs reach mean {s_rm['CRPS']:.3f} (cov {s_to['coverage95']:.2f}); discrepancy common {ch.disc_split['common_sd_log']:.2f} independent {ch.disc_split['independent_sd_log']:.2f}; sections on held-out days in band {cv.in_band.mean()*100:.0f}% of {len(cv)} (MAE {out['section_day_out']['MAE_m_s']:.3f} m/s, median rel err {out['section_day_out']['median_relative_error']*100:.0f}%)")
    # ---- continuity: the measured section discharges against L, and a rating through them
    cont = []
    for t in order:
        p = secs[t]
        if len(p) < 3: continue
        y = p.across_m.values; dy_ = np.diff(y); w = np.r_[dy_[0] / 2 + min(p.dist_bank_m.iloc[0], 5), (dy_[:-1] + dy_[1:]) / 2, dy_[-1] / 2 + min(p.dist_bank_m.iloc[-1], 5)]; h = p.water_depth_m.values; vb = p.a.values * (np.log(np.clip(h, np.e * z0, None) / z0) - 1)
        cont.append({"section": t, "date": p.date.iloc[0], "L_m": float(p.L_m.median()), "Q_measured_m3s": float(np.sum(vb * h * w)), "width_m": float(w.sum())})
    cont = pd.DataFrame(cont); ok = cont.Q_measured_m3s > 0; xr = np.log(np.clip(cont.L_m[ok] - primary[0], 0.02, None)); br = np.polyfit(xr, np.log(cont.Q_measured_m3s[ok]), 1); resid = np.log(cont.Q_measured_m3s[ok]) - np.polyval(br, xr); cont.to_csv(D / f"reach_{reach}_continuity.csv", index=False)
    same = cont.groupby("date").filter(lambda g: len(g) > 1).groupby("date").Q_measured_m3s.agg(lambda x: x.max() / x.min()); out["continuity"] = {"rating": f"Q = {np.exp(br[1]):.1f} (L - {primary[0]:.2f})^{br[0]:.2f}", "scatter_pct": float(resid.std(ddof=2) * 100), "same_day_max_over_min": {k: float(v) for k, v in same.items()}}
    log(f"   continuity: {out['continuity']['rating']} m3/s, scatter {resid.std(ddof=2)*100:.0f}%; same-day sections max/min {dict((k, round(v, 2)) for k, v in same.items())}")
    # ---- the season series: L(t) from the sensor through the rating; joint draws; the share below the threshold and the mean suitability E
    st2 = sm.resample("2h").mean().dropna(); Lt = np.polyval(b, st2.values); R = sample_reach(ch, Lt, secs, n=100, seed=7, threshold=THR); rows = []
    for t in order: mv, sh = R["mean_v"][t], R["share"][t]; rows.append(pd.DataFrame({"section": t, "time_utc": st2.index, "L_m": Lt, "mean_v_med": np.nanmedian(mv, 0), "mean_v_lo": np.nanquantile(mv, .05, 0), "mean_v_hi": np.nanquantile(mv, .95, 0), "share_med": np.nanmedian(sh, 0), "share_lo": np.nanquantile(sh, .05, 0), "share_hi": np.nanquantile(sh, .95, 0)}))
    pd.concat(rows).to_csv(D / f"reach_{reach}_section_series.csv", index=False); rs = R["reach_share"]
    # E with the Middle's nest-learned curve: exp(b_v v) averaged over the wetted width; per draw one b_v
    rng = np.random.default_rng(3); bvd = rng.normal(BV[0], BV[1], 100); E = np.full((100, len(Lt)), np.nan); names = order; ww = {t: width_weights(secs[t]) for t in names}; lw = R["length_weights"]
    Lsurv = (float(d.L_m.min()), float(d.L_m.max())); in_range = (Lt >= Lsurv[0]) & (Lt <= Lsurv[1])
    # E needs per-point velocities: recompute per draw from the section stats is not possible, so approximate E by the share-weighted form: E_t = sum_sections lw * mean over points exp(bv v); computed with a second joint draw at 30 stages and interpolated in L
    grid_L = np.linspace(Lt.min(), Lt.max(), 30); pts_all = pd.concat([secs[t].assign(_sec=t) for t in names], ignore_index=True); Eg = np.zeros((100, len(grid_L)))
    for j, L in enumerate(grid_L):
        smp = ch.sample(L, pts_all, n=100, seed=40 + j); v = smp["vbar"]; wet = smp["h"] > np.e * z0
        acc = np.zeros(100); tot = np.zeros(100)
        for t in names:
            ix = np.where(pts_all._sec.values == t)[0]; w = ww[t][None, :] * wet[:, ix]; Sv = np.exp(bvd[:, None] * v[:, ix]); acc += lw[t] * np.where(w.sum(1) > 0, (Sv * w).sum(1) / np.maximum(w.sum(1), 1e-12), np.nan); tot += lw[t]
        Eg[:, j] = acc / tot
    for k in range(100): E[k] = np.interp(Lt, grid_L, Eg[k])
    reach_ser = pd.DataFrame({"time_utc": st2.index, "sensor_stage_m": st2.values, "L_m": Lt, "in_surveyed_range": in_range, "share_med": np.nanmedian(rs, 0), "share_lo": np.nanquantile(rs, .05, 0), "share_hi": np.nanquantile(rs, .95, 0), "E_med": np.nanmedian(E, 0), "E_lo": np.nanquantile(E, .05, 0), "E_hi": np.nanquantile(E, .95, 0)}); reach_ser.to_csv(D / f"reach_{reach}_reach_series.csv", index=False)
    out["season"] = {}
    for yr in sorted(set(st2.index.year)):
        r_ = reach_ser[reach_ser.time_utc.dt.year == yr]; out["season"][str(yr)] = {"share_median": float(np.nanmedian(r_.share_med)), "share_median_in_range": float(np.nanmedian(r_[r_.in_surveyed_range].share_med)) if r_.in_surveyed_range.any() else np.nan, "fraction_hours_in_range": float(r_.in_surveyed_range.mean()), "band_halfwidth_points_in_range": float(np.nanmedian((r_[r_.in_surveyed_range].share_hi - r_[r_.in_surveyed_range].share_lo) / 2) * 100) if r_.in_surveyed_range.any() else np.nan, "E_median": float(np.nanmedian(r_.E_med))}
    log(f"   season: {json.dumps(out['season'])}")
    json.dump(out, open(D / f"reach_{reach}_chain_summary.json", "w"), indent=1, default=str)
    # ---- figure
    fig = plt.figure(figsize=(18, 11)); gs = fig.add_gridspec(3, 3, hspace=0.55, wspace=0.3); years = sorted(set(st2.index.year)); yr_main = years[-1]
    ax = fig.add_subplot(gs[0, 0:2]); ax.plot(st2.index, Lt, color=INK, lw=0.5)
    for dy in sorted(d.date.unique()): ax.axvline(pd.Timestamp(dy, tz="UTC"), color=BLUE, lw=0.7)
    ax.axhspan(*Lsurv, color=GRID, alpha=0.6, lw=0); ax.set_ylabel("sections' water level L (m, RTK datum)"); ax.set_title(f"{reach}: the sections' level from the sensor through the rating L = {b[1]:.1f} + {b[0]:.2f} s (scatter {sdL*100:.0f} cm); blue: survey days; grey: surveyed range", fontsize=8.5); ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %y")); ax.grid(True); ax.set_axisbelow(True)
    ax = fig.add_subplot(gs[0, 2]); ax.scatter(day.s, day.L, s=26, color=BLUE, edgecolor=INK, lw=0.3, zorder=3)
    for dy, r_ in day.iterrows(): ax.annotate(dy[5:], (r_.s, r_.L), fontsize=6, xytext=(3, 3), textcoords="offset points")
    xx = np.linspace(sm.min(), sm.max(), 20); ax.plot(xx, np.polyval(b, xx), color=INK, lw=0.8, ls="--"); ax.set_xlabel("sensor stage (m)"); ax.set_ylabel("sections' level L (m)"); ax.set_title(f"The sensor-to-section rating on the tied days\n(the sections move {b[0]:.1f} m per m of sensor stage)", fontsize=8); ax.grid(True); ax.set_axisbelow(True)
    ax = fig.add_subplot(gs[1, 0:2]); s_ = pd.concat(rows); s_ = s_[s_.time_utc.dt.year == yr_main]; piv = s_.pivot(index="section", columns="time_utc", values="mean_v_med").reindex(order); tt = piv.columns; xe = mdates.date2num(tt.tz_convert(None).to_pydatetime()); xe = np.r_[xe, xe[-1] + (xe[-1] - xe[-2])]
    im = ax.pcolormesh(xe, np.arange(len(order) + 1), piv.values, cmap="Blues", vmin=0, vmax=0.6, shading="flat"); ax.set_yticks(np.arange(len(order)) + 0.5); ax.set_yticklabels([f"{t} ({secs[t].along_m.median():.0f} m)" for t in order], fontsize=5.5); ax.invert_yaxis(); ax.xaxis_date(); ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b")); fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02).set_label("section mean velocity (m/s)"); ax.set_title(f"{yr_main}: section mean velocity every two hours (rows: sections, upstream at the top)", fontsize=8.5)
    ax = fig.add_subplot(gs[1, 2])
    for h0v, grp in lg.groupby("h0"): grp = grp.sort_values("gamma"); ax.plot(grp.gamma, grp.CRPS, marker="o", ms=3, lw=0.8, label=f"h0 = {h0v:.2f}")
    ax.axhline(crps_none, color=ORANGE, lw=0.8, ls="--"); ax.axhspan(best.CRPS, best.CRPS + 0.01, color=GRID, alpha=0.7, lw=0); ax.set_xlabel("gamma"); ax.set_ylabel("CRPS, day held out"); ax.set_title(f"The law in L: best {best.CRPS:.3f}, no law {crps_none:.3f}; set of {len(ridge)}\ntransect-out CRPS {s_to['CRPS']:.3f} vs reach mean {s_rm['CRPS']:.3f}", fontsize=8); ax.legend(fontsize=5.5, ncol=2); ax.grid(True); ax.set_axisbelow(True)
    for k, yr in enumerate(years[-2:]):
        ax = fig.add_subplot(gs[2, k]); r_ = reach_ser[reach_ser.time_utc.dt.year == yr]; ax.fill_between(r_.time_utc, r_.share_lo * 100, r_.share_hi * 100, color=BLUE, alpha=0.2, lw=0); ax.plot(r_.time_utc, r_.share_med * 100, color=BLUE, lw=0.8, label=f"share below {THR} m/s")
        ax.fill_between(r_.time_utc, r_.E_lo * 100, r_.E_hi * 100, color=DARK, alpha=0.12, lw=0); ax.plot(r_.time_utc, r_.E_med * 100, color=DARK, lw=0.8, ls="--", label="mean suitability E (Middle's nest curve)"); oo = r_[~r_.in_surveyed_range]; ax.scatter(oo.time_utc, np.full(len(oo), 1.5), s=2, marker="|", color=INK2, lw=0.5)
        ax.set_ylim(0, 100); ax.set_ylabel("%"); ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b")); sy = out["season"][str(yr)]; ax.set_title(f"{yr}: share below {THR} m/s (median {sy['share_median']*100:.0f}%; ±{sy['band_halfwidth_points_in_range']:.0f} points in range) and E (median {sy['E_median']:.2f}); ticks: outside the surveyed range", fontsize=7.5); ax.grid(True); ax.set_axisbelow(True)
        if k == 0: ax.legend(fontsize=6.5, loc="upper left")
    if len(years) == 1: fig.add_subplot(gs[2, 1]).axis("off")
    ax = fig.add_subplot(gs[2, 2]); ax.scatter(cont.L_m, cont.Q_measured_m3s, s=18, marker="s", color=INK, zorder=4); xx = np.linspace(cont.L_m.min() - 0.05, cont.L_m.max() + 0.05, 30); ax.plot(xx, np.exp(np.polyval(br, np.log(np.clip(xx - primary[0], 0.02, None)))), color=INK, lw=0.8, ls="--")
    for _, r_ in cont.iterrows(): ax.annotate(r_.date[5:], (r_.L_m, r_.Q_measured_m3s), fontsize=5.5, xytext=(3, 2), textcoords="offset points")
    ax.set_xlabel("sections' level L on the survey day (m)"); ax.set_ylabel("measured discharge (m³/s)"); ax.set_title(f"Continuity: {out['continuity']['rating']}, {resid.std(ddof=2)*100:.0f}% scatter\nsame-day sections max/min {', '.join(f'{v:.2f}' for v in same.values)}", fontsize=8); ax.grid(True); ax.set_axisbelow(True)
    fig.suptitle(f"Step 9b: the {reach} reach through the chain with its sensor", fontsize=10, y=0.95); fig.savefig(F / f"15_{reach}_sensor_chain.png", bbox_inches="tight"); plt.close(fig); log(f"   saved 15_{reach}")
    return out
res = {r: run(r) for r in (sys.argv[1:] or ["Lower", "Upper"])}; print(json.dumps({r: {k: v for k, v in o.items() if k in ("law", "transect_out", "section_day_out", "continuity", "season")} for r, o in res.items()}, indent=1, default=str))
