"""Step 9a. The Upper and Lower reaches with their own stage sensors (MaxBotix ARB013 on the Upper, ARB112 on the Lower; Lower files added to
data/ on 15 and 18 September 2026, the Upper record replaced on 18 November 2026). Sensor cleaning the Middle's way (implausible readings out, spikes against a 6 h rolling median out,
10-min grid, 3 h median; the daily wobble is measured and reported, not removed), stage at each reading's time, the RTK water surface tied
to the sensor per session (one slope), the bed on the sensor datum, features at a reference stage, and the day-out grid over (h0, gamma).
Bank features from bank set C at each point's survey-day level (banklines_3dep.StageFeatures, entry 30); its lines go to data_clean/banks_c/.
Outputs: data_clean/reach_{r}_stage_10min.csv, reach_{r}_points_sensor.csv, reach_{r}_tie.csv, reach_{r}_law_grid.csv, reach_{r}_sensor_qc.json, figures/14_other_reaches_sensor_prep.png"""
import json, sys, time, warnings
import numpy as np, pandas as pd
from pathlib import Path
RAW = Path(__file__).resolve().parent / "data_raw"
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import matplotlib.dates as mdates
from gpcore import *; from uncertainty import *
warnings.filterwarnings("ignore"); t0 = time.time(); log = lambda *a: (print(f"[{time.time()-t0:5.0f}s]", *a), sys.stdout.flush())
FT = 0.3048; TZ = "America/Detroit"
CFG = {"Upper": {"files": [str(RAW / "26_depthSensorData_Upper.csv")], "offset_mm": 3048.0, "site": "ARB013", "lat_lon": (42.51365577, -83.68990820), "note": "sensor on the main stem upstream of the surveyed stretch (record replaced 18 Nov 2026: ARB013, offset tag 3048 mm)"},
       "Lower": {"files": [str(RAW / "25_depthSensorData_Lower.csv"), str(RAW / "26_depthSensorData_Lower.csv")], "offset_mm": 6705.6, "site": "ARB112", "lat_lon": None, "along_sensor_m": -258.0, "note": "sensor 258 m above the first section, below Flat Rock dam; the offset tag 5943.6 before 8 June 2026 is a tag error, the range is continuous"}}
def sensor_along_m(reach, lat_lon, pts):
    """The sensor's position along the reach centreline (m, same origin as the points' along_m), from its coordinates."""
    import geopandas as gpd
    from shapely.geometry import Point
    from shapely.ops import linemerge
    from pyproj import Transformer
    cl = gpd.read_file(RAW.parent / "data_clean" / "banklines_survey" / f"centerline_union_{reach}.geojson"); u = cl.geometry.union_all() if hasattr(cl.geometry, "union_all") else cl.geometry.unary_union; line = u if u.geom_type == "LineString" else linemerge(u)
    x, y = Transformer.from_crs(4326, cl.crs, always_xy=True).transform(lat_lon[1], lat_lon[0]); p0 = pts.iloc[0]
    origin = line.project(Point(p0.Easting, p0.Northing)) - p0.along_m / FT   # the points' along_m origin, in the line's units (ft)
    return float((line.project(Point(x, y)) - origin) * FT), float(line.distance(Point(x, y)) * FT)
def clean_stage(files, offset_mm):
    d = pd.concat([pd.read_csv(f, skiprows=3) for f in files]); d = d[d["_field"] == "value"].copy(); d["time_utc"] = pd.to_datetime(d["_time"], utc=True); d["stage_raw_m"] = (offset_mm - d["_value"].astype(float)) / 1000
    d = d.sort_values("time_utc").drop_duplicates("time_utc").set_index("time_utc"); n0 = len(d); ok = d.stage_raw_m.between(0.05, 3.0); d = d[ok]; n1 = len(d)
    med = d.stage_raw_m.rolling("6h", center=True).median(); spike = (d.stage_raw_m - med).abs() > 0.15; d = d[~spike]; n2 = len(d)
    g = d.stage_raw_m.resample("10min").mean().to_frame("stage_10_m"); g["stage_m"] = g.stage_10_m.rolling("3h", center=True).median(); g["stage_m"] = g.stage_m.interpolate(limit=36, limit_area="inside")   # gaps up to 6 h bridged
    hr = g.stage_m.dropna(); hour = hr.index.hour; wob = hr.groupby(hour).mean() - hr.mean()
    return g, {"rows": n0, "implausible_removed": int(n0 - n1), "spikes_removed": int(n1 - n2), "coverage_of_10min_grid": float(g.stage_m.notna().mean()), "daily_wobble_peak_to_trough_cm": float((wob.max() - wob.min()) * 100), "period": [str(g.index.min().date()), str(g.index.max().date())]}
def stage_at(sm, t):
    if pd.isna(t): return np.nan
    i = sm.index.searchsorted(t); lo, hi = max(i - 1, 0), min(i, len(sm) - 1); gap = min(abs((sm.index[lo] - t).total_seconds()), abs((sm.index[hi] - t).total_seconds())) / 3600
    if gap > 2: return np.nan
    w = 0.5 if sm.index[hi] == sm.index[lo] else (t - sm.index[lo]) / (sm.index[hi] - sm.index[lo]); return float(sm.iloc[lo] * (1 - w) + sm.iloc[hi] * w)
raw = pd.read_csv(str(RAW / "26_ecohydroData.csv"), dtype=str); raw.columns = [c.strip() for c in raw.columns]; raw["date"] = pd.to_datetime(raw["Collection Date"].str.strip(), format="mixed").dt.strftime("%Y-%m-%d")
raw["time_local"] = pd.to_datetime(raw["Collection Date"].str.strip() + " " + raw["Collection Time"].fillna("12:00").str.strip(), format="mixed", errors="coerce"); raw["time_utc"] = raw.time_local.dt.tz_localize(TZ, ambiguous="NaT", nonexistent="NaT").dt.tz_convert("UTC")
raw["location_id"] = raw.date + "|" + raw["Transect ID"].astype(str).str.strip() + "|" + raw["Measurement Point ID"].astype(str).str.strip(); tloc = raw.groupby("location_id").time_utc.first()
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, DARK, MUTED = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#0d366b", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig = plt.figure(figsize=(18, 10)); gs = fig.add_gridspec(2, 4, wspace=0.3, hspace=0.45, width_ratios=[2, 1, 1, 1])
for i, (reach, cfg) in enumerate((r, c) for r, c in CFG.items() if r in (sys.argv[1:] or CFG)):   # python 37_other_reaches_prep.py [Upper] [Lower]
    g, qc = clean_stage(cfg["files"], cfg["offset_mm"]); g.to_csv(D / f"reach_{reach}_stage_10min.csv"); sm = g.stage_m.dropna(); qc.update({"site": cfg["site"], "offset_mm": cfg["offset_mm"], "note": cfg["note"], "along_sensor_m": cfg.get("along_sensor_m")})
    pts = pd.read_csv(D / f"reach_{reach}_points.csv")
    if cfg.get("lat_lon"):
        cfg["along_sensor_m"], off_line = sensor_along_m(reach, cfg["lat_lon"], pts); qc.update({"along_sensor_m": cfg["along_sensor_m"], "sensor_offset_from_centreline_m": off_line, "sensor_lat_lon": list(cfg["lat_lon"])})
    pts["time_utc"] = pts.location_id.map(tloc); pts["stage_m"] = [stage_at(sm, t) for t in pts.time_utc]; pts["stage_day_median_m"] = pts.date.map(sm.groupby(sm.index.date).median().rename(index=str))
    pts["stage_m"] = pts.stage_m.fillna(pts.stage_day_median_m); qc["readings_without_stage_at_time"] = int(pts.stage_m.isna().sum())
    # ---- the tie: RTK water surface against the sensor, one offset per session (day) and one slope, on the transects with a good fix
    T = pd.read_csv(D / f"reach_{reach}_transects.csv"); T["stage_m"] = [float(pts[pts.transect == t].stage_m.median()) for t in T.transect]; good = T[T.rtk_good & T.stage_m.notna()].copy(); good["diff"] = good.wse_rtk_m - good.stage_m; good["dx"] = good.along_m - cfg["along_sensor_m"]
    days = sorted(good.date.unique()); X = np.column_stack([(good.date.values[:, None] == np.array(days)[None, :]).astype(float), -good.dx.values]); b, *_ = np.linalg.lstsq(X, good["diff"].values, rcond=None); S = float(b[-1]); res = good["diff"].values - X @ b
    off = dict(zip(days, b[:-1])); qc["tie"] = {"slope_cm_per_100m": S * 1e4, "residual_sd_cm": float(res.std(ddof=1) * 100), "n_transects": int(len(good)), "session_offsets_m": {d: float(v) for d, v in off.items()}, "offset_spread_sd_m": float(np.std(list(off.values()), ddof=1)), "meaning": "WSE_rtk - stage_sensor = offset(day) - slope * (along - along_sensor); a constant offset across days means the sensor level transfers to the sections; its spread is the tie's error"}
    T["tie_resid_m"] = T.transect.map(dict(zip(good.transect, res))); T["session_offset_m"] = T.date.map(off)
    # the bed on the sensor datum: the RTK bed minus the day's offset where the fix is good, else stage at the visit minus the rod depth
    pts["session_offset_m"] = pts.date.map(off); pts["bed_sensor_m"] = np.where(pts.rtk_good & pts.session_offset_m.notna(), pts.elev_rtk_m - pts.session_offset_m + S * (pts.along_m - cfg["along_sensor_m"]), pts.stage_m - pts.water_depth_m)
    # a check on the same thing the other way: stage at the visit minus depth, against the RTK bed on the sensor datum
    chk = pts[pts.rtk_good & pts.session_offset_m.notna()]; qc["bed_check_rtk_minus_rod_cm"] = {"median": float(np.median((chk.elev_rtk_m - chk.session_offset_m + S * (chk.along_m - cfg["along_sensor_m"])) - (chk.stage_m - chk.water_depth_m)) * 100), "sd": float(np.std((chk.elev_rtk_m - chk.session_offset_m + S * (chk.along_m - cfg["along_sensor_m"])) - (chk.stage_m - chk.water_depth_m)) * 100)}
    pts["bed_m"] = pts.bed_sensor_m; ref = float(pts.stage_m.median()); qc["reference_stage_m"] = ref
    pts["h_ref_m"] = np.clip(ref - S * (pts.along_m - cfg["along_sensor_m"]) - pts.bed_m, 0, None); gsec = pts.groupby("transect"); pts["rel_depth"] = (pts.h_ref_m / gsec.h_ref_m.transform("max").replace(0, np.nan)).fillna(0); pts["section_mean_depth_m"] = gsec.h_ref_m.transform("mean")
    pts["slope_cm_per_100m"] = S * 1e4; pts["along_sensor_m"] = cfg["along_sensor_m"]; pts["local_stage_m"] = pts.stage_m - S * (pts.along_m - cfg["along_sensor_m"])   # the water level at the section on the sensor datum
    pts.to_csv(D / f"reach_{reach}_points_sensor.csv", index=False); T.to_csv(D / f"reach_{reach}_tie.csv", index=False)
    # ---- bank set C (entry 30): needs the tie just written; the bank features at each point's survey-day level, and the lines on disk
    import banklines_3dep as B
    SF = B.StageFeatures(reach); k_ = B.tie(reach); L_rec = k_["a"] + k_["b"] * sm.values; SF.write(D / "banks_c", float(L_rec.max()))
    pts = SF.training(pts); pts.to_csv(D / f"reach_{reach}_points_sensor.csv", index=False); qc["bank_set_C"] = {"record_max_level_m": float(L_rec.max()), "dist_bank_m_median": float(pts.dist_bank_m.median())}
    # ---- the day-out grid over (h0, gamma) with the sensor stage (as step 3b for the Middle): h0 below the lowest stage the river reached
    d = pts[pts.a > 0].copy(); smin = float(sm.min()); rows = []
    for h0 in np.round(np.arange(max(smin - 0.45, -0.3), smin + 0.01, 0.1), 2):
        for gam in [0.5, 0.75, 1.0, 1.25]:
            dd = with_target(d, ("stage", gam, h0)); Y, MU, SD = [], [], []
            for day in sorted(dd.date.unique()):
                tr, te = dd[dd.date != day], dd[dd.date == day]
                if len(tr) < 20: continue
                m = Model(n_restarts=1).fit(tr); mu, sd = m.predict(te); Y.append(te.y.values); MU.append(mu); SD.append(np.sqrt(sd ** 2 + te.sd_y.values ** 2))
            sc = scores(np.concatenate(Y), np.concatenate(MU), np.concatenate(SD)); rows.append({"h0": float(h0), "gamma": gam, **{k: v for k, v in sc.items()}})
    lg = pd.DataFrame(rows).sort_values("CRPS"); lg.to_csv(D / f"reach_{reach}_law_grid.csv", index=False)
    dn = with_target(d, ("none",)); Y, MU, SD = [], [], []
    for day in sorted(dn.date.unique()):
        tr, te = dn[dn.date != day], dn[dn.date == day]
        if len(tr) < 20: continue
        m = Model(n_restarts=1).fit(tr); mu, sd = m.predict(te); Y.append(te.y.values); MU.append(mu); SD.append(np.sqrt(sd ** 2 + te.sd_y.values ** 2))
    crps_none = scores(np.concatenate(Y), np.concatenate(MU), np.concatenate(SD))["CRPS"]; best = lg.iloc[0]; qc["law"] = {"best": [float(best.h0), float(best.gamma)], "CRPS_best": float(best.CRPS), "CRPS_none": float(crps_none), "set_within_0.01": [[float(r.h0), float(r.gamma)] for _, r in lg[lg.CRPS <= best.CRPS + 0.01].iterrows()], "stage_range_survey_days_m": [float(d.stage_m.min()), float(d.stage_m.max())], "stage_range_record_m": [float(sm.min()), float(sm.max())]}
    json.dump(qc, open(D / f"reach_{reach}_sensor_qc.json", "w"), indent=1, default=str)
    log(f"{reach} ({cfg['site']}): cleaned {qc['rows']} rows ({qc['implausible_removed']} implausible, {qc['spikes_removed']} spikes), daily wobble {qc['daily_wobble_peak_to_trough_cm']:.1f} cm; tie on {len(good)} transects: slope {S*1e4:.1f} cm/100 m, residual {res.std(ddof=1)*100:.1f} cm, session offsets sd {qc['tie']['offset_spread_sd_m']*100:.0f} cm; bed check RTK - rod {qc['bed_check_rtk_minus_rod_cm']['median']:.0f} ± {qc['bed_check_rtk_minus_rod_cm']['sd']:.0f} cm; law: best {qc['law']['best']} CRPS {best.CRPS:.3f} vs none {crps_none:.3f}; set {qc['law']['set_within_0.01']}; survey stages {d.stage_m.min():.2f} to {d.stage_m.max():.2f} m, record {sm.min():.2f} to {sm.max():.2f} m")
    log("   session offsets (m): " + ", ".join(f"{k}: {v:.2f}" for k, v in off.items()))
    # ---- figure row
    ax = fig.add_subplot(gs[i, 0]); ax.plot(sm.index, sm.values, color=INK, lw=0.5)
    for dday in sorted(pts.date.unique()): ax.axvline(pd.Timestamp(dday, tz="UTC"), color=BLUE, lw=0.7)
    ax.set_ylabel("stage (m, sensor datum)"); ax.set_title(f"{reach}: the cleaned sensor record ({cfg['site']}; {qc['implausible_removed']} implausible and {qc['spikes_removed']} spikes removed; daily wobble {qc['daily_wobble_peak_to_trough_cm']:.1f} cm); blue: survey days", fontsize=8.5); ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %y")); ax.grid(True); ax.set_axisbelow(True)
    ax = fig.add_subplot(gs[i, 1]); cmap = plt.get_cmap("tab10"); dcol = {dd_: cmap(k % 10) for k, dd_ in enumerate(sorted(T.date.unique()))}
    for _, t in good.iterrows(): ax.scatter(t.dx, t["diff"], s=22, color=dcol[t.date], edgecolor=INK, lw=0.3); ax.annotate(t.date[5:], (t.dx, t["diff"]), fontsize=5.5, xytext=(3, 2), textcoords="offset points")
    xx = np.linspace(good.dx.min() - 50, good.dx.max() + 50, 10)
    for dday, o in off.items(): ax.plot(xx, o - S * xx, color=dcol[dday], lw=0.6, alpha=0.7)
    ax.set_xlabel("along the reach from the sensor (m)"); ax.set_ylabel("RTK water surface − sensor stage (m)"); ax.set_title(f"The tie: one offset per session, one slope ({S*1e4:.1f} cm/100 m)\nresidual {res.std(ddof=1)*100:.1f} cm; offsets spread {qc['tie']['offset_spread_sd_m']*100:.0f} cm across days", fontsize=8); ax.grid(True); ax.set_axisbelow(True)
    ax = fig.add_subplot(gs[i, 2]); dl = d.groupby("date").agg(st=("stage_m", "median"), la=("a", lambda v: np.log(v).median()), lo=("a", lambda v: np.log(v).quantile(.25)), hi=("a", lambda v: np.log(v).quantile(.75)), n=("a", "size"))
    for dday, r in dl.iterrows(): ax.errorbar(r.st, r.la, yerr=[[r.la - r.lo], [r.hi - r.la]], fmt="o", color=dcol.get(dday, INK), ms=4, elinewidth=0.6); ax.annotate(f"{dday[5:]} ({int(r.n)})", (r.st, r.la), fontsize=5.5, xytext=(3, 2), textcoords="offset points")
    xx = np.linspace(d.stage_m.min() - 0.02, d.stage_m.max() + 0.02, 30); ax.plot(xx, np.log(np.clip(xx - best.h0, 0.02, None)) * best.gamma + (dl.la - np.log(np.clip(dl.st - best.h0, 0.02, None)) * best.gamma).mean(), color=INK, lw=0.8, ls="--")
    ax.set_xlabel("sensor stage on the survey day (m)"); ax.set_ylabel("log a, day median (bars: quartiles)"); ax.set_title(f"Amplitude against the sensor stage by day\n(dashed: the best law, h0 {best.h0:.2f} m, gamma {best.gamma:.2f})", fontsize=8); ax.grid(True); ax.set_axisbelow(True)
    ax = fig.add_subplot(gs[i, 3])
    for h0v, grp in lg.groupby("h0"): grp = grp.sort_values("gamma"); ax.plot(grp.gamma, grp.CRPS, marker="o", ms=3, lw=0.8, label=f"h0 = {h0v:.2f} m")
    ax.axhline(crps_none, color=ORANGE, lw=0.8, ls="--"); ax.text(0.52, crps_none + 0.003, "no law", color=ORANGE, fontsize=7); ax.axhspan(best.CRPS, best.CRPS + 0.01, color=GRID, alpha=0.7, lw=0)
    ax.set_xlabel("gamma"); ax.set_ylabel("CRPS of log a, each day held out"); ax.set_title(f"The law grid with the sensor: best {best.CRPS:.3f} against no law {crps_none:.3f}\n(grey: within 0.01 of the best)", fontsize=8); ax.legend(fontsize=6); ax.grid(True); ax.set_axisbelow(True)
fig.suptitle("Step 9a: the Upper and Lower reaches with their own stage sensors: cleaned records, the RTK tie, the day levels, the law grid", fontsize=10, y=0.96); fig.savefig(F / "14_other_reaches_sensor_prep.png", bbox_inches="tight"); plt.close(fig); log("saved 14")
