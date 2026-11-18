"""Step 8a. The other two reaches (Upper, Lower), prepared the way the Middle was, but without a stage sensor.
Readings -> one amplitude a per location (log profile, one z0 per reach). RTK -> one water-surface elevation per transect (the median of
bed + depth over its points; the surface across a transect is flat, so the scatter grades the session's GNSS fix), bed = that surface
minus the rod depth, the reach slope from the good transects, a local stage on each good day and its relation to the nearest gauge
(Upper: 04172000 near Hamburg, 6 km downstream; Lower: 04174500 at Ann Arbor, 42 km upstream, three dams between). Transects whose
session had no usable fix get their surface from the gauge rating, extrapolated along the slope, and are marked untied when that
extrapolation is too long to trust. Single-reading 'BT' transects (one reading at 0.4 of the depth, elevations on another datum) are set aside.
Outputs: data_clean/reach_{r}_points.csv, reach_{r}_transects.csv, reach_{r}_days.csv, reach_{r}_qc.json, figures/08a_other_reaches_inputs.png"""
import json, sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd, geopandas as gpd
from shapely.geometry import Point
from shapely.ops import linemerge
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore"); t0 = time.time(); log = lambda *a: (print(f"[{time.time()-t0:5.0f}s]", *a), sys.stdout.flush())
ROOT = Path(__file__).resolve().parents[1]; HERE = Path(__file__).resolve().parent; D = HERE / "data_clean"; F = HERE / "figures"; RAW = HERE / "data_raw"; MAP = RAW
FT = 0.3048; ZMIN, VMIN = 0.03, 0.005; SUBS = ["cobble", "gravel", "sand", "silt", "vegetation", "woody debris"]
GOOD_SCATTER_M = 0.15; TIE_MAX_M = 0.20   # a transect's surface is 'good' when its points agree to this; a rating-derived surface is 'tied' when its uncertainty is below this
# Upper: the southern ~1 km only (Northing below 365,000 ft); the August transects further upstream were exploratory and are left out
REACHES = {"Upper": {"band": (340000, 365000), "gauge": "04172000", "gauge_name": "Huron River near Hamburg (04172000), 6 km downstream"},
           "Lower": {"band": (-1e9, 255000), "gauge": "04174500", "gauge_name": "Huron River at Ann Arbor (04174500), 42 km upstream, three dams between"}}
raw = pd.read_csv(RAW / "26_ecohydroData.csv", dtype=str); raw.columns = [c.strip() for c in raw.columns]
for c in ["Total water depth", "Depth", "Velocity", "Northing", "Easting", "Elevation"]: raw[c] = pd.to_numeric(raw[c], errors="coerce")
for c in ["Transect ID", "Measurement Point ID", "Type", "Point Notes"]: raw[c] = raw[c].astype(str).str.strip().replace({"nan": np.nan})
raw["date"] = pd.to_datetime(raw["Collection Date"].str.strip(), format="mixed").dt.strftime("%Y-%m-%d")
gauges = {"04172000": pd.read_csv(D / "usgs_04172000_daily.csv").set_index("date").discharge_cfs, "04174500": pd.read_csv(D / "usgs_daily.csv").set_index("date").discharge_cfs}
def geometry(reach):
    cl = gpd.read_file(MAP / f"centerline_union_{reach}.geojson"); u = cl.geometry.union_all() if hasattr(cl.geometry, "union_all") else cl.geometry.unary_union; line = u if u.geom_type == "LineString" else linemerge(u)
    bk = gpd.read_file(MAP / f"bank_lines_{reach}.geojson"); return line, bk[bk.bank == "left"].geometry.iloc[0], bk[bk.bank == "right"].geometry.iloc[0]
def robust_sd(x): x = np.asarray(x, float); x = x[~np.isnan(x)]; return float(1.4826 * np.median(np.abs(x - np.median(x)))) if len(x) >= 3 else np.nan
out_all = {}
for reach, cfg in REACHES.items():
    qc = {"reach": reach, "gauge": cfg["gauge_name"]}
    r = raw[(raw.Northing > cfg["band"][0]) & (raw.Northing < cfg["band"][1])].copy(); qc["rows_by_type"] = r.Type.value_counts(dropna=False).to_dict()
    tr = r[r.Type == "Transect"].copy(); isbt = tr["Transect ID"].str.contains("BT", na=False); qc["BT_transects_set_aside"] = {"transects": sorted(tr[isbt]["Transect ID"].unique().tolist()), "rows": int(isbt.sum()), "why": "one reading per point at 0.4 of the depth, elevations 7 m above the reach's surface (another datum or no fix)"}; tr = tr[~isbt]
    tr["height_m"] = tr.Depth * FT; tr["water_depth_m"] = tr["Total water depth"] * FT; tr["velocity_ms"] = tr.Velocity * FT; tr["elev_rtk_m"] = tr.Elevation * FT
    n0 = len(tr); tr = tr[tr.velocity_ms.notna() & tr.height_m.notna() & tr.water_depth_m.notna()]; n1 = len(tr); tr = tr[tr.velocity_ms >= 0]; n2 = len(tr); tr = tr[~(tr.height_m > tr.water_depth_m + 0.02)]; n3 = len(tr)
    qc["transect_filter"] = {"rows": n0, "with velocity, height and depth": n1, "velocity >= 0": n2, "height not above the water surface": n3}
    tr["location_id"] = tr.date + "|" + tr["Transect ID"].astype(str) + "|" + tr["Measurement Point ID"].astype(str); tr["transect"] = tr["Transect ID"]
    # ---- profiles: one z0 for the reach from profiles with three or more usable readings, then a per location
    use = tr[(tr.height_m >= ZMIN) & (tr.velocity_ms > VMIN)]; nper = use.groupby("location_id").size(); fitset = use[use.location_id.isin(nper[nper >= 3].index)]
    def sse_for(z0):
        tot = 0.0
        for _, p in fitset.groupby("location_id"):
            w = np.log(p.height_m.values / z0); v = p.velocity_ms.values; a = (v * w).sum() / (w * w).sum(); tot += ((v - a * w) ** 2).sum()
        return tot
    z0_grid = np.logspace(-4, -1, 61); sse = np.array([sse_for(z) for z in z0_grid]); z0 = float(z0_grid[sse.argmin()]); flat = z0_grid[sse <= sse.min() * 1.05]
    qc["z0_m"] = z0; qc["z0_within_5pct_of_min_m"] = [float(flat.min()), float(flat.max())]; qc["profiles_in_z0_fit"] = int(fitset.location_id.nunique())
    rows = []
    for lid, p in tr.groupby("location_id"):
        ok = p[(p.height_m >= ZMIN) & (p.velocity_ms > VMIN)]; h = float(p.water_depth_m.iloc[0]); rec = {"location_id": lid, "n_readings": len(p), "n_used": len(ok), "water_depth_m": h}
        if len(ok) >= 1:
            w = np.log(ok.height_m.values / z0); v = ok.velocity_ms.values; a = float((v * w).sum() / (w * w).sum()); implied = v / w
            rec.update({"a": a, "a_sd": float(implied.std(ddof=1)) if len(ok) >= 2 else np.nan, "fit_r2": float(1 - ((v - a * w) ** 2).sum() / max(((v - v.mean()) ** 2).sum(), 1e-12)) if len(ok) >= 2 else np.nan, "v_mean": a * (np.log(h / z0) - 1) if h > z0 * np.e else 0.0})
        else: rec.update({"a": 0.0, "a_sd": np.nan, "fit_r2": np.nan, "v_mean": 0.0})   # every reading at or below the velocity floor: standing water
        for c in ["transect", "date", "Easting", "Northing", "elev_rtk_m"]: rec[c] = p[c].iloc[0]
        rec["notes"] = " | ".join(p["Point Notes"].dropna().unique()); rows.append(rec)
    loc = pd.DataFrame(rows); loc["a_sd"] = loc.a_sd.fillna(loc.a_sd.median()); loc["year"] = 2026
    qc["locations"] = {"n": int(len(loc)), "a_zero_standing": int((loc.a <= 0).sum()), "a_median": float(loc.a[loc.a > 0].median()), "a_sd_median": float(loc.a_sd.median()), "fit_r2_median": float(loc.fit_r2.median())}
    # ---- geometry: along the centreline from the most upstream location, signed across, banks
    line, left, right = geometry(reach); pts = [Point(e, n) for e, n in zip(loc.Easting, loc.Northing)]; d_along = np.array([line.project(p) for p in pts]); origin = d_along.min()
    def across(p, d):
        q = line.interpolate(d); d1, d2 = max(d - 5, 0), min(d + 5, line.length); t = np.array(line.interpolate(d2).coords[0]) - np.array(line.interpolate(d1).coords[0]); t /= np.linalg.norm(t); off = np.array([p.x - q.x, p.y - q.y]); return float(t[0] * off[1] - t[1] * off[0]) * FT
    loc["along_m"] = (d_along - origin) * FT; loc["across_m"] = [across(p, d) for p, d in zip(pts, d_along)]; o = line.interpolate(origin); loc["x_m"] = (loc.Easting - o.x) * FT; loc["y_m"] = (loc.Northing - o.y) * FT
    loc["dist_left_bank_m"] = [left.distance(p) * FT for p in pts]; loc["dist_right_bank_m"] = [right.distance(p) * FT for p in pts]; loc["dist_bank_m"] = loc[["dist_left_bank_m", "dist_right_bank_m"]].min(axis=1)
    loc["width_m"] = loc.dist_left_bank_m + loc.dist_right_bank_m; loc["rel_across"] = (loc.dist_right_bank_m - loc.dist_left_bank_m) / loc.width_m
    qc["geometry"] = {"reach_length_m": float(loc.along_m.max()), "max_offset_from_centreline_m": float(loc.across_m.abs().max()), "corr_rel_across_vs_across": float(np.corrcoef(loc.rel_across, loc.across_m)[0, 1]), "width_m_median": float(loc.width_m.median())}
    # ---- the water surface, one per transect: median of (RTK bed + rod depth); its scatter grades the fix
    loc["wse_pt_m"] = loc.elev_rtk_m + loc.water_depth_m
    T = loc.groupby("transect").agg(date=("date", "first"), n=("a", "size"), along_m=("along_m", "median"), wse_rtk_m=("wse_pt_m", "median"), wse_scatter_m=("wse_pt_m", robust_sd), depth_mean_m=("water_depth_m", "mean"), width_m=("width_m", "median")).reset_index()
    T["rtk_good"] = (T.wse_scatter_m < GOOD_SCATTER_M) & (T.n >= 3); good = T[T.rtk_good]
    # slope from the good transects: wse = c_day - S along, days with two or more good transects
    gd = good.groupby("date").filter(lambda g: len(g) >= 2); days_s = sorted(gd.date.unique()); X = np.column_stack([(gd.date.values[:, None] == np.array(days_s)[None, :]).astype(float), -gd.along_m.values]); b, res_ss, *_ = np.linalg.lstsq(X, gd.wse_rtk_m.values, rcond=None)
    S = float(b[-1]); dof = max(len(gd) - X.shape[1], 1); resid = gd.wse_rtk_m.values - X @ b; s2 = (resid ** 2).sum() / dof; S_se = float(np.sqrt(s2 * np.linalg.inv(X.T @ X)[-1, -1]))
    along_ref = float(good.along_m.median()); qc["slope"] = {"cm_per_100m": S * 1e4, "se_cm_per_100m": S_se * 1e4, "from_good_transects": int(len(gd)), "days": days_s, "resid_sd_cm": float(np.sqrt(s2) * 100), "along_ref_m": along_ref}
    # local stage per good day at the reference position, and its relation to the gauge
    good = good.assign(stage_ref_m=good.wse_rtk_m + S * (good.along_m - along_ref)); dd = good.groupby("date").agg(stage_local_m=("stage_ref_m", "mean"), stage_spread_m=("stage_ref_m", lambda s: s.max() - s.min()), n_good_transects=("transect", "size"), along_median_m=("along_m", "median")).reset_index()
    Q = gauges[cfg["gauge"]]; dd["Q_cfs"] = dd.date.map(Q); ok = dd.Q_cfs.notna(); br = np.polyfit(np.log(dd.Q_cfs[ok]), dd.stage_local_m[ok], 1) if ok.sum() >= 3 else np.array([np.nan, np.nan])
    dd["stage_from_rating_m"] = np.polyval(br, np.log(dd.Q_cfs)); rr = dd.stage_local_m - dd.stage_from_rating_m; rating_sd = float(rr[ok].std(ddof=2)) if ok.sum() >= 3 else np.nan
    qc["rating"] = {"stage_vs_lnQ_slope_m": float(br[0]), "intercept_m": float(br[1]), "residual_sd_cm": rating_sd * 100, "n_days": int(ok.sum()), "Q_range_cfs": [float(dd.Q_cfs.min()), float(dd.Q_cfs.max())], "stage_range_m": [float(dd.stage_local_m.min()), float(dd.stage_local_m.max())], "corr_stage_lnQ": float(np.corrcoef(np.log(dd.Q_cfs[ok]), dd.stage_local_m[ok])[0, 1]) if ok.sum() >= 3 else np.nan}
    # every transect gets a surface: RTK when good, else from the rating at the day's flow carried along the slope; the tie uncertainty says whether that can be trusted
    T["Q_cfs"] = T.date.map(Q); T["wse_rating_m"] = np.polyval(br, np.log(T.Q_cfs)) - S * (T.along_m - along_ref); T["tie_sd_m"] = np.sqrt(rating_sd ** 2 + (S_se * (T.along_m - along_ref)) ** 2)
    T["wse_m"] = np.where(T.rtk_good, T.wse_rtk_m, T.wse_rating_m); T["wse_source"] = np.where(T.rtk_good, "rtk", "rating"); T["tied"] = T.rtk_good | (T.tie_sd_m < TIE_MAX_M); T["wse_sd_m"] = np.where(T.rtk_good, 0.03, T.tie_sd_m)
    T["rtk_minus_rating_m"] = T.wse_rtk_m - T.wse_rating_m
    qc["transects"] = {"n": int(len(T)), "rtk_good": int(T.rtk_good.sum()), "tied_total": int(T.tied.sum()), "untied": sorted(T[~T.tied].transect.tolist()), "days_all": sorted(T.date.unique().tolist()), "days_good": sorted(good.date.unique().tolist())}
    loc = loc.merge(T[["transect", "wse_m", "wse_source", "wse_sd_m", "tied", "rtk_good", "wse_scatter_m"]], on="transect"); loc["bed_m"] = loc.wse_m - loc.water_depth_m; loc["bed_rtk_resid_m"] = loc.elev_rtk_m - loc.bed_m
    loc["stage_local_m"] = loc.wse_m   # the surface at the point on its day (NAVD88 m)
    # ---- features within the section from the survey-day depths (no datum needed), notes
    g = loc.groupby("transect"); loc["rel_depth"] = loc.water_depth_m / g.water_depth_m.transform("max").replace(0, np.nan); loc["rel_depth"] = loc.rel_depth.fillna(0); loc["section_mean_depth_m"] = g.water_depth_m.transform("mean"); loc["section_width_m"] = g.width_m.transform("median")
    loc = loc.sort_values(["transect", "across_m"]).reset_index(drop=True); g = loc.groupby("transect")
    loc["lateral_bed_slope"] = g.apply(lambda d: pd.Series(np.gradient(d.bed_m.values, d.across_m.values) if (len(d) > 2 and d.across_m.nunique() == len(d)) else np.zeros(len(d)), index=d.index), include_groups=False).reset_index(level=0, drop=True); loc["lateral_bed_slope_abs"] = loc.lateral_bed_slope.abs().fillna(0)
    notes = loc.notes.fillna("").str.lower(); sub = pd.Series("not noted", index=loc.index)
    for word, lab in [("cobble", "cobble"), ("gravel", "gravel"), ("pebble", "gravel"), ("sand", "sand"), ("silt", "silt"), ("clay", "silt"), ("muck", "silt"), ("veg", "vegetation"), ("grass", "vegetation"), ("wood", "woody debris"), ("debris", "woody debris")]: sub[notes.str.contains(word) & (sub == "not noted")] = lab
    loc["substrate"] = sub; loc["wake_flag"] = notes.str.contains("debris|wood|eddy|behind|slow|log|branch").astype(int)
    for s_ in SUBS: loc[f"sub_{s_}"] = (loc.substrate == s_).astype(float)
    qc["substrate_counts"] = loc.substrate.value_counts().to_dict(); qc["wake_flagged"] = int(loc.wake_flag.sum())
    loc["Q_cfs"] = loc.date.map(Q); loc["Q_m3s"] = loc.Q_cfs * 0.0283168
    loc.to_csv(D / f"reach_{reach}_points.csv", index=False); T.to_csv(D / f"reach_{reach}_transects.csv", index=False); dd.to_csv(D / f"reach_{reach}_days.csv", index=False); json.dump(qc, open(D / f"reach_{reach}_qc.json", "w"), indent=2, default=str)
    out_all[reach] = dict(loc=loc, T=T, dd=dd, qc=qc, z0_grid=z0_grid, sse=sse, z0=z0, line=line, left=left, right=right, S=S, br=br, along_ref=along_ref)
    log(f"{reach}: {len(loc)} locations, {len(T)} transects ({int(T.rtk_good.sum())} with a good fix, {int(T.tied.sum())} tied), {T.date.nunique()} days, reach {loc.along_m.max():.0f} m; z0 {z0*1000:.1f} mm ({flat.min()*1000:.1f} to {flat.max()*1000:.1f}); a median {qc['locations']['a_median']:.3f}")
    log(f"   slope {S*1e4:.1f} ± {S_se*1e4:.1f} cm/100 m from {len(gd)} good transects on {len(days_s)} days (resid {np.sqrt(s2)*100:.1f} cm); local stage vs ln Q on {int(ok.sum())} good days: r {qc['rating']['corr_stage_lnQ']:.2f}, slope {br[0]:.3f} m per ln(cfs), residual {rating_sd*100:.1f} cm, Q {dd.Q_cfs.min():.0f} to {dd.Q_cfs.max():.0f} cfs; untied: {qc['transects']['untied']}")
# ---- figure
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, DARK, MUTED = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#0d366b", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig = plt.figure(figsize=(19, 10)); gs = fig.add_gridspec(2, 5, wspace=0.38, hspace=0.45, width_ratios=[1.5, 1, 1.2, 1, 1])
for i, reach in enumerate(REACHES):
    o = out_all[reach]; loc, T, dd = o["loc"], o["T"], o["dd"]; cmap = plt.get_cmap("tab10"); days = sorted(loc.date.unique()); dcol = {d: cmap(k % 10) for k, d in enumerate(days)}
    ax = fig.add_subplot(gs[i, 0]); i0 = loc.along_m.idxmin(); ox, oy = o["line"].interpolate(o["line"].project(Point(loc.Easting[i0], loc.Northing[i0]))).coords[0]
    for gml in (o["left"], o["right"]): xs, ys = np.array(gml.coords).T; ax.plot((xs - ox) * FT, (ys - oy) * FT, color="#8a877f", lw=0.6)
    xs, ys = np.array(o["line"].coords).T; ax.plot((xs - ox) * FT, (ys - oy) * FT, color=AXIS, lw=0.5, ls=":")
    sc = ax.scatter(loc.x_m, loc.y_m, c=loc.a, cmap="Blues", s=9, vmin=0, vmax=np.nanquantile(loc.a, .98), edgecolor=INK, lw=0.15); fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02).set_label("amplitude a (m/s)")
    for _, t in T[~T.tied].iterrows(): p = loc[loc.transect == t.transect]; ax.scatter(p.x_m.mean(), p.y_m.mean(), s=90, facecolor="none", edgecolor="red", lw=0.7)
    pad = 40; ax.set_xlim(loc.x_m.min() - pad, loc.x_m.max() + pad); ax.set_ylim(loc.y_m.min() - pad, loc.y_m.max() + pad); ax.set_aspect("equal"); ax.set_title(f"{reach}: {len(loc)} locations, {len(T)} transects, {len(days)} days, {loc.along_m.max():.0f} m of channel\n(red rings: transects with no vertical tie to the others)", fontsize=8.5); ax.grid(True); ax.set_axisbelow(True); ax.tick_params(labelsize=6); ax.set_xlabel("east (m)"); ax.set_ylabel("north (m)")
    ax = fig.add_subplot(gs[i, 1]); ax.plot(o["z0_grid"] * 1000, o["sse"], color=INK, lw=1.1); ax.axvline(o["z0"] * 1000, color=BLUE, lw=1); ax.axvline(11, color=ORANGE, lw=0.8, ls="--"); ax.set_xscale("log"); ax.set_xlabel("z0 (mm)"); ax.set_ylabel("sum of squared residuals"); ax.set_title(f"Roughness scale: z0 = {o['z0']*1000:.1f} mm (orange: the Middle's 11 mm)\nprofiles with 3+ readings: {o['qc']['profiles_in_z0_fit']}, median fit R² {o['qc']['locations']['fit_r2_median']:.2f}", fontsize=8); ax.grid(True); ax.set_axisbelow(True)
    ax = fig.add_subplot(gs[i, 2])
    for _, t in T.iterrows():
        c = dcol[t.date]; ax.errorbar(t.along_m, t.wse_rtk_m, yerr=t.wse_scatter_m, fmt="o" if t.rtk_good else "o", mfc=c if t.rtk_good else "white", mec=c, ecolor=c, ms=4, elinewidth=0.6, capsize=0)
        if not t.rtk_good: ax.scatter(t.along_m, t.wse_rating_m, marker="_", s=60, color=c, lw=1.2)
    xx = np.linspace(0, loc.along_m.max(), 50)
    for d in dd.date: ax.plot(xx, dd.set_index("date").stage_local_m[d] - o["S"] * (xx - o["along_ref"]), color=dcol[d], lw=0.5, alpha=0.7)
    lo, hi = np.nanpercentile(T.wse_rating_m, [0, 100]); ax.set_ylim(lo - 1.5, hi + 1.5); ax.set_xlabel("along the reach (m, downstream)"); ax.set_ylabel("water surface (m NAVD88)")
    ax.set_title(f"One water surface per transect from the RTK (dot: median of bed + depth; bar: the points' scatter)\nfilled: good fix; open: no usable fix, dash = surface from the gauge rating; lines: good days along the slope {o['S']*1e4:.1f} cm/100 m", fontsize=7.5); ax.grid(True); ax.set_axisbelow(True)
    ax = fig.add_subplot(gs[i, 3]); rt = o["qc"]["rating"]; ax.scatter(dd.Q_cfs * 0.0283168, dd.stage_local_m, s=26, color=[dcol[d] for d in dd.date], edgecolor=INK, lw=0.3, zorder=3)
    qq = np.linspace(min(dd.Q_cfs.min(), T.Q_cfs.min()) * 0.8, max(dd.Q_cfs.max(), T.Q_cfs.max()) * 1.2, 50); ax.plot(qq * 0.0283168, np.polyval(o["br"], np.log(qq)), color=INK, lw=0.8); ax.fill_between(qq * 0.0283168, np.polyval(o["br"], np.log(qq)) - 2 * rt["residual_sd_cm"] / 100, np.polyval(o["br"], np.log(qq)) + 2 * rt["residual_sd_cm"] / 100, color=MUTED, alpha=0.2, lw=0)
    for _, rr_ in dd.iterrows(): ax.annotate(rr_.date[5:], (rr_.Q_cfs * 0.0283168, rr_.stage_local_m), fontsize=5.5, xytext=(3, 3), textcoords="offset points", color=INK2)
    bad_days = sorted(set(T[~T.rtk_good].date) - set(dd.date))
    for d in bad_days: ax.axvline(T[T.date == d].Q_cfs.iloc[0] * 0.0283168, color=dcol[d], lw=0.6, ls=":")
    ax.set_xscale("log"); ax.set_xlabel(f"gauge discharge on the day (m³/s), {REACHES[reach]['gauge']}"); ax.set_ylabel(f"local stage at {o['along_ref']:.0f} m along (m NAVD88)"); ax.set_title(f"Local stage from the good days against the gauge: r {rt['corr_stage_lnQ']:.2f}, residual {rt['residual_sd_cm']:.1f} cm, n {rt['n_days']}\n(dotted: flows on the days with no usable fix; the rating has to reach them)", fontsize=7.5); ax.grid(True, which="both"); ax.set_axisbelow(True)
    ax = fig.add_subplot(gs[i, 4])
    for d in days: p = loc[loc.date == d]; ax.errorbar(p.water_depth_m, p.a, yerr=1.96 * p.a_sd, fmt="o", color=dcol[d], ms=2.5, elinewidth=0.4, alpha=0.75, label=f"{d[5:]} ({p.Q_cfs.iloc[0]*0.0283168:.1f} m³/s)")
    ax.set_xlabel("water depth at the location on its day (m)"); ax.set_ylabel("amplitude a (m/s)"); ax.set_title("a against depth by day (bars: 95% from the readings' spread)", fontsize=8); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=5.5, ncol=2)
fig.suptitle("Step 8a: the Upper and Lower reaches prepared the same way, with the RTK water surface and the nearest gauge in place of the stage sensor", fontsize=10, y=0.96); fig.savefig(F / "08a_other_reaches_inputs.png", bbox_inches="tight"); plt.close(fig); log("saved 08a")
