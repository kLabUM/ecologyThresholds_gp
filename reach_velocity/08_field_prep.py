"""Step 0e. Prepare the hand-held readings as model inputs, and lay everything out for quality control.
1. Log profile: v(z) = a * ln(z / z0), one roughness scale z0 for the reach (grid search on all profiles with three or more
   readings), one amplitude a per location by least squares. Two targets per location: the velocity at 1 ft above the bed
   (v_1ft) and the depth-averaged velocity (v_mean = a * (ln(h/z0) - 1)); the spread of the per-reading implied values
   gives the uncertainty. Readings below 0.03 m or at or below 0.005 m/s are not used in the fit.
2. Reach coordinates from the centreline: distance along (m, zero at the ADCP, positive downstream) and offset across
   (m, positive toward the left bank looking downstream).
3. Checks: duplicated positions, points far from the centreline, inverted profiles, locations shallower than 1 ft,
   transfer error where a reading near 1 ft exists, coordinate consistency between State Plane and lat/lon.
Outputs: data_clean/locations.csv (one row per location: all model inputs and targets), readings_transect_fitted.csv,
field_qc.json, figures/00l_field_inputs.png, 00m_logprofile.png"""
import json
from pathlib import Path
import numpy as np, pandas as pd, geopandas as gpd
from shapely.geometry import Point
from shapely.ops import linemerge
from pyproj import Transformer
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; HERE = Path(__file__).resolve().parent; D = HERE / "data_clean"; F = HERE / "figures"; RAW = HERE / "data_raw"
FT = 0.3048; Z_REF = 0.3048; ADCP_E, ADCP_N = 13244857.0, 319995.5; ZMIN, VMIN = 0.03, 0.005
tr = pd.read_csv(D / "readings_transect.csv", parse_dates=["time_utc"]); qc = {}
# ---------------------------------------------------------------- 1. log profile
use = tr[(tr.height_m >= ZMIN) & (tr.velocity_ms > VMIN)].copy(); nper = use.groupby("location_id").size()
fitset = use[use.location_id.isin(nper[nper >= 3].index)]
def sse_for(z0):
    tot = 0.0
    for _, p in fitset.groupby("location_id"):
        w = np.log(p.height_m.values / z0); v = p.velocity_ms.values; a = (v * w).sum() / (w * w).sum(); tot += ((v - a * w) ** 2).sum()
    return tot
z0_grid = np.logspace(-4, -1, 61); sse = np.array([sse_for(z) for z in z0_grid]); z0 = float(z0_grid[sse.argmin()])
flat = z0_grid[sse <= sse.min() * 1.05]; qc["z0_m"] = z0; qc["z0_within_5pct_of_min_m"] = [float(flat.min()), float(flat.max())]; qc["profiles_in_z0_fit"] = int(fitset.location_id.nunique())
rows = []; fitted = []
for loc, p in tr.groupby("location_id"):
    ok = p[(p.height_m >= ZMIN) & (p.velocity_ms > VMIN)]; h = float(p.water_depth_m.iloc[0]); r = {"location_id": loc, "n_readings": len(p), "n_used": len(ok), "water_depth_m": h, "h_min_m": float(p.height_m.min()), "h_max_m": float(p.height_m.max())}
    if len(ok) == 0:   # all readings at or below the floor: a still location
        r.update({"a": 0.0, "v_1ft": 0.0, "v_1ft_sd": np.nan, "v_mean": 0.0, "v_mean_sd": np.nan, "fit_r2": np.nan, "method": "all readings at or below 0.005 m/s"})
    else:
        w = np.log(ok.height_m.values / z0); v = ok.velocity_ms.values; a = float((v * w).sum() / (w * w).sum()); implied = v * np.log(Z_REF / z0) / w
        r2 = float(1 - ((v - a * w) ** 2).sum() / ((v - v.mean()) ** 2).sum()) if len(ok) >= 3 and v.std() > 0 else np.nan
        r.update({"a": a, "v_1ft": a * np.log(Z_REF / z0), "v_1ft_sd": float(implied.std(ddof=1)) if len(ok) >= 2 else np.nan, "v_mean": a * (np.log(h / z0) - 1) if h > z0 * np.e else 0.0, "fit_r2": r2, "method": "log profile, pooled z0"})
        r["v_mean_sd"] = r["v_1ft_sd"] * r["v_mean"] / r["v_1ft"] if r["v_1ft"] > 0 and not np.isnan(r["v_1ft_sd"]) else np.nan
        for _, q in ok.assign(implied_v_1ft=implied).iterrows(): fitted.append(q)
    near = p[(p.height_m - Z_REF).abs() <= 0.03]; r["v_measured_near_1ft"] = float(near.velocity_ms.mean()) if len(near) else np.nan
    if len(near) and len(ok) > len(near):   # converted from the OTHER readings only, for the transfer check
        oth = ok[(ok.height_m - Z_REF).abs() > 0.03]; w2 = np.log(oth.height_m.values / z0); a2 = float((oth.velocity_ms.values * w2).sum() / (w2 * w2).sum()) if len(oth) else np.nan; r["v_1ft_from_others"] = a2 * np.log(Z_REF / z0) if len(oth) else np.nan
    else: r["v_1ft_from_others"] = np.nan
    first = p.iloc[0]; r.update({"year": int(first.year), "date": first.date, "time_utc": p.time_utc.min(), "transect": first["Transect ID"], "point": first["Measurement Point ID"], "Easting": first.Easting, "Northing": first.Northing, "Latitude": first.Latitude, "Longitude": first.Longitude, "x_m": first.x_m, "y_m": first.y_m, "stage_m": float(p.stage_m.median()), "stage_raw_m": float(p.stage_raw_m.median()), "notes": first["Point Notes"]})
    rows.append(r)
loc = pd.DataFrame(rows); pd.DataFrame(fitted).to_csv(D / "readings_transect_fitted.csv", index=False)
tc = loc.dropna(subset=["v_1ft_from_others", "v_measured_near_1ft"]); qc["transfer_check"] = {"n_profiles": int(len(tc)), "MAE_m_s": float((tc.v_1ft_from_others - tc.v_measured_near_1ft).abs().mean()), "bias_m_s": float((tc.v_1ft_from_others - tc.v_measured_near_1ft).mean()), "r": float(np.corrcoef(tc.v_1ft_from_others, tc.v_measured_near_1ft)[0, 1])}
loc["v_1ft_sd"] = loc.v_1ft_sd.fillna(qc["transfer_check"]["MAE_m_s"])   # single-reading locations get the pooled transfer error
loc["shallower_than_1ft"] = loc.water_depth_m < Z_REF; loc["inverted_profile"] = loc.a < 0; loc["extrapolated_up"] = (loc.h_max_m < 0.15) & (loc.n_used > 0)
# ---------------------------------------------------------------- 2. reach coordinates
cl = gpd.read_file(RAW / "centerline_union_Middle.geojson"); cl = cl.to_crs("ESRI:102690") if (cl.crs is None or cl.crs.is_geographic) else cl
u = cl.geometry.union_all() if hasattr(cl.geometry, "union_all") else cl.geometry.unary_union
line = u if u.geom_type == "LineString" else linemerge(u)
if line.geom_type != "LineString": line = max(line.geoms, key=lambda g: g.length); qc["centreline"] = "MultiLineString: longest part used"
def proj(e, n):
    pt = Point(e, n); d = line.project(pt); q = line.interpolate(d); d1, d2 = max(d - 5, 0), min(d + 5, line.length); t = np.array(line.interpolate(d2).coords[0]) - np.array(line.interpolate(d1).coords[0]); t /= np.linalg.norm(t)
    off = np.array([e - q.x, n - q.y]); return d, float(t[0] * off[1] - t[1] * off[0]), float(np.linalg.norm(off))
P = np.array([proj(e, n) for e, n in zip(loc.Easting, loc.Northing)]); d_adcp, _, _ = proj(ADCP_E, ADCP_N)
along = (P[:, 0] - d_adcp) * FT; downstream_sign = -1.0 if np.corrcoef(along, loc.y_m)[0, 1] > 0 else 1.0   # downstream = decreasing Northing (the gauge at Ann Arbor lies south-east)
loc["along_m"] = along * downstream_sign; loc["across_m"] = P[:, 1] * FT * downstream_sign; loc["dist_to_centreline_m"] = P[:, 2] * FT
qc["centreline_orientation"] = "native direction is downstream" if downstream_sign > 0 else "native direction is upstream; along and across flipped"
# ---------------------------------------------------------------- 3. checks
tf = Transformer.from_crs("EPSG:4269", "ESRI:102690", always_xy=True); ex, ny = tf.transform(loc.Longitude.values, loc.Latitude.values); loc["latlon_vs_stateplane_m"] = np.hypot(ex - loc.Easting, ny - loc.Northing) * FT
dup = loc.groupby(["Easting", "Northing"]).size(); qc["duplicate_positions"] = int((dup > 1).sum())
qc["locations"] = int(len(loc)); qc["shallower_than_1ft"] = int(loc.shallower_than_1ft.sum()); qc["inverted_profiles"] = int(loc.inverted_profile.sum()); qc["single_reading_locations"] = int((loc.n_used == 1).sum()); qc["still_locations"] = int((loc.n_used == 0).sum())
qc["far_from_centreline_over_40m"] = int((loc.dist_to_centreline_m > 40).sum()); qc["latlon_vs_stateplane_over_20m"] = int((loc.latlon_vs_stateplane_m > 20).sum()); qc["latlon_missing"] = int(loc.Latitude.isna().sum())
qc["fit_r2_quartiles"] = [float(x) for x in loc.fit_r2.quantile([.25, .5, .75])]; qc["v_1ft_range_m_s"] = [float(loc.v_1ft.min()), float(loc.v_1ft.max())]; qc["v_1ft_sd_median_m_s"] = float(loc.v_1ft_sd.median())
cols = ["location_id", "year", "date", "time_utc", "transect", "point", "along_m", "across_m", "dist_to_centreline_m", "x_m", "y_m", "Easting", "Northing", "Latitude", "Longitude", "latlon_vs_stateplane_m", "stage_m", "stage_raw_m", "water_depth_m", "n_readings", "n_used", "h_min_m", "h_max_m", "a", "fit_r2", "v_1ft", "v_1ft_sd", "v_mean", "v_mean_sd", "v_measured_near_1ft", "v_1ft_from_others", "shallower_than_1ft", "inverted_profile", "extrapolated_up", "method", "notes"]
loc = loc[cols].sort_values(["date", "transect", "along_m", "across_m"]).reset_index(drop=True); loc.to_csv(D / "locations.csv", index=False); json.dump(qc, open(D / "field_qc.json", "w"), indent=2, default=str); print(json.dumps(qc, indent=2, default=str))
# ---------------------------------------------------------------- figures
INK, INK2, AXIS, GRID, BLUE, ORANGE, C25, C26 = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#2a78d6", "#1a9e74"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig = plt.figure(figsize=(14, 12)); gs = fig.add_gridspec(5, 6, hspace=0.75, wspace=0.35)
ax = fig.add_subplot(gs[:, :2]); bk = gpd.read_file(RAW / "bank_lines_Middle.geojson"); bk = bk.to_crs("ESRI:102690") if (bk.crs is None or bk.crs.is_geographic) else bk
for geom in bk.geometry:
    for g_ in (geom.geoms if hasattr(geom, "geoms") else [geom]): xs, ys = np.array(g_.coords).T; ax.plot((xs - ADCP_E) * FT, (ys - ADCP_N) * FT, color="#b9b8b0", lw=0.8)
xs, ys = np.array(line.coords).T; ax.plot((xs - ADCP_E) * FT, (ys - ADCP_N) * FT, color="#6f6e68", lw=0.9, ls="--")
sc = ax.scatter(loc.x_m, loc.y_m, c=loc.v_1ft, cmap="Blues", vmin=0, vmax=0.6, s=16, edgecolor=INK, lw=0.2, zorder=4); ax.scatter(loc.x_m[loc.shallower_than_1ft], loc.y_m[loc.shallower_than_1ft], facecolor="none", edgecolor=ORANGE, s=40, lw=0.8, zorder=5, label="shallower than 1 ft")
ax.scatter([0], [0], marker="*", s=150, color=ORANGE, edgecolor="white", lw=0.6, zorder=6); fig.colorbar(sc, ax=ax, shrink=0.4, pad=0.02).set_label("velocity at 1 ft above the bed (m/s)")
pad = 60; ax.set_xlim(loc.x_m.min() - pad, loc.x_m.max() + pad); ax.set_ylim(loc.y_m.min() - pad, loc.y_m.max() + pad); ax.set_aspect("equal"); ax.set_xlabel("east of the ADCP (m)"); ax.set_ylabel("north of the ADCP (m)"); ax.set_title("Every location, coloured by the 1 ft velocity", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7, loc="lower right")
trs = loc.groupby("transect").agg(date=("date", "first"), stage=("stage_m", "median"), n=("v_1ft", "size"), along=("along_m", "median")).sort_values(["date", "along"]).reset_index()
for i, t in trs.iterrows():
    if i >= 20: break
    ax = fig.add_subplot(gs[i // 4, 2 + i % 4]); p = loc[loc.transect == t.transect].sort_values("across_m"); col = C25 if p.year.iloc[0] == 2025 else C26
    ax.errorbar(p.across_m, p.v_1ft, yerr=1.96 * p.v_1ft_sd, fmt="o-", color=col, ms=3, lw=0.8, elinewidth=0.6, capsize=1.5); ax.scatter(p.across_m[p.shallower_than_1ft], p.v_1ft[p.shallower_than_1ft], facecolor="none", edgecolor=ORANGE, s=40, lw=0.8, zorder=5)
    ax.axhline(0.2, color=INK2, lw=0.5, ls=":"); ax.set_title(f"{t.transect}  {t.date}\nstage {t.stage:.2f} m, {int(t.along):+d} m along, n={int(t.n)}", fontsize=6.5); ax.set_ylim(-0.05, 0.9); ax.grid(True); ax.set_axisbelow(True); ax.tick_params(labelsize=6)
    if i % 4 == 0: ax.set_ylabel("v at 1 ft (m/s)", fontsize=7)
    if i // 4 == 4: ax.set_xlabel("across (m), + = left bank", fontsize=7)
fig.suptitle("The field inputs: one velocity per location at 1 ft above the bed (bars: 95% from the spread of the readings; circles: water shallower than 1 ft; dotted: 0.2 m/s)", fontsize=9, y=0.93)
fig.savefig(F / "00l_field_inputs.png", bbox_inches="tight"); plt.close(fig)
fig, axes = plt.subplots(1, 4, figsize=(14, 3.4), gridspec_kw={"wspace": 0.35})
ax = axes[0]; ax.plot(z0_grid * 1000, sse, color=INK, lw=1.2); ax.axvline(z0 * 1000, color=BLUE, lw=1); ax.set_xscale("log"); ax.set_xlabel("z0 (mm)"); ax.set_ylabel("sum of squared residuals"); ax.set_title(f"Roughness scale: z0 = {z0*1000:.1f} mm (flat between {flat.min()*1000:.1f} and {flat.max()*1000:.1f})", fontsize=8); ax.grid(True); ax.set_axisbelow(True)
ax = axes[1]; ff = pd.read_csv(D / "readings_transect_fitted.csv")
for lid, p in ff.groupby("location_id"):
    if p.location_id.iloc[0] in loc.location_id.values[:0]: pass
ex = loc[loc.n_used >= 6].sample(12, random_state=3).location_id
for lid in ex:
    p = ff[ff.location_id == lid].sort_values("height_m"); a = loc.set_index("location_id").a[lid]; z = np.linspace(z0 * 1.01, max(p.height_m.max(), Z_REF) * 1.05, 50); ax.plot(a * np.log(z / z0), z, color=BLUE, lw=0.6, alpha=0.6); ax.scatter(p.velocity_ms, p.height_m, s=6, color=INK, alpha=0.6, lw=0)
ax.axhline(Z_REF, color=INK, lw=0.6, ls=(0, (4, 3))); ax.set_xlabel("velocity (m/s)"); ax.set_ylabel("height above bed (m)"); ax.set_title("Twelve profiles and their fitted log curves", fontsize=8); ax.grid(True); ax.set_axisbelow(True)
ax = axes[2]; ax.plot([0, 0.7], [0, 0.7], color=AXIS, lw=0.8); ax.scatter(tc.v_measured_near_1ft, tc.v_1ft_from_others, s=12, color=BLUE, edgecolor="white", lw=0.4)
ax.set_xlabel("measured within 3 cm of 1 ft (m/s)"); ax.set_ylabel("moved to 1 ft from the other readings (m/s)"); ax.set_title(f"Transfer check: {len(tc)} profiles, MAE {qc['transfer_check']['MAE_m_s']*1000:.0f} mm/s, bias {qc['transfer_check']['bias_m_s']*1000:+.0f}", fontsize=8); ax.grid(True); ax.set_axisbelow(True)
ax = axes[3]; ax.hist(loc.fit_r2.dropna(), bins=np.linspace(-0.5, 1, 31), color=BLUE, edgecolor="white", lw=0.4); ax.set_xlabel("R2 of the log fit (locations with 3+ readings)"); ax.set_ylabel("locations"); ax.set_title(f"Fit quality: median R2 {loc.fit_r2.median():.2f}; {int(loc.inverted_profile.sum())} inverted profiles", fontsize=8); ax.grid(True); ax.set_axisbelow(True)
fig.savefig(F / "00m_logprofile.png", bbox_inches="tight"); plt.close(fig); print("saved 00l, 00m")
print(loc[["location_id", "along_m", "across_m", "water_depth_m", "n_used", "v_1ft", "v_1ft_sd", "v_mean", "fit_r2", "stage_m"]].describe().round(3).to_string())
