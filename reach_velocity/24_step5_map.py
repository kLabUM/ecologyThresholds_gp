"""Step 5. From the sections to the reach: the bed between transects (GP on the RTK cross-sections and the bank points),
the channel grid, the conveyance field queried on the grid, velocity and the share below a threshold at chosen stages,
and where the bands are wide. Numbers still come from the sections; this is the map with its uncertainty.
Outputs: grid.csv, step5_reach_totals.csv, figures/05_step5_map.png"""
import json, sys, time, warnings
import numpy as np, pandas as pd, geopandas as gpd
from shapely.geometry import Point
from shapely.ops import linemerge
from scipy import linalg
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel as C, Matern, WhiteKernel
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from gpcore import *; from uncertainty import *
warnings.filterwarnings("ignore"); t0 = time.time(); log = lambda *a: (print(f"[{time.time()-t0:5.0f}s]", *a), sys.stdout.flush())
RAW = HERE / "data_raw"; FT = 0.3048; ADCP_E, ADCP_N = 13244857.0, 319995.5; STAGE_REF = 0.90; THR = 0.20; STAGES = [0.75, 0.90, 1.10]
df = load_points(); df["bed_m"] = df.bed_rtk_sensor_datum_m.fillna(df.bed_elev_m)
# ---- geometry: centreline, banks, channel polygon
cl = gpd.read_file(RAW / "centerline_union_Middle.geojson"); cl = cl.to_crs("ESRI:102690") if (cl.crs is None or cl.crs.is_geographic) else cl
u = cl.geometry.union_all() if hasattr(cl.geometry, "union_all") else cl.geometry.unary_union; line = u if u.geom_type == "LineString" else linemerge(u)
bk = gpd.read_file(RAW / "bank_lines_Middle.geojson"); bk = bk.to_crs("ESRI:102690") if (bk.crs is None or bk.crs.is_geographic) else bk; left, right = bk[bk.bank == "left"].geometry.iloc[0], bk[bk.bank == "right"].geometry.iloc[0]
ext = gpd.read_file((RAW / "channel_extent_Middle.geojson") if (RAW / "channel_extent_Middle.geojson").exists() else (ROOT / "mapping" / "channel_extent_Middle.geojson")); ext = ext.to_crs("ESRI:102690") if (ext.crs is None or ext.crs.is_geographic) else ext; poly = ext.geometry.iloc[0]
d_adcp = line.project(Point(ADCP_E, ADCP_N))
def proj(e, n):
    pt = Point(e, n); d = line.project(pt); q = line.interpolate(d); d1, d2 = max(d - 5, 0), min(d + 5, line.length); t = np.array(line.interpolate(d2).coords[0]) - np.array(line.interpolate(d1).coords[0]); t /= np.linalg.norm(t); off = np.array([e - q.x, n - q.y])
    return (d - d_adcp) * FT, float(t[0] * off[1] - t[1] * off[0]) * FT
# ---- bed GP on the surveyed points plus the bank points
banks = pd.read_csv(D / "banks.csv").dropna(subset=["bank_sensor_datum_m"]); bp = np.array([proj(e, n) for e, n in zip(banks.Easting.astype(float), banks.Northing.astype(float))])
bed_pts = pd.DataFrame({"along": np.r_[df.along_m.values, bp[:, 0]], "across": np.r_[df.across_m.values, bp[:, 1]], "z": np.r_[df.bed_m.values, banks.bank_sensor_datum_m.values], "sd": np.r_[np.full(len(df), 0.06), np.full(len(banks), 0.10)]})
Xb = bed_pts[["along", "across"]].values; mxb, sxb = Xb.mean(0), Xb.std(0); zb_m, zb_s = bed_pts.z.mean(), bed_pts.z.std()
kb = C(1.0, (1e-2, 1e2)) * Matern([1.0, 1.0], [(0.02, 20.0), (0.02, 20.0)], nu=1.5) + WhiteKernel(0.05, (1e-4, 1.0))
gpb = GaussianProcessRegressor(kb, alpha=(bed_pts.sd.values / zb_s) ** 2, n_restarts_optimizer=3, random_state=0).fit((Xb - mxb) / sxb, (bed_pts.z.values - zb_m) / zb_s)
lsb = gpb.kernel_.k1.k2.length_scale * sxb; log(f"bed GP: along scale {lsb[0]:.0f} m, across scale {lsb[1]:.1f} m, nugget sd {np.sqrt(gpb.kernel_.k2.noise_level)*zb_s*100:.0f} cm, on {len(bed_pts)} points")
# ---- the channel grid (2 m raster inside the channel polygon)
minx, miny, maxx, maxy = poly.bounds; step = 2.0 / FT; gx, gy = np.meshgrid(np.arange(minx, maxx, step), np.arange(miny, maxy, step)); gx, gy = gx.ravel(), gy.ravel()
inside = np.array([poly.contains(Point(x, y)) for x, y in zip(gx, gy)]); gx, gy = gx[inside], gy[inside]; P = np.array([proj(x, y) for x, y in zip(gx, gy)])
grid = pd.DataFrame({"x_m": (gx - ADCP_E) * FT, "y_m": (gy - ADCP_N) * FT, "along_m": P[:, 0], "across_m": P[:, 1]}); grid = grid[(grid.along_m > df.along_m.min() - 15) & (grid.along_m < df.along_m.max() + 15)].reset_index(drop=True)
pts = [Point(x / FT + ADCP_E, y / FT + ADCP_N) for x, y in zip(grid.x_m, grid.y_m)]; grid["dist_left"] = [left.distance(p) * FT for p in pts]; grid["dist_right"] = [right.distance(p) * FT for p in pts]
grid["dist_bank_m"] = grid[["dist_left", "dist_right"]].min(axis=1); grid["rel_across"] = (grid.dist_right - grid.dist_left) / (grid.dist_left + grid.dist_right)
mb, sb = gpb.predict((grid[["along_m", "across_m"]].values - mxb) / sxb, return_std=True); grid["bed_m"] = mb * zb_s + zb_m; grid["bed_sd_m"] = sb * zb_s
grid["h_ref_m"] = np.clip(STAGE_REF - grid.bed_m, 0, None); abin = (grid.along_m // 10) * 10; grid["rel_depth"] = grid.h_ref_m / grid.groupby(abin).h_ref_m.transform("max").replace(0, np.nan); grid["rel_depth"] = grid.rel_depth.fillna(0)
sec = df.groupby("transect").agg(along=("along_m", "median"), smd=("section_mean_depth_m", "first")).sort_values("along"); grid["section_mean_depth_m"] = np.interp(grid.along_m, sec.along, sec.smd)
grid["lateral_bed_slope_abs"] = float(df.lateral_bed_slope_abs.median())
for c in NOTE_COLS: grid[c] = float(df[c].mean())   # the grid carries the population share of each note, not the 'nothing noted' class, so the map's mean matches the surveyed points
log(f"grid: {len(grid)} cells of 2 m inside the channel polygon, along {grid.along_m.min():.0f} to {grid.along_m.max():.0f} m")
# ---- the conveyance field on the grid
ch = Chain(df); m = ch.models[ch.primary]; muA, sdA = m.predict(grid); grid["logA"], grid["logA_sd"] = muA, sdA
# nearest transect distance for each cell (how far from data)
grid["dist_to_transect_m"] = [float(np.min(np.abs(sec.along.values - a))) for a in grid.along_m]
grid.to_csv(D / "grid.csv", index=False); log(f"field on the grid: median sd of log A {np.median(sdA):.2f}; on transects {np.median(sdA[grid.dist_to_transect_m < 3]):.2f}, between (> 15 m) {np.median(sdA[grid.dist_to_transect_m > 15]):.2f}")
# ---- reach totals from joint draws on a thinned grid (every other row/col ~ 4 m), and from the sections weighted by length
thin = grid.iloc[::3].reset_index(drop=True); totals = []
for st in STAGES + [1.30]:
    smp = ch.sample(st, thin, n=200, seed=5); wet = smp["h"] > 0.02; area_share = np.array([np.mean(smp["vbar"][k][wet[k]] < THR) if wet[k].any() else np.nan for k in range(200)]); wet_area = wet.mean(axis=1) * len(grid) * 4.0
    totals.append({"stage_m": st, "source": "map (joint draws, 4 m grid)", "share_below_thr_median": float(np.nanmedian(area_share)), "lo": float(np.nanquantile(area_share, .05)), "hi": float(np.nanquantile(area_share, .95)), "wet_area_m2": float(np.median(wet_area))})
    # sections: the joint reach sampler (one GP draw over every section, common stage, z0 and day-level discrepancy; width-weighted shares, length-weighted reach total)
    secs = {t_: df[df.transect == t_].assign(_x=df.x_m.round(2), _y=df.y_m.round(2)).drop_duplicates(["_x", "_y"]).sort_values("across_m").reset_index(drop=True) for t_ in sec.index}
    rs = sample_reach(ch, [st], secs, n=200, seed=100)["reach_share"][:, 0]
    totals.append({"stage_m": st, "source": "sections, length-weighted", "share_below_thr_median": float(np.nanmedian(rs)), "lo": float(np.nanquantile(rs, .05)), "hi": float(np.nanquantile(rs, .95)), "wet_area_m2": np.nan}); log(f"stage {st:.2f}: share below {THR} m/s, map {totals[-2]['share_below_thr_median']:.2f} [{totals[-2]['lo']:.2f}, {totals[-2]['hi']:.2f}], sections {totals[-1]['share_below_thr_median']:.2f} [{totals[-1]['lo']:.2f}, {totals[-1]['hi']:.2f}]")
tot = pd.DataFrame(totals); tot.to_csv(D / "step5_reach_totals.csv", index=False)
# ---- figure
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, DARK, MUTED = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#0d366b", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig = plt.figure(figsize=(16, 11)); gs = fig.add_gridspec(1, 5, wspace=0.35)
def bankplot(ax):
    for g in (left, right): xs, ys = np.array(g.coords).T; ax.plot((xs - ADCP_E) * FT, (ys - ADCP_N) * FT, color="#8a877f", lw=0.7)
panels = [("bed_m", "viridis_r", "bed elevation (m, sensor datum)", None, None), ("logA", "Blues", "conveyance A", None, None), ("logA_sd", "magma_r", "sd of log A (how uncertain the field is)", 0.2, 0.9)]
for j, (col, cmap, lab, vmin, vmax) in enumerate(panels):
    ax = fig.add_subplot(gs[0, j]); vals = np.exp(grid[col]) if col == "logA" else grid[col]; sc = ax.scatter(grid.x_m, grid.y_m, c=vals, cmap=cmap, s=2.5, lw=0, vmin=vmin, vmax=vmax); bankplot(ax); ax.scatter(df.x_m, df.y_m, s=1.5, color=INK, alpha=0.5, lw=0)
    fig.colorbar(sc, ax=ax, shrink=0.5, pad=0.02).set_label(lab); ax.set_aspect("equal"); ax.set_title(["The bed, from the surveyed cross-sections and bank points", "The conveyance field on the grid", "Where the field is a guess: the band between transects"][j], fontsize=8.5); ax.tick_params(labelsize=6); ax.set_xlabel("east of the ADCP (m)")
    if j == 0: ax.set_ylabel("north of the ADCP (m)")
smp = ch.sample(STAGE_REF, grid, n=1, sources=(), seed=0); vb = smp["vbar"][0]; wet = smp["h"][0] > 0.02
ax = fig.add_subplot(gs[0, 3]); sc = ax.scatter(grid.x_m[wet], grid.y_m[wet], c=vb[wet], cmap="Blues", s=2.5, lw=0, vmin=0, vmax=0.6); ax.scatter(grid.x_m[~wet], grid.y_m[~wet], s=2.5, color="#e6e2d8", lw=0); bankplot(ax); fig.colorbar(sc, ax=ax, shrink=0.5, pad=0.02).set_label(f"depth-averaged velocity at stage {STAGE_REF} m (m/s)"); ax.set_aspect("equal"); ax.set_title(f"Velocity at {STAGE_REF} m of stage (dry cells in beige)", fontsize=8.5); ax.tick_params(labelsize=6); ax.set_xlabel("east of the ADCP (m)")
ax = fig.add_subplot(gs[0, 4]); below = np.where(wet, vb < THR, np.nan); ax.scatter(grid.x_m[wet], grid.y_m[wet], c=np.where(vb[wet] < THR, 1, 0), cmap="coolwarm_r", s=2.5, lw=0, vmin=0, vmax=1); ax.scatter(grid.x_m[~wet], grid.y_m[~wet], s=2.5, color="#e6e2d8", lw=0); bankplot(ax); ax.set_aspect("equal"); ax.set_title(f"Below {THR} m/s at {STAGE_REF} m (blue) or not (red)", fontsize=8.5); ax.tick_params(labelsize=6); ax.set_xlabel("east of the ADCP (m)")
fig.suptitle("Step 5: from the sections to the reach. The map is the interpolation; the numbers come from the sections", fontsize=10, y=0.93); fig.savefig(F / "05_step5_map.png", bbox_inches="tight"); plt.close(fig); log("saved 05")
print(tot.round(3).to_string(index=False))
