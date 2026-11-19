"""Candidate input features for the conveyance field, built from what is measured, drawn, and checked against A.
Coordinates: along, across. Hydraulic and geometric: relative position across the channel, distance to the nearest bank,
local width, bed elevation, depth at a reference stage, relative depth within the section, lateral bed slope, section mean
depth and width (pool/riffle), channel curvature and inner/outer bank, substrate class and a wake flag from the notes.
Bank lines: set C (logbook entry 30) at each point's survey-day water level, from banklines_3dep.StageFeatures; centreline from 44_.
Outputs: data_clean/features.csv, data_clean/banks_c/C_Middle_*.geojson, C_outline_Middle.geojson, C_Middle_reference.geojson, figures/00v_features.png"""
import json, re
from pathlib import Path
import numpy as np, pandas as pd, geopandas as gpd
from shapely.geometry import Point
from shapely.ops import linemerge
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"; F = HERE / "figures"; RAW = HERE / "data_raw"
FT = 0.3048; ADCP_E, ADCP_N = 13244857.0, 319995.5; STAGE_REF = 0.90
loc = pd.read_csv(D / "locations.csv"); loc = loc[loc.a > 0].copy(); loc["A"] = loc.a / np.sqrt(loc.water_depth_m)
loc["bed_m"] = loc.bed_rtk_sensor_datum_m.fillna(loc.bed_elev_m)
# ---- banks: set C (logbook entry 30) at each point's own survey-day water level: distance to each bank line, local width, relative position.
# The lines on disk (C at each survey day, at the highest stage of the record, and that outline) are written here.
import banklines_3dep as B
SF = B.StageFeatures("Middle"); SF.write(D / "banks_c", float(pd.read_csv(D / "series_10min.csv").stage_m.max()))
SF.lines(STAGE_REF).to_file(D / "banks_c" / "C_Middle_reference.geojson", driver="GeoJSON")   # C at the reference stage 0.90 m: the lines drawn on the overview and paper maps
loc = SF.training(loc)
for side in ("left", "right"):
    loc[f"dist_{side}_bank_m"] = np.nan   # filled below from the day's lines
for day, ix in loc.groupby("date").indices.items():
    bk = SF.lines(SF.day_stage[day]); P = [Point(e, n) for e, n in zip(loc.Easting.values[ix], loc.Northing.values[ix])]
    for side in ("left", "right"): loc.iloc[ix, loc.columns.get_loc(f"dist_{side}_bank_m")] = [bk[bk.bank == side].geometry.iloc[0].distance(p) * FT for p in P]
sign_check = np.corrcoef(loc.rel_across, loc.across_m)[0, 1]; print(f"relative position vs across_m: r = {sign_check:.2f} (should be strongly positive if 'left' means the same bank in both)")
# ---- depth at a reference stage, relative depth within the section, lateral bed slope, section summaries
loc["h_ref_m"] = np.clip(STAGE_REF - loc.bed_m, 0, None)
g = loc.groupby("transect"); loc["rel_depth"] = loc.h_ref_m / g.h_ref_m.transform("max"); loc["section_mean_depth_m"] = g.h_ref_m.transform("mean"); loc["section_width_m"] = g.width_m.transform("median")
loc = loc.sort_values(["transect", "across_m"]); loc["lateral_bed_slope"] = g.apply(lambda d: pd.Series(np.gradient(d.bed_m.values, d.across_m.values) if len(d) > 2 else np.zeros(len(d)), index=d.index), include_groups=False).reset_index(level=0, drop=True)
loc["lateral_bed_slope_abs"] = loc.lateral_bed_slope.abs()
# ---- curvature of the centreline at each location, and inner/outer bank
cl = gpd.read_file(D / "banklines_survey" / "centerline_union_Middle.geojson"); cl = cl.to_crs("ESRI:102690") if (cl.crs is None or cl.crs.is_geographic) else cl
u = cl.geometry.union_all() if hasattr(cl.geometry, "union_all") else cl.geometry.unary_union; line = u if u.geom_type == "LineString" else linemerge(u)
def curvature(e, n, ds=25 / FT):
    d = line.project(Point(e, n)); p0, p1, p2 = [np.array(line.interpolate(max(min(d + k, line.length), 0)).coords[0]) for k in (-ds, 0, ds)]
    t1, t2 = p1 - p0, p2 - p1; ang = np.arctan2(t1[0] * t2[1] - t1[1] * t2[0], t1 @ t2); return ang / (2 * ds * FT)   # rad per m, positive = turning left
loc["curvature_1_per_m"] = [curvature(e, n) for e, n in zip(loc.Easting, loc.Northing)]
loc["outer_bank_score"] = -np.sign(loc.curvature_1_per_m) * loc.rel_across * (np.abs(loc.curvature_1_per_m) * 1000)   # positive toward the outer bank of a bend, scaled by how sharp the bend is
# ---- notes: substrate class and wake flag
notes = loc.notes.fillna("").str.lower(); sub = pd.Series("not noted", index=loc.index)
for word, lab in [("cobble", "cobble"), ("gravel", "gravel"), ("sand", "sand"), ("silt", "silt"), ("clay", "silt"), ("veg", "vegetation"), ("grass", "vegetation"), ("wood", "woody debris"), ("debris", "woody debris")]: sub[notes.str.contains(word) & (sub == "not noted")] = lab
loc["substrate"] = sub; loc["wake_flag"] = notes.str.contains("debris|wood|eddy|behind|slow|log|branch").astype(int)
feats = ["along_m", "across_m", "rel_across", "dist_bank_m", "width_m", "bed_m", "h_ref_m", "rel_depth", "lateral_bed_slope_abs", "section_mean_depth_m", "section_width_m", "curvature_1_per_m", "outer_bank_score", "wake_flag"]
loc[["location_id", "transect", "date", "year", "A", "a_sd"] + feats + ["substrate"]].to_csv(D / "features.csv", index=False)
rows = [{"feature": f, "spearman_r_with_A": spearmanr(loc[f], loc.A)[0], "spearman_r_with_logA": spearmanr(loc[f], np.log(loc.A))[0]} for f in feats]
sub_med = loc.groupby("substrate").A.agg(["median", "size"]); print(pd.DataFrame(rows).round(2).to_string(index=False)); print("\nA by substrate class:\n", sub_med.round(3).to_string())
json.dump({"features": feats, "spearman": {r["feature"]: round(r["spearman_r_with_A"], 3) for r in rows}, "A_by_substrate": sub_med["median"].round(3).to_dict()}, open(D / "features_summary.json", "w"), indent=2)
# ---- figure
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, MUTED = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig = plt.figure(figsize=(16, 12)); gs = fig.add_gridspec(3, 6, height_ratios=[1.5, 1, 1], hspace=0.5, wspace=0.45)
maps = [("rel_across", "relative position across\n(-1 right bank, +1 left bank)", "RdBu"), ("dist_bank_m", "distance to nearest bank (m)", "viridis"), ("h_ref_m", f"depth at stage {STAGE_REF} m (m)", "viridis"), ("rel_depth", "relative depth in the section", "viridis"), ("outer_bank_score", "outer-bank score\n(bend direction × side × sharpness)", "RdBu"), ("lateral_bed_slope_abs", "lateral bed slope |dz/dy|", "magma")]
for j, (col, lab, cmap) in enumerate(maps):
    ax = fig.add_subplot(gs[0, j]); sc = ax.scatter(loc.x_m, loc.y_m, c=loc[col], cmap=cmap, s=9, lw=0); ax.plot([0], [0], "*", color=ORANGE, ms=9); fig.colorbar(sc, ax=ax, shrink=0.7, pad=0.02); ax.set_aspect("equal"); ax.set_title(lab, fontsize=8); ax.tick_params(labelsize=6); ax.grid(True); ax.set_axisbelow(True)
scat = [("rel_across", "relative position across"), ("dist_bank_m", "distance to nearest bank (m)"), ("h_ref_m", f"depth at stage {STAGE_REF} m (m)"), ("rel_depth", "relative depth in the section"), ("outer_bank_score", "outer-bank score"), ("lateral_bed_slope_abs", "lateral bed slope"), ("section_mean_depth_m", "section mean depth at 0.90 m (pool vs riffle)"), ("width_m", "local width (m)"), ("curvature_1_per_m", "centreline curvature (1/m)"), ("along_m", "along the reach (m)"), ("wake_flag", "wake flag from the notes (0/1)"), ("substrate", "substrate class from the notes")]
for k, (col, lab) in enumerate(scat):
    ax = fig.add_subplot(gs[1 + k // 6, k % 6])
    if col == "substrate":
        order = ["cobble", "gravel", "sand", "silt", "vegetation", "woody debris", "not noted"]; data = [loc.A[loc.substrate == o].values for o in order]; ax.boxplot(data, tick_labels=[f"{o}\n({len(d)})" for o, d in zip(order, data)], widths=0.6, medianprops=dict(color=BLUE)); ax.tick_params(axis="x", labelsize=5.5, rotation=30)
    else:
        jitter = np.random.default_rng(0).uniform(-0.05, 0.05, len(loc)) if col == "wake_flag" else 0; ax.scatter(loc[col] + jitter, loc.A, s=6, color=BLUE, alpha=0.5, lw=0); r = spearmanr(loc[col], loc.A)[0]; ax.set_title(f"{lab}\nSpearman r = {r:+.2f}", fontsize=7.5)
    if col == "substrate": ax.set_title(lab, fontsize=7.5)
    ax.set_ylabel("A", fontsize=7); ax.tick_params(labelsize=6); ax.grid(True); ax.set_axisbelow(True)
fig.suptitle("Candidate input features for the conveyance field, built from measured things, and how the conveyance A responds to each", fontsize=10, y=0.94)
fig.savefig(F / "00v_features.png", bbox_inches="tight"); plt.close(fig); print("saved 00v")
