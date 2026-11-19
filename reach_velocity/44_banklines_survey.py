"""Step 10b. The origin of the old bank lines (set A in logbook entries 27 and 28), made self-contained.

The old lines (bank_lines, channel_extent, centerline_union for each reach) were first copied into data_raw from the
main branch's mapping/ folder; those copies are gone (logbook entry 30), and this script now writes the lines. The code that
wrote them is dataSorting.ipynb on the main branch at commit 2346281 ("Modelling visual tweaks", 18 September 2026), cells 1, 2, 6, 7, 10, 11, 12 and 13.
This script carries that code over unchanged in what it computes (the plots, the CSV and GeoJSON round trips and the
upstream-depth columns, which do not affect the lines, are left out) and checks, by fingerprint, that it reproduces the copies.

Method (the notebook's own description, cell 13): bank points noted in the field ('bank height', 'edge of bank') are
placed on their transects with a straight-line fit of offset against position through the transect's in-water points;
each bank is a curve of offset from the main-stem centreline against distance along it; where no bank point was
surveyed the edge is the outermost survey point plus the margin measured where both exist; the curves are smoothed
over 300 ft, kept 2 ft outside the survey, and pushed out until every survey point is inside the polygon. The in-water
points are the survey visits measured at 1 ft above the bed in water at least 1 ft deep (cell 7, measured-only), plus
the ADCP's 5-hour medians at its own position (cell 6), trimmed to the northing limits of cell 12.

Inputs:  data_raw/25_ecohydroData.csv, 26_ecohydroData.csv, 25_ADCPData.csv, Huron_RiverSystem.geojson
Outputs: data_clean/banklines_survey/bank_lines_{r}.geojson, channel_extent_{r}.geojson, centerline_union_{r}.geojson,
         in_water_points_{r}.csv (type and position only), summary.json (including the check against the September copies)"""
import io, json, tempfile
from pathlib import Path
import numpy as np, pandas as pd, geopandas as gpd
from shapely.geometry import Point, box, LineString, Polygon
from shapely.ops import linemerge

HERE = Path(__file__).resolve().parent; RAW = HERE / "data_raw"; OUT = HERE / "data_clean" / "banklines_survey"; OUT.mkdir(parents=True, exist_ok=True)
CRS = "ESRI:102690"
SECTIONS = ["Upper", "Middle", "Lower"]

# ---------------------------------------------------------------- settings (as in the notebook) ---
NORTHING_UPPER, NORTHING_MIDDLE, NORTHING_LOWER = 340000, (280000, 323000), 255000   # cell 2: sections by northing
ADCP_E, ADCP_N = 13244857, 319995.5                                                     # cell 6
TARGET_HEIGHT_FT, FIT_MIN_HEIGHT_FT = 1.0, 0.5                                           # cell 7
Z_P_BOUNDS_FT, Z_P_STEPS, MIN_LEVELS = (1e-4, 0.30), 60, 2
NORTHING_MAX_FT = {"Upper": 365000, "Middle": 321000, "Lower": 220000}                 # cell 12, applied after the centreline
BANK_NOTES = ["bank height", "edge of bank"]                                            # cell 13
MIN_WET_POINTS, CALIB_LIMIT_FT, OUTSIDE_TOL_FT = 4, 150.0, 1.0
GRID_FT, BIN_FT, WINDOW_FT, BANK_REACH_FT, SMOOTH_FT, SAFETY_FT = 50.0, 100.0, 200.0, 400.0, 300.0, 2.0
MIN_WIDTH_FT, MAX_GAP_FT, MAX_PASSES = 10.0, 800.0, 12
_EPS_FT = 1.0
SEPTEMBER_FINGERPRINTS = {"Upper": {"bank_lines": "8d297be14bb8cfbb", "centerline_union": "903090e1317b9743", "channel_extent": "fe704475ab24da89"},
                          "Middle": {"bank_lines": "6e479317e200e0a2", "centerline_union": "9b1fa5fb8a329ba5", "channel_extent": "7c50b8a6cf2410df"},
                          "Lower": {"bank_lines": "0d13deb6d75532e9", "centerline_union": "2bddaa2f4b27100d", "channel_extent": "c2b999cde8f1c214"}}


def fingerprint(path):
    import hashlib
    h = hashlib.sha256()
    for g in gpd.read_file(path).to_crs(CRS).geometry:
        h.update(np.round(np.asarray(g.exterior.coords if g.geom_type == "Polygon" else g.coords), 1).tobytes())
    return h.hexdigest()[:16]


# ---------------------------------------------------------------- cells 1, 2, 6: the survey rows by section ---
def survey_by_section():
    complete = pd.concat([pd.read_csv(RAW / "25_ecohydroData.csv"), pd.read_csv(RAW / "26_ecohydroData.csv")], ignore_index=True)
    complete.loc[complete.Northing > NORTHING_UPPER, "Section"] = "Upper"
    complete.loc[(complete.Northing <= NORTHING_MIDDLE[1]) & (complete.Northing >= NORTHING_MIDDLE[0]), "Section"] = "Middle"
    complete.loc[complete.Northing < NORTHING_LOWER, "Section"] = "Lower"
    complete["Section"] = complete["Section"].replace("", np.nan).fillna("Middle")
    secs = {s: complete[complete.Section == s].copy() for s in SECTIONS}
    a = pd.read_csv(RAW / "25_ADCPData.csv", skiprows=3); a["_time"] = pd.to_datetime(a["_time"]); a = a.set_index("_time")
    agg = a["_value"].resample("300min").median().dropna(); agg = agg[np.abs(agg - agg.mean()) <= 2.0 * agg.std()] / 0.3048
    adcp = pd.DataFrame({"Collection Date": agg.index, "Velocity": agg.values, "Depth": 1, "Northing": ADCP_N, "Easting": ADCP_E, "Type": "Sontek", "Section": "Middle"})
    secs["Middle"]["Collection Date"] = pd.to_datetime(secs["Middle"]["Collection Date"], format="mixed", dayfirst=False)
    secs["Middle"] = pd.concat([secs["Middle"], adcp], ignore_index=True)
    return complete, secs


# ---------------------------------------------------------------- cell 7: one row per visit measured at 1 ft ---
def _shape(z, z_p):
    return np.log(np.log(z / z_p))


def _fit_profile_scale(df):
    q = df[(df["Velocity"] > 0) & (df["Depth"] >= FIT_MIN_HEIGHT_FT)].copy()
    q = q[q.groupby(q["_pt"])["Depth"].transform("nunique") >= MIN_LEVELS]
    if q["_pt"].nunique() < 10:
        return np.nan
    lv = np.log(q["Velocity"].values); grp = q["_pt"].values
    lv_c = lv - pd.Series(lv).groupby(grp).transform("mean").values; z = q["Depth"].values
    lo, hi = Z_P_BOUNDS_FT; hi = min(hi, 0.9 * z.min())
    if not (hi > lo):
        return np.nan
    scores = []
    for z_p in np.logspace(np.log10(lo), np.log10(hi), Z_P_STEPS):
        F = _shape(z, z_p); F_c = F - pd.Series(F).groupby(grp).transform("mean").values
        scores.append((np.sum((lv_c - F_c) ** 2), z_p))
    return min(scores, key=lambda t: t[0])[1] if scores else np.nan


def _visit_id(df):
    collected = pd.to_datetime(df["Collection Date"], errors="coerce"); survey_day = collected.dt.strftime("%Y-%m-%d")
    mpid = df["Measurement Point ID"]
    has_mpid = mpid.notna() & mpid.astype(str).str.strip().ne("") & mpid.astype(str).ne("nan")
    fixed = ~has_mpid | (df["Type"].astype(str).eq("Sontek") if "Type" in df.columns else pd.Series(False, index=df.index))
    key = survey_day + "|" + mpid.astype(str)
    return key.where(~fixed, collected.astype(str) + "|ts|" + df.index.astype(str))


def in_water_points(raw):
    """Cell 7 with MEASURED_ONLY for all three sections: the visits measured at 1 ft, in water at least 1 ft deep."""
    work = raw.dropna(subset=["Velocity", "Depth"]).copy(); work = work[work["Depth"] > 0]; work = work[work["Velocity"] >= 0]
    work["_pt"] = _visit_id(work)
    if not np.isfinite(_fit_profile_scale(work)):          # the notebook skipped a section whose profile could not be fitted
        return None
    rows = []
    for _, g in work.groupby("_pt"):
        water_vals = pd.to_numeric(g["Total water depth"], errors="coerce").to_numpy(float) if "Total water depth" in g.columns else np.array([])
        water_vals = water_vals[np.isfinite(water_vals)]; water = float(np.median(water_vals)) if water_vals.size else np.nan
        if np.isfinite(water) and water < TARGET_HEIGHT_FT:
            continue
        at_target = g[np.isclose(g["Depth"], TARGET_HEIGHT_FT)]
        if len(at_target):
            rows.append(at_target.iloc[0].copy())
    out = pd.DataFrame(rows).drop(columns=["_pt"]).reset_index(drop=True)
    return out[np.isclose(pd.to_numeric(out["Depth"], errors="coerce"), TARGET_HEIGHT_FT)].copy()


# ---------------------------------------------------------------- cells 10, 11: main stem and distances ---
def build_main_stem(line_union, points, tol=1.0):
    merged = linemerge(line_union) if line_union.geom_type == "MultiLineString" else line_union
    if merged.geom_type == "LineString":
        return merged
    parts = list(merged.geoms); counts = [0] * len(parts)
    for pt in points:
        if pt is None or pt.is_empty:
            continue
        counts[int(np.argmin([pt.distance(p) for p in parts]))] += 1
    close = lambda a, b: abs(a[0] - b[0]) <= tol and abs(a[1] - b[1]) <= tol
    seed = int(np.argmax(counts)); used = {seed}; chain = [list(parts[seed].coords)]

    def find_next(node):
        best_key, best_j = None, None
        for j, p in enumerate(parts):
            if j in used:
                continue
            c = list(p.coords)
            if close(c[0], node) or close(c[-1], node):
                key = (counts[j], p.length)
                if best_key is None or key > best_key:
                    best_key, best_j = key, j
        return best_j
    while True:
        node = chain[-1][-1]; j = find_next(node)
        if j is None: break
        used.add(j); c = list(parts[j].coords)
        if close(c[-1], node): c.reverse()
        chain.append(c)
    while True:
        node = chain[0][0]; j = find_next(node)
        if j is None: break
        used.add(j); c = list(parts[j].coords)
        if close(c[0], node): c.reverse()
        chain.insert(0, c)
    stem = []
    for c in chain:
        stem.extend(c[1:] if stem else c)
    return LineString(stem)


def distances_to_line(points, line):
    eps = 1.0; L = line.length; signed = np.empty(len(points)); along = np.empty(len(points))
    for i, p in enumerate(points):
        if p is None or p.is_empty:
            signed[i] = along[i] = np.nan; continue
        d = line.project(p); p1 = line.interpolate(max(d - eps, 0)); p2 = line.interpolate(min(d + eps, L))
        cross = (p2.x - p1.x) * (p.y - p1.y) - (p2.y - p1.y) * (p.x - p1.x)
        signed[i] = -p.distance(line) if cross > 0 else p.distance(line); along[i] = d
    return signed, along


# ---------------------------------------------------------------- cell 13: the bank lines ---
def _to_map_frame(line, station, offset):
    station = np.atleast_1d(np.asarray(station, float)); offset = np.atleast_1d(np.asarray(offset, float))
    L = line.length; clipped = np.clip(station, 0.0, L); east, north = np.empty(clipped.shape), np.empty(clipped.shape)
    for i, sv in enumerate(clipped.ravel()):
        here = line.interpolate(float(sv)); back = line.interpolate(max(float(sv) - _EPS_FT, 0.0)); ahead = line.interpolate(min(float(sv) + _EPS_FT, L))
        dx, dy = ahead.x - back.x, ahead.y - back.y; d = np.hypot(dx, dy)
        if d < 1e-12: dx, dy, d = 1.0, 0.0, 1.0
        tx, ty = dx / d, dy / d; nx, ny = ty, -tx                     # right of travel
        j = np.unravel_index(i, clipped.shape); east[j] = here.x + nx * offset.ravel()[i]; north[j] = here.y + ny * offset.ravel()[i]
    return east, north


_e, _n = _to_map_frame(LineString([(0.0, 0.0), (100.0, 0.0)]), [50.0], [10.0])
assert abs(_e[0] - 50.0) < 1e-6 and _n[0] < -9.99, "reach-frame sign convention is mirrored"


def bank_lines(section, centerline, wet, complete):
    wet = wet.dropna(subset=["Easting", "Northing", "dist_along_centerline", "distance_to_centerline"]).copy()
    wet["day"] = pd.to_datetime(wet["Collection Date"], format="mixed").dt.date.astype(str)
    s_water, n_water = wet.dist_along_centerline.values.astype(float), wet.distance_to_centerline.values.astype(float)
    tidy = complete["Point Notes"].astype(str).str.lower().str.replace(r"\s+", " ", regex=True).str.strip()
    is_bank = np.zeros(len(complete), bool)
    for w in BANK_NOTES:
        is_bank |= tidy.str.contains(w, regex=False).values
    bank = complete[(complete.Section == section) & is_bank].dropna(subset=["Easting", "Northing", "Elevation", "Transect ID"]).copy()
    bank = bank.drop_duplicates(subset=["Transect ID", "Easting", "Northing"])
    bank["day"] = pd.to_datetime(bank["Collection Date"], format="mixed").dt.date.astype(str)
    placed, rejected = [], []
    pool = wet[wet.distance_to_centerline.abs() <= CALIB_LIMIT_FT]
    for (tid, day), points in bank.groupby(["Transect ID", "day"]):
        t = pool[(pool["Transect ID"] == tid) & (pool.day == day)].drop_duplicates(subset=["Easting", "Northing"])
        if len(t) < MIN_WET_POINTS:
            continue
        xy = t[["Easting", "Northing"]].values.astype(float); middle = xy.mean(0); direction = np.linalg.svd(xy - middle)[2][0]
        u = (xy - middle) @ direction
        if np.ptp(u) < 1.0:
            continue
        n_t = t.distance_to_centerline.values; slope, intercept = np.polyfit(u, n_t, 1); left_end, right_end = n_t.min(), n_t.max()
        proj = (points[["Easting", "Northing"]].values.astype(float) - middle) @ direction
        for (_, pt), off in zip(points.iterrows(), slope * proj + intercept):
            past_left, past_right = left_end - off, off - right_end
            if max(past_left, past_right) < -OUTSIDE_TOL_FT:
                rejected.append(tid); continue
            placed.append(dict(side="left" if past_left > past_right else "right", station=float(t.dist_along_centerline.median()), offset=float(off), elev=float(pt["Elevation"])))
    P = pd.DataFrame(placed)
    grid = np.arange(s_water.min(), s_water.max() + GRID_FT, GRID_FT)
    cuts = np.arange(s_water.min(), s_water.max() + BIN_FT, BIN_FT); which = np.clip(np.digitize(s_water, cuts) - 1, 0, len(cuts) - 1)
    b_s, b_lo, b_hi = [], [], []
    for b in np.unique(which):
        m = which == b; b_s.append(s_water[m].mean()); b_lo.append(n_water[m].min()); b_hi.append(n_water[m].max())
    survey_edge = {"left": np.interp(grid, b_s, b_lo), "right": np.interp(grid, b_s, b_hi)}
    has_data = np.min(np.abs(grid[:, None] - s_water[None, :]), axis=1) <= WINDOW_FT
    curves, margins = {}, {}; k = int(round(SMOOTH_FT / GRID_FT))
    for side, outward in (("left", -1.0), ("right", 1.0)):
        q = P[P.side == side].sort_values("station"); edge = survey_edge[side]
        if len(q) >= 2:
            near = np.min(np.abs(grid[:, None] - q.station.values[None, :]), axis=1) <= BANK_REACH_FT
            measured = np.where(near, np.interp(grid, q.station.values, q.offset.values), np.nan); ok = np.isfinite(measured)
            margin = max(float(np.median(outward * (measured[ok] - edge[ok]))), 0.0)
        else:
            measured, margin = np.full(grid.shape, np.nan), 10.0
        margins[side] = margin
        curve = np.where(np.isfinite(measured), measured, edge + outward * margin)
        curve = np.convolve(curve, np.ones(k), "same") / np.convolve(np.ones_like(curve), np.ones(k), "same")
        limit = edge + outward * SAFETY_FT; curves[side] = np.where(outward * (limit - curve) > 0, limit, curve)
    narrow = (curves["right"] - curves["left"]) < MIN_WIDTH_FT; mid = 0.5 * (curves["left"] + curves["right"])
    curves["left"] = np.where(narrow, mid - MIN_WIDTH_FT / 2, curves["left"]); curves["right"] = np.where(narrow, mid + MIN_WIDTH_FT / 2, curves["right"])
    gaps = np.flatnonzero(np.diff(grid[has_data]) > MAX_GAP_FT) + 1
    blocks = [b for b in np.split(np.flatnonzero(has_data), gaps) if len(b) >= 2]
    pts = gpd.GeoSeries(gpd.points_from_xy(wet.Easting, wet.Northing), crs=CRS)

    def build_shapes():
        edges, shapes = {}, []
        for side in ("left", "right"):
            edges[side] = [LineString(np.column_stack(_to_map_frame(centerline, grid[idx], curves[side][idx]))) for idx in blocks]
        for j in range(len(blocks)):
            poly = Polygon(np.vstack([np.array(edges["left"][j].coords), np.array(edges["right"][j].coords)[::-1]]))
            shapes.append(poly if poly.is_valid else poly.buffer(0))
        return edges, shapes
    for _ in range(MAX_PASSES):
        edges, polys = build_shapes(); extent = gpd.GeoSeries(polys, crs=CRS).union_all(); outside = ~pts.within(extent.buffer(0.01))
        if not outside.any():
            break
        push = pts[outside].distance(extent) + 0.5; need = {"left": np.zeros(grid.shape), "right": np.zeros(grid.shape)}
        for ps, pn, dist in zip(s_water[outside.values], n_water[outside.values], push.values):
            i = int(np.argmin(np.abs(grid - ps))); side = "left" if abs(pn - curves["left"][i]) <= abs(pn - curves["right"][i]) else "right"
            near = np.abs(grid - ps) <= WINDOW_FT; need[side][near] = np.maximum(need[side][near], dist)
        curves["left"] -= need["left"]; curves["right"] += need["right"]
    kept = np.concatenate(blocks); width = (curves["right"] - curves["left"])[kept]
    ext = gpd.GeoDataFrame([dict(section=section, stretch=j + 1, area_acres=round(p.area / 43560.0, 2), bank_elev_min_ft=round(float(P.elev.min()), 2),
                                 bank_elev_max_ft=round(float(P.elev.max()), 2), note="streambed at maximum stage; surveyed bank where available, otherwise the survey edge plus the measured margin",
                                 geometry=p) for j, p in enumerate(polys)], crs=CRS)
    lines = gpd.GeoDataFrame([dict(section=section, bank=side, stretch=j + 1, geometry=g) for side in ("left", "right") for j, g in enumerate(edges[side])], crs=CRS)
    summ = {"bank_points": int(len(bank)), "bank_points_placed": int(len(P)), "bank_points_rejected": len(rejected), "margin_ft": margins,
            "median_width_ft": round(float(np.median(width)), 1), "acres": round(sum(p.area for p in polys) / 43560.0, 2), "points_outside": int(outside.sum()), "in_water_points": int(len(wet))}
    return lines, ext, summ


if __name__ == "__main__":
    complete, secs = survey_by_section()
    network = gpd.read_file(RAW / "Huron_RiverSystem.geojson").to_crs(CRS)
    summary = {}
    for s in SECTIONS:
        wet = in_water_points(secs[s])
        # the notebook wrote these rows to mapping/data_{s}.geojson (cell 8) and read them back for the distances (cell 11),
        # then through sorted_data/*.csv (cells 11, 12); both round trips are repeated so dates and empty points behave the same
        tmp = Path(tempfile.mkdtemp()) / f"data_{s}.geojson"
        gpd.GeoDataFrame(wet, geometry=[Point(xy) for xy in zip(wet.Easting, wet.Northing)], crs=CRS).to_file(tmp, driver="GeoJSON")
        g = gpd.read_file(tmp); wet = pd.DataFrame(g.drop(columns="geometry"))
        minx, miny, maxx, maxy = g.total_bounds                    # the centreline is built before the northing trim (cell 11 before cell 12)
        stem = build_main_stem(network[network.geometry.intersects(box(minx - 500, miny - 500, maxx + 500, maxy + 500))].geometry.union_all(), g.geometry.values)
        wet["distance_to_centerline"], wet["dist_along_centerline"] = distances_to_line(g.geometry.values, stem)
        buf = io.StringIO(); wet.to_csv(buf, index=False); buf.seek(0); wet = pd.read_csv(buf)
        wet = wet[wet.Northing <= NORTHING_MAX_FT[s]]
        wet[["Type", "Easting", "Northing"]].to_csv(OUT / f"in_water_points_{s}.csv", index=False)   # type and position only: they set the 3DEP transects' extent (banklines_3dep.transect_points)
        lines, ext, summ = bank_lines(s, stem, wet, complete)
        gpd.GeoDataFrame(geometry=[stem], crs=CRS).to_file(OUT / f"centerline_union_{s}.geojson", driver="GeoJSON")
        lines.to_file(OUT / f"bank_lines_{s}.geojson", driver="GeoJSON"); ext.to_file(OUT / f"channel_extent_{s}.geojson", driver="GeoJSON")
        # the check: does this reproduce the files the analysis used from September to November 2026 (copied then from the main branch's
        # mapping/ folder, now removed from data_raw so that every line in the folder has a script that writes it)? Fingerprints of those
        # copies: sha256 of the coordinates in ft rounded to 0.1 ft, first 16 hex digits.
        summ["reproduces_september_copies"] = {k: fingerprint(OUT / f"{k}_{s}.geojson") == SEPTEMBER_FINGERPRINTS[s][k] for k in SEPTEMBER_FINGERPRINTS[s]}
        summary[s] = summ; print(s, json.dumps(summ))
    json.dump(summary, open(OUT / "summary.json", "w"), indent=1)
