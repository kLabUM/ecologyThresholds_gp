"""Bank lines for the Upper, Middle and Lower reaches from USGS 3DEP lidar (hydro-flattened water).

The bank lines and channel polygons in mapping/ (bank_lines_{Section}.geojson,
channel_extent_{Section}.geojson) were scraped from an older code chunk and are not reliable.
This script rebuilds them from an external, reproducible source, the way the Huron full-river
model did it (SWMMAnywhereHuron, twin_kayla/scripts/build_bed.py step 1.7 'Top width' and
build_profile.py '3DEP').

Method (per reach)
------------------
1. Transects perpendicular to mapping/centerline_union_{Section}.geojson every TRANSECT_EVERY_M,
   over the surveyed stretch plus PAD_M at each end, sampled every TRANSECT_STEP_M out to
   +/- TRANSECT_HALF_M. Elevations come from the USGS 3DEP ImageServer getSamples (best
   available resolution, 1 m lidar here, nearest neighbour). Requests are cached and resumable.
2. Open water in 3DEP lidar DEMs is hydro-flattened: across a river the water is a run of cells
   at one elevation. On each transect the lowest cell within CENTRE_SEARCH_M of the centreline
   is the water; the water is the contiguous run of cells within FLAT_TOL_M of it, and its two
   ends are the water's edges. A run that reaches the end of the transect (a pond) or is
   shorter than MIN_FLAT_CELLS (no flat water) is set aside and filled from its neighbours.
3. The edge offsets are cleaned along the reach (rolling median over EDGE_SMOOTH transects),
   turned into a left and a right bank line, and the polygon between them. 'Left' and 'right'
   follow the existing mapping/bank_lines files, so code that reads them does not change sign.

What the lines are: the wetted edge on the lidar day (USGS flights are flown at low to
moderate flow), not bankfull. The check figure compares them with the old lines and with the
surveyed transects, whose end points were taken at the water's edge.

Outputs (the old files in mapping/ are not touched)
-------
banklines/bank_lines_{Section}_3dep.geojson       left and right bank lines (ESRI:102690, US ft)
banklines/channel_extent_{Section}_3dep.geojson   the polygon between them
banklines/edges_{Section}.csv                     per transect: edges, width, lidar-day water level, note
banklines/summary.json                            widths and the checks against the survey, old vs new
banklines/banklines_3dep_check.png                the visual check, all three reaches
banklines/raw/3dep_transects_{Section}.csv        cached 3DEP samples (delete, or --refresh, to resample)

Run from the repository root (needs internet the first time; about 90 requests, a few minutes):
    python banklines/build_banklines_3dep.py
    python banklines/build_banklines_3dep.py --sites Middle --refresh
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import linemerge

# ---------------------------------------------------------------- assumptions ---
SITES = ("Upper", "Middle", "Lower")
TRANSECT_EVERY_M = 10.0       # along the reach (the full-river model used 250 m; these reaches are under 1 km)
TRANSECT_HALF_M = 80.0        # half length of each transect (widest reach about 40 m, centreline off by up to about 25 m)
TRANSECT_STEP_M = 1.0         # the 3DEP 1 m grid
PAD_M = 50.0                  # transects extend this far beyond the surveyed stretch at each end
CENTRE_SEARCH_M = 25.0        # lowest cell searched within this of the centreline
FLAT_TOL_M = 0.10             # water = contiguous cells within this of the lowest cell (as in twin_kayla)
MIN_FLAT_CELLS = 3            # fewer flat cells = no flat water seen
EDGE_SMOOTH = 5               # rolling median over this many transects (50 m)
SURVEY_TYPES = ("Transect", "Other", "Thalweg")   # survey rows that mark the surveyed stretch
EPQS = "https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/getSamples"
BATCH = 500                   # points per getSamples request
CRS_FT, CRS_GEO = "ESRI:102690", "EPSG:4326"   # the survey's StatePlane Michigan South (US ft)
FT = 0.3048

ROOT = Path(__file__).resolve().parents[1]
MAP, OUT = ROOT / "mapping", ROOT / "banklines"
RAW = OUT / "raw"


# ------------------------------------------------------------------ geometry ---
def centreline(site: str):
    cl = gpd.read_file(MAP / f"centerline_union_{site}.geojson").to_crs(CRS_FT)
    u = cl.geometry.union_all() if hasattr(cl.geometry, "union_all") else cl.geometry.unary_union
    return u if u.geom_type == "LineString" else linemerge(u)


def survey(site: str) -> gpd.GeoDataFrame:
    d = gpd.read_file(MAP / f"data_{site}.geojson").to_crs(CRS_FT)
    d = d[d["Type"].isin(SURVEY_TYPES)].copy()
    d["x"], d["y"] = d.geometry.x, d.geometry.y
    return d


def transects(line, xs, ys) -> pd.DataFrame:
    """Points across the centreline every TRANSECT_EVERY_M over the surveyed stretch (+ PAD_M)."""
    d = np.array([line.project(Point(x, y)) for x, y in zip(xs, ys)])
    s0, s1 = max(d.min() - PAD_M / FT, 0), min(d.max() + PAD_M / FT, line.length)
    stations = np.arange(s0, s1 + 1e-9, TRANSECT_EVERY_M / FT)
    offs = np.arange(-TRANSECT_HALF_M, TRANSECT_HALF_M + 1e-9, TRANSECT_STEP_M)
    rows = []
    for i, s in enumerate(stations):
        p = line.interpolate(s)
        a, b = line.interpolate(max(s - 10 / FT, 0)), line.interpolate(min(s + 10 / FT, line.length))
        tx, ty = b.x - a.x, b.y - a.y
        nrm = np.hypot(tx, ty) or 1.0
        nx, ny = -ty / nrm, tx / nrm                      # unit normal, left of the digitised direction
        for o in offs:
            rows.append((i, (s - s0) * FT, o, p.x + o / FT * nx, p.y + o / FT * ny))
    t = pd.DataFrame(rows, columns=["transect", "station_m", "offset_m", "x_ft", "y_ft"])
    ll = gpd.GeoSeries(gpd.points_from_xy(t.x_ft, t.y_ft), crs=CRS_FT).to_crs(CRS_GEO)
    t["lon"], t["lat"] = ll.x.round(7), ll.y.round(7)
    return t


# --------------------------------------------------------------------- 3DEP ---
def sample_3dep(pts: pd.DataFrame, cache: Path) -> np.ndarray:
    """3DEP elevation at every point. Each finished request is appended to `cache` at once, so
    an interrupted run resumes; a complete cache is reused with no network access."""
    import requests
    n_batches = int(np.ceil(len(pts) / BATCH))
    done = {}
    if cache.exists():
        prev = pd.read_csv(cache)
        if len(prev) and prev["n_points"].iloc[0] == len(pts):      # same point list
            done = {int(k): g for k, g in prev.groupby("batch")}
        else:
            cache.unlink()
    todo = [k for k in range(n_batches) if k not in done]
    print(f"  3DEP at {len(pts)} points: {len(done)} of {n_batches} requests cached, {len(todo)} to go")
    for j, k in enumerate(todo):
        chunk = pts.iloc[k * BATCH:(k + 1) * BATCH]
        geom = {"points": chunk[["lon", "lat"]].values.tolist(), "spatialReference": {"wkid": 4326}}
        params = {"geometry": json.dumps(geom), "geometryType": "esriGeometryMultipoint",
                  "returnFirstValueOnly": "true", "interpolation": "RSP_NearestNeighbor",
                  "outFields": "", "f": "json"}
        body = None
        for attempt in range(6):                           # the service returns 5xx under load
            try:
                r = requests.post(EPQS, data=params, timeout=180)
                r.raise_for_status()
                body = r.json()
                if "error" in body:
                    raise RuntimeError(body["error"])
                break
            except Exception as ex:
                if attempt == 5:
                    raise SystemExit(f"3DEP request {k + 1}/{n_batches} failed after 6 tries ({ex}). "
                                     f"{len(done)} requests are saved; run the script again to continue.")
                time.sleep(15 * (attempt + 1))
        vals, res = np.full(len(chunk), np.nan), np.full(len(chunk), np.nan)
        for smp in body.get("samples", []):
            try:
                v = float(smp["value"])
            except (TypeError, ValueError):
                v = np.nan
            vals[int(smp["locationId"])] = v if v > -1000 else np.nan   # NoData is a large negative
            res[int(smp["locationId"])] = smp.get("resolution", np.nan)
        rows = pd.DataFrame({"batch": k, "n_points": len(pts), "i": np.arange(len(chunk)),
                             "elev_m": vals, "resolution_m": res})
        rows.to_csv(cache, mode="a", header=not cache.exists(), index=False)
        done[k] = rows
        if (j + 1) % 10 == 0 or j + 1 == len(todo):
            print(f"    {len(done)}/{n_batches} requests done")
        time.sleep(0.2)
    out = np.full(len(pts), np.nan)
    for k, g in done.items():
        out[k * BATCH + g["i"].to_numpy()] = g["elev_m"].to_numpy()
    return out


def flat_edges(t: pd.DataFrame) -> pd.DataFrame:
    """Per transect: the water's two edge offsets (m from the centreline) and its level."""
    out = []
    for i, g in t.groupby("transect"):
        g = g.sort_values("offset_m")
        z, o = g["elev_m"].to_numpy(float), g["offset_m"].to_numpy(float)
        rec = {"transect": i, "station_m": float(g.station_m.iloc[0]),
               "edge_neg_m": np.nan, "edge_pos_m": np.nan, "water_m": np.nan}
        centre = np.where(np.abs(o) <= CENTRE_SEARCH_M)[0]
        if not len(centre) or np.all(np.isnan(z[centre])):
            out.append({**rec, "note": "no data"})
            continue
        k = centre[np.nanargmin(z[centre])]
        flat = np.abs(z - z[k]) <= FLAT_TOL_M
        lo, hi = k, k
        while lo > 0 and flat[lo - 1]:
            lo -= 1
        while hi < len(z) - 1 and flat[hi + 1]:
            hi += 1
        if lo == 0 or hi == len(z) - 1:
            note = "pond (flat beyond the transect)"
        elif hi - lo + 1 < MIN_FLAT_CELLS:
            note = "no flat water (not hydro-flattened)"
        else:
            note = "3DEP flat water"
            rec.update(edge_neg_m=o[lo] - TRANSECT_STEP_M / 2, edge_pos_m=o[hi] + TRANSECT_STEP_M / 2,
                       water_m=float(z[k]))           # the edge is half a cell beyond the last flat cell
        out.append({**rec, "note": note})
    e = pd.DataFrame(out).sort_values("transect").reset_index(drop=True)
    ok = e["note"] == "3DEP flat water"
    for c in ("edge_neg_m", "edge_pos_m"):
        v = e[c].where(ok).rolling(EDGE_SMOOTH, center=True, min_periods=2).median()
        e[c + "_smooth"] = v.interpolate(limit_direction="both")
    e["width_m"] = e["edge_pos_m_smooth"] - e["edge_neg_m_smooth"]
    return e


def build_banks(site, line, t, e, old):
    """Bank lines and the polygon between them, with 'left' on the side the old file calls left."""
    ends = t.groupby("transect")
    c = t[t.offset_m == 0].set_index("transect")
    nx = ((ends.x_ft.last() - ends.x_ft.first()) / (2 * TRANSECT_HALF_M / FT)).to_numpy()
    ny = ((ends.y_ft.last() - ends.y_ft.first()) / (2 * TRANSECT_HALF_M / FT)).to_numpy()

    def edge_line(col):
        off = e[col].to_numpy() / FT
        return LineString(np.column_stack([c.x_ft.to_numpy() + off * nx, c.y_ft.to_numpy() + off * ny]))

    pos, neg = edge_line("edge_pos_m_smooth"), edge_line("edge_neg_m_smooth")
    ol = old[old.bank == "left"].geometry.iloc[0]
    m = ol.interpolate(0.5, normalized=True)
    s = line.project(m); p = line.interpolate(s)
    a, b = line.interpolate(max(s - 30, 0)), line.interpolate(min(s + 30, line.length))
    left_is_pos = ((b.x - a.x) * (m.y - p.y) - (b.y - a.y) * (m.x - p.x)) > 0
    left, right = (pos, neg) if left_is_pos else (neg, pos)
    src = "USGS 3DEP 1 m hydro-flattened water edge (banklines/build_banklines_3dep.py)"
    banks = gpd.GeoDataFrame({"section": [site, site], "bank": ["left", "right"], "stretch": [1, 1], "source": [src, src]},
                             geometry=[left, right], crs=CRS_FT)
    poly = Polygon(list(left.coords) + list(right.coords)[::-1]).buffer(0)
    ext = gpd.GeoDataFrame({"section": [site], "stretch": [1], "area_acres": [round(poly.area / 43560.0, 2)],
                            "note": ["3DEP hydro-flattened water on the lidar day; edges every 10 m, 50 m running median"]},
                           geometry=[poly], crs=CRS_FT)
    return banks, ext


# -------------------------------------------------------------------- checks ---
def transect_ends(d: gpd.GeoDataFrame) -> pd.DataFrame:
    """The two outermost points of every surveyed transect: the crew's water edges."""
    rows = []
    for (_, tid), g in d[d["Type"] == "Transect"].groupby(["Collection Date", "Transect ID"]):
        if len(g) < 3:
            continue
        xy = g[["x", "y"]].to_numpy()
        dd = np.hypot(*(xy[:, None, :] - xy[None, :, :]).transpose(2, 0, 1))
        i, j = np.unravel_index(dd.argmax(), dd.shape)
        rows += [xy[i], xy[j]]
    return pd.DataFrame(rows, columns=["x", "y"])


def checks(d, ends, banks, ext) -> dict:
    lines, poly = banks.geometry.tolist(), ext.geometry.iloc[0]
    dist = np.array([min(l.distance(Point(x, y)) for l in lines) * FT for x, y in zip(ends.x, ends.y)])
    outside = np.array([not poly.buffer(1 / FT).contains(Point(x, y)) for x, y in zip(d.x, d.y)])
    return {"transect_ends": int(len(ends)), "transect_end_to_bank_median_m": round(float(np.median(dist)), 2),
            "transect_end_to_bank_p90_m": round(float(np.quantile(dist, 0.9)), 2),
            "survey_points_outside_polygon": int(outside.sum()), "survey_points": int(len(d)),
            "area_acres": float(ext.area_acres.iloc[0])}


# ---------------------------------------------------------------------- figure ---
def draw(results, path):
    n = len(results)
    fig = plt.figure(figsize=(6.4 * n, 11))
    gs = fig.add_gridspec(2, n, height_ratios=[3.2, 1], hspace=0.25, wspace=0.15)
    BLUE, ORANGE, GREY, INK = "#2563eb", "#ea580c", "#9ca3af", "#1f2937"
    for j, (site, r) in enumerate(results.items()):
        d, ends, banks, ext, old, line, e, sm = r["d"], r["ends"], r["banks"], r["ext"], r["old"], r["line"], r["e"], r["summary"]
        x0, y0 = d.x.median(), d.y.median()
        X = lambda v: (np.asarray(v) - x0) * FT
        Y = lambda v: (np.asarray(v) - y0) * FT
        ax = fig.add_subplot(gs[0, j]); ax.set_aspect("equal")
        px, py = ext.geometry.iloc[0].exterior.xy; ax.fill(X(px), Y(py), color="#dbeafe", lw=0, zorder=0)
        lx, ly = line.xy; ax.plot(X(lx), Y(ly), color=GREY, lw=0.8, ls=":", zorder=1)
        for g in old.geometry:
            gx, gy = g.xy; ax.plot(X(gx), Y(gy), color=ORANGE, lw=1.3, ls=(0, (4, 2)), zorder=3)
        for g in banks.geometry:
            gx, gy = g.xy; ax.plot(X(gx), Y(gy), color=BLUE, lw=1.7, zorder=4)
        ax.scatter(X(d.x), Y(d.y), s=5, color=INK, lw=0, zorder=5)
        ax.scatter(X(ends.x), Y(ends.y), s=22, facecolor="white", edgecolor=INK, lw=0.8, zorder=6)
        b = ext.geometry.iloc[0].bounds
        ax.set_xlim(min(X(d.x).min(), X(b[0])) - 20, max(X(d.x).max(), X(b[2])) + 20)
        ax.set_ylim(min(Y(d.y).min(), Y(b[1])) - 20, max(Y(d.y).max(), Y(b[3])) + 20)
        ax.set_xlabel("east (m)"); ax.set_ylabel("north (m)" if j == 0 else "")
        n3, no = sm["new"], sm["old"]
        ax.set_title(f"{site}: water width {sm['width_median_m']:.0f} m (median), {sm['flat_water_transects']} of {sm['transects']} transects flat\n"
                     f"transect ends to bank line (median): 3DEP {n3['transect_end_to_bank_median_m']:.1f} m, old {no['transect_end_to_bank_median_m']:.1f} m\n"
                     f"survey points outside the polygon: 3DEP {n3['survey_points_outside_polygon']}, old {no['survey_points_outside_polygon']} (of {n3['survey_points']})",
                     fontsize=9.5, color=INK, loc="left")
        a2 = fig.add_subplot(gs[1, j]); ok = e["note"] == "3DEP flat water"
        a2.plot(e.station_m[ok], (e.edge_pos_m - e.edge_neg_m)[ok], ".", color=GREY, ms=3, label="raw flat-water width")
        a2.plot(e.station_m, e.width_m, color=BLUE, lw=1.6, label="smoothed (50 m running median)")
        if (~ok).any():
            a2.plot(e.station_m[~ok], np.full((~ok).sum(), 0.5), "|", color=ORANGE, ms=9, label="no flat water or pond (filled)")
        a2.set_ylim(bottom=0); a2.grid(alpha=0.3)
        a2.set_xlabel("distance along the reach (m)"); a2.set_ylabel("water width (m)" if j == 0 else "")
        if j == 0:
            a2.legend(fontsize=8, frameon=False)
    fig.legend(handles=[Line2D([], [], color=BLUE, lw=1.7, label="3DEP bank lines (new)"),
                        Line2D([], [], color=ORANGE, lw=1.3, ls=(0, (4, 2)), label="old bank lines"),
                        Line2D([], [], color=GREY, lw=0.8, ls=":", label="centreline"),
                        Line2D([], [], marker="o", ls="", color=INK, ms=3, label="surveyed points"),
                        Line2D([], [], marker="o", ls="", mfc="white", mec=INK, ms=6, label="transect ends (water's edge)")],
               loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=5, frameon=False, fontsize=10)
    fig.suptitle("Bank lines from USGS 3DEP hydro-flattened water (wetted edge on the lidar day) against the old bank lines",
                 y=0.985, fontsize=12)
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", nargs="+", default=list(SITES), choices=SITES)
    ap.add_argument("--refresh", action="store_true", help="resample 3DEP")
    a = ap.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    results, summary = {}, {}
    for site in a.sites:
        print(f"{site}")
        line, d = centreline(site), survey(site)
        old = gpd.read_file(MAP / f"bank_lines_{site}.geojson").to_crs(CRS_FT)
        old_ext = gpd.read_file(MAP / f"channel_extent_{site}.geojson").to_crs(CRS_FT)
        t = transects(line, d.x.to_numpy(), d.y.to_numpy())
        cache = RAW / f"3dep_transects_{site}.csv"
        if a.refresh:
            cache.unlink(missing_ok=True)
        t["elev_m"] = sample_3dep(t, cache)
        e = flat_edges(t)
        e.to_csv(OUT / f"edges_{site}.csv", index=False)
        banks, ext = build_banks(site, line, t, e, old)
        banks.to_file(OUT / f"bank_lines_{site}_3dep.geojson", driver="GeoJSON")
        ext.to_file(OUT / f"channel_extent_{site}_3dep.geojson", driver="GeoJSON")
        ends = transect_ends(d)
        ok = e["note"] == "3DEP flat water"
        sm = {"transects": int(len(e)), "flat_water_transects": int(ok.sum()), "notes": e["note"].value_counts().to_dict(),
              "width_median_m": round(float(e.width_m.median()), 1),
              "width_p10_p90_m": [round(float(e.width_m.quantile(.1)), 1), round(float(e.width_m.quantile(.9)), 1)],
              "lidar_water_level_median_m": round(float(e.water_m.median()), 2),
              "new": checks(d, ends, banks, ext), "old": checks(d, ends, old, old_ext)}
        summary[site] = sm
        print(json.dumps(sm, indent=1))
        results[site] = {"d": d, "ends": ends, "banks": banks, "ext": ext, "old": old, "line": line, "e": e, "summary": sm}
    json.dump(summary, open(OUT / "summary.json", "w"), indent=1)
    draw(results, OUT / "banklines_3dep_check.png")
    print(f"Wrote banklines/: bank lines, polygons, edges, summary.json and banklines_3dep_check.png")


if __name__ == "__main__":
    main()
