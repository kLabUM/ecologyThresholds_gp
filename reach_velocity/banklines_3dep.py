"""Bank lines for the three reaches: the old lines (A), a fixed survey + 3DEP hybrid (B), and stage-dependent lines with
their uncertainty (C). The rule for choosing between them is logbook entry 27; 43_bankline_comparison.py applies it.

3DEP part (after the Huron full-river model, SWMMAnywhereHuron twin_kayla/scripts/build_bed.py step 1.7): transects every
10 m across the reach, sampled every 1 m out to +/- 80 m from the USGS 3DEP ImageServer (cached in
data_raw/3dep_bank_transects_{reach}.csv; sample_3dep() refetches if the cache is deleted). Open water in 3DEP lidar is
hydro-flattened, so on each transect the water is the run of cells within 0.10 m of the lowest cell near the centreline,
and its ends are the lidar-day water edges.

Survey part: each surveyed transect's two end points, moved 1 m outward (the crew stops at a wadeable depth), are the
water edges on that survey day; their level is the sensor's water level carried to the section (water_level()).

Offsets are metres from the 3DEP station centre along the station's normal ('neg' and 'pos' sides); 'left' and 'right'
follow the old files' convention (side_map()).
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import linemerge

HERE = Path(__file__).resolve().parent
RAW, D = HERE / "data_raw", HERE / "data_clean"
FT = 0.3048
CRS_FT = "ESRI:102690"

# ---------------------------------------------------------------- assumptions (logbook entry 27) ---
FLAT_TOL_M = 0.10            # water = contiguous cells within this of the lowest cell (as in twin_kayla)
CENTRE_SEARCH_M = 25.0       # lowest cell searched within this of the centreline
MIN_FLAT_CELLS = 3           # fewer flat cells = no flat water seen
EDGE_SMOOTH = 5              # running median over this many 10 m stations (50 m)
SURVEY_MARGIN_M = 1.0        # surveyed transect ends moved this far outward (decision 3)
CANOPY_RATIO = 0.8           # lidar flat width below this x the surveyed width: use the surveyed edges (decision 5)
LIDAR_SD_M = 0.10            # 3DEP vertical error
GAP_FRAC = 0.25              # interpolation uncertainty: this fraction of the survey-to-lidar edge gap
SUBSTITUTE_SD_M = 2.0        # uncertainty where the surveyed edge stands in for the lidar
EDGE_SD_MAX_M = 6.0          # cap on one edge's sd
CORR_M = 50.0                # along-reach correlation length of the edge errors
N_DRAWS = 20                 # draws of C (decision 6)
TIE = {  # sensor stage s (m) -> water level (m NAVD88) at a station: level = a + b s - slope * (along - along_ref)
    "Middle": {"a": 253.55, "b": 1.0, "slope": 0.050 / 100, "along_ref": 0.0, "sd": 0.02},   # entry 10: zero at 253.55 m, 5.0 cm / 100 m, along from the ADCP
}
STATION_FIELDS = ["transect", "station_m", "offset_m", "x_ft", "y_ft", "lon", "lat", "elev_m"]


# ---------------------------------------------------------------- inputs ---
SURVEY_LINES = D / "banklines_survey"   # written by 44_banklines_survey.py (centreline, set A, in-water point positions)


def centreline(reach):
    cl = gpd.read_file(SURVEY_LINES / f"centerline_union_{reach}.geojson").to_crs(CRS_FT)
    u = cl.geometry.union_all() if hasattr(cl.geometry, "union_all") else cl.geometry.unary_union
    return u if u.geom_type == "LineString" else linemerge(u)


def sensor_rating(reach):
    """Upper and Lower: the sections' water level L at the reference position (RTK surface, slope-corrected) against the sensor stage s on
    the tied days, L = a + b s, with the fit's scatter (entries 21 and 26). Inputs from 37_other_reaches_prep.py; 38_ uses the same fit."""
    pts = pd.read_csv(D / f"reach_{reach}_points_sensor.csv"); T = pd.read_csv(D / f"reach_{reach}_tie.csv")
    S = float(pts.slope_cm_per_100m.iloc[0]) / 1e4; along_ref = float(pts.along_m.median())
    good = T[T.rtk_good & T.stage_m.notna()].copy(); good["L"] = good.wse_rtk_m + S * (good.along_m - along_ref)
    day = good.groupby("date").agg(s=("stage_m", "median"), L=("L", "mean"), n=("transect", "size"))
    b = np.polyfit(day.s, day.L, 1); day["L_fit"] = np.polyval(b, day.s); sd = float((day.L - day.L_fit).std(ddof=2))
    return {"a": float(b[1]), "b": float(b[0]), "slope": S, "along_ref": along_ref, "sd": sd, "day": day, "coef": b}


def tie(reach):
    """The sensor-to-section tie: Middle from entry 10, Upper and Lower from their sensor-to-section rating (entries 21, 26)."""
    if reach in TIE:
        return TIE[reach]
    k = sensor_rating(reach)
    return {c: k[c] for c in ("a", "b", "slope", "along_ref", "sd")}


def survey_points(reach):
    """Every surveyed location with its date, transect, position, along_m and the sensor stage at the reading."""
    if reach == "Middle":
        p = pd.read_csv(D / "locations.csv")[["location_id", "date", "transect", "Easting", "Northing", "along_m", "stage_m", "a"]]
    else:
        p = pd.read_csv(D / f"reach_{reach}_points_sensor.csv")[["location_id", "date", "transect", "Easting", "Northing", "along_m", "stage_m", "a"]]
    return p.dropna(subset=["Easting", "Northing"]).reset_index(drop=True)


def load_transects(reach, check=True):
    """The 3DEP transects with their elevations (data_raw cache, made by fetch_transects). Their positions are checked against
    transect_points(): a cache made on other positions is refused, so the elevations always belong to this folder's geometry."""
    t = pd.read_csv(RAW / f"3dep_bank_transects_{reach}.csv")[STATION_FIELDS]
    if check:
        err = check_transects(reach, t)
        if err > 0.01:
            raise SystemExit(f"{reach}: the cached 3DEP transects sit up to {err:.2f} m from transect_points(); refetch with "
                             f"'python banklines_3dep.py fetch {reach}' (needs internet)")
    return t


# ---------------------------------------------------------------- the 3DEP transects (geometry and fetch) ---
TRANSECT_EVERY_M, TRANSECT_HALF_M, TRANSECT_STEP_M, PAD_M = 10.0, 80.0, 1.0, 50.0   # every 10 m, +/- 80 m, 1 m step, 50 m beyond the survey
SURVEY_TYPES = ("Transect", "Other", "Thalweg")


def transect_points(reach):
    """The points where 3DEP is sampled: transects across the centreline (44_) every 10 m over the stretch of the in-water survey
    points (44_, the types that mark the surveyed stretch), padded 50 m at each end; normals to the left of the digitised direction.
    Carried over from the fetch first run outside this folder (Code for Paper banklines/build_banklines_3dep.py, 18 Nov 2026)."""
    line = centreline(reach); w = pd.read_csv(SURVEY_LINES / f"in_water_points_{reach}.csv"); w = w[w.Type.isin(SURVEY_TYPES)]
    d = np.array([line.project(Point(x, y)) for x, y in zip(w.Easting, w.Northing)])
    s0, s1 = max(d.min() - PAD_M / FT, 0), min(d.max() + PAD_M / FT, line.length)
    offs = np.arange(-TRANSECT_HALF_M, TRANSECT_HALF_M + 1e-9, TRANSECT_STEP_M); rows = []
    for i, s_ in enumerate(np.arange(s0, s1 + 1e-9, TRANSECT_EVERY_M / FT)):
        p = line.interpolate(s_); a, b = line.interpolate(max(s_ - 10 / FT, 0)), line.interpolate(min(s_ + 10 / FT, line.length))
        tx, ty = b.x - a.x, b.y - a.y; nrm = np.hypot(tx, ty) or 1.0; nx, ny = -ty / nrm, tx / nrm
        rows += [(i, (s_ - s0) * FT, o, p.x + o / FT * nx, p.y + o / FT * ny) for o in offs]
    t = pd.DataFrame(rows, columns=["transect", "station_m", "offset_m", "x_ft", "y_ft"])
    ll = gpd.GeoSeries(gpd.points_from_xy(t.x_ft, t.y_ft), crs=CRS_FT).to_crs(4326); t["lon"], t["lat"] = ll.x.round(7), ll.y.round(7)
    return t


def check_transects(reach, t=None):
    """Largest distance (m) between the cached transect points and transect_points(); inf if the counts differ."""
    t = pd.read_csv(RAW / f"3dep_bank_transects_{reach}.csv") if t is None else t; g = transect_points(reach)
    if len(t) != len(g):
        return np.inf
    a, b = t.sort_values(["transect", "offset_m"]), g.sort_values(["transect", "offset_m"])
    return float(np.max(np.hypot(a.x_ft.to_numpy() - b.x_ft.to_numpy(), a.y_ft.to_numpy() - b.y_ft.to_numpy())) * FT)


def fetch_transects(reach):
    """Sample 3DEP at transect_points() and write the cache data_raw/3dep_bank_transects_{reach}.csv (needs internet)."""
    t = transect_points(reach); t["elev_m"] = sample_3dep(t); t.to_csv(RAW / f"3dep_bank_transects_{reach}.csv", index=False); return t


def stations(reach, t):
    """One row per 10 m station: centre point, unit normal (pos side), along_m in the reach's own frame."""
    g = t.groupby("transect"); half = (t.offset_m.max() - t.offset_m.min()) / 2
    c = t[t.offset_m.abs() < 1e-6].set_index("transect")
    st = pd.DataFrame({"station_m": c.station_m, "x": c.x_ft, "y": c.y_ft,
                       "nx": (g.x_ft.last() - g.x_ft.first()) / (2 * half / FT), "ny": (g.y_ft.last() - g.y_ft.first()) / (2 * half / FT)})
    line, p = centreline(reach), survey_points(reach)
    origin = np.mean([line.project(Point(e, n)) * FT - a for e, n, a in zip(p.Easting, p.Northing, p.along_m)])
    st["along_m"] = [line.project(Point(x, y)) * FT - origin for x, y in zip(st.x, st.y)]
    return st.reset_index()


def water_level(reach, stage, along_m):
    k = tie(reach)
    return k["a"] + k["b"] * np.asarray(stage) - k["slope"] * (np.asarray(along_m) - k["along_ref"])


# ---------------------------------------------------------------- 3DEP edges ---
def flat_edges(t):
    out = []
    for i, g in t.groupby("transect"):
        g = g.sort_values("offset_m"); z, o = g.elev_m.to_numpy(float), g.offset_m.to_numpy(float)
        rec = {"transect": i, "edge_neg_m": np.nan, "edge_pos_m": np.nan, "water_m": np.nan, "note": ""}
        centre = np.where(np.abs(o) <= CENTRE_SEARCH_M)[0]
        if not len(centre) or np.all(np.isnan(z[centre])):
            out.append({**rec, "note": "no data"}); continue
        k = centre[np.nanargmin(z[centre])]; flat = np.abs(z - z[k]) <= FLAT_TOL_M; lo = hi = k
        while lo > 0 and flat[lo - 1]: lo -= 1
        while hi < len(z) - 1 and flat[hi + 1]: hi += 1
        if lo == 0 or hi == len(z) - 1: rec["note"] = "pond"
        elif hi - lo + 1 < MIN_FLAT_CELLS: rec["note"] = "no flat water"
        else: rec.update(edge_neg_m=o[lo] - 0.5, edge_pos_m=o[hi] + 0.5, water_m=float(z[k]), note="flat water")
        out.append(rec)
    e = pd.DataFrame(out).sort_values("transect").reset_index(drop=True); ok = e.note == "flat water"
    for c in ("edge_neg_m", "edge_pos_m", "water_m"):
        e[c + "_s"] = e[c].where(ok).rolling(EDGE_SMOOTH, center=True, min_periods=2).median().interpolate(limit_direction="both")
    e["flat_ok"] = ok
    return e


# ---------------------------------------------------------------- survey edges ---
def survey_anchors(reach, st):
    """Each surveyed transect's outermost points as edges (1 m outward), on the station frame, with the day's level."""
    p = survey_points(reach); cl = LineString(st[["x", "y"]].to_numpy())
    rows = []
    for (date, tr), g in p.groupby(["date", "transect"]):
        if len(g) < 3:
            continue
        s_pt = np.array([cl.project(Point(e, n)) for e, n in zip(g.Easting, g.Northing)])
        i = np.clip(np.searchsorted(np.r_[0, np.cumsum(np.hypot(np.diff(st.x), np.diff(st.y)))], s_pt) - 1, 0, len(st) - 1)
        off = ((g.Easting.to_numpy() - st.x.to_numpy()[i]) * st.nx.to_numpy()[i] + (g.Northing.to_numpy() - st.y.to_numpy()[i]) * st.ny.to_numpy()[i]) * FT
        k = int(np.round(np.median(i))); stage = float(g.stage_m.median())
        lev = float(water_level(reach, stage, st.along_m.iat[k]))
        rows.append({"date": date, "transect": tr, "station": k, "station_m": float(st.station_m.iat[k]), "stage_m": stage, "level_m": lev,
                     "neg_m": off.min() - SURVEY_MARGIN_M if off.min() < 0 else np.nan, "pos_m": off.max() + SURVEY_MARGIN_M if off.max() > 0 else np.nan})
    return pd.DataFrame(rows).sort_values("station").reset_index(drop=True)


def along_interp(st, an, col):
    a = an.dropna(subset=[col]).groupby("station")[col].mean()
    return np.interp(np.arange(len(st)), a.index.to_numpy(float), a.to_numpy(float))


# ---------------------------------------------------------------- the three sets ---
def set_A(reach):
    """The old lines, as written by 44_banklines_survey.py (the main branch's dataSorting.ipynb at commit 2346281)."""
    return gpd.read_file(SURVEY_LINES / f"bank_lines_{reach}.geojson").to_crs(CRS_FT)


def offsets_to_lines(reach, st, neg, pos, old):
    x, y, nx, ny = st.x.to_numpy(), st.y.to_numpy(), st.nx.to_numpy(), st.ny.to_numpy()
    lneg = LineString(np.column_stack([x + neg / FT * nx, y + neg / FT * ny]))
    lpos = LineString(np.column_stack([x + pos / FT * nx, y + pos / FT * ny]))
    left_is_pos = side_map(reach, st, old)
    left, right = (lpos, lneg) if left_is_pos else (lneg, lpos)
    return gpd.GeoDataFrame({"section": [reach] * 2, "bank": ["left", "right"], "stretch": [1, 1]}, geometry=[left, right], crs=CRS_FT)


def side_map(reach, st, old):
    """True if the old file's 'left' bank is on the stations' pos side."""
    ol = old[old.bank == "left"].geometry.iloc[0]; m = ol.interpolate(0.5, normalized=True)
    i = int(np.argmin(np.hypot(st.x - m.x, st.y - m.y)))
    return ((m.x - st.x.iat[i]) * st.nx.iat[i] + (m.y - st.y.iat[i]) * st.ny.iat[i]) > 0


def push_out(st, neg, pos, reach, p=None):
    """Move edges outward wherever a surveyed point (all, or those given) would sit outside them (1 m margin)."""
    p = survey_points(reach) if p is None else p; cl = LineString(st[["x", "y"]].to_numpy()); cum = np.r_[0, np.cumsum(np.hypot(np.diff(st.x), np.diff(st.y)))]
    for e, n in zip(p.Easting, p.Northing):
        i = int(np.clip(np.searchsorted(cum, cl.project(Point(e, n))) - 1, 0, len(st) - 1))
        o = ((e - st.x.iat[i]) * st.nx.iat[i] + (n - st.y.iat[i]) * st.ny.iat[i]) * FT
        sl = slice(max(i - 1, 0), i + 2)
        neg[sl] = np.minimum(neg[sl], o - SURVEY_MARGIN_M); pos[sl] = np.maximum(pos[sl], o + SURVEY_MARGIN_M)
    return neg, pos


def set_B(reach, st, e, an):
    """3DEP edges, shifted along the reach so that they pass through the surveyed edges, then pushed out to the points."""
    out = {}
    for side in ("neg", "pos"):
        lid = e[f"edge_{side}_m_s"].to_numpy()
        a = an.dropna(subset=[f"{side}_m"]).groupby("station")[f"{side}_m"].mean()
        delta = a.to_numpy() - lid[a.index.to_numpy()]
        out[side] = lid + np.interp(np.arange(len(st)), a.index.to_numpy(float), delta)
    neg, pos = push_out(st, out["neg"].copy(), out["pos"].copy(), reach)
    return neg, pos


class StageBanks:
    """Set C: the edge at any water level, with its uncertainty, per station and side."""

    def __init__(self, reach, t, st, e, an):
        self.reach, self.st, self.e = reach, st, e
        self.prof = {i: (g.sort_values("offset_m").offset_m.to_numpy(float), g.sort_values("offset_m").elev_m.to_numpy(float)) for i, g in t.groupby("transect")}
        self.W = e.water_m_s.to_numpy()
        self.es = {s: along_interp(st, an, f"{s}_m") for s in ("neg", "pos")}
        self.hs = along_interp(st, an, "level_m")
        self.eW = {s: e[f"edge_{s}_m_s"].to_numpy() for s in ("neg", "pos")}
        survey_w = self.es["pos"] - self.es["neg"]; lidar_w = self.eW["pos"] - self.eW["neg"]
        self.valid = e.flat_ok.to_numpy() & (lidar_w >= CANOPY_RATIO * survey_w)

    def _walk(self, i, side, h):
        o, z = self.prof[i]; start = self.eW[side][i]; sgn = -1 if side == "neg" else 1
        idx = np.where((sgn * o >= sgn * start) & (z > h))[0]
        if not len(idx):
            return o[0] if side == "neg" else o[-1], EDGE_SD_MAX_M
        j = idx[np.argmin(sgn * o[idx])]
        jj = j - sgn if 0 <= j - sgn < len(o) else j
        slope = abs(z[j] - z[jj]) / max(abs(o[j] - o[jj]), 1e-6)
        return o[j] - sgn * 0.5, min(LIDAR_SD_M / max(slope, 1e-3), EDGE_SD_MAX_M)

    def edges(self, h):
        """Best-estimate offsets and sds at water level h (array over stations)."""
        h = np.broadcast_to(np.asarray(h, float), (len(self.st),))
        out = {}
        for side in ("neg", "pos"):
            sgn = -1 if side == "neg" else 1; E = np.empty(len(self.st)); S = np.empty(len(self.st))
            for i in range(len(self.st)):
                es, hs, eW, W = self.es[side][i], self.hs[i], self.eW[side][i], self.W[i]
                if not self.valid[i]:
                    E[i], S[i] = es, SUBSTITUTE_SD_M
                elif h[i] >= W:
                    E[i], S[i] = self._walk(i, side, h[i])
                else:
                    if sgn * eW <= sgn * es or abs(W - hs) < 0.02:      # lidar edge not outside the surveyed edge: no information on the slope
                        E[i], S[i] = es, SUBSTITUTE_SD_M
                    else:
                        dedh = (eW - es) / (W - hs); E[i] = es + dedh * (h[i] - hs)
                        S[i] = min(np.hypot(LIDAR_SD_M * abs(dedh), GAP_FRAC * abs(eW - es)), EDGE_SD_MAX_M)
                if h[i] >= hs - 0.02:                                   # at or above the survey level the surveyed edge is inside
                    E[i] = es if sgn * es > sgn * E[i] else E[i]
                E[i] = sgn * max(sgn * E[i], 0.5)                       # never across the centreline
            E = pd.Series(E).rolling(EDGE_SMOOTH, center=True, min_periods=1).median().to_numpy()
            out[side] = (E, S)
        return out

    def draw(self, h, z_neg, z_pos):
        ed = self.edges(h)
        neg = ed["neg"][0] + z_neg * ed["neg"][1]; pos = ed["pos"][0] + z_pos * ed["pos"][1]
        return np.minimum(neg, -0.5), np.maximum(pos, 0.5)


def correlated_normals(n, rng, step_m=10.0):
    r = np.exp(-step_m / CORR_M); z = np.empty(n); z[0] = rng.standard_normal()
    for i in range(1, n):
        z[i] = r * z[i - 1] + np.sqrt(1 - r * r) * rng.standard_normal()
    return z


# ---------------------------------------------------------------- features from bank lines ---
def bank_features(pts, banks):
    """dist_left/right, dist_bank_m, width_m, rel_across (-1 right bank, +1 left bank), as in 19_ and 27_."""
    left = banks[banks.bank == "left"].geometry.iloc[0]; right = banks[banks.bank == "right"].geometry.iloc[0]
    P = [Point(e, n) for e, n in zip(pts.Easting, pts.Northing)]
    dl = np.array([left.distance(p) for p in P]) * FT; dr = np.array([right.distance(p) for p in P]) * FT
    return pd.DataFrame({"dist_bank_m": np.minimum(dl, dr), "width_m": dl + dr, "rel_across": (dr - dl) / (dl + dr)}, index=pts.index)


def polygon(banks):
    left = banks[banks.bank == "left"].geometry.iloc[0]; right = banks[banks.bank == "right"].geometry.iloc[0]
    return Polygon(list(left.coords) + list(right.coords)[::-1]).buffer(0)


# ---------------------------------------------------------------- 3DEP fetch (only if the cache is deleted) ---
EPQS = "https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/getSamples"


def sample_3dep(pts, batch=500):
    """3DEP elevations at pts[lon, lat] (needs internet). The cached tables in data_raw were made with this call."""
    import requests
    out = np.full(len(pts), np.nan)
    for k in range(int(np.ceil(len(pts) / batch))):
        chunk = pts.iloc[k * batch:(k + 1) * batch]
        params = {"geometry": json.dumps({"points": chunk[["lon", "lat"]].values.tolist(), "spatialReference": {"wkid": 4326}}),
                  "geometryType": "esriGeometryMultipoint", "returnFirstValueOnly": "true", "interpolation": "RSP_NearestNeighbor", "outFields": "", "f": "json"}
        for attempt in range(6):
            try:
                r = requests.post(EPQS, data=params, timeout=180); r.raise_for_status(); body = r.json()
                if "error" in body: raise RuntimeError(body["error"])
                break
            except Exception:
                if attempt == 5: raise
                time.sleep(15 * (attempt + 1))
        for s in body.get("samples", []):
            try: v = float(s["value"])
            except (TypeError, ValueError): v = np.nan
            out[k * batch + int(s["locationId"])] = v if v > -1000 else np.nan
    return out


# ---------------------------------------------------------------- set C in the model (logbook entry 30) ---
LEVEL_STEP_M = 0.05     # the chain's stage grid for the bank features (decision in entry 30)
C_SEED = 30


class StageFeatures:
    """Set C as the model uses it: bank lines, distance to the bank and relative position at any stage, best estimate or draw j.

    stage is in the chain's own unit: the sensor stage for the Middle, the sections' level L at the reference position (m NAVD88)
    for the Upper and Lower. At a stage, the edges are C's edges at the water level that stage gives at each station, moved out so that
    every point surveyed on a day with the same or a lower water level is inside (1 m margin); the same rule makes the training
    features (each point at its own survey day's level). Draw j: the edge errors (lidar, gap, substitute) as correlated normals fixed
    per draw, and for the survey days also the day's tie error; away from survey days the tie error enters through the chain's stage
    source, which moves the stage at which the features are read."""

    def __init__(self, reach, n_draws=N_DRAWS, seed=C_SEED):
        self.reach, self.n_draws = reach, n_draws
        t = load_transects(reach); self.st = stations(reach, t); e = flat_edges(t); an = survey_anchors(reach, self.st)
        self.C = StageBanks(reach, t, self.st, e, an); self.old = set_A(reach); self.k = tie(reach)
        self.along = self.st.along_m.to_numpy(); self.p = survey_points(reach)
        sday = self.p.groupby("date").stage_m.median()
        self.day_stage = sday.map(self.from_sensor).to_dict() if reach != "Middle" else sday.to_dict()   # in the chain's unit
        rng = np.random.default_rng(seed); n = len(self.st)
        self.zn = [correlated_normals(n, rng) for _ in range(n_draws)]; self.zp = [correlated_normals(n, rng) for _ in range(n_draws)]
        self.ztie = [{d: float(rng.standard_normal()) for d in sorted(self.day_stage)} for _ in range(n_draws)]
        self._cache = {}

    def from_sensor(self, s):
        return s if self.reach == "Middle" else self.k["a"] + self.k["b"] * s

    def level(self, stage):
        if self.reach == "Middle":
            return water_level(self.reach, stage, self.along)
        return stage - self.k["slope"] * (self.along - self.k["along_ref"])

    def offsets(self, stage, j=None, day=None):
        key = (round(float(stage), 4), j, day)
        if key in self._cache:
            return self._cache[key]
        lev = self.level(stage) + (self.ztie[j][day] * self.k["sd"] if (j is not None and day is not None) else 0.0)
        if j is None:
            ed = self.C.edges(lev); neg, pos = ed["neg"][0].copy(), ed["pos"][0].copy()
        else:
            neg, pos = self.C.draw(lev, self.zn[j], self.zp[j])
        lower = [d for d, s in self.day_stage.items() if s <= float(stage) + 1e-9]
        neg, pos = push_out(self.st, neg, pos, self.reach, self.p[self.p.date.isin(lower)])
        self._cache[key] = (neg, pos); return neg, pos

    def lines(self, stage, j=None, day=None):
        neg, pos = self.offsets(stage, j, day); return offsets_to_lines(self.reach, self.st, neg, pos, self.old)

    def features(self, E, N, stage, j=None, day=None):
        """dist_bank_m, rel_across, width_m at positions (ft, state plane) for one stage."""
        import shapely
        bk = self.lines(stage, j, day); left = bk[bk.bank == "left"].geometry.iloc[0]; right = bk[bk.bank == "right"].geometry.iloc[0]
        P = shapely.points(np.asarray(E, float), np.asarray(N, float)); dl = shapely.distance(P, left) * FT; dr = shapely.distance(P, right) * FT
        return np.minimum(dl, dr), (dr - dl) / (dl + dr), dl + dr

    def training(self, df, j=None):
        """df with the bank features at each point's own survey day (C best, or draw j)."""
        d = df.copy(); dist, rel, wid = np.empty(len(d)), np.empty(len(d)), np.empty(len(d))
        for day, ix in d.groupby("date").indices.items():
            r = self.features(d.Easting.values[ix], d.Northing.values[ix], self.day_stage[day], j, day); dist[ix], rel[ix], wid[ix] = r
        d["dist_bank_m"], d["rel_across"], d["width_m"] = dist, rel, wid; return d

    def grid(self, lo, hi):
        """The stage grid covering [lo, hi] at LEVEL_STEP_M."""
        return np.round(np.arange(np.floor(lo / LEVEL_STEP_M) * LEVEL_STEP_M, np.ceil(hi / LEVEL_STEP_M) * LEVEL_STEP_M + 1e-9, LEVEL_STEP_M), 4)

    def top_survey_stage(self):
        """The highest survey-day stage: the map's outline is C there, so the grid stays inside what the survey supports (entry 30)."""
        return float(max(self.day_stage.values()))

    def outline(self, stage):
        """The channel polygon of the best estimate at a stage (the map uses the highest survey-day stage)."""
        return polygon(self.lines(stage))

    def write(self, out_dir):
        """The lines on disk: C at each survey day and the outline at the highest survey-day stage."""
        out_dir.mkdir(parents=True, exist_ok=True)
        for day, s in self.day_stage.items():
            self.lines(s).to_file(out_dir / f"C_{self.reach}_{day}.geojson", driver="GeoJSON")
        top = self.top_survey_stage()
        gpd.GeoDataFrame({"reach": [self.reach], "stage": [top]}, geometry=[self.outline(top)], crs=CRS_FT).to_file(out_dir / f"C_outline_{self.reach}.geojson", driver="GeoJSON")


if __name__ == "__main__":
    import sys
    if len(sys.argv) >= 3 and sys.argv[1] == "fetch":
        fetch_transects(sys.argv[2])
    elif len(sys.argv) >= 2 and sys.argv[1] == "check":
        for r in ("Upper", "Middle", "Lower"):
            print(r, "largest distance from transect_points (m):", check_transects(r))
