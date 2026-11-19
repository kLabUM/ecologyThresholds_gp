"""Step 0a. Load the raw files into clean tables, with every choice written down.
Inputs (data_raw/, copies of the repo data/ files): 25_ecohydroData.csv, 26_ecohydroData.csv (hand-held surveys), 25_ADCPData.csv (SonTek at one point),
25_/26_depthSensorData_Middle.csv (MaxBotix ultrasonic stage). Outputs (data_clean/): readings_transect.csv, nests.csv,
adcp.csv, stage_10min.csv, usgs_daily.csv, summary.json. Middle reach only."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[1]; HERE = Path(__file__).resolve().parent; OUT = HERE / "data_clean"; OUT.mkdir(exist_ok=True)
RAW = HERE / "data_raw" if (HERE / "data_raw").exists() else ROOT / "data"   # a copy of the raw files lives next to the logbook
import provenance; provenance.check()   # every line or polygon file in the folder has a script that writes it (entry 30); stops here otherwise
FT = 0.3048
NORTHING_BAND = (280000, 323000)          # the field lead's rule for 'Middle' (State Plane MI South, ft)
STAGE_DATUM_MM = 7313.2                   # sensor mount height above the bed datum, mm (offset tag in the 2025/2026 files)
STAGE_PLAUSIBLE_M = (0.3, 3.0)            # outside this the reading is a sensor fault (obstruction episode Aug 2025 read 6 m)
ADCP_E, ADCP_N = 13244857.0, 319995.5     # SonTek location (from dataSorting.ipynb)
TZ = "America/Detroit"
notes = {}

# ---------------------------------------------------------------- hand-held surveys
frames = []
for yr, f in [(2025, "25_ecohydroData.csv"), (2026, "26_ecohydroData.csv")]:
    r = pd.read_csv(RAW / f, dtype=str); r["year"] = yr; frames.append(r)
s = pd.concat(frames, ignore_index=True)
for c in ["Total water depth", "Depth", "Actual depth", "Velocity", "Northing", "Easting", "Latitude", "Longitude"]:
    s[c] = pd.to_numeric(s.get(c), errors="coerce")
for c in ["Transect ID", "Measurement Point ID", "Type", "Section"]: s[c] = s[c].astype(str).str.strip().replace({"nan": np.nan})   # ' AT11' carried a leading space
# named transect-label corrections, confirmed by the author (the field lead, 21 Sep 2026): the crew sometimes left the previous day's
# transect ID in the form. These points carry AT07's names and coordinates and sit on AT07's line, 26 m from the rest of AT06.
RELABEL = {("7/14/2025", "AT0710"): "AT07", ("7/14/2025", "AT0711"): "AT07", ("7/14/2025", "AT0714"): "AT07"}
for (d, pt), t in RELABEL.items(): s.loc[(s["Collection Date"].str.strip() == d) & (s["Measurement Point ID"] == pt), "Transect ID"] = t
s["time_local"] = pd.to_datetime(s["Collection Date"].str.strip() + " " + s["Collection Time"].fillna("12:00").str.strip(), format="mixed", errors="coerce")
s["time_utc"] = s["time_local"].dt.tz_localize(TZ, ambiguous="NaT", nonexistent="NaT").dt.tz_convert("UTC")
s["date"] = s["time_local"].dt.date.astype(str)
notes["surveys_raw_rows"] = {"2025": int((s.year == 2025).sum()), "2026": int((s.year == 2026).sum())}
mid = s[(s.Northing > NORTHING_BAND[0]) & (s.Northing < NORTHING_BAND[1])].copy()
notes["middle_rows_by_type"] = mid.Type.value_counts(dropna=False).to_dict()
# height above the bed: 2025 'Depth' is a rod reading, 'Actual depth' (= 0.4 x Depth) is the height in ft; 2026 'Depth' is the height in ft
mid["height_ft"] = np.where(mid["Actual depth"].notna(), mid["Actual depth"], mid["Depth"])
mid["height_m"] = mid.height_ft * FT; mid["water_depth_m"] = mid["Total water depth"] * FT; mid["velocity_ms"] = mid.Velocity * FT
mid["x_m"] = (mid.Easting - ADCP_E) * FT; mid["y_m"] = (mid.Northing - ADCP_N) * FT      # local frame in metres, origin at the ADCP
mid["location_id"] = mid.date + "|" + mid["Transect ID"].astype(str) + "|" + mid["Measurement Point ID"].astype(str)
keep = ["year", "date", "time_local", "time_utc", "Type", "Transect ID", "Measurement Point ID", "location_id", "Easting", "Northing", "Latitude", "Longitude", "x_m", "y_m",
        "water_depth_m", "height_ft", "height_m", "velocity_ms", "Point Notes", "SMB Presence", "Investigators"]
tr = mid[mid.Type == "Transect"][keep].copy()
n0 = len(tr); tr = tr[tr.velocity_ms.notna() & tr.height_m.notna()]; n1 = len(tr); tr = tr[tr.velocity_ms >= 0]; n2 = len(tr)
tr = tr[~(tr.height_m > tr.water_depth_m + 0.02)]; n3 = len(tr)
# named entry errors (each one listed with its reason, so the exclusion is visible)
ENTRY_ERRORS = {("2025-07-14", "AT06", "AT0605", 0.08): "3.60 ft/s at 0.08 ft while the other 11 readings of the profile are 0.5 to 0.8 ft/s; a stray digit"}
bad = pd.Series([(d, t, m, round(h, 2)) in ENTRY_ERRORS for d, t, m, h in zip(tr.date, tr["Transect ID"], tr["Measurement Point ID"], tr.height_ft)], index=tr.index)
tr = tr[~bad]; n4 = len(tr)
notes["transect_filter"] = {"rows": n0, "with velocity and height": n1, "velocity >= 0": n2, "height not above the water surface": n3, "named entry errors removed": n3 - n4}
tr = tr.sort_values(["time_utc", "location_id", "height_m"]).reset_index(drop=True); tr.to_csv(OUT / "readings_transect.csv", index=False)
nests = mid[mid.Type == "Nest"][keep].copy(); nests.to_csv(OUT / "nests.csv", index=False)
notes["dropped_types"] = {k: int(v) for k, v in mid.Type.value_counts().items() if k not in ("Transect", "Nest")}

# ---------------------------------------------------------------- stage sensor (10-min ultrasonic range to the water surface)
st = []
for f in ["25_depthSensorData_Middle.csv", "26_depthSensorData_Middle.csv"]:
    d = pd.read_csv(RAW / f, skiprows=3); d = d[d["_field"] == "value"]; st.append(pd.DataFrame({"time_utc": pd.to_datetime(d["_time"], utc=True), "range_mm": d["_value"].astype(float), "offset_tag_mm": d["offset"].astype(float)}))
st = pd.concat(st).sort_values("time_utc").drop_duplicates("time_utc").reset_index(drop=True)
st["stage_m"] = (STAGE_DATUM_MM - st.range_mm) / 1000.0
ok = st.stage_m.between(*STAGE_PLAUSIBLE_M); notes["stage_rows"] = {"total": int(len(st)), "implausible_removed": int((~ok).sum())}
st = st[ok].set_index("time_utc")
med = st.stage_m.rolling("6h", center=True).median(); spike = (st.stage_m - med).abs() > 0.15
notes["stage_rows"]["spikes_removed"] = int(spike.sum()); st = st[~spike]
st10 = st.stage_m.resample("10min").mean().to_frame("stage_m"); st10["stage_3h_smooth"] = st10.stage_m.rolling("3h", center=True).mean()
st10.to_csv(OUT / "stage_10min.csv")
notes["stage_period"] = {"2025": [str(st10[st10.index.year == 2025].dropna().index.min().date()), str(st10[st10.index.year == 2025].dropna().index.max().date())],
                         "2026": [str(st10[st10.index.year == 2026].dropna().index.min().date()), str(st10[st10.index.year == 2026].dropna().index.max().date())]}

# stage at the time of each reading, interpolated from the corrected 10-min grid (series_10min.csv from 06_align_10min.py:
# daily sensor pattern removed against the night line, median filtered) when it exists, else from the raw 3 h smooth; both
# are kept as columns. Run order: 00 -> 04 -> 05 -> 06 -> 00 again (this block) -> 01. Flag readings with no stage within 2 h.
sm_raw = st10.stage_3h_smooth.dropna(); grid_path = OUT / "series_10min.csv"
if grid_path.exists():   # 06_align_10min.py: stage corrected in 05 (median filter, daily wobble removed against the night line, 3 h median) plus a 2 h median; ADCP on the same grid, same filters
    grid = pd.read_csv(grid_path, parse_dates=["time_utc"]).set_index("time_utc"); grid.index = pd.to_datetime(grid.index, utc=True); sm = grid.stage_m.dropna().sort_index(); notes["stage_used_for_readings"] = "stage_m from series_10min.csv (05 + 06)"
else: grid = None; sm = sm_raw; notes["stage_used_for_readings"] = "raw 3 h smooth (run 04, 05, 06, then 00_load.py again)"
def stage_at(t, sm=sm):
    if pd.isna(t): return np.nan, np.nan
    i = sm.index.searchsorted(t); lo, hi = max(i - 1, 0), min(i, len(sm) - 1)
    gap = min(abs((sm.index[lo] - t).total_seconds()), abs((sm.index[hi] - t).total_seconds())) / 3600
    if gap > 2: return np.nan, gap
    w = 0.5 if sm.index[hi] == sm.index[lo] else (t - sm.index[lo]) / (sm.index[hi] - sm.index[lo]); return float(sm.iloc[lo] * (1 - w) + sm.iloc[hi] * w), gap
for df in (tr, nests):
    vals = [stage_at(t) for t in df.time_utc]; df["stage_m"] = [v[0] for v in vals]; df["stage_gap_h"] = [v[1] for v in vals]
    df["stage_raw_m"] = [stage_at(t, sm_raw)[0] for t in df.time_utc]
tr.to_csv(OUT / "readings_transect.csv", index=False); nests.to_csv(OUT / "nests.csv", index=False)
notes["readings_without_stage"] = int(tr.stage_m.isna().sum())

# ---------------------------------------------------------------- ADCP (SonTek beam velocity, m/s, one location)
if grid is not None:   # smoothed ADCP on the common 10-min grid (06_align_10min.py), with the stage from the same grid row
    g_ = grid.dropna(subset=["adcp_ms"]); ad = pd.DataFrame({"time_utc": g_.index, "velocity_ms": g_.adcp_ms.values, "velocity_raw_ms": g_.adcp_raw_ms.values, "stage_m": g_.stage_m.values, "stage_raw_m": g_.stage_raw_m.values}).reset_index(drop=True)
    notes["adcp_source"] = "series_10min.csv (smoothed, 10-min grid)"
else:
    a = pd.read_csv(RAW / "25_ADCPData.csv", skiprows=3); a = a[a["_field"] == "value"]
    ad = pd.DataFrame({"time_utc": pd.to_datetime(a["_time"], utc=True), "velocity_ms": a["_value"].astype(float)}).sort_values("time_utc").drop_duplicates("time_utc").reset_index(drop=True)
    vals = [stage_at(t) for t in ad.time_utc]; ad["stage_m"] = [v[0] for v in vals]; ad["stage_raw_m"] = [stage_at(t, sm_raw)[0] for t in ad.time_utc]; notes["adcp_source"] = "raw file"
ad = ad.dropna(subset=["stage_m"]).reset_index(drop=True)
ad["x_m"], ad["y_m"] = 0.0, 0.0; ad.to_csv(OUT / "adcp.csv", index=False)
notes["adcp"] = {"rows": int(len(ad)), "period": [str(ad.time_utc.min().date()), str(ad.time_utc.max().date())], "median_step_min": float(ad.time_utc.diff().dt.total_seconds().median() / 60), "velocity_range_ms": [float(ad.velocity_ms.min()), float(ad.velocity_ms.max())]}

# ---------------------------------------------------------------- USGS daily discharge at Ann Arbor (04174500), if the download is present
for site, out_name in [("04174500", "usgs_daily.csv"), ("04172000", "usgs_04172000_daily.csv"), ("04174040", "usgs_04174040_daily.csv")]:   # Ann Arbor (reference for the Middle), Hamburg and Zeeb Road (other reaches); downloaded 21 Sep 2026 from waterservices.usgs.gov, parameter 00060, daily
    usgs = RAW / f"usgs_{site}_dv.json"
    if usgs.exists():
        j = json.load(open(usgs)); v = j["value"]["timeSeries"][0]["values"][0]["value"]
        q = pd.DataFrame({"date": [x["dateTime"][:10] for x in v], "discharge_cfs": [float(x["value"]) for x in v]}); q["discharge_m3s"] = q.discharge_cfs * 0.0283168; q.to_csv(OUT / out_name, index=False); notes[f"usgs_{site}"] = {"rows": int(len(q))}

# ---------------------------------------------------------------- summary
loc = tr.groupby("location_id").agg(year=("year", "first"), date=("date", "first"), n=("velocity_ms", "size"), water_depth_m=("water_depth_m", "first"), stage_m=("stage_m", "first"))
summ = {"notes": notes, "by_year": {}}
for yr, g in tr.groupby("year"):
    L = loc[loc.year == yr]
    summ["by_year"][str(yr)] = {"survey_days": int(g.date.nunique()), "transects": int(g["Transect ID"].nunique()), "locations": int(len(L)), "readings": int(len(g)),
                                "readings_per_location": [int(L.n.min()), float(L.n.median()), int(L.n.max())], "height_m_range": [float(g.height_m.min()), float(g.height_m.max())],
                                "water_depth_m_median": float(g.water_depth_m.median()), "velocity_ms_median_iqr": [float(g.velocity_ms.median()), float(g.velocity_ms.quantile(.25)), float(g.velocity_ms.quantile(.75))],
                                "stage_m_range_on_survey_days": [float(np.nanmin(L.stage_m)), float(np.nanmax(L.stage_m))], "days": sorted(g.date.unique().tolist())}
json.dump(summ, open(OUT / "summary.json", "w"), indent=2, default=str)
print(json.dumps(summ, indent=2, default=str))
