# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "openpyxl>=3.1",
#   "pandas>=2.0",
#   "requests>=2.31",
# ]
# ///
"""
Fetch daily precipitation (mm) for each main-stem site from Open-Meteo's
archive API, using site coordinates from siteLatLong.xlsx. Writes a long-form
CSV that ridgeline_publication.py reads to draw the rainfall bars.

Run once per date-range change:
  uv run --script fetch_rainfall.py
  (or:  pip install openpyxl pandas requests; python fetch_rainfall.py)

Inputs:
  SITE_XLSX  — Excel with columns Site, Lat, Long (one row per site).
  START/STOP — must match ridgeline_publication.py.

Output:
  rainfall_by_site.csv  columns: site_id, date, precip_mm
"""

from __future__ import annotations

import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

# ────────────────────────────────────────────────────────────────────────────
# CONFIG — keep in sync with ridgeline_publication.py
# ────────────────────────────────────────────────────────────────────────────

SITE_XLSX = Path("siteLatLong.xlsx")            # adjust if you keep it elsewhere
OUT_CSV = Path("rainfall_by_site.csv")

MAIN_STEM: set[str] = {
    "ARB034", "ARB026", "ARB013", "ARB012", "ARB027", "ARB025", "ARB031",
    "ARB059", "ARB109", "ARB077", "ARB112", "ARB118", "ARB120",
}

START: datetime = datetime(2026, 4, 1, tzinfo=timezone.utc)
STOP: datetime = datetime(2026, 9, 1, tzinfo=timezone.utc)

# Open-Meteo archive API is free, no key required.
# Docs: https://open-meteo.com/en/docs/historical-weather-api
ENDPOINT = "https://archive-api.open-meteo.com/v1/archive"
REQUEST_TIMEOUT_S = 30
PAUSE_BETWEEN_S = 0.15          # be polite to the free API

# ────────────────────────────────────────────────────────────────────────────

def load_sites() -> pd.DataFrame:
    if not SITE_XLSX.exists():
        raise FileNotFoundError(
            f"Could not find {SITE_XLSX.resolve()}. Put the xlsx next to this "
            f"script, or edit SITE_XLSX at the top of the file."
        )
    df = pd.read_excel(SITE_XLSX)
    # Normalize column names, tolerating variants like 'Longitude'.
    lower = {c.lower().strip(): c for c in df.columns}
    def pick(*aliases):
        for a in aliases:
            if a in lower: return lower[a]
        raise KeyError(f"Missing a column matching any of {aliases!r}. Got {list(df.columns)!r}.")
    site_col = pick("site", "site_id", "id")
    lat_col  = pick("lat", "latitude")
    lon_col  = pick("long", "lon", "lng", "longitude")
    df = df[[site_col, lat_col, lon_col]].rename(
        columns={site_col: "site_id", lat_col: "lat", lon_col: "lon"}
    )
    df["site_id"] = df["site_id"].astype(str).str.strip()
    df = df.dropna(subset=["site_id", "lat", "lon"]).drop_duplicates(subset=["site_id"])
    return df

def fetch_one(site_id: str, lat: float, lon: float) -> pd.DataFrame:
    params = {
        "latitude": lat,
        "longitude": lon,
        "start_date": START.date().isoformat(),
        "end_date":   STOP.date().isoformat(),
        "daily": "precipitation_sum",
        "timezone": "UTC",
    }
    r = requests.get(ENDPOINT, params=params, timeout=REQUEST_TIMEOUT_S)
    r.raise_for_status()
    j = r.json()
    daily = j.get("daily") or {}
    dates  = daily.get("time") or []
    precip = daily.get("precipitation_sum") or []
    if not dates:
        return pd.DataFrame(columns=["site_id", "date", "precip_mm"])
    return pd.DataFrame({
        "site_id": site_id,
        "date": pd.to_datetime(dates, utc=True),
        "precip_mm": pd.to_numeric(pd.Series(precip), errors="coerce"),
    })

def main() -> None:
    sites = load_sites()
    wanted = sites[sites["site_id"].isin(MAIN_STEM)].copy()
    missing = MAIN_STEM - set(wanted["site_id"])
    if missing:
        print(f"WARNING — no coordinates in xlsx for: {sorted(missing)}", file=sys.stderr)
    if wanted.empty:
        raise RuntimeError("None of the main-stem sites were found in the xlsx.")

    frames = []
    for i, row in enumerate(wanted.itertuples(index=False), 1):
        sid = row.site_id
        print(f"[{i:>2}/{len(wanted)}] {sid}  ({row.lat:.5f}, {row.lon:.5f})", flush=True)
        try:
            frames.append(fetch_one(sid, float(row.lat), float(row.lon)))
        except requests.HTTPError as e:
            print(f"    HTTP {e.response.status_code}: {e.response.text[:200]}", file=sys.stderr)
        except Exception as e:
            print(f"    failed: {e}", file=sys.stderr)
        time.sleep(PAUSE_BETWEEN_S)

    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
        columns=["site_id", "date", "precip_mm"]
    )
    # Write date as ISO YYYY-MM-DD for easy parsing downstream.
    out["date"] = pd.to_datetime(out["date"]).dt.strftime("%Y-%m-%d")
    out.to_csv(OUT_CSV, index=False)
    print(
        f"\nSaved → {OUT_CSV.resolve()}   "
        f"({len(out):,} rows across {out['site_id'].nunique()} sites)"
    )

if __name__ == "__main__":
    main()
