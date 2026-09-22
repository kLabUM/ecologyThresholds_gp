"""Transect-label check. A transect is a cross-section, so its points should share one position along the reach.
For every (date, transect) group: the spread of along_m, points more than 10 m from the group's median, points that sit
closer to another transect's line than to their own, point names whose prefix does not match the transect, and
coordinates duplicated across dates. Output: data_clean/transect_id_check.csv and a printed summary."""
import re
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"
loc = pd.read_csv(D / "locations.csv")
grp = loc.groupby(["date", "transect"]).along_m.agg(["median", "min", "max", "size"]).reset_index(); grp["spread_m"] = grp["max"] - grp["min"]
print("along-reach spread of each (date, transect) group, m:\n", grp.sort_values("spread_m", ascending=False).head(8).round(1).to_string(index=False))
loc["group_median_along"] = loc.groupby(["date", "transect"]).along_m.transform("median"); loc["off_own_line_m"] = (loc.along_m - loc.group_median_along).abs()
lines = grp.set_index("transect")["median"]
def nearest_other(r):
    others = lines.drop(r.transect, errors="ignore"); d = (others - r.along_m).abs(); return others.index[d.argmin()], float(d.min())
no = loc.apply(nearest_other, axis=1, result_type="expand"); loc["nearest_other_transect"], loc["dist_to_other_line_m"] = no[0], no[1]
loc["name_prefix_mismatch"] = [(str(p)[:4] != str(t)[:4]) if re.match(r"^AT\d{4}$", str(p)) else False for p, t in zip(loc.point, loc.transect)]
loc["coords_duplicated_across_dates"] = loc.duplicated(["Easting", "Northing"], keep=False)
sus = loc[(loc.off_own_line_m > 10) | loc.name_prefix_mismatch | loc.coords_duplicated_across_dates].copy()
sus["verdict"] = np.where((sus.off_own_line_m > 10) & (sus.dist_to_other_line_m < 5), "on another transect's line: " + sus.nearest_other_transect, np.where(sus.name_prefix_mismatch, "name says another transect", "coordinates seen on another date"))
cols = ["date", "transect", "point", "along_m", "group_median_along", "off_own_line_m", "nearest_other_transect", "dist_to_other_line_m", "name_prefix_mismatch", "coords_duplicated_across_dates", "verdict"]
sus[cols].round(1).to_csv(D / "transect_id_check.csv", index=False); print(f"\nsuspect locations: {len(sus)}\n", sus[cols].round(1).to_string(index=False))
