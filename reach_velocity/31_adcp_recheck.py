"""The ADCP re-read after the reviewer pass (22 September 2026). The instrument sits under about 0.8 m of water at 0.75 m of stage (the
nearest surveyed profile, 1.7 m away, has its bed 0.04 m below the sensor datum), not the 0.2 m the fitted h0 of 0.55 m implied. So its
record is a near-bed value: at 0.73 to 0.77 m of stage it reads what the neighbouring log profile gives 3 cm above the bed, a third of
that profile's column mean. Its stage response is then re-fitted with the surveyed depth under both readings. Output: adcp_recheck.json"""
import json, numpy as np, pandas as pd
from gpcore import *
Z0 = json.load(open(D / "field_qc.json"))["z0_m"]; H0 = 0.45
ad = pd.read_csv(D / "adcp.csv"); loc = pd.read_csv(D / "locations.csv"); near = loc.assign(d=np.hypot(loc.x_m, loc.y_m)).nsmallest(3, "d")
bed_adcp = float(near.bed_rtk_sensor_datum_m.iloc[0]); x = np.log(np.clip(ad.stage_m - H0, 0.02, None)); v = ad.velocity_ms.values; ok = v > 0.005
h = ad.stage_m - bed_adcp; a_col = v / (np.log(np.clip(h, np.e * Z0, None) / Z0) - 1)
band = (ad.stage_m > 0.73) & (ad.stage_m < 0.77); v_band = float(v[band].mean()); a_near = float(near.a.iloc[0]); h_near = float(near.water_depth_m.iloc[0])
out = {"nearest_profile": near.location_id.iloc[0], "distance_m": float(near.d.iloc[0]), "bed_at_adcp_sensor_datum_m": bed_adcp, "water_depth_over_adcp_at_0.75_m": 0.75 - bed_adcp, "h0_previously_fitted_m": 0.55,
       "adcp_mean_at_0.73_to_0.77_m_s": v_band, "neighbour_profile_amplitude_a": a_near, "neighbour_column_mean_m_s": a_near * (np.log(h_near / Z0) - 1), "log_law_height_where_v_equals_adcp_m": float(Z0 * np.exp(v_band / a_near)),
       "adcp_exponent_on_stage_minus_0.45": {"point reading (v proportional to a)": float(np.polyfit(x[ok], np.log(v[ok]), 1)[0]), "column reading with the surveyed depth": float(np.polyfit(x[ok], np.log(a_col[ok]), 1)[0])},
       "field_A_at_adcp_over_record": {"note": "the field's A at the ADCP position against the record's implied A under the point reading at 3 cm; the ratio is what an instrument model has to explain"},
       "reading": "Vb is a near-bed quantity (a bottom cell, or one slanted beam's component, which would be about 0.42 of the along-stream speed); the SonTek model and configuration are still to be confirmed by the author"}
json.dump(out, open(D / "adcp_recheck.json", "w"), indent=1); print(json.dumps(out, indent=1))
