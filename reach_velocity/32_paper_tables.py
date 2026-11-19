"""The paper's tables, from the saved tables in data_clean/: the sections (Table 1), the parameters and settings of the chain (Table 2),
the validation scores (Table 3), the ablations (Table 4). Outputs paper/TABLES.md and paper/table{1..4}_*.csv."""
import json, numpy as np, pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"; P = HERE / "paper"
from uncertainty import PRIMARY, LAW_RIDGE; H0, GAM = PRIMARY   # the stage-law set, read from the saved law grid (entry 30)
loc = pd.read_csv(D / "locations.csv"); loc = loc[loc.a > 0]; fe = pd.read_csv(D / "features.csv")
# ---- Table 1: the sections
cn = pd.read_csv(D / "step6_continuity.csv"); cv = pd.read_csv(D / "step6_section_coverage.csv"); npts = loc.groupby("transect").agg(n=("a", "size"), along_m=("along_m", "median")).reset_index().rename(columns={"transect": "section"})
t1 = cn.merge(cv[["section", "mean_v_measured", "share_measured", "mean_v_med", "mean_v_lo", "mean_v_hi"]], on="section").merge(npts, on="section").sort_values("along_m")
t1 = t1[["section", "date", "along_m", "n", "width_m", "stage_m", "mean_v_measured", "share_measured", "Q_measured_m3s", "mean_v_med", "mean_v_lo", "mean_v_hi"]]
t1.columns = ["section", "survey day", "along (m)", "locations", "width (m)", "stage (m)", "mean velocity measured (m/s)", "share below 0.2 m/s measured", "discharge measured (m3/s)", "mean velocity, day held out (m/s)", "90% lo", "90% hi"]
t1.round(3).to_csv(P / "table1_sections.csv", index=False)
# ---- Table 2: parameters and settings
fq = json.load(open(D / "field_qc.json")); cal = json.load(open(D / "uncertainty_calibration.json")); s3 = json.load(open(D / "step3_summary.json")); sc = s3["scales_best_all_days"]; rtk = pd.read_csv(D / "rtk_session_offsets.csv")
sdy = np.clip(fe.a_sd / fe.A / 1.0, 0.02, 1.0)   # a_sd/a: the per-point noise in log (A = a/g, same relative sd)
law_set = ", ".join(f"({h0}, {g})" for h0, g in LAW_RIDGE)   # the rule applied to the saved law grid (uncertainty.law_set, entry 30)
rows = [
 ("Readings to one number", "vertical profile", "v(z) = a ln(z / z0); readings at z >= 3 cm and v > 5 mm/s", "fixed form"),
 ("", "roughness scale z0", f"{fq['z0_m']*1000:.1f} mm for the reach (flat between {fq['z0_within_5pct_of_min_m'][0]*1000:.0f} and {fq['z0_within_5pct_of_min_m'][1]*1000:.0f} mm), from {fq['profiles_in_z0_fit']} profiles with three or more readings", "fitted, pooled grid search"),
 ("", "amplitude a per location", f"{len(loc)} locations with a > 0; per-location sd from the spread of the readings' implied a (median {np.median(fe.a_sd / fe.A):.2f} relative)", "fitted per location"),
 ("Stage record", "sensor", "ultrasonic range every 10 min; daytime sensor pattern removed against the night line, median filtered; 2 h means for the series", "processing"),
 ("", "stage error in the chain", "1.5 cm, one draw per time shared by every section", "assumed"),
 ("Stage law", "form", "a = A (stage - h0)^gamma, one law for the reach", "fixed form"),
 ("", "primary pair", f"h0 = {H0} m, gamma = {GAM} (best day-out CRPS)", "chosen on a grid"),
 ("", "law set drawn in the chain", f"every grid pair within 0.01 CRPS of the best with h0 below 0.66 m: {law_set}", "uniform over the set"),
 ("Spatial model", "target", "log A, A = a / (stage - h0)^gamma", ""),
 ("", "inputs", "along the reach, relative position across, distance to the nearer bank, relative depth in the section, section mean depth, lateral bed slope", "measured"),
 ("", "mean", "weighted least squares on an intercept, six substrate classes from the notes and a wake flag", "fitted"),
 ("", "kernel", "ARD Matern 3/2 plus a white nugget; per-point noise fixed from the readings' spread", "fitted by marginal likelihood, 2 to 3 restarts"),
 ("", "fitted scales", f"amplitude sd {sc['amp_sd']:.2f}, nugget sd {sc['nugget_sd']:.2f} (log); length scales: along {sc['along_m']:.0f} m, relative across {sc['rel_across']:.1f}, bank distance {sc['dist_bank_m']:.0f} m, relative depth {sc['rel_depth']:.2f}, section mean depth {sc['section_mean_depth_m']:.2f} m, lateral slope {sc['lateral_bed_slope_abs']:.2f}", "fitted"),
 ("Uncertainty chain", "field", "one joint Gaussian process draw over every point of every section", "posterior"),
 ("", "discrepancy", f"inflation {cal['inflation']:.3f} from held-out days (95% coverage before inflation {cal['held_out_coverage95_before']:.3f}, n {cal['n']}); total {cal['disc_sd_log']:.3f} in log, split into a reach-level part {cal['discrepancy_split']['common_sd_log']:.3f} (one draw per reach) and a point-level part {cal['discrepancy_split']['independent_sd_log']:.3f}", "calibrated"),
 ("", "bed", f"RTK cross-sections on the sensor datum through per-session offsets (one water-surface slope of 5 cm per 100 m, residual 1.3 cm, {len(rtk)} sessions); 6 cm per point in the chain", "measured; error assumed"),
 ("", "roughness", "z0 drawn between 7 and 16 mm, one per draw; a rescaled through the readings' typical height (0.25 m)", "assumed range"),
 ("", "draws", "100 per two-hourly stage for the series (3,106 stages), 300 for the rating curves", "setting"),
 ("Map (illustration only)", "bed between sections", "Gaussian process on the RTK sections and bank points (scales 117 m along, 29 m across, nugget 5 cm); 2 m grid inside the mapped bank lines", "fitted"),
 ("Ecological number", "definition", "share of the wetted width below 0.2 m/s on the depth-averaged velocity (alternative: 0.10 m above the bed); width-weighted within a section, length-weighted over the 20 sections; wetted width ends at the outermost surveyed points", "choice, to be tied to the ecology"),
]
t2 = pd.DataFrame(rows, columns=["block", "item", "value", "how set"]); t2.to_csv(P / "table2_parameters.csv", index=False)
# ---- Table 3: scores
m2 = pd.read_csv(D / "step2_metrics.csv"); to = m2[m2.scheme == "transect_out"].set_index("model"); rf = m2[m2.scheme == "random_10fold"].set_index("model"); m3 = pd.read_csv(D / "step3_metrics.csv"); do = m3[m3.scheme == "day_out"].set_index("law"); lg = pd.read_csv(D / "step3b_law_grid.csv").set_index(["h0", "gamma"])
rh = pd.read_csv(D / "step6_readings_heldout_adopted_law.csv"); mae = np.mean(np.abs(rh.velocity_ms - rh.v_pred)); mad = np.mean(np.abs(rh.velocity_ms - rh.velocity_ms.median())); cov = np.mean((rh.velocity_ms >= rh.v_lo) & (rh.velocity_ms <= rh.v_hi)); below = rh.velocity_ms < 0.2; hit = np.mean((rh.v_pred < 0.2) == below); base = max(below.mean(), 1 - below.mean())
ad = pd.read_csv(D / "step3_adcp.csv").dropna(subset=["velocity_ms", "stage_m"]); g_ = np.clip(ad.stage_m - H0, 0.02, None) ** GAM; c_ = np.median(ad.velocity_ms / g_); mae_adcp = np.mean(np.abs(ad.velocity_ms - c_ * g_)) * 1000; adc = json.load(open(D / "adcp_recheck.json"))
cs = pd.read_csv(D / "step6_continuity_common_stage.csv"); sp = cs.groupby("stage_m").Q_med.agg(lambda x: np.std(np.log(x), ddof=1)); rr = json.load(open(D / "step6_reach_rating.json")); same = cn.groupby("date").filter(lambda g: len(g) > 1).groupby("date").Q_measured_m3s.agg(lambda x: x.max() / x.min())
s7 = json.load(open(D / "step7_summary.json")); tot = pd.read_csv(D / "step5_reach_totals.csv")
rows = [
 ("Section held out (20 folds, 257 points)", "CRPS of log A: six features + notes / bank distance alone / reach mean", f"{to.CRPS_log['full features + notes']:.3f} / {to.CRPS_log['generic shape (bank distance)']:.3f} / {to.CRPS_log['reach mean']:.3f}"),
 ("", "MAE of log A: six features + notes / reach mean", f"{to.MAE_log['full features + notes']:.3f} / {to.MAE_log['reach mean']:.3f}"),
 ("", "95% band coverage, six features + notes", f"{to.coverage95['full features + notes']:.2f}"),
 ("", "random ten-fold CRPS for the same three models (the cautionary comparison)", f"{rf.CRPS_log['full features + notes']:.3f} / {rf.CRPS_log['generic shape (bank distance)']:.3f} / {rf.CRPS_log['reach mean']:.3f}"),
 ("Day held out (11 folds, 257 points)", f"CRPS of log a: stage law ({H0}, {GAM}) / local depth^0.5 / no law", f"{lg.CRPS[PRIMARY]:.3f} / {do.CRPS['depth^0.5']:.3f} / {do.CRPS['none']:.3f}"),
 ("", f"CRPS across the law set ({len(LAW_RIDGE)} pairs)", f"{lg.CRPS.min():.3f} to {lg.CRPS.loc[LAW_RIDGE].max():.3f}"),
 ("", "day level: sd of the held-out day means before / after the law (log)", "0.34 / 0.17"),
 ("Readings on held-out days (1,208 readings)", "MAE / mean absolute deviation of the readings (m/s)", f"{mae:.3f} / {mad:.3f}"),
 ("", "95% band coverage", f"{cov:.2f}"),
 ("", "below or above 0.2 m/s called right / base rate", f"{hit*100:.0f}% / {base*100:.0f}%"),
 ("Sections on held-out days (20)", "measured section mean velocity inside the 90% band", f"{cv.mean_v_in_band.mean()*100:.0f}% (MAE {np.mean(np.abs(cv.mean_v_measured - cv.mean_v_med)):.3f} m/s, median relative error {np.median(np.abs(cv.mean_v_measured - cv.mean_v_med) / cv.mean_v_measured)*100:.0f}%)"),
 ("", "measured share below 0.2 m/s inside the 90% band", f"{cv.share_in_band.mean()*100:.0f}% (MAE {np.mean(np.abs(cv.share_measured - cv.share_med))*100:.0f} points)"),
 ("ADCP through 2025 (shape check, level fitted)", f"MAE of the law c (stage - {H0})^{GAM} against the 10-min record", f"{mae_adcp:.0f} mm/s"),
 ("", "sd of the log ratio over the season", f"{json.load(open(D / 'step3b_adcp_shape.json'))['adcp_shape'][f'reach stage law (stage-{H0:.2f})^{GAM:.2f} | ADCP read as point at fixed height']['log_ratio_sd_over_season']:.3f}"),
 ("", "its own exponent on (stage - 0.45) at the surveyed depth: point / column reading", f"{adc['adcp_exponent_on_stage_minus_0.45']['point reading (v proportional to a)']:.2f} / {adc['adcp_exponent_on_stage_minus_0.45']['column reading with the surveyed depth']:.2f}"),
 ("Continuity (measured)", "20 section discharges from 11 days on one rating", f"Q = {rr['coef']:.0f} (stage - 0.45)^{rr['exponent']:.2f} m3/s, {rr['scatter_pct']:.0f}% scatter"),
 ("", "same-day sections, max/min of the measured discharge", ", ".join(f"{v:.2f}" for v in same.values)),
 ("Continuity (modelled, common stage)", "sd of log Q across the 20 sections at 0.70 / 0.90 / 1.30 m", f"{sp[0.70]:.3f} / {sp[0.90]:.3f} / {sp[1.30]:.3f}"),
 ("Reach share below 0.2 m/s at 0.90 m", "sections (width- and length-weighted) / map with margins to the bank line", f"{float(tot[(tot.stage_m == 0.9) & tot.source.str.startswith('sections')].share_below_thr_median.iloc[0]):.3f} / {float(tot[(tot.stage_m == 0.9) & tot.source.str.startswith('map')].share_below_thr_median.iloc[0]):.3f}"),
 ("", "90% band half-width inside the surveyed range (points)", f"{s7['season']['2025']['band_halfwidth_points_in_range']:.1f}"),
 ("Season medians", "2025 record / hours inside the surveyed range", f"{s7['season']['2025']['median_share']*100:.0f}% / {s7['season']['2025']['median_share_in_range_hours']*100:.0f}% ({(1-s7['season']['2025']['fraction_of_hours_outside_surveyed_range'])*100:.0f}% of hours)"),
 ("", "2026 record / hours inside the surveyed range", f"{s7['season']['2026']['median_share']*100:.0f}% / {s7['season']['2026']['median_share_in_range_hours']*100:.0f}% ({(1-s7['season']['2026']['fraction_of_hours_outside_surveyed_range'])*100:.0f}% of hours)"),
 ("", "May 2026 (extrapolated: 24% of its hours above the range)", f"{s7['season']['May_2026']['share_median']*100:.0f}% [{s7['season']['May_2026']['lo']*100:.0f}, {s7['season']['May_2026']['hi']*100:.0f}]"),
 ("Nests (59 in wet cells)", "AUC of nests against the wet channel: field velocity / bank distance / place factor / depth", f"{s7['nests']['AUC_available_vs_used']:.2f} / {s7['nests']['baselines']['bank distance']['AUC']:.2f} / {s7['nests']['baselines']['conveyance A (no stage)']['AUC']:.2f} / {s7['nests']['baselines']['depth']['AUC']:.2f}"),
 ("", "field velocity within bank-distance strata (+-2 m)", f"{s7['nests']['baselines']['velocity within bank-distance strata (±2 m)']['AUC']:.2f} (p = {s7['nests']['baselines']['velocity within bank-distance strata (±2 m)']['ks_p']:.2f})"),
]
t3 = pd.DataFrame(rows, columns=["test", "quantity", "value"]); t3.to_csv(P / "table3_scores.csv", index=False)
# ---- Table 4: ablations
fa = pd.read_csv(D / "step2_feature_ablation.csv"); ka = pd.read_csv(D / "step2_kernel_ablation.csv")
t4 = pd.concat([fa.assign(kind="features").rename(columns={"features": "variant"})[["kind", "variant", "CRPS", "sections_better_than_full", "wilcoxon_p_vs_full"]], ka.assign(kind="kernel", variant=ka.kernel + ", " + ka.restarts.astype(str) + " restarts")[["kind", "variant", "CRPS"]]], ignore_index=True); t4.round(4).to_csv(P / "table4_ablations.csv", index=False)
# ---- markdown
def md(df): 
    cols = list(df.columns); out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows(): out.append("| " + " | ".join("" if (isinstance(v, float) and np.isnan(v)) else (f"{v:.3f}" if isinstance(v, float) else str(v)) for v in r.values) + " |")
    return "\n".join(out)
txt = f"""# Tables for the paper (Middle reach with the stage sensor and the ADCP; 22 September 2026)

Built by `32_paper_tables.py` from `data_clean/`. CSVs sit beside this file.

## Table 1. The surveyed sections

Twenty sections, eleven days. Velocity and share are width-weighted over the section's locations; the discharge is the sum of v h w over them; the last three columns are the chain's prediction for the section with its own day held out.

{md(t1.round(3))}

## Table 2. Parameters and settings of the chain

{md(t2)}

## Table 3. Validation scores

{md(t3)}

## Table 4. Ablations (section held out, stage law)

{md(t4.round(3))}
"""
open(P / "TABLES.md", "w").write(txt); print(txt[:3000])
