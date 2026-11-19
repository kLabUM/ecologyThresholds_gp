"""Step 7. Deliverable 2: the ecological number. The share of the reach's wetted width below a velocity threshold over
time (length-weighted over the sections, with bands from the section draws), the ecological rating curve (share against
stage) for the depth-averaged and the near-bed reference, the sensitivity to the threshold, and the check against the
nest sites held aside since day one. Outputs: step7_reach_series.csv, step7_sensitivity.csv, step7_nests.csv,
step7_summary.json, figures/07_step7_ecology.png"""
import json, sys, time, warnings
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import matplotlib.dates as mdates
from gpcore import *; from uncertainty import *
warnings.filterwarnings("ignore"); t0 = time.time(); log = lambda *a: (print(f"[{time.time()-t0:5.0f}s]", *a), sys.stdout.flush())
THR = 0.20; Z_BED = 0.10
df = load_points(); df["bed_m"] = df.bed_rtk_sensor_datum_m.fillna(df.bed_elev_m); ch = Chain(df, banks=banks_for("Middle"))
sec_pts = {t: g.assign(_x=g.x_m.round(2), _y=g.y_m.round(2)).drop_duplicates(["_x", "_y"]).sort_values("across_m").reset_index(drop=True) for t, g in df.groupby("transect")}
order = sorted(sec_pts, key=lambda t: sec_pts[t].along_m.median()); al = np.array([sec_pts[t].along_m.median() for t in order]); edges = np.r_[al[0] - (al[1] - al[0]) / 2, (al[1:] + al[:-1]) / 2, al[-1] + (al[-1] - al[-2]) / 2]; wts = np.diff(edges); wts = wts / wts.sum()
# ---- reach series from the section draws (length-weighted per draw)
z = np.load(D / "step6_section_draws.npz", allow_pickle=True); times = pd.to_datetime(z["times"]); stages = z["stages"]; share_draws = z["reach_share"]; order = list(z["sections"]); wts = z["length_weights"]
SURVEYED = (0.74, 1.09); in_range = (stages >= SURVEYED[0]) & (stages <= SURVEYED[1])
reach = pd.DataFrame({"time_utc": times, "stage_m": stages, "share_med": np.nanmedian(share_draws, 0), "share_lo": np.nanquantile(share_draws, .05, 0), "share_hi": np.nanquantile(share_draws, .95, 0), "in_surveyed_range": in_range}); reach.to_csv(D / "step7_reach_series.csv", index=False)
season = {}
for yr in [2025, 2026]:
    r = reach[reach.time_utc.dt.year == yr]; ri = r[r.in_surveyed_range]; season[str(yr)] = {"median_share": float(np.nanmedian(r.share_med)), "median_share_in_range_hours": float(np.nanmedian(ri.share_med)), "fraction_of_hours_outside_surveyed_range": float(1 - r.in_surveyed_range.mean()), "fraction_of_hours_above_range": float((r.stage_m > SURVEYED[1]).mean()), "band_halfwidth_points_in_range": float(np.nanmedian((ri.share_hi - ri.share_lo) / 2) * 100), "band_halfwidth_points_above_range": float(np.nanmedian((r[r.stage_m > SURVEYED[1]].share_hi - r[r.stage_m > SURVEYED[1]].share_lo) / 2) * 100) if (r.stage_m > SURVEYED[1]).any() else np.nan}
may26 = reach[(reach.time_utc >= "2026-05-01") & (reach.time_utc < "2026-06-01")]; season["May_2026"] = {"share_median": float(np.nanmedian(may26.share_med)), "lo": float(np.nanmedian(may26.share_lo)), "hi": float(np.nanmedian(may26.share_hi)), "fraction_above_range": float((may26.stage_m > SURVEYED[1]).mean())}
log(f"reach share below {THR} m/s: {json.dumps(season)}")
# ---- ecological rating curves: share against stage, two references, and the threshold sweep
grid_st = np.linspace(0.65, 1.35, 29); rows = []; secs_all = {t: sec_pts[t] for t in order}
for ref, zr in [("depth-averaged", None), (f"near-bed at {Z_BED:.2f} m", Z_BED)]:
    for thr in [0.10, 0.15, 0.20, 0.30, 0.40]:
        R = sample_reach(ch, grid_st, secs_all, n=300, seed=11, threshold=thr, z_ref=zr); acc = R["reach_share"]
        for j, st in enumerate(grid_st): rows.append({"reference": ref, "threshold_m_s": thr, "stage_m": st, "share_med": float(np.nanmedian(acc[:, j])), "share_lo": float(np.nanquantile(acc[:, j], .05)), "share_hi": float(np.nanquantile(acc[:, j], .95))})
    log(f"{ref}: done")
# the alternative-law scenario above the surveyed range: the flattest member of the law set alone, against the primary alone
scen = {}
for lab, key in [(f"primary {ch.primary}", ch.primary), ("flattest in the set (0.65, 0.5)", (0.65, 0.5))]:
    ch.primary, keep = key, ch.primary; R = sample_reach(ch, grid_st, secs_all, n=200, seed=12, threshold=THR, sources=("A_field", "disc", "bed", "stage", "z0")); ch.primary = keep
    scen[lab] = {f"{st:.2f}": float(np.nanmedian(R["reach_share"][:, j])) for j, st in enumerate(grid_st) if abs(st - round(st * 20) / 20) < 1e-6}
json.dump(scen, open(D / "step7_law_scenarios.json", "w"), indent=1); log(f"law scenarios (share below {THR}): " + "; ".join(f"{k}: 1.30 m {v['1.30']:.2f}" for k, v in scen.items()))
sens = pd.DataFrame(rows); sens.to_csv(D / "step7_sensitivity.csv", index=False)
# ---- nests: predicted velocity at each nest against what was available, at the nest's own stage
nests = pd.read_csv(D / "nests.csv"); nests = nests.drop_duplicates("location_id").copy(); grid = pd.read_csv(D / "grid.csv")
st10 = pd.read_csv(D / "series_10min.csv", parse_dates=["time_utc"]).set_index("time_utc").stage_m; st10.index = pd.to_datetime(st10.index, utc=True)
def stage_for(t):
    t = pd.Timestamp(t, tz="UTC") if pd.Timestamp(t).tzinfo is None else pd.Timestamp(t); i = st10.index.searchsorted(t); i = min(max(i, 0), len(st10) - 1); return float(st10.iloc[i]) if abs((st10.index[i] - t).total_seconds()) < 6 * 3600 else np.nan
nests["stage_use"] = [stage_for(t) if pd.notna(t) else np.nan for t in nests.time_utc]
day_med = st10.groupby(st10.index.date).median(); day_med.index = pd.to_datetime(day_med.index).strftime("%Y-%m-%d")
nests.loc[nests.stage_use.isna(), "stage_use"] = nests.loc[nests.stage_use.isna(), "date"].map(day_med)   # no time recorded: the day's median stage
nests.loc[nests.stage_use.isna() & (nests.date == "2025-05-27"), "stage_use"] = float(st10["2025-05-30"].median())   # the sensor started 30 May; three days earlier, same base flow at the gauge
assert nests.stage_use.notna().all(), nests[nests.stage_use.isna()]
from scipy.spatial import cKDTree
tree = cKDTree(grid[["x_m", "y_m"]].values); dist, idx = tree.query(nests[["x_m", "y_m"]].values); nests["grid_dist_m"] = dist; nests = nests[nests.grid_dist_m < 6].copy(); idx = idx[nests.index.map(lambda i: True)] if False else tree.query(nests[["x_m", "y_m"]].values)[1]
gcell = grid.iloc[idx].reset_index(drop=True); nests = nests.reset_index(drop=True); GLV = GridLevels(); nest_dist = []; pct, pct_near, vpred, vpred_h, hcell = [], [], [], [], []
near_along = np.array([np.min(np.abs(gcell.along_m.values - a)) < 30 for a in grid.along_m])   # the restricted null: cells within 30 m along of some nest
def v_of(logA, h, st, z=None):
    a = np.exp(logA) * np.clip(st - ch.primary[0], 0.02, None) ** ch.primary[1]
    if z is None: return np.where(h > np.e * Z0, a * (np.log(np.clip(h, np.e * Z0, None) / Z0) - 1), np.nan)   # depth-averaged, nan where dry
    return np.where(h > z, a * np.log(z / Z0), np.nan)   # at a height above the bed
for k, r in nests.iterrows():
    st = r.stage_use; cell = gcell.iloc[k]; gA, gD = GLV.at(st); cA = gA[idx[k]]; nest_dist.append(float(gD[idx[k]])); h = max(st - cell.bed_m, 0); hcell.append(h); v = float(v_of(cA, h, st)); vpred.append(v)   # field and bank distance at the nest's own stage (bank set C)
    vpred_h.append(float(v_of(cA, h, st, z=r.height_m)) if pd.notna(r.height_m) else np.nan)
    hg = np.clip(st - grid.bed_m.values, 0, None); vg = v_of(gA, hg, st); ok = ~np.isnan(vg)
    pct.append(float(np.mean(vg[ok] < v)) if not np.isnan(v) else np.nan); pct_near.append(float(np.mean(vg[ok & near_along] < v)) if not np.isnan(v) else np.nan)
nests["h_cell_m"] = hcell; nests["v_pred_ms"] = vpred; nests["v_pred_at_height_ms"] = vpred_h; nests["percentile_in_reach"] = pct; nests["percentile_near"] = pct_near
nests[["location_id", "date", "x_m", "y_m", "stage_use", "h_cell_m", "water_depth_m", "height_m", "velocity_ms", "v_pred_ms", "v_pred_at_height_ms", "percentile_in_reach", "percentile_near", "grid_dist_m"]].to_csv(D / "step7_nests.csv", index=False)
dry = nests.v_pred_ms.isna(); log(f"nests in cells dry at their stage (bed above the water in the 2 m grid): {int(dry.sum())} of {len(nests)}; these are set aside from the percentile test")
st_ref = float(np.nanmedian(nests.stage_use)); hg = np.clip(st_ref - grid.bed_m.values, 0, None); gA_ref, gD_ref = GLV.at(st_ref); ag = np.exp(gA_ref) * np.clip(st_ref - ch.primary[0], 0.02, None) ** ch.primary[1]; v_avail = np.where(hg > np.e * Z0, ag * (np.log(np.clip(hg, np.e * Z0, None) / Z0) - 1), np.nan); v_avail = v_avail[~np.isnan(v_avail)]
from scipy.stats import kstest
q = nests.dropna(subset=["percentile_in_reach"]); ks = kstest(q.percentile_in_reach, "uniform"); ks_near = kstest(q.percentile_near, "uniform")
# baselines: the same rank test with other scores (lower = 'better' habitat by that score), and the model's velocity within bank-distance strata (±2 m)
pct_b = {"bank distance": [], "conveyance A (no stage)": [], "depth": []}; pct_strat = []
for k, r in nests.iterrows():
    if pd.isna(r.percentile_in_reach): continue
    st = r.stage_use; cell = gcell.iloc[k]; gA, gD = GLV.at(st); cA, cD = gA[idx[k]], gD[idx[k]]; hg = np.clip(st - grid.bed_m.values, 0, None); vg = v_of(gA, hg, st); ok = ~np.isnan(vg); hc = max(st - cell.bed_m, 0)
    pct_b["bank distance"].append(float(np.mean(gD[ok] < cD))); pct_b["conveyance A (no stage)"].append(float(np.mean(gA[ok] < cA))); pct_b["depth"].append(float(np.mean(hg[ok] < hc)))
    strat = ok & (np.abs(gD - cD) <= 2.0); pct_strat.append(float(np.mean(vg[strat] < float(v_of(cA, hc, st)))))
baselines = {k: {"AUC": float(1 - np.mean(v)), "ks_p": float(kstest(v, "uniform").pvalue)} for k, v in pct_b.items()}; baselines["velocity within bank-distance strata (±2 m)"] = {"AUC": float(1 - np.mean(pct_strat)), "ks_p": float(kstest(pct_strat, "uniform").pvalue)}
day_means = q.groupby("date").percentile_in_reach.mean(); from scipy.stats import ttest_1samp, wilcoxon
cluster = {"n_days": int(len(day_means)), "mean_of_day_means": float(day_means.mean()), "t_p_vs_0.5": float(ttest_1samp(day_means, 0.5).pvalue), "wilcoxon_p": float(wilcoxon(day_means - 0.5).pvalue) if len(day_means) >= 6 else np.nan}
nest_res = {"n_nests": int(len(nests)), "n_in_wet_cells": int(len(q)), "median_percentile_in_reach": float(q.percentile_in_reach.median()), "AUC_available_vs_used": float(1 - q.percentile_in_reach.mean()), "ks_p_vs_uniform": float(ks.pvalue),
            "AUC_restricted_to_30m_of_nests": float(1 - q.percentile_near.mean()), "ks_p_restricted": float(ks_near.pvalue), "share_of_nests_below_thr_predicted": float(np.mean(q.v_pred_ms < THR)), "share_of_reach_below_thr_at_median_nest_stage": float(np.mean(v_avail < THR)),
            "measured_nest_velocity_median": float(np.nanmedian(nests.velocity_ms)), "predicted_at_nests_depth_averaged_median": float(q.v_pred_ms.median()), "predicted_at_nests_at_reading_height_median": float(np.nanmedian(nests.v_pred_at_height_ms)), "nest_reading_height_median_m": float(np.nanmedian(nests.height_m)),
            "baselines": baselines, "day_cluster_inference": cluster, "nest_bank_distance_median_m": float(np.median(nest_dist)), "wet_cell_bank_distance_median_m": float(np.median(gD_ref[hg > np.e * Z0]))}
log(f"nests: {nest_res}"); json.dump({"season": season, "nests": nest_res}, open(D / "step7_summary.json", "w"), indent=2)
# ---- figure
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, DARK, MUTED = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#0d366b", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig = plt.figure(figsize=(16, 10)); gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.1], hspace=0.45, wspace=0.3)
for k, yr in enumerate([2025, 2026]):
    ax = fig.add_subplot(gs[0, k]); r = reach[reach.time_utc.dt.year == yr]; ax.fill_between(r.time_utc, r.share_lo * 100, r.share_hi * 100, color=BLUE, alpha=0.2, lw=0); ax.plot(r.time_utc, r.share_med * 100, color=BLUE, lw=0.9)
    ax2 = ax.twinx(); ax2.plot(r.time_utc, r.stage_m, color=ORANGE, lw=0.6, alpha=0.8); ax2.set_ylabel("stage (m)", color=ORANGE); ax2.set_ylim(0.5, 1.8)
    for d in df[df.year == yr].date.unique(): ax.axvline(pd.Timestamp(d, tz="UTC"), color=INK, lw=0.5, ls=":", alpha=0.6)
    oo = r[~r.in_surveyed_range]; ax.scatter(oo.time_utc, np.full(len(oo), 2), s=2, color=MUTED, lw=0, marker="|")
    ax.set_ylim(0, 100); ax.set_ylabel(f"share of the reach's wetted width below {THR} m/s (%)"); ax.set_title(f"{yr}: the ecological share every two hours, 90% band (dotted: survey days; grey ticks at the bottom: stage outside the surveyed range)", fontsize=8.5); ax.grid(True); ax.set_axisbelow(True); ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
ax = fig.add_subplot(gs[0, 2])
for ref, col in [("depth-averaged", BLUE), (f"near-bed at {Z_BED:.2f} m", DARK)]:
    s = sens[(sens.reference == ref) & (sens.threshold_m_s == THR)]; ax.fill_between(s.stage_m, s.share_lo * 100, s.share_hi * 100, color=col, alpha=0.15, lw=0); ax.plot(s.stage_m, s.share_med * 100, color=col, lw=1.5, label=ref)
ax.axvspan(0.74, 1.09, color=GRID, alpha=0.5, lw=0); ax.set_xlabel("stage (m)"); ax.set_ylabel(f"share below {THR} m/s (%)"); ax.set_title("The ecological rating curve: two references, 90% bands", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7); ax.set_ylim(0, 100)
ax = fig.add_subplot(gs[1, 0]); cmap = plt.get_cmap("Blues")
for i, thr in enumerate([0.10, 0.15, 0.20, 0.30, 0.40]):
    s = sens[(sens.reference == "depth-averaged") & (sens.threshold_m_s == thr)]; ax.plot(s.stage_m, s.share_med * 100, color=cmap(0.3 + 0.7 * i / 4), lw=1.4, label=f"{thr:.2f} m/s")
ax.axvspan(0.74, 1.09, color=GRID, alpha=0.5, lw=0); ax.set_xlabel("stage (m)"); ax.set_ylabel("share of wetted width below the threshold (%)"); ax.set_title("Sensitivity to the threshold (depth-averaged reference)", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7, title="threshold"); ax.set_ylim(0, 100)
ax = fig.add_subplot(gs[1, 1]); bins = np.linspace(0, 0.8, 33); ax.hist(v_avail, bins=bins, density=True, color=MUTED, alpha=0.6, label=f"available: every wet cell at stage {st_ref:.2f} m"); ax.hist(q.v_pred_ms, bins=bins, density=True, color=GREEN, alpha=0.7, label=f"predicted at the {len(q)} nests, each at its own stage")
ax.axvline(THR, color=INK, lw=0.6, ls="--"); ax.set_xlabel("depth-averaged velocity (m/s)"); ax.set_ylabel("density"); ax.set_title(f"Nests sit in the slow margin: AUC {nest_res['AUC_available_vs_used']:.2f}; bank distance alone {baselines['bank distance']['AUC']:.2f}; velocity within bank-distance strata {baselines['velocity within bank-distance strata (±2 m)']['AUC']:.2f} (p {baselines['velocity within bank-distance strata (±2 m)']['ks_p']:.2f})", fontsize=7.5); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7)
ax = fig.add_subplot(gs[1, 2]); qq = nests.dropna(subset=["velocity_ms", "v_pred_at_height_ms"]); ax.scatter(qq.velocity_ms, qq.v_pred_at_height_ms, s=14, color=GREEN, edgecolor="white", lw=0.4); ax.plot([-0.05, 0.6], [-0.05, 0.6], color=AXIS, lw=0.8); ax.axhline(THR, color=INK, lw=0.5, ls="--"); ax.axvline(THR, color=INK, lw=0.5, ls="--")
ax.set_xlabel("measured at the nest (m/s, one reading at the nest's height)"); ax.set_ylabel("predicted at the same height in the nest's 2 m cell (m/s)"); ax.set_title(f"Like for like at the reading height ({nest_res['nest_reading_height_median_m']:.2f} m): measured median {nest_res['measured_nest_velocity_median']:.2f}, predicted {nest_res['predicted_at_nests_at_reading_height_median']:.2f}", fontsize=8); ax.grid(True); ax.set_axisbelow(True); ax.set_xlim(-0.05, 0.6); ax.set_ylim(-0.05, 0.6)
fig.suptitle("Step 7: the ecological number, with its uncertainty, and the nests", fontsize=10, y=0.95); fig.savefig(F / "07_step7_ecology.png", bbox_inches="tight"); plt.close(fig); log("saved 07")
