"""Step 10c. The switch to bank set C (logbook entry 30): what the lines look like and what moved.
Before: the outputs of commit 4298715 (bank set A), kept in data_clean/banks_c/before_set_A/. After: the same tables written by the
rerun with set C. Panels: (a to c) A against C at three stages on each reach, with C's outline at the record's highest stage; (d) the
ecological rating curve (share below 0.2 m/s, depth-averaged) before and after; (e) the mean suitability E against stage; (f) where the
width of the band comes from at one section (AT10, 0.90 m), with the new bank-line source.
Outputs: data_clean/banks_c/before_after.csv, figures/19_bank_set_c_switch.png"""
import json
from pathlib import Path
import numpy as np, pandas as pd, geopandas as gpd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE = Path(__file__).resolve().parent; D = HERE / "data_clean"; F = HERE / "figures"; BEF = D / "banks_c" / "before_set_A"; FT = 0.3048
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, MUTED = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})


def jl(p): return json.load(open(p)) if Path(p).exists() else None
import banklines_3dep as B
READY = {r for r in ("Upper", "Middle", "Lower") if B.check_transects(r) <= 0.01}   # a reach whose 3DEP cache does not match the geometry has no set C yet


# ---------------------------------------------------------------- the numbers before and after
rows = []
def add(item, before, after, unit=""): rows.append({"item": item, "set_A_before": before, "set_C_after": after, "unit": unit})
vb, va = pd.read_csv(BEF / "step6_validation_table.csv").set_index("test").value, pd.read_csv(D / "step6_validation_table.csv").set_index("test").value
for k in vb.index:
    if k in va.index: add(f"Middle: {k}", vb[k], va[k])
for lab, path, key in [("Middle: reach share below 0.2 m/s at 0.90 m (sections)", "step5_reach_totals.csv", 0.9)]:
    b, a = [pd.read_csv(p / path) for p in (BEF, D)]
    for nm, t in (("before", b), ("after", a)):
        r = t[(t.stage_m == key) & t.source.str.startswith("sections")].iloc[0]; t_ = f"{r.share_below_thr_median:.3f} [{r.lo:.3f}, {r.hi:.3f}]"
        if nm == "before": bb = t_
        else: add(lab, bb, t_)
nb, na = jl(BEF / "nest_bayes.json"), jl(D / "nest_bayes.json")
add("Middle: mean suitability E at 0.90 m", "{median:.3f} [{lo:.3f}, {hi:.3f}]".format(**nb["reach_mean_suitability"]["at_0.90_m"]), "{median:.3f} [{lo:.3f}, {hi:.3f}]".format(**na["reach_mean_suitability"]["at_0.90_m"]))
add("Middle: velocity effect on nest use", "{:.2f} +/- {:.2f}".format(*nb["velocity_coef_per_m_s"]["posterior"]), "{:.2f} +/- {:.2f}".format(*na["velocity_coef_per_m_s"]["posterior"]), "per m/s")
for yr in ("2025", "2026"): add(f"Middle: E median over the {yr} record", round(nb["reach_mean_suitability"][yr]["median"], 3), round(na["reach_mean_suitability"][yr]["median"], 3))
ab, aa = jl(BEF / "adcp_single_beam.json"), jl(D / "adcp_single_beam.json"); add("Middle: ADCP level ratio (record / field)", round(ab["level_ratio_record_over_field_median"], 3), round(aa["level_ratio_record_over_field_median"], 3))
cb, ca = jl(BEF / "uncertainty_calibration.json"), jl(D / "uncertainty_calibration.json"); add("Middle: chain calibration inflation", round(cb["inflation"], 3), round(ca["inflation"], 3))
for r in ("Upper", "Lower"):
    sb, sa = jl(BEF / f"reach_{r}_chain_summary.json"), jl(D / f"reach_{r}_chain_summary.json")
    if sa is None or r not in READY: continue
    add(f"{r}: section held out CRPS, model / reach mean", f"{sb['transect_out']['CRPS_model']:.3f} / {sb['transect_out']['CRPS_reach_mean']:.3f}", f"{sa['transect_out']['CRPS_model']:.3f} / {sa['transect_out']['CRPS_reach_mean']:.3f}")
    add(f"{r}: sections on held-out days inside the band", f"{sb['section_day_out']['in_band']:.2f} (MAE {sb['section_day_out']['MAE_m_s']:.3f} m/s)", f"{sa['section_day_out']['in_band']:.2f} (MAE {sa['section_day_out']['MAE_m_s']:.3f} m/s)")
    for yr in sa["season"]:
        b_, a_ = sb["season"].get(yr, {}), sa["season"][yr]
        add(f"{r}: {yr} share below 0.2 m/s, median", round(b_.get("share_median", np.nan), 3), round(a_["share_median"], 3)); add(f"{r}: {yr} E median", round(b_.get("E_median", np.nan), 3), round(a_["E_median"], 3))
ba = pd.DataFrame(rows); ba.to_csv(D / "banks_c" / "before_after.csv", index=False); print(ba.to_string(index=False))

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(17, 11)); gs = fig.add_gridspec(2, 3, height_ratios=[1.35, 1], hspace=0.3, wspace=0.25)
LEVELS = {"Middle": ("sensor stage", [0.75, 0.90, 1.10])}
for r in ("Upper", "Lower"):
    k = jl(D / f"reach_{r}_chain_summary.json")
    if k is not None and r in READY:
        Lr = [v[1] for v in k["sensor_to_section_rating"]["days"].values()]; LEVELS[r] = ("section level L (m NAVD88)", list(np.round(np.quantile(Lr, [0.1, 0.5, 0.9]), 2)))
CS = ["#93c5fd", "#3b82f6", "#1e3a8a"]
for j, r in enumerate(("Upper", "Middle", "Lower")):
    ax = fig.add_subplot(gs[0, j])
    if r not in LEVELS: ax.text(0.5, 0.5, f"{r}: not rerun yet", ha="center", transform=ax.transAxes); ax.set_axis_off(); continue
    SF = B.StageFeatures(r); A = B.set_A(r); p = SF.p; x0, y0 = p.Easting.median(), p.Northing.median()
    out = gpd.read_file(D / "banks_c" / f"C_outline_{r}.geojson").to_crs(B.CRS_FT).geometry.iloc[0]; xs, ys = np.array(out.exterior.coords).T; ax.fill((xs - x0) * FT, (ys - y0) * FT, color=GRID, alpha=0.6, lw=0)
    for g in A.geometry: xs, ys = np.array(g.coords).T; ax.plot((xs - x0) * FT, (ys - y0) * FT, color=ORANGE, lw=1.0)
    for c, lv in zip(CS, LEVELS[r][1]):
        for g in SF.lines(lv).geometry: xs, ys = np.array(g.coords).T; ax.plot((xs - x0) * FT, (ys - y0) * FT, color=c, lw=0.8)
    ax.scatter((p.Easting - x0) * FT, (p.Northing - y0) * FT, s=2, color=INK, lw=0, alpha=0.6)
    pad = 40; ax.set_xlim((p.Easting.min() - x0) * FT - pad, (p.Easting.max() - x0) * FT + pad); ax.set_ylim((p.Northing.min() - y0) * FT - pad, (p.Northing.max() - y0) * FT + pad); ax.set_aspect("equal"); ax.grid(True); ax.set_axisbelow(True); ax.tick_params(labelsize=6)
    ax.set_xlabel("east (m)"); ax.set_ylabel("north (m)") if j == 0 else None
    ax.legend(handles=[Line2D([], [], color=ORANGE, label="A (old lines, one for all stages)")] + [Line2D([], [], color=c, label=f"C at {LEVELS[r][0]} {lv:g}") for c, lv in zip(CS, LEVELS[r][1])] + [Line2D([], [], color=GRID, lw=6, label="C at the highest survey-day stage (the map's outline)")], fontsize=6, loc="lower left")
    ax.set_title(f"({'abc'[j]}) {r}: bank lines A and C", fontsize=9)
ax = fig.add_subplot(gs[1, 0]); sb, sa = [pd.read_csv(p_ / "step7_sensitivity.csv") for p_ in (BEF, D)]
for t, c, lab in ((sb, ORANGE, "before (A)"), (sa, BLUE, "after (C)")):
    q = t[(t.reference == "depth-averaged") & np.isclose(t.threshold_m_s, 0.2)]; ax.fill_between(q.stage_m, q.share_lo * 100, q.share_hi * 100, color=c, alpha=0.15, lw=0); ax.plot(q.stage_m, q.share_med * 100, color=c, lw=1.4, label=lab)
ax.axvspan(0.74, 1.09, color=GRID, alpha=0.5, lw=0); ax.set_xlabel("sensor stage (m)"); ax.set_ylabel("share of the wetted width below 0.2 m/s (%)"); ax.set_title("(d) Middle: the ecological rating curve, 90% bands (grey: surveyed range)", fontsize=9); ax.legend(fontsize=7); ax.grid(True); ax.set_axisbelow(True); ax.set_ylim(0, 100)
ax = fig.add_subplot(gs[1, 1]); eb, ea = [pd.read_csv(p_ / "nest_bayes_curve.csv") for p_ in (BEF, D)]
for t, c, lab in ((eb, ORANGE, "before (A)"), (ea, BLUE, "after (C)")): ax.fill_between(t.stage_m, t.E_lo, t.E_hi, color=c, alpha=0.15, lw=0); ax.plot(t.stage_m, t.E_med, color=c, lw=1.4, label=lab)
ax.axvspan(0.74, 1.09, color=GRID, alpha=0.5, lw=0); ax.set_xlabel("sensor stage (m)"); ax.set_ylabel("mean suitability E (0 to 1)"); ax.set_ylim(0, 1); ax.set_title("(e) Middle: the mean suitability against stage, 90% bands", fontsize=9); ax.legend(fontsize=7); ax.grid(True); ax.set_axisbelow(True)
ax = fig.add_subplot(gs[1, 2]); ab_, aa_ = [pd.read_csv(p_ / "step4_attribution.csv") for p_ in (BEF, D)]; src = ["A_field", "disc", "bed", "stage", "law", "z0", "bank"]; labs = ["A field", "discrepancy", "bed", "stage", "law", "z0", "bank lines"]
for k, (t, c, lab) in enumerate(((ab_, ORANGE, "before (A)"), (aa_, BLUE, "after (C)"))):
    q = t[(t.section == "AT10") & np.isclose(t.stage_m, 0.90)].set_index("source"); ax.bar(np.arange(len(src)) + (k - 0.5) * 0.38, [q.mean_v_iqr90.get(s_, 0.0) for s_ in src], 0.38, color=c, label=lab)
ax.set_xticks(range(len(src))); ax.set_xticklabels(labs, fontsize=7); ax.set_ylabel("90% interval width of the section mean velocity (m/s)"); ax.set_title("(f) Where the band comes from: section AT10 at 0.90 m, one source at a time", fontsize=9); ax.legend(fontsize=7); ax.grid(True, axis="y"); ax.set_axisbelow(True)
fig.suptitle("Entry 30: the model on bank set C. Stage-dependent bank lines (a to c) and what moved at the Middle (d to f)", fontsize=10.5, y=0.94)
fig.savefig(F / "19_bank_set_c_switch.png", bbox_inches="tight"); plt.close(fig); print("saved 19")
