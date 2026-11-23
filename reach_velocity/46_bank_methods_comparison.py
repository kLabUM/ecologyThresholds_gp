"""Step 10d. The bank-line methods explored, side by side (logbook entries 27 to 31): one figure and one table to show what each
method is and what it does to the results.

The methods: A, the old lines (44_, the original notebook code); B, the fixed hybrid (surveyed edges with 3DEP between them, 43_);
C moving with stage (the edges at the sensor's water level, entry 30), once with the map's outline at the record's highest stage and once
at the highest survey-day stage; and C fixed at each reach's median survey-day level,
the version the model uses (entry 31). Fixed C is C's method (the surveyed edges first, the 3DEP lidar where the survey has no edge,
interpolated between them by water level) read at one level and pushed out so that every surveyed point is inside.

Inputs: data_clean/bank_methods/<method>/ (the tables each earlier run wrote, kept for this comparison only), the live data_clean tables
(the version in use), data_clean/bankcmp/ (43_: set B and the section scores under the entry-27 rule, scores_points_fixed_<reach>.csv for
the fixed versions).
Outputs: data_clean/bank_methods/comparison_table.csv and .md, figures/20_bank_methods_comparison.png"""
import json
from pathlib import Path
import numpy as np, pandas as pd, geopandas as gpd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import banklines_3dep as B

HERE = Path(__file__).resolve().parent; D = HERE / "data_clean"; F = HERE / "figures"; M = D / "bank_methods"; CMP = D / "bankcmp"; FT = 0.3048
REACHES = ("Upper", "Middle", "Lower")
INK, INK2, AXIS, GRID = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6,
                     "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})

# key, label, folder with its tables (None: not run through the model), colour
METHODS = [("A", "A: old lines", M / "set_A", "#eb6834"),
           ("B", "B: fixed hybrid", None, "#1a9e74"),
           ("C_stage_record", "C moving with stage, outline at record high", M / "C_stage_record_outline", "#c4b5fd"),
           ("C_stage", "C moving with stage", M / "C_stage", "#7c3aed"),
           ("C_fixed_median", "C fixed at the median survey level (in use)", D, "#1d4ed8")]
LAB = {k: l for k, l, _, _ in METHODS}; SHORT = {"A": "A: old lines", "B": "B: fixed hybrid", "C_stage_record": "C moving, record-high outline", "C_stage": "C moving with stage", "C_fixed_median": "C fixed at the median (in use)"}; COL = {k: c for k, _, _, c in METHODS}; DIR = {k: p for k, _, p, _ in METHODS}


def jl(p):
    return json.load(open(p)) if Path(p).exists() else None


def band(d, f="{:.3f}"):
    return f"{f.format(d['median'])} [{f.format(d['lo'])}, {f.format(d['hi'])}]"


# ---------------------------------------------------------------- the entry-27 rule's scores (sections held out, all methods)
S = pd.read_csv(CMP / "scores_sections.csv").replace({"set": {"C": "C_stage"}})
fx = [pd.read_csv(CMP / f"scores_points_fixed_{r}.csv") for r in REACHES if (CMP / f"scores_points_fixed_{r}.csv").exists()]
if fx:
    P = pd.concat(fx).replace({"set": {"Cfix_median": "C_fixed_median"}})
    S = pd.concat([S, P.groupby(["reach", "set", "transect"]).agg(crps=("crps", "mean"), in90=("in90", "mean"), n=("crps", "size")).reset_index()])
per = S.groupby(["reach", "set"]).agg(crps=("crps", "mean"), in90=("in90", "mean"))
W = S.pivot_table(index=["reach", "transect"], columns="set", values="crps")


def paired(a, b, scope=None):
    w = W if scope is None else W.loc[scope]; d = (w[a] - w[b]).dropna()
    return d.mean(), d.std(ddof=1) / np.sqrt(len(d)), len(d)


# ---------------------------------------------------------------- surveyed points outside each method's lines
def outside(pts, banks):
    """Surveyed points outside the lines, with 43_'s 0.5 m tolerance."""
    from shapely.geometry import Point
    poly = B.polygon(banks).buffer(0.5 / FT)
    return int(sum(not poly.contains(Point(x, y)) for x, y in zip(pts.Easting, pts.Northing)))


build = jl(CMP / "summary_build.json"); out = {r: dict(build[r]["points_outside"]) for r in REACHES}
SFs = {r: B.StageFeatures(r) for r in REACHES}   # the version in use (fixed at the median)
for r in REACHES:
    out[r]["C_fixed_median"] = outside(SFs[r].p, SFs[r].lines(SFs[r].fixed))

# ---------------------------------------------------------------- the table
rows = []
def add(item, **v):
    rows.append({"item": item, **{LAB[k]: v.get(k, "") for k in LAB}})

add("What it is",
    A="the original notebook's lines, one set for all stages", B="surveyed edges, 3DEP edges between transects, one set",
    C_stage_record="edges at the sensor's water level each time step; map outline at the record's highest stage",
    C_stage="edges at the sensor's water level each time step; map outline at the highest survey-day stage",
    C_fixed_median="C's edges at one level per reach, the median of the survey-day levels, all points inside")
lev = {r: SFs[r].fixed for r in REACHES}
add("Level of the fixed lines (Upper / Middle / Lower)",
    C_fixed_median=f"{lev['Upper']:.2f} m NAVD88 / {lev['Middle']:.3f} m sensor / {lev['Lower']:.2f} m NAVD88")
add("Surveyed points outside the lines (Upper / Middle / Lower)",
    A=" / ".join(str(out[r]["A"]) for r in REACHES), B=" / ".join(str(out[r]["B"]) for r in REACHES),
    C_stage="0 / 0 / 0 (each point on its own day)", C_stage_record="0 / 0 / 0 (each point on its own day)",
    C_fixed_median=" / ".join(str(out[r]["C_fixed_median"]) for r in REACHES))
for r in REACHES:
    add(f"{r}: section held out, CRPS of ln A (share of points in the 90% band)",
        **{k: f"{per.loc[(r, k), 'crps']:.3f} ({per.loc[(r, k), 'in90']:.2f})" for k in LAB if (r, k) in per.index})
add("All sections: CRPS minus A's, mean +/- se (sections better than A)",
    **{k: "{:+.3f} +/- {:.3f}".format(*paired(k, "A")[:2]) + f" ({int(((W[k] - W['A']).dropna() < 0).sum())} of {paired(k, 'A')[2]})" for k in LAB if k in W and k != "A"})
add("Run through the whole model?", A="yes", B="no (scored under the entry-27 rule only)", C_stage_record="yes (Middle)", C_stage="yes", C_fixed_median="yes")


def mid(k, f):
    p = DIR[k]
    return None if p is None else f(p)


vt = {k: mid(k, lambda p: pd.read_csv(p / "step6_validation_table.csv").set_index("test").value if (p / "step6_validation_table.csv").exists() else None) for k in LAB}
def vrow(test, lab):
    add(lab, **{k: str(v[test]) for k, v in vt.items() if v is not None and test in v.index})
vrow("readings on held-out days: MAE / their mean absolute deviation (m/s)", "Middle: readings on held-out days, MAE / spread (m/s)")
def share(p):
    t = pd.read_csv(p / "step5_reach_totals.csv"); r = t[(t.stage_m == 0.9) & t.source.str.startswith("sections")].iloc[0]
    return f"{r.share_below_thr_median:.3f} [{r.lo:.3f}, {r.hi:.3f}]"
add("Middle: share of the width below 0.2 m/s at 0.90 m", **{k: mid(k, share) for k in LAB if DIR[k] is not None})
nb = {k: mid(k, lambda p: jl(p / "nest_bayes.json")) for k in LAB}
add("Middle: mean suitability E at 0.90 m", **{k: band(v["reach_mean_suitability"]["at_0.90_m"]) for k, v in nb.items() if v})
add("Middle: E median over 2025 / 2026", **{k: f"{v['reach_mean_suitability']['2025']['median']:.3f} / {v['reach_mean_suitability']['2026']['median']:.3f}" for k, v in nb.items() if v})
add("Middle: velocity effect on nest use (per m/s)", **{k: "{:.2f} +/- {:.2f}".format(*v["velocity_coef_per_m_s"]["posterior"]) for k, v in nb.items() if v})
add("Middle: share of draws with velocity helping nest use", **{k: (lambda x: "below 0.1%" if x < 0.001 else f"{x:.1%}")(v['v50_m_s_posterior']['share_of_draws_with_positive_velocity_effect']) for k, v in nb.items() if v})
for r in ("Upper", "Lower"):
    ss = {k: mid(k, lambda p: jl(p / f"reach_{r}_chain_summary.json")) for k in LAB}
    add(f"{r}: E median per season", **{k: ", ".join(f"{y} {a['E_median']:.3f}" for y, a in v["season"].items()) for k, v in ss.items() if v})
T = pd.DataFrame(rows); T.to_csv(M / "comparison_table.csv", index=False)
with open(M / "comparison_table.md", "w") as fh:
    fh.write("# Bank-line methods explored (logbook entries 27 to 31)\n\nScores: each surveyed section held out in turn (the rule of entry 27; lower CRPS is better). "
             "Ecological numbers: the Middle reach at the sensor stage 0.90 m, 90% bands in brackets. Blank: not run for that method. "
             "Figure: figures/20_bank_methods_comparison.png. Written by 46_bank_methods_comparison.py.\n\n")
    fh.write("| | " + " | ".join(T.columns[1:]) + " |\n|" + "---|" * len(T.columns) + "\n")
    for _, r_ in T.iterrows(): fh.write("| " + " | ".join(str(x) if pd.notna(x) else "" for x in r_.values) + " |\n")
print(T.to_string(index=False))

# ---------------------------------------------------------------- the figure
fig = plt.figure(figsize=(17, 10.5)); gs = fig.add_gridspec(2, 3, height_ratios=[1.3, 1], hspace=0.5, wspace=0.42)
for j, r in enumerate(REACHES):
    ax = fig.add_subplot(gs[0, j]); sf = SFs[r]; p = sf.p; x0, y0 = p.Easting.median(), p.Northing.median()
    mv = B.StageFeatures(r, fixed=None); lv = sorted(mv.day_stage.values()); lo, hi = lv[0], lv[-1]
    def draw(bk, c, lw, ls="-"):
        for g in bk.geometry:
            xs, ys = np.array(g.coords).T; ax.plot((xs - x0) * FT, (ys - y0) * FT, color=c, lw=lw, ls=ls)
    draw(B.set_A(r), COL["A"], 1.0); draw(gpd.read_file(CMP / f"banks_B_{r}.geojson").to_crs(B.CRS_FT), COL["B"], 0.9, "--")
    draw(mv.lines(lo), COL["C_stage"], 0.6, ":"); draw(mv.lines(hi), COL["C_stage"], 0.6, ":")
    draw(sf.lines(sf.fixed), COL["C_fixed_median"], 1.5)
    ax.scatter((p.Easting - x0) * FT, (p.Northing - y0) * FT, s=2, color=INK, lw=0, alpha=0.6)
    pad = 40; ax.set_xlim((p.Easting.min() - x0) * FT - pad, (p.Easting.max() - x0) * FT + pad); ax.set_ylim((p.Northing.min() - y0) * FT - pad, (p.Northing.max() - y0) * FT + pad)
    ax.set_aspect("equal"); ax.grid(True); ax.set_axisbelow(True); ax.tick_params(labelsize=6); ax.set_xlabel("east (m)")
    if j == 0: ax.set_ylabel("north (m)")
    u = "m sensor stage" if r == "Middle" else "m NAVD88"
    ax.set_title(f"({'abc'[j]}) {r}: survey days {lo:.2f} to {hi:.2f}, median {sf.fixed:.2f} {u}", fontsize=9)

fig.legend(handles=[Line2D([], [], color=COL["A"], label="A: old lines"), Line2D([], [], color=COL["B"], ls="--", label="B: fixed hybrid (surveyed edges, 3DEP between)"),
                    Line2D([], [], color=COL["C_stage"], ls=":", label="C moving with stage, drawn at the lowest and highest survey day"),
                    Line2D([], [], color=COL["C_fixed_median"], lw=1.5, label="C fixed at the median survey-day level (in use)"), Line2D([], [], color=INK, marker="o", ls="", ms=3, label="surveyed points")],
           loc="upper center", bbox_to_anchor=(0.5, 0.488), ncol=5, fontsize=7.5)
ax = fig.add_subplot(gs[1, 0]); ks = [k for k in LAB if k in W.columns]; x = np.arange(len(REACHES)); w = 0.8 / len(ks)
for i, k in enumerate(ks):
    v = [per.loc[(r, k), "crps"] if (r, k) in per.index else np.nan for r in REACHES]; ax.bar(x + (i - (len(ks) - 1) / 2) * w, v, w, color=COL[k], label=SHORT[k])
ax.set_xticks(x); ax.set_xticklabels(REACHES); ax.set_ylabel("CRPS of ln A, each section held out (lower is better)"); ax.set_ylim(0, None)
ax.legend(fontsize=6, loc="upper right"); ax.grid(True, axis="y"); ax.set_axisbelow(True); ax.set_title("(d) Velocity model, each section held out", fontsize=9)

ax = fig.add_subplot(gs[1, 1])
for k in LAB:
    p = DIR[k]
    if p is None or not (p / "nest_bayes_curve.csv").exists(): continue
    t = pd.read_csv(p / "nest_bayes_curve.csv"); main = k in ("A", "C_fixed_median")
    if main: ax.fill_between(t.stage_m, t.E_lo, t.E_hi, color=COL[k], alpha=0.12, lw=0)
    ax.plot(t.stage_m, t.E_med, color=COL[k], lw=1.6 if main else 1.0, label=SHORT[k])
ax.axvspan(0.74, 1.09, color=GRID, alpha=0.5, lw=0); ax.set_ylim(0, 1); ax.set_xlabel("sensor stage (m)"); ax.set_ylabel("mean suitability E (0 to 1)")
ax.legend(fontsize=6); ax.grid(True); ax.set_axisbelow(True); ax.set_title("(e) Middle: mean suitability E against stage\n(90% bands for A and the version in use; grey: surveyed range)", fontsize=9)

ax = fig.add_subplot(gs[1, 2]); ks2 = [k for k in LAB if nb.get(k)]
for i, k in enumerate(ks2):
    e = nb[k]["reach_mean_suitability"]["at_0.90_m"]; ax.errorbar([e["median"]], [i], xerr=[[e["median"] - e["lo"]], [e["hi"] - e["median"]]], fmt="o", color=COL[k], capsize=3)
    b_, s_ = nb[k]["velocity_coef_per_m_s"]["posterior"]; ax.text(0.98, i + 0.25, f"velocity effect {b_:.1f} +/- {s_:.1f} per m/s", fontsize=6.5, ha="right", va="center", color=INK2, transform=ax.get_yaxis_transform())
ax.set_yticks(range(len(ks2))); ax.set_yticklabels([SHORT[k] for k in ks2], fontsize=7); ax.set_xlim(0, 1.35); ax.set_ylim(-0.5, len(ks2) - 0.3); ax.set_xlabel("mean suitability E at 0.90 m, median and 90% band")
ax.grid(True, axis="x"); ax.set_axisbelow(True); ax.set_title("(f) Middle: E at 0.90 m and the nest model's velocity effect", fontsize=9, loc="right")
fig.suptitle("Bank-line methods explored (entries 27 to 31). The model uses C fixed at each reach's median survey-day level: the surveyed edges first, the 3DEP lidar where the survey has none", fontsize=10.5, y=0.97)
fig.savefig(F / "20_bank_methods_comparison.png", bbox_inches="tight"); plt.close(fig); print("saved 20")
