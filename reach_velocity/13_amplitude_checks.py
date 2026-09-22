"""Step 0f, part 2. Review material and checks for the amplitude method.
(a) One panel per location: readings, fitted log curve, water surface, 1 ft line, the derived numbers. Pages of 48 in
    figures/profiles/ and one PDF, profiles_all.pdf, for paging through.
(b) Checks that interrogate the method: z0 by substrate (from the point notes), by year and by crew; amplitude from a
    least-squares fit in linear space against one in log space; residuals of the log law against relative height z/h;
    the locations that fit worse than a flat line and where they are; bed elevation along each transect and at the one
    revisited pair; and the ADCP's depth scaling. Figure 00o_amplitude_checks.png, table amplitude_checks.json."""
import json, re
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"; F = HERE / "figures"; (F / "profiles").mkdir(exist_ok=True)
z0 = json.load(open(D / "field_qc.json"))["z0_m"]; ZMIN, VMIN = 0.03, 0.005
loc = pd.read_csv(D / "locations.csv").set_index("location_id"); tr = pd.read_csv(D / "readings_transect.csv")
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, RED = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#c0392b"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 7, "axes.edgecolor": AXIS, "axes.linewidth": 0.5, "grid.color": GRID, "grid.linewidth": 0.4, "legend.frameon": False, "savefig.dpi": 150, "axes.facecolor": "#fcfcfb"})
# ---------------------------------------------------------------- (a) per-location pages
order = loc.sort_values(["date", "transect", "across_m"]).index.tolist(); per_page = 48; pages = [order[i:i + per_page] for i in range(0, len(order), per_page)]
pdf = PdfPages(F / "profiles_all.pdf")
for pi, ids in enumerate(pages):
    fig, axes = plt.subplots(6, 8, figsize=(16, 12), gridspec_kw={"hspace": 0.6, "wspace": 0.35}); axes = axes.ravel()
    for ax, lid in zip(axes, ids):
        p = tr[tr.location_id == lid].sort_values("height_m"); r = loc.loc[lid]; ok = p[(p.height_m >= ZMIN) & (p.velocity_ms > VMIN)]
        ax.axhspan(0, r.water_depth_m, color="#e8f0fa", lw=0); ax.axhline(r.water_depth_m, color=BLUE, lw=0.6); ax.axhline(0.3048, color=INK, lw=0.5, ls=(0, (3, 2)))
        if r.a > 0:
            z = np.linspace(z0 * 1.02, max(r.water_depth_m, p.height_m.max()) * 1.05, 60); ax.fill_betweenx(z, (r.a - 1.96 * r.a_sd) * np.log(z / z0), (r.a + 1.96 * r.a_sd) * np.log(z / z0), color=BLUE, alpha=0.15, lw=0); ax.plot(r.a * np.log(z / z0), z, color=BLUE, lw=1.0)
        ax.scatter(p.velocity_ms, p.height_m, s=9, color=INK, zorder=4); ax.scatter(p.velocity_ms[~p.index.isin(ok.index)], p.height_m[~p.index.isin(ok.index)], s=14, facecolor="none", edgecolor=RED, lw=0.7, zorder=5)
        flag = ("R2<0 " if (r.fit_r2 < 0) else "") + ("1 reading " if r.n_used == 1 else "") + ("shallow " if r.shallower_than_1ft else "")
        ax.set_title(f"{r.transect} {r.point}  {r.date}\na={r.a:.3f}±{r.a_sd:.3f}  R²={r.fit_r2:.2f}  h={r.water_depth_m:.2f} m  v̄={r.v_mean:.2f}\n{flag}", fontsize=6, color=RED if flag else INK)
        ax.set_xlim(-0.02, max(0.5, p.velocity_ms.max() * 1.15)); ax.set_ylim(0, max(0.5, r.water_depth_m * 1.15)); ax.grid(True); ax.set_axisbelow(True); ax.tick_params(labelsize=5.5)
    for ax in axes[len(ids):]: ax.axis("off")
    fig.suptitle(f"Every location, page {pi+1} of {len(pages)}: readings (dots; red ring = not used), fitted log curve with 95% band, water surface (blue), 1 ft (dashed). v̄ = depth-averaged velocity from the fit (m/s)", fontsize=9, y=0.995)
    fig.savefig(F / "profiles" / f"page_{pi+1:02d}.png", bbox_inches="tight"); pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)
pdf.close(); print(f"per-location pages: {len(pages)}")
# ---------------------------------------------------------------- (b) checks
use = tr[(tr.height_m >= ZMIN) & (tr.velocity_ms > VMIN)].copy(); nper = use.groupby("location_id").size(); fitset = use[use.location_id.isin(nper[nper >= 3].index)].copy()
def fit_z0(df):
    grid = np.logspace(-4, -1, 61); sse = []
    for z in grid:
        tot = 0.0
        for _, p in df.groupby("location_id"): w = np.log(p.height_m.values / z); v = p.velocity_ms.values; a = (v * w).sum() / (w * w).sum(); tot += ((v - a * w) ** 2).sum()
        sse.append(tot)
    sse = np.array(sse); return float(grid[sse.argmin()]), float(grid[sse <= sse.min() * 1.05].min()), float(grid[sse <= sse.min() * 1.05].max()), int(df.location_id.nunique())
notes = loc.notes.fillna("").str.lower(); sub = pd.Series("not noted", index=loc.index)
for word, lab in [("cobble", "cobble"), ("gravel", "gravel"), ("sand", "sand"), ("silt", "silt"), ("clay", "silt"), ("veg", "vegetation"), ("grass", "vegetation"), ("wood", "woody debris"), ("debris", "woody debris")]: sub[notes.str.contains(word) & (sub == "not noted")] = lab
fitset["substrate"] = fitset.location_id.map(sub); fitset["year"] = fitset.location_id.map(loc.year); crew_raw = fitset.Investigators.fillna("?").str.lower().str.replace(" ", ""); crew_code = {c: f"crew {i+1}" for i, c in enumerate(sorted(crew_raw.unique()))}; fitset["crew"] = crew_raw.map(crew_code)   # crews are coded, not named, in every output
checks = {"z0_by_substrate": {}, "z0_by_year": {}, "z0_by_crew": {}}
for k, v in fitset.groupby("substrate"):
    if v.location_id.nunique() >= 8: checks["z0_by_substrate"][k] = fit_z0(v)
for k, v in fitset.groupby("year"): checks["z0_by_year"][str(k)] = fit_z0(v)
for k, v in fitset.groupby("crew"):
    if v.location_id.nunique() >= 8: checks["z0_by_crew"][k[:30]] = fit_z0(v)
# linear-space vs log-space amplitude
lin, lg = {}, {}
for lid, p in use.groupby("location_id"):
    w = np.log(p.height_m.values / z0); v = p.velocity_ms.values; lin[lid] = (v * w).sum() / (w * w).sum(); lg[lid] = np.exp(np.mean(np.log(v) - np.log(w)))
cmpdf = pd.DataFrame({"lin": lin, "log": lg}); cmpdf["ratio"] = cmpdf.log / cmpdf.lin; checks["log_vs_linear_amplitude"] = {"median_ratio": float(cmpdf.ratio.median()), "IQR": [float(cmpdf.ratio.quantile(.25)), float(cmpdf.ratio.quantile(.75))]}
# residual vs relative height
use["a"] = use.location_id.map(loc.a); use["h"] = use.location_id.map(loc.water_depth_m); use["pred"] = use.a * np.log(use.height_m / z0); use["rel"] = (use.velocity_ms - use.pred) / use.a; use["zh"] = use.height_m / use.h
bins = pd.cut(use.zh, [0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.85, 1.01]); prof = use.groupby(bins, observed=True).rel.agg(["median", "size"]); checks["residual_by_relative_height_in_units_of_a"] = {str(k): [float(v["median"]), int(v["size"])] for k, v in prof.iterrows()}
# bad fits: where
bad = loc[loc.fit_r2 < 0]; checks["fits_worse_than_flat"] = {"n": int(len(bad)), "median_|across|_m": float(bad.across_m.abs().median()), "all_median_|across|_m": float(loc.across_m.abs().median()), "median_depth_m": float(bad.water_depth_m.median()), "median_a": float(bad.a.median()), "all_median_a": float(loc.a.median()), "share_with_flow_notes": float(bad.notes.fillna("").str.contains("debris|wood|eddy|behind|slow", case=False).mean())}
# bed elevation: revisit pair and along-transect roughness
pair = loc[loc.point.isin(["AT0710", "AT0711"])][["point", "date", "stage_m", "water_depth_m", "bed_elev_m"]]; checks["revisit_pair_bed_elevation"] = pair.round(3).to_dict("records")
step = loc.sort_values(["transect", "across_m"]).groupby("transect").bed_elev_m.diff().abs(); checks["bed_elev_step_between_neighbours_m"] = {"median": float(step.median()), "90th": float(step.quantile(.9))}
json.dump(checks, open(D / "amplitude_checks.json", "w"), indent=2, default=str); print(json.dumps(checks, indent=2, default=str))
# ---------------------------------------------------------------- figure of checks
plt.rcParams.update({"font.size": 8}); fig = plt.figure(figsize=(14, 9)); gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.3)
ax = fig.add_subplot(gs[0, 0]); items = [(f"substrate: {k} ({v[3]})", v) for k, v in checks["z0_by_substrate"].items()] + [(f"year {k} ({v[3]})", v) for k, v in checks["z0_by_year"].items()] + [(f"crew {k[:14]} ({v[3]})", v) for k, v in checks["z0_by_crew"].items()]
for i, (lab, v) in enumerate(items): ax.plot([v[1] * 1000, v[2] * 1000], [i, i], color=BLUE, lw=3, alpha=0.4); ax.plot(v[0] * 1000, i, "o", color=BLUE)
ax.axvline(z0 * 1000, color=INK, lw=0.8, ls="--"); ax.set_yticks(range(len(items))); ax.set_yticklabels([l for l, _ in items], fontsize=7); ax.set_xscale("log"); ax.set_xlabel("fitted z0 (mm), dot = best, bar = within 5% of the minimum"); ax.set_title("Does the roughness scale differ by group? (n locations)", fontsize=9); ax.grid(True); ax.set_axisbelow(True)
ax = fig.add_subplot(gs[0, 1]); ax.plot([0, 0.2], [0, 0.2], color=AXIS, lw=0.8); ax.scatter(cmpdf.lin, cmpdf.log, s=8, color=BLUE, alpha=0.6, lw=0); ax.set_xlabel("amplitude, least squares in linear space"); ax.set_ylabel("amplitude, least squares in log space"); ax.set_title(f"Fitting choice: median ratio {checks['log_vs_linear_amplitude']['median_ratio']:.2f}", fontsize=9); ax.grid(True); ax.set_axisbelow(True)
ax = fig.add_subplot(gs[0, 2]); mids = [iv.mid for iv in prof.index]; ax.plot(mids, prof["median"], "o-", color=BLUE); ax.axhline(0, color=INK, lw=0.6)
for m_, (_, row) in zip(mids, prof.iterrows()): ax.annotate(f"n={int(row['size'])}", (m_, row["median"]), xytext=(0, 6), textcoords="offset points", ha="center", fontsize=6.5, color=INK2)
ax.set_xlabel("relative height above the bed, z / h"); ax.set_ylabel("median residual, in units of a"); ax.set_title("Where the log law misses: residual by relative height", fontsize=9); ax.grid(True); ax.set_axisbelow(True)
ax = fig.add_subplot(gs[1, 0]); ax.scatter(loc.across_m, loc.bed_elev_m, s=6, color=AXIS, lw=0, label="all locations"); ax.scatter(bad.across_m, bad.bed_elev_m, s=14, color=RED, lw=0, label=f"fit worse than flat ({len(bad)})"); ax.set_xlabel("across (m)"); ax.set_ylabel("bed elevation (m, stage datum)"); ax.set_title("The poorly fitted profiles: where in the channel", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7)
ax = fig.add_subplot(gs[1, 1])
for t, p in loc.sort_values("across_m").groupby("transect"): ax.plot(p.across_m, p.bed_elev_m, "-", lw=0.7, alpha=0.7, color=BLUE if p.year.iloc[0] == 2025 else GREEN)
for _, r in pair.iterrows(): ax.plot(loc[loc.point == r.point].across_m.iloc[0], r.bed_elev_m, "s", color=ORANGE, ms=5)
ax.set_xlabel("across (m)"); ax.set_ylabel("bed elevation (m, stage datum)"); ax.set_title("Bed elevation along each transect (orange squares: the revisited pair, two days)", fontsize=9); ax.grid(True); ax.set_axisbelow(True)
ax = fig.add_subplot(gs[1, 2]); ad = pd.read_csv(D / "adcp_amplitude_2h.csv"); m = ad[ad.local_depth_m > 0.02]; info = json.load(open(D / "amplitude_inputs.json"))["adcp_depth_scaling"]
ax.scatter(m.local_depth_m, m.velocity_ms, s=4, color=INK, alpha=0.3, lw=0); hh = np.linspace(m.local_depth_m.min(), m.local_depth_m.max(), 50)
for b, col, lab in [(0.5, GREEN, "1/2"), (2 / 3, BLUE, "2/3"), (info["exponent_free"], ORANGE, f"free {info['exponent_free']:.2f}")]:
    c = np.exp(np.mean(np.log(m.velocity_ms) - b * np.log(m.local_depth_m))); ax.plot(hh, c * hh ** b, color=col, lw=1.2, label=f"exponent {lab}")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("local depth at the ADCP (m)"); ax.set_ylabel("ADCP velocity (m/s)"); ax.set_title("ADCP: velocity against local depth, log axes", fontsize=9); ax.grid(True, which="both"); ax.set_axisbelow(True); ax.legend(fontsize=7)
fig.suptitle("Interrogating the amplitude method", fontsize=9, y=0.95); fig.savefig(F / "00o_amplitude_checks.png", bbox_inches="tight"); plt.close(fig); print("saved 00o")
