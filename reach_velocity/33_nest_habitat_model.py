"""The ecological number from the nests instead of a fixed threshold (the PI, 22 September 2026). Three things:
(1) nest-calibrated thresholds: the quantiles of the modelled cell velocity at the 59 nests (each on its day), with day-cluster bootstrap
    intervals, and the reach share below each of them against stage with the chain's band;
(2) a used-against-available selection model (logistic) on velocity, bank distance and depth, with day-cluster bootstrap intervals, to
    see what the nests support once bank distance is in the model;
(3) the probability that the reach share exceeds a given level, straight from the draws.
Outputs: nest_habitat_model.json, nest_thresholds_curves.csv, figures/10_nest_habitat.png"""
import json, numpy as np, pandas as pd
from scipy.spatial import cKDTree
from sklearn.linear_model import LogisticRegression
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from gpcore import *; from uncertainty import *
rng = np.random.default_rng(0); Z0 = json.load(open(D / "field_qc.json"))["z0_m"]; H0, GAM = PRIMARY   # the primary stage law (uncertainty.law_set, entry 30)
df = load_points(); df["bed_m"] = df.bed_rtk_sensor_datum_m.fillna(df.bed_elev_m); ch = Chain(df, banks=banks_for("Middle"))
grid = pd.read_csv(D / "grid.csv"); nests = pd.read_csv(D / "step7_nests.csv"); nests = nests.dropna(subset=["v_pred_ms"]).reset_index(drop=True); tree = cKDTree(grid[["x_m", "y_m"]].values); cell = grid.iloc[tree.query(nests[["x_m", "y_m"]].values)[1]].reset_index(drop=True); cix = tree.query(nests[["x_m", "y_m"]].values)[1]; GLV = GridLevels()   # the field and bank distance at each nest's own stage (bank set C, entry 30)
def vbar(logA, bed, st):
    h = np.clip(st - bed, 0, None); a = np.exp(logA) * np.clip(st - H0, 0.02, None) ** GAM; return np.where(h > np.e * Z0, a * (np.log(np.clip(h, np.e * Z0, None) / Z0) - 1), np.nan), h
# ---- (1) nest-calibrated thresholds with day-cluster bootstrap
days = nests.date.values; udays = np.unique(days); qs = [0.5, 0.75, 0.9]
def quant(v): return {f"q{int(q*100)}": float(np.quantile(v, q)) for q in qs}
point = quant(nests.v_pred_ms); boots = {k: [] for k in point}
for b in range(2000):
    pick = rng.choice(udays, len(udays), replace=True); v = np.concatenate([nests.v_pred_ms[days == d].values for d in pick]); qb = quant(v)
    for k in boots: boots[k].append(qb[k])
thr = {k: {"value": point[k], "ci90": [float(np.quantile(boots[k], .05)), float(np.quantile(boots[k], .95))]} for k in point}
# the reach share below each nest-calibrated threshold, against stage, from the joint chain
sec_pts = {t: g.assign(_x=g.x_m.round(2), _y=g.y_m.round(2)).drop_duplicates(["_x", "_y"]).sort_values("across_m").reset_index(drop=True) for t, g in df.groupby("transect")}
grid_st = np.linspace(0.65, 1.35, 15); rows = []; PR = {}
for k, t_ in thr.items():
    R = sample_reach(ch, grid_st, sec_pts, n=300, seed=5, threshold=t_["value"]); rs = R["reach_share"]
    for j, st in enumerate(grid_st): rows.append({"threshold_name": k, "threshold_m_s": t_["value"], "stage_m": st, "share_med": float(np.nanmedian(rs[:, j])), "share_lo": float(np.nanquantile(rs[:, j], .05)), "share_hi": float(np.nanquantile(rs[:, j], .95)), "P_share_above_0.3": float(np.nanmean(rs[:, j] > 0.3)), "P_share_above_0.5": float(np.nanmean(rs[:, j] > 0.5))})
    PR[k] = rs
curves = pd.DataFrame(rows); curves.to_csv(D / "nest_thresholds_curves.csv", index=False)
# ---- (2) used against available: logistic selection model with velocity, bank distance and depth; availability sampled at each nest's own stage
X_used = []; X_avail = []; cl_used = []; cl_avail = []
for i, r in nests.iterrows():
    st = r.stage_use; gA, gD = GLV.at(st); vg, hg = vbar(gA, grid.bed_m.values, st); ok = ~np.isnan(vg); c = cell.iloc[i]; vc, hc = vbar(gA[cix[i]], c.bed_m, st)
    X_used.append([float(vc), gD[cix[i]], float(hc)]); cl_used.append(r.date); pick = rng.choice(np.where(ok)[0], 20, replace=False); X_avail += [[vg[j], gD[j], hg[j]] for j in pick]; cl_avail += [r.date] * 20
X = np.array(X_used + X_avail); y = np.r_[np.ones(len(X_used)), np.zeros(len(X_avail))]; cl = np.array(cl_used + cl_avail); names = ["velocity (m/s)", "bank distance (m)", "depth (m)"]
mu_, sd_ = X.mean(0), X.std(0); Xs = (X - mu_) / sd_
def fit(Xs_, y_, cols): m = LogisticRegression(C=10.0, max_iter=1000).fit(Xs_[:, cols], y_); return m.coef_[0]
specs = {"velocity alone": [0], "bank distance alone": [1], "velocity + bank distance": [0, 1], "velocity + bank distance + depth": [0, 1, 2]}; sel = {}
for name, cols in specs.items():
    b0 = fit(Xs, y, cols); B = []
    for b in range(500):
        pick = rng.choice(udays, len(udays), replace=True); idx = np.concatenate([np.where(cl == d)[0] for d in pick]); B.append(fit(Xs[idx], y[idx], cols))
    B = np.array(B); sel[name] = {names[c]: {"coef_per_sd": float(b0[k]), "ci90": [float(np.quantile(B[:, k], .05)), float(np.quantile(B[:, k], .95))], "per_unit": float(b0[k] / sd_[c])} for k, c in enumerate(cols)}
out = {"n_nests": int(len(nests)), "nest_calibrated_thresholds_m_s": thr, "share_of_nests_below_0.2": float(np.mean(nests.v_pred_ms < 0.2)), "selection_model": sel, "note": "coefficients are log-odds of use per standard deviation of the covariate (sd: " + ", ".join(f"{n} {s:.2f}" for n, s in zip(names, sd_)) + "); 20 available cells per nest at the nest's own stage; day-cluster bootstrap, 9 days"}
json.dump(out, open(D / "nest_habitat_model.json", "w"), indent=1); print(json.dumps(out, indent=1))
# ---- figure
INK, AXIS, GRID, BLUE, GREEN, DARK, MUTED, ORANGE = "#0b0b0b", "#c3c2b7", "#e1e0d9", "#2a78d6", "#1a9e74", "#0d366b", "#898781", "#eb6834"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig, axes = plt.subplots(1, 3, figsize=(16, 4.4))
ax = axes[0]; vs = np.sort(nests.v_pred_ms); ax.step(vs, np.arange(1, len(vs) + 1) / len(vs), color=GREEN, lw=1.4, label="modelled cell velocity at the nests")
for k, c in zip(["q50", "q75", "q90"], [BLUE, DARK, INK]): ax.axvline(thr[k]["value"], color=c, lw=0.9, ls="--"); ax.axvspan(*thr[k]["ci90"], color=c, alpha=0.12, lw=0); ax.text(thr[k]["value"] + 0.01, 0.08 + 0.1 * ["q50", "q75", "q90"].index(k), f"{k}: {thr[k]['value']:.2f} [{thr[k]['ci90'][0]:.2f}, {thr[k]['ci90'][1]:.2f}]", fontsize=7, color=c)
ax.axvline(0.2, color=ORANGE, lw=0.9); ax.text(0.205, 0.95, "0.2 m/s", color=ORANGE, fontsize=7); ax.set_xlabel("depth-averaged velocity in the nest's cell on its day (m/s)"); ax.set_ylabel("fraction of nests below"); ax.set_title("Nest-calibrated thresholds: the velocity quantiles of used cells\n(bands: 90% day-cluster bootstrap over 9 days)", fontsize=8.5); ax.grid(True); ax.set_axisbelow(True); ax.set_xlim(0, 0.7)
ax = axes[1]
for k, c in zip(["q50", "q75", "q90"], [BLUE, DARK, INK]): q = curves[curves.threshold_name == k]; ax.fill_between(q.stage_m, q.share_lo * 100, q.share_hi * 100, color=c, alpha=0.12, lw=0); ax.plot(q.stage_m, q.share_med * 100, color=c, lw=1.3, label=f"below {thr[k]['value']:.2f} m/s (nest {k})")
ax.axvspan(0.74, 1.09, color=GRID, alpha=0.5, lw=0); ax.set_xlabel("stage (m)"); ax.set_ylabel("share of the reach's wetted width (%)"); ax.set_title("The reach share below each nest-calibrated threshold, 90% bands", fontsize=8.5); ax.legend(fontsize=7); ax.grid(True); ax.set_axisbelow(True); ax.set_ylim(0, 100)
ax = axes[2]; specs_show = ["velocity alone", "bank distance alone", "velocity + bank distance", "velocity + bank distance + depth"]; yy = 0; ticks = []
for name in specs_show:
    for cov, c in [("velocity (m/s)", GREEN), ("bank distance (m)", MUTED), ("depth (m)", BLUE)]:
        if cov in sel[name]: v_ = sel[name][cov]; ax.errorbar(v_["coef_per_sd"], yy, xerr=[[v_["coef_per_sd"] - v_["ci90"][0]], [v_["ci90"][1] - v_["coef_per_sd"]]], fmt="o", color=c, ms=5, elinewidth=1); ticks.append((yy, f"{name}: {cov.split(' (')[0]}")); yy += 1
    yy += 0.5
ax.axvline(0, color=INK, lw=0.6); ax.set_yticks([t for t, _ in ticks]); ax.set_yticklabels([l for _, l in ticks], fontsize=7); ax.invert_yaxis(); ax.set_xlabel("log-odds of a nest per sd of the covariate (90% day-cluster bootstrap)"); ax.set_title("Used against available: what the nests select on", fontsize=8.5); ax.grid(True, axis="x"); ax.set_axisbelow(True)
fig.suptitle("The ecological number from the nests: calibrated thresholds, their reach shares, and the selection model", fontsize=10, y=1.03); fig.savefig(F / "10_nest_habitat.png", bbox_inches="tight"); plt.close(fig); print("saved 10")
