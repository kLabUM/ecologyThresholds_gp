"""Step 2. The conveyance field. GP on log A over the surveyed points with the agreed inputs, fixed per-point noise plus a
learned nugget, a mean adjustment for substrate and wake notes. Baselines: the reach mean; the generic cross-channel shape
(A against distance to bank only); coordinates only (along, relative across). Validation: leave-one-transect-out (the
spatial transfer that matters), random 10-fold (the cautionary comparison), and the single day 25 June 2026 alone.
Outputs: step2_metrics.csv, step2_predictions.csv, step2_kernels.json, figures/02_step2_conveyance_gp.png"""
import json, sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import norm
from sklearn.model_selection import KFold
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel as C, Matern, WhiteKernel
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"; F = HERE / "figures"
t0 = time.time(); log = lambda *a: (print(f"[{time.time()-t0:5.0f}s]", *a), sys.stdout.flush())
BETA, SD_H, N_RESTARTS = 0.5, 0.06, 3
fe = pd.read_csv(D / "features.csv"); lc = pd.read_csv(D / "locations.csv")[["location_id", "a", "water_depth_m", "x_m", "y_m", "across_m"]]
df = fe.drop(columns=["across_m"]).merge(lc, on="location_id"); df = df[df.a > 0].reset_index(drop=True)
df["A"] = df.a / df.water_depth_m ** BETA; df["y"] = np.log(df.A); df["sd_y"] = np.clip(np.sqrt((df.a_sd / df.a) ** 2 + (BETA * SD_H / df.water_depth_m) ** 2), 0.02, 1.0)
SUBS = ["cobble", "gravel", "sand", "silt", "vegetation", "woody debris"]
for s_ in SUBS: df[f"sub_{s_}"] = (df.substrate == s_).astype(float)
NOTE_COLS = [f"sub_{s_}" for s_ in SUBS] + ["wake_flag"]
FEATS = {"coords": ["along_m", "rel_across"], "coords+bank": ["along_m", "rel_across", "dist_bank_m"], "full": ["along_m", "rel_across", "dist_bank_m", "rel_depth", "section_mean_depth_m", "lateral_bed_slope_abs"]}
log(f"{len(df)} points; median sd of log A {df.sd_y.median():.2f}; substrate noted at {int((df.substrate != 'not noted').sum())}, wake flag at {int(df.wake_flag.sum())}")
def wls(Xb, y, w):
    Wt = np.sqrt(w)[:, None]; b, *_ = np.linalg.lstsq(Xb * Wt, y * Wt[:, 0], rcond=None); return b
class Model:
    """mean = WLS on [1, notes] (optional); GP with ARD Matern 3/2 on the features (optional); fixed per-point noise + nugget"""
    def __init__(self, feats=None, notes=False, name=""): self.feats, self.notes, self.name = feats, notes, name
    def design(self, d): return np.column_stack([np.ones(len(d))] + ([d[c].values for c in NOTE_COLS] if self.notes else []))
    def fit(self, d):
        w = 1 / d.sd_y.values ** 2; self.b = wls(self.design(d), d.y.values, w); r = d.y.values - self.design(d) @ self.b
        if self.feats is None: self.rsd = float(np.sqrt(np.average(r ** 2, weights=w))); return self
        X = d[self.feats].values.astype(float); self.mx, self.sx = X.mean(0), X.std(0) + 1e-9; self.ys = r.std() + 1e-9
        k = C(1.0, (1e-2, 1e2)) * Matern(np.ones(len(self.feats)), [(0.05, 50.0)] * len(self.feats), nu=1.5) + WhiteKernel(0.1, (1e-4, 4.0))
        self.gp = GaussianProcessRegressor(k, alpha=(d.sd_y.values / self.ys) ** 2, n_restarts_optimizer=N_RESTARTS, random_state=0).fit((X - self.mx) / self.sx, r / self.ys); return self
    def predict(self, d):
        mu = self.design(d) @ self.b
        if self.feats is None: return mu, np.full(len(d), self.rsd)
        m, sd = self.gp.predict((d[self.feats].values.astype(float) - self.mx) / self.sx, return_std=True); return mu + m * self.ys, sd * self.ys
    def scales(self):
        if self.feats is None: return {}
        k = self.gp.kernel_; ls = np.atleast_1d(k.k1.k2.length_scale) * self.sx; return {"amp_sd": float(np.sqrt(k.k1.k1.constant_value) * self.ys), "nugget_sd": float(np.sqrt(k.k2.noise_level) * self.ys), **{f: float(l) for f, l in zip(self.feats, ls)}}
MODELS = {"reach mean": Model(None, False), "reach mean + notes": Model(None, True), "generic shape (bank distance)": Model(["dist_bank_m"], False), "coords": Model(FEATS["coords"], False), "coords + bank": Model(FEATS["coords+bank"], False), "full features": Model(FEATS["full"], False), "full features + notes": Model(FEATS["full"], True)}
def scores(y, mu, sd):
    z = (y - mu) / sd; crps = sd * (z * (2 * norm.cdf(z) - 1) + 2 * norm.pdf(z) - 1 / np.sqrt(np.pi))
    return {"MAE_log": float(np.mean(np.abs(y - mu))), "CRPS_log": float(np.mean(crps)), "coverage95": float(np.mean(np.abs(z) <= 1.96)), "mean_sd": float(np.mean(sd)), "MAE_A": float(np.mean(np.abs(np.exp(y) - np.exp(mu)))), "n": int(len(y))}
schemes = {"transect_out": [(np.where(df.transect != t)[0], np.where(df.transect == t)[0]) for t in sorted(df.transect.unique())], "random_10fold": list(KFold(10, shuffle=True, random_state=0).split(df))}
day = df[df.date == "2026-06-25"].reset_index(drop=True); schemes_day = {"june25_transect_out": [(np.where(day.transect != t)[0], np.where(day.transect == t)[0]) for t in sorted(day.transect.unique())]}
rows, preds, kern = [], [], {}
for sname, folds in list(schemes.items()) + list(schemes_day.items()):
    data = day if sname.startswith("june25") else df
    for mname, mdl in MODELS.items():
        if sname.startswith("june25") and "along" in (mdl.feats or []) and mname != "full features": pass
        MU, SD = np.full(len(data), np.nan), np.full(len(data), np.nan)
        for tr, te in folds:
            m = Model(mdl.feats, mdl.notes).fit(data.iloc[tr]); mu, sd = m.predict(data.iloc[te]); MU[te], SD[te] = mu, sd
        sd_tot = np.sqrt(SD ** 2 + data.sd_y.values ** 2)   # the held-out value carries its own measurement noise
        sc = scores(data.y.values, MU, sd_tot); rows.append({"scheme": sname, "model": mname, **sc}); log(f"{sname:20s} {mname:30s} MAE_log {sc['MAE_log']:.3f} CRPS_log {sc['CRPS_log']:.3f} cov {sc['coverage95']:.2f} MAE_A {sc['MAE_A']:.3f}")
        for i in range(len(data)): preds.append((sname, mname, data.location_id.iat[i], data.transect.iat[i], float(data.y.iat[i]), float(MU[i]), float(sd_tot[i])))
met = pd.DataFrame(rows); met.to_csv(D / "step2_metrics.csv", index=False); pd.DataFrame(preds, columns=["scheme", "model", "location_id", "transect", "y", "mu", "sd"]).to_csv(D / "step2_predictions.csv", index=False)
for mname in ["coords", "full features", "full features + notes"]: m = Model(MODELS[mname].feats, MODELS[mname].notes).fit(df); kern[mname] = {"scales": m.scales(), "kernel": str(m.gp.kernel_), "mean_coefs": dict(zip(["intercept"] + (NOTE_COLS if m.notes else []), [float(b) for b in m.b]))}; log(f"full fit {mname}: {m.scales()}")
mday = Model(FEATS["full"], False).fit(day); kern["june25 full features"] = {"scales": mday.scales(), "kernel": str(mday.gp.kernel_)}; log(f"june25 fit: {mday.scales()}")
json.dump(kern, open(D / "step2_kernels.json", "w"), indent=2)
# ---------------------------------------------------------------- figure
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, DARK, MUTED = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#0d366b", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig = plt.figure(figsize=(15, 10)); gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.3)
ax = fig.add_subplot(gs[0, 0]); order = list(MODELS); x = np.arange(len(order)); w = 0.27
for k, (sname, col) in enumerate([("transect_out", BLUE), ("random_10fold", MUTED), ("june25_transect_out", GREEN)]):
    v = [met[(met.scheme == sname) & (met.model == mn)].CRPS_log.iloc[0] for mn in order]; ax.bar(x + (k - 1) * w, v, w, color=col, label={"transect_out": "transect held out (all days)", "random_10fold": "random folds", "june25_transect_out": "25 June 2026 only, transect held out"}[sname])
ax.set_xticks(x); ax.set_xticklabels(order, rotation=25, ha="right", fontsize=7); ax.set_ylabel("CRPS of log A (lower is better)"); ax.set_title("Held-out score by model and scheme", fontsize=9); ax.grid(True, axis="y"); ax.set_axisbelow(True); ax.legend(fontsize=7)
ax = fig.add_subplot(gs[0, 1]); sc_ = kern["full features + notes"]["scales"]; names = FEATS["full"]; vals = [sc_[f] for f in names]; ranges = [float(df[f].max() - df[f].min()) for f in names]
ax.barh(range(len(names)), [v / r for v, r in zip(vals, ranges)], color=BLUE); ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=7); ax.set_xscale("log"); ax.axvline(1, color=INK, lw=0.6, ls="--"); ax.set_xlabel("fitted length scale / range of the feature (log axis; > 1 means the feature barely matters)"); ax.set_title(f"ARD length scales, full fit (amplitude sd {sc_['amp_sd']:.2f}, nugget sd {sc_['nugget_sd']:.2f} in log A)", fontsize=9); ax.grid(True, axis="x"); ax.set_axisbelow(True)
ax = fig.add_subplot(gs[0, 2]); pp = pd.read_csv(D / "step2_predictions.csv"); q = pp[(pp.scheme == "transect_out") & (pp.model == "full features + notes")]
ax.errorbar(np.exp(q.y), np.exp(q.mu), yerr=[np.exp(q.mu) - np.exp(q.mu - 1.96 * q.sd), np.exp(q.mu + 1.96 * q.sd) - np.exp(q.mu)], fmt="o", color=BLUE, ms=3, elinewidth=0.4, alpha=0.6); ax.plot([0, 0.35], [0, 0.35], color=AXIS, lw=0.8)
sc2 = met[(met.scheme == "transect_out") & (met.model == "full features + notes")].iloc[0]; ax.set_xlabel("measured A (transect held out)"); ax.set_ylabel("predicted A with 95% band"); ax.set_title(f"Transect held out, full features + notes: MAE {sc2.MAE_log*100:.0f}% in log, coverage {sc2.coverage95:.2f}", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.set_xlim(0, 0.35); ax.set_ylim(0, 0.35)
for k, t in enumerate(["AT10", "0625-T02", "AT0401"]):
    ax = fig.add_subplot(gs[1, k]); p = df[df.transect == t].sort_values("across_m"); qq = q.set_index("location_id").loc[p.location_id]
    ax.errorbar(p.across_m, p.A, yerr=1.96 * p.A * p.sd_y, fmt="s", color=INK, ms=4, elinewidth=0.6, capsize=1.5, label="measured, with its own error"); ax.fill_between(p.across_m, np.exp(qq.mu - 1.96 * qq.sd), np.exp(qq.mu + 1.96 * qq.sd), color=BLUE, alpha=0.15, lw=0); ax.plot(p.across_m, np.exp(qq.mu), color=BLUE, lw=1.4, label="predicted with this transect held out")
    qm = pp[(pp.scheme == "transect_out") & (pp.model == "reach mean")].set_index("location_id").loc[p.location_id]; ax.plot(p.across_m, np.exp(qm.mu), color=MUTED, lw=1.0, ls="--", label="reach mean")
    ax.set_xlabel("across (m)"); ax.set_ylabel("conveyance A"); ax.set_title(f"{t} ({p.date.iloc[0]}), predicted from the other transects", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.set_ylim(0, 0.4)
    if k == 0: ax.legend(fontsize=7)
fig.suptitle("Step 2: the conveyance field, interpolated between transects", fontsize=10, y=0.96); fig.savefig(F / "02_step2_conveyance_gp.png", bbox_inches="tight"); plt.close(fig); log("saved 02")
print(met.round(3).to_string(index=False))
