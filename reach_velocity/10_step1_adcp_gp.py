"""Step 1. One-dimensional GP: velocity at the ADCP against stage (2025), done properly.
Training and evaluation grid: 2-hourly means of the smoothed series (series_10min.csv). Hyperparameters by marginal
likelihood, refitted inside every fold. Cross-validation: (a) contiguous 10-day blocks, leave one out; (b) stage bands:
hold out the highest 15% of stages, then the lowest 15% (extrapolation); (c) random 10-fold, shown only to see how much it
flatters. Models: GP with RBF, Matern 3/2, Matern 5/2, Matern 3/2 + linear kernel; baselines: linear and quadratic
regression on stage with a Gaussian residual. Scores: MAE, RMSE, CRPS, 95% coverage, mean predictive sd, NLPD, and R2
for reference. Then: thinning check (training on daily medians), residuals over time, hysteresis, and the 2026 prediction
from the 2026 stage record. Outputs: data_clean/step1_cv_metrics.csv, step1_cv_predictions.csv, step1_fit_curves.csv,
step1_2026_prediction.csv, step1_summary.json, figures/01_step1_adcp_gp.png"""
import json, sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import norm
from sklearn.model_selection import KFold
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel as C, RBF, Matern, DotProduct, WhiteKernel
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import matplotlib.dates as mdates
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"; F = HERE / "figures"
t0 = time.time(); log = lambda *a: (print(f"[{time.time()-t0:5.0f}s]", *a), sys.stdout.flush())
g = pd.read_csv(D / "series_10min.csv", parse_dates=["time_utc"]).set_index("time_utc"); g.index = pd.to_datetime(g.index, utc=True)
s = g.dropna(subset=["adcp_ms", "stage_m"]); h = s[["stage_m", "adcp_ms"]].resample("2h").mean().dropna()
X, y = h.stage_m.values[:, None], h.adcp_ms.values; tdays = ((h.index - h.index[0]).total_seconds() / 86400).values
log(f"training grid: {len(h)} two-hourly points, {h.index[0].date()} to {h.index[-1].date()}, stage {X.min():.2f} to {X.max():.2f} m")
def kernels():
    return {"RBF": C(1.0, (1e-2, 1e2)) * RBF(1.0, (0.05, 50.0)) + WhiteKernel(1e-2, (1e-6, 1.0)),
            "Matern32": C(1.0, (1e-2, 1e2)) * Matern(1.0, (0.05, 50.0), nu=1.5) + WhiteKernel(1e-2, (1e-6, 1.0)),
            "Matern52": C(1.0, (1e-2, 1e2)) * Matern(1.0, (0.05, 50.0), nu=2.5) + WhiteKernel(1e-2, (1e-6, 1.0)),
            "Matern32+linear": C(1.0, (1e-2, 1e2)) * Matern(1.0, (0.05, 50.0), nu=1.5) + C(1.0, (1e-3, 1e2)) * DotProduct(1.0, (1e-3, 1e2)) + WhiteKernel(1e-2, (1e-6, 1.0))}
class GPfit:
    def __init__(self, name): self.name = name
    def fit(self, X, y):
        self.mx, self.sx, self.my, self.sy = X.mean(), X.std(), y.mean(), y.std()
        self.gp = GaussianProcessRegressor(kernels()[self.name], n_restarts_optimizer=2, random_state=0).fit((X - self.mx) / self.sx, (y - self.my) / self.sy); return self
    def predict(self, X):
        mu, sd = self.gp.predict((X - self.mx) / self.sx, return_std=True); return mu * self.sy + self.my, sd * self.sy
    def kernel(self): return str(self.gp.kernel_)
class Poly:
    def __init__(self, deg): self.deg = deg; self.name = {1: "linear", 2: "quadratic"}[deg]
    def fit(self, X, y): self.b = np.polyfit(X.ravel(), y, self.deg); self.sd = float(np.std(y - np.polyval(self.b, X.ravel()))); return self
    def predict(self, X): return np.polyval(self.b, X.ravel()), np.full(len(X), self.sd)
def scores(y, mu, sd):
    z = (y - mu) / sd; crps = sd * (z * (2 * norm.cdf(z) - 1) + 2 * norm.pdf(z) - 1 / np.sqrt(np.pi))
    return {"MAE": float(np.mean(np.abs(y - mu))), "RMSE": float(np.sqrt(np.mean((y - mu) ** 2))), "CRPS": float(np.mean(crps)), "coverage95": float(np.mean(np.abs(z) <= 1.96)), "mean_sd": float(np.mean(sd)),
            "NLPD": float(np.mean(0.5 * np.log(2 * np.pi * sd ** 2) + 0.5 * z ** 2)), "R2": float(1 - np.sum((y - mu) ** 2) / np.sum((y - y.mean()) ** 2)), "n": int(len(y))}
block = (tdays // 10).astype(int); q85, q15 = np.quantile(X, 0.85), np.quantile(X, 0.15)
schemes = {"blocked_10day": [(np.where(block != b)[0], np.where(block == b)[0]) for b in np.unique(block)],
           "stage_top15_heldout": [(np.where(X.ravel() <= q85)[0], np.where(X.ravel() > q85)[0])], "stage_bottom15_heldout": [(np.where(X.ravel() >= q15)[0], np.where(X.ravel() < q15)[0])],
           "random_10fold": list(KFold(10, shuffle=True, random_state=0).split(X))}
models = ["linear", "quadratic", "RBF", "Matern32", "Matern52", "Matern32+linear"]
rows, preds, kern = [], [], []
for sname, folds in schemes.items():
    for mname in models:
        MU, SD, IDX = np.full(len(y), np.nan), np.full(len(y), np.nan), np.zeros(len(y), bool)
        for fi, (tr, te) in enumerate(folds):
            m = (Poly(1) if mname == "linear" else Poly(2) if mname == "quadratic" else GPfit(mname)).fit(X[tr], y[tr]); mu, sd = m.predict(X[te]); MU[te], SD[te], IDX[te] = mu, sd, True
            if hasattr(m, "kernel"): kern.append({"scheme": sname, "model": mname, "fold": fi, "kernel": m.kernel()})
        sc = scores(y[IDX], MU[IDX], SD[IDX]); rows.append({"scheme": sname, "model": mname, **sc}); log(f"{sname:24s} {mname:16s} MAE {sc['MAE']*1000:5.1f} mm/s  CRPS {sc['CRPS']*1000:5.1f}  cov {sc['coverage95']:.2f}  sd {sc['mean_sd']*1000:4.1f}  R2 {sc['R2']:6.2f}")
        for i in np.where(IDX)[0]: preds.append((sname, mname, str(h.index[i]), float(X[i, 0]), float(y[i]), float(MU[i]), float(SD[i])))
met = pd.DataFrame(rows); met.to_csv(D / "step1_cv_metrics.csv", index=False); pd.DataFrame(preds, columns=["scheme", "model", "time_utc", "stage_m", "y_ms", "mu_ms", "sd_ms"]).to_csv(D / "step1_cv_predictions.csv", index=False)
best = met[(met.scheme == "blocked_10day") & met.model.str.contains("Matern|RBF")].sort_values("CRPS").model.iloc[0]; log(f"best GP kernel under blocked CV by CRPS: {best}")
# thinning: train on daily medians of the training blocks, test on the 2-hourly points of the held-out block
MU, SD = np.full(len(y), np.nan), np.full(len(y), np.nan); daily = h.resample("1D").median().dropna(); tday_d = ((daily.index - h.index[0]).total_seconds() / 86400).values
for b in np.unique(block):
    trd = np.where((tday_d // 10).astype(int) != b)[0]; te = np.where(block == b)[0]
    m = GPfit(best).fit(daily.stage_m.values[trd][:, None], daily.adcp_ms.values[trd]); MU[te], SD[te] = m.predict(X[te])
thin = scores(y, MU, SD); log(f"thinning check, {best} trained on daily medians ({len(daily)} points): MAE {thin['MAE']*1000:.1f} mm/s, CRPS {thin['CRPS']*1000:.1f}, cov {thin['coverage95']:.2f}")
# full fits, curves, residuals, hysteresis, 2026
full = {name: GPfit(name).fit(X, y) for name in ["Matern32", best]} if best != "Matern32" else {"Matern32": GPfit("Matern32").fit(X, y)}
grid = np.linspace(0.55, 1.75, 241)[:, None]; curves = pd.DataFrame({"stage_m": grid.ravel()})
for name, m in full.items(): mu, sd = m.predict(grid); curves[f"{name}_mu"], curves[f"{name}_sd"] = mu, sd
curves.to_csv(D / "step1_fit_curves.csv", index=False)
bp = pd.DataFrame(preds, columns=["scheme", "model", "time_utc", "stage_m", "y_ms", "mu_ms", "sd_ms"]); bp = bp[(bp.scheme == "blocked_10day") & (bp.model == best)].copy(); bp["time_utc"] = pd.to_datetime(bp.time_utc, utc=True); bp["resid"] = bp.y_ms - bp.mu_ms
dstage = h.stage_m.diff(3); rising = (dstage.reindex(bp.time_utc).values > 0.01); falling = (dstage.reindex(bp.time_utc).values < -0.01)
hyst = {"rising_mean_resid_mm_s": float(bp.resid[rising].mean() * 1000), "falling_mean_resid_mm_s": float(bp.resid[falling].mean() * 1000), "n_rising": int(rising.sum()), "n_falling": int(falling.sum())}
monthly = bp.groupby(bp.time_utc.dt.month).resid.mean() * 1000; hyst["monthly_mean_resid_mm_s"] = {int(k): float(v) for k, v in monthly.items()}
s26 = g[(g.index.year == 2026)].stage_m.resample("2h").mean().dropna(); p26 = pd.DataFrame({"stage_m": s26.values}, index=s26.index)
for name, m in full.items(): mu, sd = m.predict(s26.values[:, None]); p26[f"{name}_mu"], p26[f"{name}_sd"] = mu, sd
p26.index.name = "time_utc"; p26.to_csv(D / "step1_2026_prediction.csv")
summ = {"n_train_2h": int(len(h)), "stage_range_train_m": [float(X.min()), float(X.max())], "best_kernel_blocked_crps": best, "full_fit_kernels": {k: m.kernel() for k, m in full.items()},
        "thinning_daily_medians": thin, "hysteresis": hyst, "share_2026_stage_above_training_max": float((s26 > X.max()).mean()), "share_2026_stage_below_training_min": float((s26 < X.min()).mean())}
json.dump(summ, open(D / "step1_summary.json", "w"), indent=2); pd.DataFrame(kern).to_csv(D / "step1_fold_kernels.csv", index=False)
# figure
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, MUTED = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig = plt.figure(figsize=(13, 9)); gs = fig.add_gridspec(3, 2, height_ratios=[1.3, 1, 1.1], hspace=0.5, wspace=0.22)
ax = fig.add_subplot(gs[0, 0]); ax.axvspan(X.min(), X.max(), color=GRID, alpha=0.5, lw=0); ax.scatter(X, y, s=5, color=INK, alpha=0.35, lw=0, label="2-hourly means, 2025")
for name, col, ls in [("Matern32", ORANGE, "--"), (best, BLUE, "-")]:
    if name not in curves.columns.str.replace("_mu", "").values: continue
    ax.fill_between(curves.stage_m, curves[f"{name}_mu"] - 1.96 * curves[f"{name}_sd"], curves[f"{name}_mu"] + 1.96 * curves[f"{name}_sd"], color=col, alpha=0.12, lw=0); ax.plot(curves.stage_m, curves[f"{name}_mu"], color=col, lw=1.4, ls=ls, label=f"GP {name}, 95% band")
    if name == best and best == "Matern32": break
ax.set_xlim(0.55, 1.75); ax.set_ylim(-0.05, 0.4); ax.set_xlabel("stage (m)"); ax.set_ylabel("velocity at the ADCP (m/s)"); ax.set_title("The fit, and what each kernel does beyond the data (grey: the 2025 stage range)", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7, loc="upper left")
ax = fig.add_subplot(gs[0, 1]); mb = met[met.scheme != "random_10fold"]; order = models; x = np.arange(len(order)); w = 0.26
for k, (sname, col) in enumerate([("blocked_10day", BLUE), ("stage_top15_heldout", ORANGE), ("stage_bottom15_heldout", GREEN)]):
    v = [mb[(mb.scheme == sname) & (mb.model == mn)].CRPS.iloc[0] * 1000 for mn in order]; ax.bar(x + (k - 1) * w, v, w, color=col, label=sname.replace("_", " "))
ax.set_xticks(x); ax.set_xticklabels(order, rotation=20, ha="right", fontsize=7); ax.set_ylabel("CRPS (mm/s), lower is better"); ax.set_title("Held-out score by model and validation scheme", fontsize=9); ax.grid(True, axis="y"); ax.set_axisbelow(True); ax.legend(fontsize=7)
ax = fig.add_subplot(gs[1, :]); ax.plot(bp.time_utc, bp.resid * 1000, color=BLUE, lw=0.7); ax.axhline(0, color=INK, lw=0.5)
ax2 = ax.twinx(); ax2.plot(h.index, h.stage_m, color=MUTED, lw=0.7, alpha=0.8); ax2.set_ylabel("stage (m)", color=MUTED)
ax.set_ylabel("held-out residual (mm/s)"); ax.set_title(f"Residuals over the season, {best}, each 10-day block predicted from the others (grey: stage). Rising limbs {hyst['rising_mean_resid_mm_s']:+.1f} mm/s, falling {hyst['falling_mean_resid_mm_s']:+.1f} mm/s", fontsize=9); ax.grid(True); ax.set_axisbelow(True)
ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
ax = fig.add_subplot(gs[2, :]); ax.axhspan(X.min(), X.max(), color=GRID, alpha=0.5, lw=0); ax.plot(p26.index, p26.stage_m, color=MUTED, lw=0.7); ax.set_ylabel("2026 stage (m)", color=MUTED); ax.set_ylim(0.5, 1.75)
ax2 = ax.twinx()
for name, col, ls in [("Matern32", ORANGE, "--"), (best, BLUE, "-")]:
    if f"{name}_mu" not in p26.columns: continue
    ax2.fill_between(p26.index, p26[f"{name}_mu"] - 1.96 * p26[f"{name}_sd"], p26[f"{name}_mu"] + 1.96 * p26[f"{name}_sd"], color=col, alpha=0.12, lw=0); ax2.plot(p26.index, p26[f"{name}_mu"], color=col, lw=1.1, ls=ls, label=f"predicted velocity, {name}")
    if name == best and best == "Matern32": break
ax2.set_ylabel("predicted velocity at the ADCP (m/s)"); ax2.set_ylim(-0.05, 0.45); ax2.legend(fontsize=7, loc="upper right")
ax.set_title(f"2026: velocity at the ADCP predicted from the 2026 stage record (no ADCP that year). Grey band: the stage range the model was trained on; {summ['share_2026_stage_above_training_max']*100:.0f}% of 2026 lies above it", fontsize=9); ax.grid(True); ax.set_axisbelow(True)
ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
fig.suptitle("Step 1: one-dimensional GP, velocity at the ADCP against stage", fontsize=9, y=0.97); fig.savefig(F / "01_step1_adcp_gp.png", bbox_inches="tight"); plt.close(fig); log("saved 01_step1")
print(met.round(4).to_string(index=False))
