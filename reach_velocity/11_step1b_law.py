"""Step 1b. The GP with a physical law as its mean. At a point in the channel, Manning's relation gives velocity rising
with the local depth as a power: v = c * (stage - h0)^b, where h0 is the stage at which that point runs dry and b is 2/3
for a wide channel. Three laws are fitted by least squares inside every fold: linear (b = 1 with a free intercept),
Manning with b fixed at 2/3, and the power law with b free. A Matern 3/2 GP is then fitted to the residuals, so the
prediction is law(stage) + GP(stage) and beyond the data it follows the law. Same validation and scores as step 1.
Outputs: data_clean/step1b_metrics.csv, step1b_laws.json, step1b_2026_prediction.csv, figures/01b_step1_law.png"""
import json, sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import norm
from sklearn.model_selection import KFold
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel as C, Matern, WhiteKernel
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import matplotlib.dates as mdates
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"; F = HERE / "figures"
t0 = time.time(); log = lambda *a: (print(f"[{time.time()-t0:5.0f}s]", *a), sys.stdout.flush())
g = pd.read_csv(D / "series_10min.csv", parse_dates=["time_utc"]).set_index("time_utc"); g.index = pd.to_datetime(g.index, utc=True)
s = g.dropna(subset=["adcp_ms", "stage_m"]); h = s[["stage_m", "adcp_ms"]].resample("2h").mean().dropna(); X, y = h.stage_m.values, h.adcp_ms.values; tdays = ((h.index - h.index[0]).total_seconds() / 86400).values
LAWS = {"linear": (lambda x, c, h0: c * (x - h0), [0.3, 0.5], ([0, -5], [10, 0.66])),
        "manning_2_3": (lambda x, c, h0: c * np.clip(x - h0, 1e-6, None) ** (2 / 3), [0.3, 0.3], ([0, -5], [10, 0.66])),
        "power_free": (lambda x, c, h0, b: c * np.clip(x - h0, 1e-6, None) ** b, [0.3, 0.3, 1.0], ([0, -5, 0.2], [10, 0.66, 3.0]))}
class LawGP:
    def __init__(self, law, gp=True): self.law = law; self.use_gp = gp
    def fit(self, X, y):
        f, p0, bounds = LAWS[self.law]; self.p, _ = curve_fit(f, X, y, p0=p0, bounds=bounds, maxfev=20000); r = y - f(X, *self.p); self.rsd = float(r.std())
        if self.use_gp:
            self.mx, self.sx = X.mean(), X.std(); self.gp = GaussianProcessRegressor(C(1.0, (1e-3, 1e2)) * Matern(1.0, (0.05, 50.0), nu=1.5) + WhiteKernel(0.5, (1e-6, 2.0)), n_restarts_optimizer=1, random_state=0).fit(((X - self.mx) / self.sx)[:, None], r / self.rsd)
        return self
    def predict(self, X):
        f = LAWS[self.law][0]; mu = f(X, *self.p)
        if self.use_gp: m, sd = self.gp.predict(((X - self.mx) / self.sx)[:, None], return_std=True); return mu + m * self.rsd, sd * self.rsd
        return mu, np.full(len(X), self.rsd)
def scores(y, mu, sd):
    z = (y - mu) / sd; crps = sd * (z * (2 * norm.cdf(z) - 1) + 2 * norm.pdf(z) - 1 / np.sqrt(np.pi))
    return {"MAE": float(np.mean(np.abs(y - mu))), "RMSE": float(np.sqrt(np.mean((y - mu) ** 2))), "CRPS": float(np.mean(crps)), "coverage95": float(np.mean(np.abs(z) <= 1.96)), "mean_sd": float(np.mean(sd)), "NLPD": float(np.mean(0.5 * np.log(2 * np.pi * sd ** 2) + 0.5 * z ** 2)), "R2": float(1 - np.sum((y - mu) ** 2) / np.sum((y - y.mean()) ** 2)), "n": int(len(y))}
block = (tdays // 10).astype(int); q85, q15 = np.quantile(X, 0.85), np.quantile(X, 0.15)
schemes = {"blocked_10day": [(np.where(block != b)[0], np.where(block == b)[0]) for b in np.unique(block)], "stage_top15_heldout": [(np.where(X <= q85)[0], np.where(X > q85)[0])], "stage_bottom15_heldout": [(np.where(X >= q15)[0], np.where(X < q15)[0])], "random_10fold": list(KFold(10, shuffle=True, random_state=0).split(X))}
models = [("linear", False), ("linear", True), ("manning_2_3", False), ("manning_2_3", True), ("power_free", False), ("power_free", True)]
rows = []; params = {}
for sname, folds in schemes.items():
    for law, use_gp in models:
        MU, SD = np.full(len(y), np.nan), np.full(len(y), np.nan); IDX = np.zeros(len(y), bool)
        for tr, te in folds:
            m = LawGP(law, use_gp).fit(X[tr], y[tr]); MU[te], SD[te] = m.predict(X[te]); IDX[te] = True
        sc = scores(y[IDX], MU[IDX], SD[IDX]); name = law + (" + GP" if use_gp else ""); rows.append({"scheme": sname, "model": name, **sc}); log(f"{sname:24s} {name:18s} MAE {sc['MAE']*1000:5.1f}  CRPS {sc['CRPS']*1000:5.1f}  cov {sc['coverage95']:.2f}")
met = pd.DataFrame(rows); met.to_csv(D / "step1b_metrics.csv", index=False)
full = {}
for law in LAWS:
    m = LawGP(law, True).fit(X, y); full[law] = m; params[law] = {"params": [float(v) for v in m.p], "residual_sd_m_s": m.rsd, "gp_kernel": str(m.gp.kernel_)}
    log(f"{law}: params {np.round(m.p, 4)} residual sd {m.rsd*1000:.1f} mm/s; GP {m.gp.kernel_}")
json.dump(params, open(D / "step1b_laws.json", "w"), indent=2)
grid = np.linspace(0.55, 1.75, 241); s26 = g[g.index.year == 2026].stage_m.resample("2h").mean().dropna(); p26 = pd.DataFrame({"stage_m": s26.values}, index=s26.index); curves = pd.DataFrame({"stage_m": grid})
for law, m in full.items():
    mu, sd = m.predict(grid); curves[f"{law}_mu"], curves[f"{law}_sd"] = mu, sd; mu, sd = m.predict(s26.values); p26[f"{law}_mu"], p26[f"{law}_sd"] = mu, sd
p26.index.name = "time_utc"; p26.to_csv(D / "step1b_2026_prediction.csv"); curves.to_csv(D / "step1b_fit_curves.csv", index=False)
old = pd.read_csv(D / "step1_fit_curves.csv"); old26 = pd.read_csv(D / "step1_2026_prediction.csv", parse_dates=["time_utc"]).set_index("time_utc")
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, MUTED, DARK = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#898781", "#0d366b"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
styles = {"linear": (GREEN, "linear law + GP"), "manning_2_3": (DARK, "Manning law (b = 2/3) + GP"), "power_free": (BLUE, f"power law + GP (fitted b = {params['power_free']['params'][2]:.2f})")}
fig = plt.figure(figsize=(13, 8.5)); gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1], hspace=0.42, wspace=0.22)
ax = fig.add_subplot(gs[0, 0]); ax.axvspan(X.min(), X.max(), color=GRID, alpha=0.5, lw=0); ax.scatter(X, y, s=5, color=INK, alpha=0.3, lw=0, label="2-hourly means, 2025")
ax.plot(old.stage_m, old["Matern32_mu"], color=MUTED, lw=1.2, ls="--", label="step 1: Matern 3/2, zero mean (falls back to the average)")
for law, (col, lab) in styles.items(): ax.fill_between(curves.stage_m, curves[f"{law}_mu"] - 1.96 * curves[f"{law}_sd"], curves[f"{law}_mu"] + 1.96 * curves[f"{law}_sd"], color=col, alpha=0.1, lw=0); ax.plot(curves.stage_m, curves[f"{law}_mu"], color=col, lw=1.4, label=lab)
ax.set_xlim(0.55, 1.75); ax.set_ylim(-0.05, 0.45); ax.set_xlabel("stage (m)"); ax.set_ylabel("velocity at the ADCP (m/s)"); ax.set_title("A law as the mean: what each one does beyond the data", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7, loc="upper left")
ax = fig.add_subplot(gs[0, 1]); mb = met[met.scheme != "random_10fold"]; order = [l + (" + GP" if g_ else "") for l, g_ in models]; x = np.arange(len(order)); w = 0.26
for k, (sname, col) in enumerate([("blocked_10day", BLUE), ("stage_top15_heldout", ORANGE), ("stage_bottom15_heldout", GREEN)]):
    v = [mb[(mb.scheme == sname) & (mb.model == mn)].CRPS.iloc[0] * 1000 for mn in order]; ax.bar(x + (k - 1) * w, v, w, color=col, label=sname.replace("_", " "))
ax.set_xticks(x); ax.set_xticklabels(order, rotation=25, ha="right", fontsize=7); ax.set_ylabel("CRPS (mm/s), lower is better"); ax.set_title("Held-out score by law, with and without the residual GP", fontsize=9); ax.grid(True, axis="y"); ax.set_axisbelow(True); ax.legend(fontsize=7)
ax = fig.add_subplot(gs[1, :]); ax.axhspan(X.min(), X.max(), color=GRID, alpha=0.5, lw=0); ax.plot(p26.index, p26.stage_m, color=MUTED, lw=0.7); ax.set_ylabel("2026 stage (m)", color=MUTED); ax.set_ylim(0.5, 1.75)
ax2 = ax.twinx(); ax2.plot(old26.index, old26["Matern32_mu"], color=MUTED, lw=1.0, ls="--", label="step 1: Matern 3/2, zero mean")
for law, (col, lab) in styles.items(): ax2.fill_between(p26.index, p26[f"{law}_mu"] - 1.96 * p26[f"{law}_sd"], p26[f"{law}_mu"] + 1.96 * p26[f"{law}_sd"], color=col, alpha=0.1, lw=0); ax2.plot(p26.index, p26[f"{law}_mu"], color=col, lw=1.1, label=lab)
ax2.set_ylabel("predicted velocity at the ADCP (m/s)"); ax2.set_ylim(-0.05, 0.5); ax2.legend(fontsize=7, loc="upper right")
ax.set_title("2026: velocity at the ADCP from the 2026 stage record, law + GP (grey band: the 2025 stage range)", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
fig.suptitle("Step 1b: the GP with a physical law as its mean", fontsize=9, y=0.97); fig.savefig(F / "01b_step1_law.png", bbox_inches="tight"); plt.close(fig); log("saved 01b")
print(met.round(4).to_string(index=False))
