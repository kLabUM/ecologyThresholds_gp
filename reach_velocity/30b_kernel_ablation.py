"""Kernel family and optimiser-restart sensitivity of the place-factor GP (Middle, stage law), every section held out.
Matern 3/2 with one, three and six restarts, Matern 5/2, and the squared exponential. Output: step2_kernel_ablation.csv (Table 4)."""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from gpcore import *; from gpcore import _wls
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel as C, Matern, WhiteKernel, RBF
df = load_points(); d = with_target(df, ("stage", 1.0, 0.45)); rows = []
class ModelK(Model):
    def __init__(self, kind, **kw): super().__init__(**kw); self.kind = kind
    def fit(self, dd):
        w = 1 / dd.sd_y.values ** 2; self.b = _wls(self.design(dd), dd.y.values, w); r = dd.y.values - self.design(dd) @ self.b
        X = dd[self.feats].values.astype(float); self.mx, self.sx = X.mean(0), X.std(0) + 1e-9; self.ys = r.std() + 1e-9
        base = RBF(np.ones(len(self.feats)), [(0.05, 50.0)] * len(self.feats)) if self.kind == "rbf" else Matern(np.ones(len(self.feats)), [(0.05, 50.0)] * len(self.feats), nu={"m32": 1.5, "m52": 2.5}[self.kind])
        k = C(1.0, (1e-2, 1e2)) * base + WhiteKernel(0.1, (1e-4, 4.0)); self.gp = GaussianProcessRegressor(k, alpha=(dd.sd_y.values / self.ys) ** 2, n_restarts_optimizer=self.n_restarts, random_state=0).fit((X - self.mx) / self.sx, r / self.ys); return self
for kind, nr in [("m32", 1), ("m32", 3), ("m32", 6), ("m52", 3), ("rbf", 3)]:
    Y, MU, SD = [], [], []
    for t in sorted(d.transect.unique()):
        te = d[d.transect == t]; m = ModelK(kind, n_restarts=nr).fit(d[d.transect != t]); mu, sd = m.predict(te); Y.append(te.y.values); MU.append(mu); SD.append(np.sqrt(sd ** 2 + te.sd_y.values ** 2))
    s = scores(np.concatenate(Y), np.concatenate(MU), np.concatenate(SD)); mf = ModelK(kind, n_restarts=nr).fit(d); ll = float(mf.gp.log_marginal_likelihood_value_)
    rows.append({"kernel": {"m32": "Matern 3/2", "m52": "Matern 5/2", "rbf": "squared exponential"}[kind], "restarts": nr, "CRPS": s["CRPS"], "MAE": s["MAE"], "coverage95": s["coverage95"], "log_marginal_likelihood_full_fit": ll}); print(rows[-1])
pd.DataFrame(rows).to_csv(D / "step2_kernel_ablation.csv", index=False)
