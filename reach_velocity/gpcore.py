"""Shared pieces for the conveyance GP: point table, the model (WLS mean on notes + ARD Matern GP with fixed per-point
noise and a learned nugget), scores."""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel as C, Matern, WhiteKernel
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"; F = HERE / "figures"
SD_H = 0.06; SUBS = ["cobble", "gravel", "sand", "silt", "vegetation", "woody debris"]; NOTE_COLS = [f"sub_{s}" for s in SUBS] + ["wake_flag"]
FULL = ["along_m", "rel_across", "dist_bank_m", "rel_depth", "section_mean_depth_m", "lateral_bed_slope_abs"]
Z0 = json.load(open(D / "field_qc.json"))["z0_m"]
def load_points():
    fe = pd.read_csv(D / "features.csv"); lc = pd.read_csv(D / "locations.csv")[["location_id", "a", "water_depth_m", "stage_m", "stage_raw_m", "x_m", "y_m", "across_m", "bed_rtk_sensor_datum_m", "bed_elev_m", "point"]]
    df = fe.drop(columns=["across_m"]).merge(lc, on="location_id"); df = df[df.a > 0].reset_index(drop=True)
    q = pd.read_csv(D / "usgs_daily.csv"); df["Q_cfs"] = df.date.map(q.set_index("date").discharge_cfs)
    for s in SUBS: df[f"sub_{s}"] = (df.substrate == s).astype(float)
    return df
def destage(df, law):
    """returns the factor g such that A = a / g; law = ('none',) | ('depth', beta) | ('Q', m) | ('stage', beta, h0)"""
    if law[0] == "none": return np.ones(len(df))
    if law[0] == "depth": return df.water_depth_m.values ** law[1]
    if law[0] == "Q": return (df.Q_cfs.values / 350.0) ** law[1]
    if law[0] == "stage": return np.clip(df.stage_m.values - law[2], 0.02, None) ** law[1]
def with_target(df, law, beta_for_noise=0.5):
    d = df.copy(); g = destage(d, law); d["g"] = g; d["A"] = d.a / g; d["y"] = np.log(d.A)
    bdep = law[1] if law[0] == "depth" else 0.0
    d["sd_y"] = np.clip(np.sqrt((d.a_sd / d.a) ** 2 + (bdep * SD_H / d.water_depth_m) ** 2), 0.02, 1.0); return d
def _wls(Xb, y, w):
    Wt = np.sqrt(w)[:, None]; b, *_ = np.linalg.lstsq(Xb * Wt, y * Wt[:, 0], rcond=None); return b
class Model:
    def __init__(self, feats=FULL, notes=True, n_restarts=3, nu=1.5): self.feats, self.notes, self.n_restarts, self.nu = feats, notes, n_restarts, nu
    def design(self, d): return np.column_stack([np.ones(len(d))] + ([d[c].values for c in NOTE_COLS] if self.notes else []))
    def fit(self, d):
        w = 1 / d.sd_y.values ** 2; self.b = _wls(self.design(d), d.y.values, w); r = d.y.values - self.design(d) @ self.b
        if self.feats is None: self.rsd = float(np.sqrt(np.average(r ** 2, weights=w))); return self
        X = d[self.feats].values.astype(float); self.mx, self.sx = X.mean(0), X.std(0) + 1e-9; self.ys = r.std() + 1e-9
        k = C(1.0, (1e-2, 1e2)) * Matern(np.ones(len(self.feats)), [(0.05, 50.0)] * len(self.feats), nu=self.nu) + WhiteKernel(0.1, (1e-4, 4.0))
        self.gp = GaussianProcessRegressor(k, alpha=(d.sd_y.values / self.ys) ** 2, n_restarts_optimizer=self.n_restarts, random_state=0).fit((X - self.mx) / self.sx, r / self.ys); return self
    def predict(self, d):
        mu = self.design(d) @ self.b
        if self.feats is None: return mu, np.full(len(d), self.rsd)
        m, sd = self.gp.predict((d[self.feats].values.astype(float) - self.mx) / self.sx, return_std=True); return mu + m * self.ys, sd * self.ys
    def predict_cov(self, d):
        mu = self.design(d) @ self.b; m, cov = self.gp.predict((d[self.feats].values.astype(float) - self.mx) / self.sx, return_cov=True); return mu + m * self.ys, cov * self.ys ** 2
    def scales(self):
        if self.feats is None: return {}
        k = self.gp.kernel_; ls = np.atleast_1d(k.k1.k2.length_scale) * self.sx; return {"amp_sd": float(np.sqrt(k.k1.k1.constant_value) * self.ys), "nugget_sd": float(np.sqrt(k.k2.noise_level) * self.ys), **{f: float(l) for f, l in zip(self.feats, ls)}}
def scores(y, mu, sd):
    z = (y - mu) / sd; crps = sd * (z * (2 * norm.cdf(z) - 1) + 2 * norm.pdf(z) - 1 / np.sqrt(np.pi))
    return {"MAE": float(np.mean(np.abs(y - mu))), "RMSE": float(np.sqrt(np.mean((y - mu) ** 2))), "CRPS": float(np.mean(crps)), "coverage95": float(np.mean(np.abs(z) <= 1.96)), "mean_sd": float(np.mean(sd)), "bias": float(np.mean(mu - y)), "n": int(len(y))}
