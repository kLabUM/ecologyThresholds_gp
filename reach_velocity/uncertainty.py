"""The uncertainty chain (entry 14). One object, one sampler, every source switchable.
stage -> local depth h = stage - bed -> amplitude a = A (stage - h0)^gamma -> velocity v(z) = a ln(z/z0), depth average
a (ln(h/z0) - 1) -> indicator below a threshold. Sources: A_field (joint GP posterior over the points), disc (the
discrepancy calibrated on held-out days), bed, stage, law (the (h0, gamma) ridge), z0, bank (the 20 draws of bank set C, entry 30).
With bank set C the bank features depend on stage: the field is computed on a 5 cm stage grid and read at each draw's stage (Field)."""
import json, numpy as np, pandas as pd
from scipy import linalg
from gpcore import *
LEVEL_STEP = 0.05   # banklines_3dep.LEVEL_STEP_M
_BANKS = {}
def banks_for(reach):
    """Bank set C for a reach (banklines_3dep.StageFeatures), built once per run."""
    import banklines_3dep as B
    if reach not in _BANKS: _BANKS[reach] = B.StageFeatures(reach)
    return _BANKS[reach]
def law_set(path=None, tol=0.01, h0_below=0.66):
    """The stage-law set by the rule of entries 16 and 22: every (h0, gamma) on the day-out grid of 22_step3b_law_refine.py within tol CRPS
    of the best (optimiser noise is 0.007), with h0 below the lowest stage the river reached (0.66 m); ordered by CRPS, the first is the
    primary law. Read from the grid so that the rule, not a hand-copied list, sets the law (entry 30)."""
    g = pd.read_csv(path or (D / "step3b_law_grid.csv")); g = g[(g.law == "stage") & (g.h0 < h0_below)].sort_values("CRPS")
    s = g[g.CRPS <= g.CRPS.iloc[0] + tol]; return [(round(float(r.h0), 2), round(float(r.gamma), 2)) for r in s.itertuples()]
LAW_RIDGE = law_set() if (D / "step3b_law_grid.csv").exists() else [(0.45, 0.75), (0.35, 1.0), (0.55, 0.5), (0.35, 0.75), (0.55, 0.75), (0.45, 1.0), (0.65, 0.5), (0.35, 1.25)]   # fallback: the rule's answer under bank set C (entry 30)
PRIMARY = LAW_RIDGE[0]
ALL_SOURCES = ("A_field", "disc", "bed", "stage", "law", "z0", "bank")   # 'bank': the 20 draws of bank set C (entry 30); used only when the chain has banks
Z0_RANGE = (0.007, 0.016); STAGE_SD = 0.015; BED_SD = 0.06; Z_MEAS = 0.25   # typical height of the readings that set a
def g_law(stage, h0, gamma): return np.clip(stage - h0, 0.02, None) ** gamma
class Chain:
    """banks: a banklines_3dep.StageFeatures (set C). With banks, the bank features of the training points are C at each point's survey
    day, and every law is also fitted once per C draw (the 'bank' source); predictions read the bank features at the stage of each draw
    on a 5 cm grid (entry 30). Without banks the features in df are used as they are."""
    def __init__(self, df, primary=None, n_restarts=2, ridge=None, banks=None):
        self.banks = banks; df = banks.training(df) if banks is not None else df; primary = PRIMARY if primary is None else primary
        self.df = df; self.primary = primary; self.models = {}; self.ridge = list(ridge) if ridge is not None else list(LAW_RIDGE)
        for (h0, gamma) in self.ridge:
            d = with_target(df, ("stage", gamma, h0)); self.models[(h0, gamma)] = Model(n_restarts=n_restarts).fit(d)
        self.bank_models = None
        if banks is not None:
            self.bank_models = []
            for j in range(banks.n_draws):
                dj = banks.training(df, j); self.bank_models.append({(h0, gamma): Model(n_restarts=n_restarts).fit(with_target(dj, ("stage", gamma, h0))) for (h0, gamma) in self.ridge})
        self.calibrate(); self.disc_split = split_discrepancy(self)   # common (day-level) and independent parts, see the note at the end of this file
    def model(self, key, j=None): return self.models[key] if j is None else self.bank_models[j][key]
    def calibrate(self):
        """day-out under the primary law: how much wider than the GP band the held-out errors are"""
        h0, gamma = self.primary; d = with_target(self.df, ("stage", gamma, h0)); days = sorted(d.date.unique()); Z, SD = [], []
        for day in days:
            tr, te = np.where(d.date != day)[0], np.where(d.date == day)[0]; m = Model(n_restarts=1).fit(d.iloc[tr]); mu, sd = m.predict(d.iloc[te]); st = np.sqrt(sd ** 2 + d.sd_y.values[te] ** 2)
            Z.append((d.y.values[te] - mu) / st); SD.append(st)
        Z, SD = np.concatenate(Z), np.concatenate(SD); q = float(np.quantile(np.abs(Z), 0.95)); self.inflate = max(1.0, q / 1.96)
        self.disc_sd = float(np.sqrt(max(self.inflate ** 2 - 1.0, 0.0) * np.mean(SD ** 2))); self.calib = {"held_out_coverage95_before": float(np.mean(np.abs(Z) <= 1.96)), "q95_abs_z": q, "inflation": self.inflate, "disc_sd_log": self.disc_sd, "n": int(len(Z))}
        return self.calib
    def draw_keys(self, rng, n, sources):
        keys = [self.ridge[i] for i in rng.integers(len(self.ridge), size=n)] if "law" in sources else [self.primary] * n
        js = list(rng.integers(self.banks.n_draws, size=n)) if (self.banks is not None and "bank" in sources) else [None] * n
        return keys, js
    def sample(self, stage, pts, n=500, sources=ALL_SOURCES, seed=0, z_ref=None):
        """pts: DataFrame with the features, NOTE_COLS, bed_m (and bed_sd_m optional; Easting, Northing when the chain has banks). Returns dict of (n, npts) arrays."""
        rng = np.random.default_rng(seed); npts = len(pts); keys, js = self.draw_keys(rng, n, sources)
        st = stage + (rng.normal(0, STAGE_SD, n) if "stage" in sources else np.zeros(n))
        F = Field(self, pts, keys, js, rng, "A_field" in sources, st.min(), st.max())
        logA = F.at(st) + ((rng.normal(0, self.disc_split["independent_sd_log"], (n, npts)) + rng.normal(0, self.disc_split["common_sd_log"], (n, 1))) if "disc" in sources else 0.0)
        bed_sd = pts.bed_sd_m.values if "bed_sd_m" in pts else np.full(npts, BED_SD); b = pts.bed_m.values[None, :] + (rng.normal(0, bed_sd, (n, npts)) if "bed" in sources else 0.0)
        z0 = np.exp(rng.uniform(np.log(Z0_RANGE[0]), np.log(Z0_RANGE[1]), n)) if "z0" in sources else np.full(n, Z0)
        h0 = np.array([k[0] for k in keys]); gam = np.array([k[1] for k in keys])
        h = np.clip(st[:, None] - b, 0.0, None); a = np.exp(logA) * g_law(st, h0, gam)[:, None] * (np.log(Z_MEAS / Z0) / np.log(Z_MEAS / z0))[:, None]; wet = h > np.e * z0[:, None]   # a was fitted with Z0; another z0 rescales it through the readings' typical height
        zr = z_ref or 0.3048; vbar = np.where(wet, a * (np.log(np.clip(h, np.e * z0[:, None], None) / z0[:, None]) - 1), 0.0); vz = np.where(h > zr, a * np.log(zr / z0[:, None]), np.nan)
        return {"a": a, "h": h, "vbar": vbar, "vz": vz}
class Field:
    """The conveyance field log A (without the discrepancy) for every draw at pts. Each draw has its law (key) and bank draw (j, None for
    C's best estimate); Z is one standard normal vector per draw, shared by every stage (common random numbers). Without banks the field
    does not depend on stage. With banks it is computed on the 5 cm stage grid for each (law, bank draw) pair over the stages that
    pair's draws need, and read between the two grid levels around each draw's stage (linear, held at the ends)."""
    def __init__(self, chain, pts, keys, js, rng, use_field, lo=None, hi=None, need=None):
        self.n, npts = len(keys), len(pts); self.Z = rng.standard_normal((self.n, npts)) if use_field else np.zeros((self.n, npts)); self.banks = chain.banks
        groups = {}
        for k, (key, j) in enumerate(zip(keys, js)): groups.setdefault((key, j), []).append(k)
        def mean_and_draws(m, P, idx):
            if not use_field: mu, _ = m.predict(P); return np.broadcast_to(mu, (len(idx), npts))
            mu, cov = m.predict_cov(P); L = linalg.cholesky(cov + 1e-8 * np.eye(npts), lower=True); return mu[None, :] + self.Z[idx] @ L.T
        if self.banks is None:
            self.G = np.zeros((self.n, npts))
            for (key, j), idx in groups.items(): self.G[idx] = mean_and_draws(chain.model(key), pts, idx)
            return
        lo = -np.inf if lo is None else lo; hi = np.inf if hi is None else hi; self.lv = self.banks.grid(lo, hi); L_ = len(self.lv)
        self.G = np.full((self.n, L_, npts), np.nan); E, N = pts.Easting.values, pts.Northing.values; feats = {}
        for (key, j), idx in groups.items():
            if need is not None: lo_k, hi_k = need(idx)
            else: lo_k, hi_k = self.lv[0], self.lv[-1]
            il = np.arange(max(int(np.searchsorted(self.lv, lo_k + 1e-9) - 1), 0), min(int(np.searchsorted(self.lv, hi_k - 1e-9)) + 1, L_ - 1) + 1)
            for i in il:
                if (i, j) not in feats:
                    dist, rel, wid = self.banks.features(E, N, self.lv[i], j); feats[(i, j)] = pts.assign(dist_bank_m=dist, rel_across=rel, width_m=wid)
                self.G[idx, i] = mean_and_draws(chain.model(key, j), feats[(i, j)], idx)
    def at(self, s, ix=None):
        """s: (n,) or (n, T) stages; returns log A as (n, npts) or (n, T, npts), columns ix only if given."""
        if self.banks is None:
            G = self.G if ix is None else self.G[:, ix]; return G if np.ndim(s) == 1 else np.broadcast_to(G[:, None, :], (self.n, np.shape(s)[1], G.shape[1]))
        G = self.G if ix is None else self.G[:, :, ix]; L_ = len(self.lv)
        if L_ == 1: return G[:, 0] if np.ndim(s) == 1 else np.broadcast_to(G[:, None, 0], (self.n, np.shape(s)[1], G.shape[2]))
        pos = np.clip((np.asarray(s, float) - self.lv[0]) / LEVEL_STEP, 0, L_ - 1); i0 = np.minimum(np.floor(pos).astype(int), L_ - 2); w = pos - i0
        r = np.arange(self.n) if np.ndim(s) == 1 else np.arange(self.n)[:, None]
        return G[r, i0] * (1 - w)[..., None] + G[r, i0 + 1] * w[..., None]
def section_stats(samples, threshold=0.2):
    """per draw: section mean of the depth-averaged velocity and the share of points below the threshold (wet points only)"""
    vb, h = samples["vbar"], samples["h"]; wet = h > 0.02; mean_v = np.array([np.mean(vb[k][wet[k]]) if wet[k].any() else np.nan for k in range(len(vb))]); share = np.array([np.mean(vb[k][wet[k]] < threshold) if wet[k].any() else np.nan for k in range(len(vb))]); return mean_v, share
def sample_series(chain, stages, pts, n=100, sources=ALL_SOURCES, seed=0, threshold=0.2):
    """vectorised over time: the draws that do not depend on stage first, then every stage through them.
    returns dict with (n, T) arrays of the section mean depth-averaged velocity and the share below the threshold, and the wet fraction"""
    rng = np.random.default_rng(seed); npts, T = len(pts), len(stages); stages = np.asarray(stages, float); keys, js = chain.draw_keys(rng, n, sources)
    S = stages[None, :] + (rng.normal(0, STAGE_SD, (n, T)) if "stage" in sources else 0.0)
    F = Field(chain, pts, keys, js, rng, "A_field" in sources, S.min(), S.max())
    disc = (rng.normal(0, chain.disc_split["independent_sd_log"], (n, npts)) + rng.normal(0, chain.disc_split["common_sd_log"], (n, 1))) if "disc" in sources else 0.0
    bed = pts.bed_m.values[None, :] + (rng.normal(0, pts.bed_sd_m.values if "bed_sd_m" in pts else BED_SD, (n, npts)) if "bed" in sources else 0.0)
    z0 = np.exp(rng.uniform(np.log(Z0_RANGE[0]), np.log(Z0_RANGE[1]), n)) if "z0" in sources else np.full(n, Z0); zfac = np.log(Z_MEAS / Z0) / np.log(Z_MEAS / z0)
    h0 = np.array([k[0] for k in keys]); gam = np.array([k[1] for k in keys])
    mean_v, share, wetfrac = np.zeros((n, T)), np.zeros((n, T)), np.zeros((n, T))
    for j in range(T):
        s = S[:, j] if np.ndim(S) == 2 else np.full(n, stages[j]); A = np.exp(F.at(s) + disc) * zfac[:, None]; h = np.clip(s[:, None] - bed, 0, None); g = np.clip(s - h0, 0.02, None) ** gam
        a = A * g[:, None]; wet = h > np.e * z0[:, None]; vb = np.where(wet, a * (np.log(np.clip(h, 1e-6, None) / z0[:, None]) - 1), 0.0)
        nw = wet.sum(1); mean_v[:, j] = np.where(nw > 0, (vb * wet).sum(1) / np.maximum(nw, 1), np.nan); share[:, j] = np.where(nw > 0, ((vb < threshold) & wet).sum(1) / np.maximum(nw, 1), np.nan); wetfrac[:, j] = nw / npts
    return {"mean_v": mean_v, "share": share, "wet_frac": wetfrac}

# ---------------------------------------------------------------- joint reach sampler (added after the reviewer pass, 22 Sep 2026)
# Reviewer finding: drawing each section separately lets errors that are common to the reach (a day-level shift left by the law, the
# stage reading, the roughness scale) average away in the length-weighted reach total. This sampler draws the whole reach at once:
# one GP draw over every point of every section, one stage noise and one z0 per draw, and the discrepancy split into a part common to
# the reach at a time (the residual day level after the law) and a part independent per point. Shares are width-weighted (each point
# represents the interval to the midpoints with its neighbours; the margins beyond the last point are not included, so above the
# surveyed stage range the share of slow water is a lower bound).
def split_discrepancy(chain):
    """common (day-level) and independent parts of the discrepancy, in log. The day-level residual after the primary law is the sd of
    the held-out day means about zero, less the sampling noise of those means; the independent part keeps the point-level calibration."""
    h0, gamma = chain.primary; d = with_target(chain.df, ("stage", gamma, h0)); days = sorted(d.date.unique()); means, sam = [], []
    for day in days:
        tr, te = np.where(d.date != day)[0], np.where(d.date == day)[0]
        if len(te) < 3: continue   # a day with one or two points says nothing about a day level
        m = Model(n_restarts=1).fit(d.iloc[tr]); mu, sd = m.predict(d.iloc[te]); r = d.y.values[te] - mu
        means.append(float(np.mean(r))); sam.append(float(np.std(r, ddof=1) / np.sqrt(len(r))))
    var_common = max(np.var(means, ddof=1) - np.mean(np.square(sam)), 0.0); common = float(np.sqrt(var_common)); indep = float(np.sqrt(max(chain.disc_sd ** 2 - var_common, 0.0)))
    return {"common_sd_log": common, "independent_sd_log": indep, "day_means": dict(zip(days, means)), "total_disc_sd_log": chain.disc_sd}
def width_weights(pts):
    y = np.asarray(pts.across_m, float); order = np.argsort(y); w = np.zeros(len(y)); ys = y[order]
    if len(ys) == 1: w[:] = 1.0; return w
    dy = np.diff(ys); wo = np.r_[dy[0] / 2, (dy[:-1] + dy[1:]) / 2, dy[-1] / 2]; w[order] = wo; return w
def sample_reach(chain, stages, sections, n=200, sources=ALL_SOURCES, seed=0, threshold=0.2, z_ref=None, disc=None, chunk=400):
    """sections: dict name -> DataFrame of that section's points (features, notes, bed_m, across_m; Easting, Northing with banks). Returns per-section
    (n, T) arrays of the width-weighted mean velocity and share below the threshold, the reach share (length-weighted, per draw), and the wet fraction."""
    rng = np.random.default_rng(seed); names = list(sections); pts = pd.concat([sections[t].assign(_sec=t) for t in names], ignore_index=True); npts, T = len(pts), len(stages); stages = np.asarray(stages, float)
    if disc is None: disc = chain.disc_split if hasattr(chain, 'disc_split') else split_discrepancy(chain)
    keys, js = chain.draw_keys(rng, n, sources)
    sn = rng.normal(0, STAGE_SD, (n, 1)) if "stage" in sources else np.zeros((n, 1))   # one stage error per draw, shared by every section and time
    F = Field(chain, pts, keys, js, rng, "A_field" in sources, stages.min() + sn.min(), stages.max() + sn.max())
    dl = (rng.normal(0, disc["independent_sd_log"], (n, npts)) + rng.normal(0, disc["common_sd_log"], (n, 1))) if "disc" in sources else np.zeros((n, npts))   # the common part is one number per draw for the whole reach
    bed = pts.bed_m.values[None, :] + (rng.normal(0, pts.bed_sd_m.values if "bed_sd_m" in pts else BED_SD, (n, npts)) if "bed" in sources else 0.0)
    z0 = np.exp(rng.uniform(np.log(Z0_RANGE[0]), np.log(Z0_RANGE[1]), n)) if "z0" in sources else np.full(n, Z0); zfac = np.log(Z_MEAS / Z0) / np.log(Z_MEAS / z0)
    h0 = np.array([k[0] for k in keys]); gam = np.array([k[1] for k in keys])
    sec_idx = {t: np.where(pts._sec.values == t)[0] for t in names}; ww = {t: width_weights(sections[t]) for t in names}
    al = np.array([float(sections[t].along_m.median()) for t in names]); o = np.argsort(al); ala = al[o]; lw = np.ones(len(names))
    if len(names) > 1: edges = np.r_[ala[0] - (ala[1] - ala[0]) / 2, (ala[1:] + ala[:-1]) / 2, ala[-1] + (ala[-1] - ala[-2]) / 2]; lw = np.zeros(len(names)); lw[o] = np.diff(edges)
    lw = lw / lw.sum()
    out = {"mean_v": {t: np.full((n, T), np.nan) for t in names}, "share": {t: np.full((n, T), np.nan) for t in names}, "wet_frac": {t: np.zeros((n, T)) for t in names}, "reach_share": np.full((n, T), np.nan), "length_weights": dict(zip(names, lw)), "disc": disc}
    if chain.banks is not None: chunk = min(chunk, 40)   # the field is read per draw and time: keep the (n, chunk, points) arrays small
    for j0 in range(0, T, chunk):
        js_ = np.arange(j0, min(j0 + chunk, T)); st = stages[js_][None, :] + sn   # (n, tj)
        g = np.clip(st - h0[:, None], 0.02, None) ** gam[:, None]; reach_acc = np.zeros((n, len(js_)))
        for t in names:
            ix = sec_idx[t]; h = np.clip(st[:, :, None] - bed[:, None, ix], 0, None); A = np.exp(F.at(st, ix) + dl[:, None, ix]) * zfac[:, None, None]; a = A * g[:, :, None]; wet = h > np.e * z0[:, None, None]
            if z_ref is None: v = np.where(wet, a * (np.log(np.clip(h, 1e-6, None) / z0[:, None, None]) - 1), 0.0)
            else: wet = h > z_ref; v = np.where(wet, a * np.log(z_ref / z0[:, None, None]), 0.0)
            w = ww[t][None, None, :] * wet; wsum = w.sum(2); ok = wsum > 0
            mv = np.where(ok, (v * w).sum(2) / np.maximum(wsum, 1e-12), np.nan); sh = np.where(ok, ((v < threshold) * w).sum(2) / np.maximum(wsum, 1e-12), np.nan)
            out["mean_v"][t][:, js_] = mv; out["share"][t][:, js_] = sh; out["wet_frac"][t][:, js_] = wet.mean(2); reach_acc += lw[names.index(t)] * np.nan_to_num(sh)
        out["reach_share"][:, js_] = reach_acc
    return out
class GridLevels:
    """The 2 m grid's field (C's best estimate, primary law) and bank distance on the 5 cm stage grid (24_), read at one stage."""
    def __init__(self, path=None):
        z = np.load(path or (D / "grid_levels.npz")); self.lv, self.logA, self.dist = z["levels"], z["logA"], z["dist_bank_m"]
    def at(self, stage):
        if len(self.lv) == 1: return self.logA[0], self.dist[0]
        pos = float(np.clip((stage - self.lv[0]) / LEVEL_STEP, 0, len(self.lv) - 1)); i = min(int(pos), len(self.lv) - 2); w = pos - i
        return self.logA[i] * (1 - w) + self.logA[i + 1] * w, self.dist[i] * (1 - w) + self.dist[i + 1] * w
