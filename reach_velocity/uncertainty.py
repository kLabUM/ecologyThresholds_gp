"""The uncertainty chain (entry 14). One object, one sampler, every source switchable.
stage -> local depth h = stage - bed -> amplitude a = A (stage - h0)^gamma -> velocity v(z) = a ln(z/z0), depth average
a (ln(h/z0) - 1) -> indicator below a threshold. Sources: A_field (joint GP posterior over the points), disc (the
discrepancy calibrated on held-out days), bed, stage, law (the (h0, gamma) ridge), z0."""
import json, numpy as np, pandas as pd
from scipy import linalg
from gpcore import *
LAW_RIDGE = [(0.45, 1.0), (0.35, 1.25), (0.55, 0.75), (0.35, 1.0), (0.45, 0.75), (0.65, 0.5), (0.55, 0.5)]   # every (h0, gamma) on the grid within 0.01 CRPS of the best (optimiser noise is 0.007) with h0 below the lowest stage the river reached (0.66 m); equivalent inside the surveyed range, they diverge above it (entries 16 and 22)
ALL_SOURCES = ("A_field", "disc", "bed", "stage", "law", "z0")
Z0_RANGE = (0.007, 0.016); STAGE_SD = 0.015; BED_SD = 0.06; Z_MEAS = 0.25   # typical height of the readings that set a
def g_law(stage, h0, gamma): return np.clip(stage - h0, 0.02, None) ** gamma
class Chain:
    def __init__(self, df, primary=(0.45, 1.0), n_restarts=2, ridge=None):
        self.df = df; self.primary = primary; self.models = {}; self.ridge = list(ridge) if ridge is not None else list(LAW_RIDGE)
        for (h0, gamma) in self.ridge:
            d = with_target(df, ("stage", gamma, h0)); self.models[(h0, gamma)] = Model(n_restarts=n_restarts).fit(d)
        self.calibrate(); self.disc_split = split_discrepancy(self)   # common (day-level) and independent parts, see the note at the end of this file
    def calibrate(self):
        """day-out under the primary law: how much wider than the GP band the held-out errors are"""
        h0, gamma = self.primary; d = with_target(self.df, ("stage", gamma, h0)); days = sorted(d.date.unique()); Z, SD = [], []
        for day in days:
            tr, te = np.where(d.date != day)[0], np.where(d.date == day)[0]; m = Model(n_restarts=1).fit(d.iloc[tr]); mu, sd = m.predict(d.iloc[te]); st = np.sqrt(sd ** 2 + d.sd_y.values[te] ** 2)
            Z.append((d.y.values[te] - mu) / st); SD.append(st)
        Z, SD = np.concatenate(Z), np.concatenate(SD); q = float(np.quantile(np.abs(Z), 0.95)); self.inflate = max(1.0, q / 1.96)
        self.disc_sd = float(np.sqrt(max(self.inflate ** 2 - 1.0, 0.0) * np.mean(SD ** 2))); self.calib = {"held_out_coverage95_before": float(np.mean(np.abs(Z) <= 1.96)), "q95_abs_z": q, "inflation": self.inflate, "disc_sd_log": self.disc_sd, "n": int(len(Z))}
        return self.calib
    def sample(self, stage, pts, n=500, sources=ALL_SOURCES, seed=0, z_ref=None):
        """pts: DataFrame with the features, NOTE_COLS, bed_m (and bed_sd_m optional). Returns dict of (n, npts) arrays."""
        rng = np.random.default_rng(seed); npts = len(pts); out = {"a": np.zeros((n, npts)), "h": np.zeros((n, npts)), "vbar": np.zeros((n, npts)), "vz": np.zeros((n, npts))}
        mus, Ls = {}, {}
        for key, m in self.models.items(): mu, cov = m.predict_cov(pts); mus[key] = mu; Ls[key] = linalg.cholesky(cov + 1e-8 * np.eye(npts), lower=True)
        bed = pts.bed_m.values; bed_sd = pts.bed_sd_m.values if "bed_sd_m" in pts else np.full(npts, BED_SD)
        for k in range(n):
            key = self.ridge[rng.integers(len(self.ridge))] if "law" in sources else self.primary; h0, gamma = key
            logA = mus[key] + (Ls[key] @ rng.standard_normal(npts) if "A_field" in sources else 0.0) + ((rng.normal(0, self.disc_split["independent_sd_log"], npts) + rng.normal(0, self.disc_split["common_sd_log"])) if "disc" in sources else 0.0)
            b = bed + (rng.normal(0, bed_sd) if "bed" in sources else 0.0); st = stage + (rng.normal(0, STAGE_SD) if "stage" in sources else 0.0)
            z0 = float(np.exp(rng.uniform(np.log(Z0_RANGE[0]), np.log(Z0_RANGE[1])))) if "z0" in sources else Z0
            h = np.clip(st - b, 0.0, None); a = np.exp(logA) * g_law(st, h0, gamma) * (np.log(Z_MEAS / Z0) / np.log(Z_MEAS / z0)); wet = h > np.e * z0   # a was fitted with Z0; another z0 rescales it through the readings' typical height
            vbar = np.where(wet, a * (np.log(np.clip(h, np.e * z0, None) / z0) - 1), 0.0); vz = np.where(h > (z_ref or 0.3048), a * np.log((z_ref or 0.3048) / z0), np.nan)
            out["a"][k], out["h"][k], out["vbar"][k], out["vz"][k] = a, h, vbar, vz
        return out
def section_stats(samples, threshold=0.2):
    """per draw: section mean of the depth-averaged velocity and the share of points below the threshold (wet points only)"""
    vb, h = samples["vbar"], samples["h"]; wet = h > 0.02; mean_v = np.array([np.mean(vb[k][wet[k]]) if wet[k].any() else np.nan for k in range(len(vb))]); share = np.array([np.mean(vb[k][wet[k]] < threshold) if wet[k].any() else np.nan for k in range(len(vb))]); return mean_v, share
def sample_series(chain, stages, pts, n=100, sources=ALL_SOURCES, seed=0, threshold=0.2):
    """vectorised over time: precompute the draws that do not depend on stage, then run every stage through them.
    returns dict with (n, T) arrays of the section mean depth-averaged velocity and the share below the threshold, and the wet fraction"""
    rng = np.random.default_rng(seed); npts, T = len(pts), len(stages); stages = np.asarray(stages, float)
    keys = [chain.ridge[i] for i in rng.integers(len(chain.ridge), size=n)] if "law" in sources else [chain.primary] * n
    mus, Ls = {}, {}
    for key in set(keys): mu, cov = chain.models[key].predict_cov(pts); mus[key] = mu; Ls[key] = linalg.cholesky(cov + 1e-8 * np.eye(npts), lower=True)
    logA = np.stack([mus[k] + (Ls[k] @ rng.standard_normal(npts) if "A_field" in sources else 0.0) + ((rng.normal(0, chain.disc_split["independent_sd_log"], npts) + rng.normal(0, chain.disc_split["common_sd_log"])) if "disc" in sources else 0.0) for k in keys])
    bed = pts.bed_m.values[None, :] + (rng.normal(0, pts.bed_sd_m.values if "bed_sd_m" in pts else BED_SD, (n, npts)) if "bed" in sources else 0.0)
    z0 = np.exp(rng.uniform(np.log(Z0_RANGE[0]), np.log(Z0_RANGE[1]), n)) if "z0" in sources else np.full(n, Z0); zfac = np.log(Z_MEAS / Z0) / np.log(Z_MEAS / z0)
    h0 = np.array([k[0] for k in keys]); gam = np.array([k[1] for k in keys]); A = np.exp(logA) * zfac[:, None]
    mean_v, share, wetfrac = np.zeros((n, T)), np.zeros((n, T)), np.zeros((n, T))
    for j, st in enumerate(stages):
        s = st + (rng.normal(0, STAGE_SD, n) if "stage" in sources else 0.0); h = np.clip(s[:, None] - bed, 0, None); g = np.clip(s - h0, 0.02, None) ** gam
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
    """sections: dict name -> DataFrame of that section's points (features, notes, bed_m, across_m). Returns per-section (n, T) arrays of the
    width-weighted mean velocity and share below the threshold, the reach share (length-weighted, per draw), and the wet fraction."""
    rng = np.random.default_rng(seed); names = list(sections); pts = pd.concat([sections[t].assign(_sec=t) for t in names], ignore_index=True); npts, T = len(pts), len(stages); stages = np.asarray(stages, float)
    if disc is None: disc = chain.disc_split if hasattr(chain, 'disc_split') else split_discrepancy(chain)
    keys = [chain.ridge[i] for i in rng.integers(len(chain.ridge), size=n)] if "law" in sources else [chain.primary] * n; mus, Ls = {}, {}
    for key in set(keys): mu, cov = chain.models[key].predict_cov(pts); mus[key] = mu; Ls[key] = linalg.cholesky(cov + 1e-8 * np.eye(npts), lower=True)
    logA = np.stack([mus[k] + (Ls[k] @ rng.standard_normal(npts) if "A_field" in sources else 0.0) for k in keys])
    if "disc" in sources: logA += rng.normal(0, disc["independent_sd_log"], (n, npts)) + rng.normal(0, disc["common_sd_log"], (n, 1))   # the common part is one number per draw for the whole reach
    bed = pts.bed_m.values[None, :] + (rng.normal(0, pts.bed_sd_m.values if "bed_sd_m" in pts else BED_SD, (n, npts)) if "bed" in sources else 0.0)
    z0 = np.exp(rng.uniform(np.log(Z0_RANGE[0]), np.log(Z0_RANGE[1]), n)) if "z0" in sources else np.full(n, Z0); zfac = np.log(Z_MEAS / Z0) / np.log(Z_MEAS / z0)
    h0 = np.array([k[0] for k in keys]); gam = np.array([k[1] for k in keys]); A = np.exp(logA) * zfac[:, None]
    sec_idx = {t: np.where(pts._sec.values == t)[0] for t in names}; ww = {t: width_weights(sections[t]) for t in names}
    al = np.array([float(sections[t].along_m.median()) for t in names]); o = np.argsort(al); ala = al[o]; lw = np.ones(len(names))
    if len(names) > 1: edges = np.r_[ala[0] - (ala[1] - ala[0]) / 2, (ala[1:] + ala[:-1]) / 2, ala[-1] + (ala[-1] - ala[-2]) / 2]; lw = np.zeros(len(names)); lw[o] = np.diff(edges)
    lw = lw / lw.sum()
    out = {"mean_v": {t: np.full((n, T), np.nan) for t in names}, "share": {t: np.full((n, T), np.nan) for t in names}, "wet_frac": {t: np.zeros((n, T)) for t in names}, "reach_share": np.full((n, T), np.nan), "length_weights": dict(zip(names, lw)), "disc": disc}
    for j0 in range(0, T, chunk):
        js = np.arange(j0, min(j0 + chunk, T)); st = stages[js][None, :] + (rng.normal(0, STAGE_SD, (n, 1)) if "stage" in sources else 0.0)   # (n, tj): one stage error per draw, shared by every section
        g = np.clip(st - h0[:, None], 0.02, None) ** gam[:, None]; reach_acc = np.zeros((n, len(js)))
        for t in names:
            ix = sec_idx[t]; h = np.clip(st[:, :, None] - bed[:, None, ix], 0, None); a = A[:, None, ix] * g[:, :, None]; wet = h > np.e * z0[:, None, None]
            if z_ref is None: v = np.where(wet, a * (np.log(np.clip(h, 1e-6, None) / z0[:, None, None]) - 1), 0.0)
            else: wet = h > z_ref; v = np.where(wet, a * np.log(z_ref / z0[:, None, None]), 0.0)
            w = ww[t][None, None, :] * wet; wsum = w.sum(2); ok = wsum > 0
            mv = np.where(ok, (v * w).sum(2) / np.maximum(wsum, 1e-12), np.nan); sh = np.where(ok, ((v < threshold) * w).sum(2) / np.maximum(wsum, 1e-12), np.nan)
            out["mean_v"][t][:, js] = mv; out["share"][t][:, js] = sh; out["wet_frac"][t][:, js] = wet.mean(2); reach_acc += lw[names.index(t)] * np.nan_to_num(sh)
        out["reach_share"][:, js] = reach_acc
    return out
