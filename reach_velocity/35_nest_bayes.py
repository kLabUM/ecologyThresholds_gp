"""A Bayesian version of the ecological number (the PI, 22 September 2026): the nests update a literature prior on how nest use falls with
velocity, and the ecological number becomes the reach's mean relative suitability, threshold-free, with a band from both the curve and
the velocity chain. Model: used-against-available logistic on velocity, bank distance and depth (unstandardised units), Gaussian priors,
Laplace posterior. The velocity coefficient's prior comes from the USFWS spawning suitability curve (Edwards et al. 1983), digitised in
data_raw/literature_priors.json and made weak (see there).
Outputs: nest_bayes.json, nest_bayes_curve.csv, nest_bayes_series.csv, figures/11_nest_bayes.png"""
import json, numpy as np, pandas as pd
from scipy.spatial import cKDTree; from scipy.optimize import minimize; from scipy import linalg
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import matplotlib.dates as mdates
from gpcore import *; from uncertainty import *
rng = np.random.default_rng(0); Z0 = json.load(open(D / "field_qc.json"))["z0_m"]; H0, GAM = PRIMARY   # the primary stage law (uncertainty.law_set, entry 30)
from scipy.optimize import minimize_scalar
LIT = json.load(open(RAW_LIT := (D.parent / "data_raw" / "literature_priors.json")))
FT = 0.3048; sv = np.array(LIT["source_1"]["spawning_velocity_SI_ft_per_s"]); vv_lit, S_lit = sv[:, 0] * FT, sv[:, 1]; av = np.array(LIT["source_1"]["adult_velocity_SI_ft_per_s_category_three_Glover_Creek"]); va_lit, Sa_lit = av[:, 0] * FT, av[:, 1]
def fit_exp(v, S):   # the exponential relative-use slope that best fits a digitised SI curve (least squares in S over 0 to 0.8 m/s)
    g = np.linspace(0, 0.8, 161); Si = np.interp(g, v, S); return float(minimize_scalar(lambda b: np.sum((np.exp(b * g) - Si) ** 2), bounds=(-30, 0), method="bounded").x)
b_spawn, b_adult = fit_exp(vv_lit, S_lit), fit_exp(va_lit, Sa_lit); sd_v = abs(b_adult - b_spawn) / 2   # weak: the steep utilisation curve sits two sd away
sd_ = np.array(LIT["source_1"]["spawning_depth_SI_ft"]); hh_lit, Sh_lit = sd_[:, 0] * FT, sd_[:, 1]; g = np.linspace(0.15, 0.9, 76); Sh = np.clip(np.interp(g, hh_lit, Sh_lit), 0.02, 1); b_depth = float(np.polyfit(g, np.log(Sh), 1)[0])   # the exponential slope over the nest depth range
PRIOR = {"velocity_per_m_s": {"mean": b_spawn, "sd": sd_v, "source": f"Edwards et al. (1983) spawning SI curve 20110 fitted by exp(b v): {b_spawn:.2f} per m/s; sd set so the adult utilisation curve ({b_adult:.2f} per m/s) is two sd away; a weak prior carrying sign and order of magnitude"},
         "bank_distance_per_m": {"mean": 0.0, "sd": 0.5, "source": "no published curve; weak, zero-centred"}, "depth_per_m": {"mean": b_depth, "sd": 5.0, "source": f"Edwards et al. (1983) spawning depth SI curve over 0.15 to 0.9 m fitted by exp(b h): {b_depth:.1f} per m; sd 5 (weak)"}, "intercept": {"mean": 0.0, "sd": 10.0}}
print("literature prior:", json.dumps({k: {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in v.items()} for k, v in PRIOR.items()}, indent=1))
df = load_points(); df["bed_m"] = df.bed_rtk_sensor_datum_m.fillna(df.bed_elev_m); ch = Chain(df, banks=banks_for("Middle"))
grid = pd.read_csv(D / "grid.csv"); nests = pd.read_csv(D / "step7_nests.csv").dropna(subset=["v_pred_ms"]).reset_index(drop=True); cell = grid.iloc[cKDTree(grid[["x_m", "y_m"]].values).query(nests[["x_m", "y_m"]].values)[1]].reset_index(drop=True); cix = cKDTree(grid[["x_m", "y_m"]].values).query(nests[["x_m", "y_m"]].values)[1]; GLV = GridLevels()   # the field and bank distance at each nest's own stage (bank set C, entry 30)
def vbar(logA, bed, st):
    h = np.clip(st - bed, 0, None); a = np.exp(logA) * np.clip(st - H0, 0.02, None) ** GAM; return np.where(h > np.e * Z0, a * (np.log(np.clip(h, np.e * Z0, None) / Z0) - 1), np.nan), h
Xu, Xa = [], []
for i, r in nests.iterrows():
    st = r.stage_use; gA, gD = GLV.at(st); vg, hg = vbar(gA, grid.bed_m.values, st); ok = ~np.isnan(vg); c = cell.iloc[i]; vc, hc = vbar(gA[cix[i]], c.bed_m, st); Xu.append([1.0, float(vc), gD[cix[i]], float(hc)])
    for j in rng.choice(np.where(ok)[0], 20, replace=False): Xa.append([1.0, vg[j], gD[j], hg[j]])
X = np.array(Xu + Xa); y = np.r_[np.ones(len(Xu)), np.zeros(len(Xa))]
pm = np.array([PRIOR["intercept"]["mean"], PRIOR["velocity_per_m_s"]["mean"], PRIOR["bank_distance_per_m"]["mean"], PRIOR["depth_per_m"]["mean"]]); ps = np.array([PRIOR["intercept"]["sd"], PRIOR["velocity_per_m_s"]["sd"], PRIOR["bank_distance_per_m"]["sd"], PRIOR["depth_per_m"]["sd"]])
def laplace(prior_sd):
    def nlp(b): z = X @ b; return float(np.sum(np.logaddexp(0, z) - y * z) + np.sum((b - pm) ** 2 / (2 * prior_sd ** 2)))
    def grad(b): p = 1 / (1 + np.exp(-(X @ b))); return X.T @ (p - y) + (b - pm) / prior_sd ** 2
    res = minimize(nlp, np.zeros(4), jac=grad, method="BFGS"); p = 1 / (1 + np.exp(-(X @ res.x))); Hm = X.T @ (X * (p * (1 - p))[:, None]) + np.diag(1 / prior_sd ** 2); return res.x, linalg.inv(Hm)
b_post, C_post = laplace(ps); b_flat, C_flat = laplace(np.array([10.0, 1e3, 1e3, 1e3]))   # flat: the data alone
samp = rng.multivariate_normal(b_post, C_post, 4000); bv = samp[:, 1]; v50 = np.log(0.5) / bv
out = {"prior": PRIOR, "n_used": int(len(Xu)), "n_available": int(len(Xa)),
       "velocity_coef_per_m_s": {"prior": [float(pm[1]), float(ps[1])], "data_alone": [float(b_flat[1]), float(np.sqrt(C_flat[1, 1]))], "posterior": [float(b_post[1]), float(np.sqrt(C_post[1, 1]))]},
       "bank_distance_coef_per_m": {"data_alone": [float(b_flat[2]), float(np.sqrt(C_flat[2, 2]))], "posterior": [float(b_post[2]), float(np.sqrt(C_post[2, 2]))]}, "depth_coef_per_m": {"data_alone": [float(b_flat[3]), float(np.sqrt(C_flat[3, 3]))], "posterior": [float(b_post[3]), float(np.sqrt(C_post[3, 3]))]},
       "v50_m_s_posterior": {"median": float(np.median(v50[bv < 0])), "ci90": [float(np.quantile(v50[bv < 0], .05)), float(np.quantile(v50[bv < 0], .95))], "share_of_draws_with_positive_velocity_effect": float(np.mean(bv >= 0))},
       "definition": "relative suitability S(v) = exp(beta_v v), the used/available odds at velocity v relative to still water, at the same bank distance and depth; the ecological number is the width- and length-weighted mean of S over the reach's wetted width"}
# ---- the reach's mean suitability against stage and through the seasons: per-point velocities per draw, one beta_v per draw
secs = {t: g.assign(_x=g.x_m.round(2), _y=g.y_m.round(2)).drop_duplicates(["_x", "_y"]).sort_values("across_m").reset_index(drop=True) for t, g in df.groupby("transect")}; names = list(secs)
pts = pd.concat([secs[t].assign(_sec=t) for t in names], ignore_index=True); npts = len(pts); ww = np.concatenate([width_weights(secs[t]) for t in names]); sec_of = pts._sec.values
al = np.array([float(secs[t].along_m.median()) for t in names]); o = np.argsort(al); ala = al[o]; edges = np.r_[ala[0] - (ala[1] - ala[0]) / 2, (ala[1:] + ala[:-1]) / 2, ala[-1] + (ala[-1] - ala[-2]) / 2]; lw = np.zeros(len(names)); lw[o] = np.diff(edges); lw /= lw.sum(); lw_pt = np.array([lw[names.index(s)] for s in sec_of])
def suit_series(stages, n, seed, chunk=40):
    """per draw and stage: the reach's wetted-width-weighted mean suitability; every chain source (law, field, discrepancy, bed, stage,
    z0, bank set C) plus one draw of the velocity coefficient. The field is read at each draw's stage on the 5 cm grid (entry 30)."""
    rg = np.random.default_rng(seed); stages = np.asarray(stages, float); keys, js = ch.draw_keys(rg, n, ALL_SOURCES); sn = rg.normal(0, STAGE_SD, (n, 1))
    F = Field(ch, pts, keys, js, rg, True, stages.min() + sn.min(), stages.max() + sn.max())
    dl = rg.normal(0, ch.disc_split["independent_sd_log"], (n, npts)) + rg.normal(0, ch.disc_split["common_sd_log"], (n, 1))
    bed = pts.bed_m.values[None, :] + rg.normal(0, BED_SD, (n, npts)); z0 = np.exp(rg.uniform(np.log(Z0_RANGE[0]), np.log(Z0_RANGE[1]), n)); zfac = np.log(Z_MEAS / Z0) / np.log(Z_MEAS / z0); h0 = np.array([k[0] for k in keys]); gam = np.array([k[1] for k in keys])
    bvd = rg.multivariate_normal(b_post, C_post, n)[:, 1]; E = np.full((n, len(stages)), np.nan)
    for j0 in range(0, len(stages), chunk):
        js_ = np.arange(j0, min(j0 + chunk, len(stages))); st = stages[js_][None, :] + sn; g = np.clip(st - h0[:, None], 0.02, None) ** gam[:, None]
        A = np.exp(F.at(st) + dl[:, None, :]) * zfac[:, None, None]
        h = np.clip(st[:, :, None] - bed[:, None, :], 0, None); wet = h > np.e * z0[:, None, None]; v = np.where(wet, A * g[:, :, None] * (np.log(np.clip(h, 1e-6, None) / z0[:, None, None]) - 1), 0.0)
        S = np.exp(bvd[:, None, None] * v); w = (ww * lw_pt)[None, None, :] * wet   # width within the section times the section's length share
        E[:, js_] = (S * w).sum(2) / np.maximum(w.sum(2), 1e-12)   # the wetted-width-weighted mean suitability over the reach
    return E
grid_st = np.linspace(0.65, 1.35, 15); Eg = suit_series(grid_st, 300, 1); curve = pd.DataFrame({"stage_m": grid_st, "E_med": np.nanmedian(Eg, 0), "E_lo": np.nanquantile(Eg, .05, 0), "E_hi": np.nanquantile(Eg, .95, 0), "P_E_above_0.5": np.nanmean(Eg > 0.5, 0), "P_E_above_0.4": np.nanmean(Eg > 0.4, 0)}); curve.to_csv(D / "nest_bayes_curve.csv", index=False)
g10 = pd.read_csv(D / "series_10min.csv", parse_dates=["time_utc"]).set_index("time_utc"); g10.index = pd.to_datetime(g10.index, utc=True); st2 = g10.stage_m.resample("2h").mean().dropna(); st2 = st2[st2.index.year >= 2025]; Es = suit_series(st2.values, 100, 2)
ser = pd.DataFrame({"time_utc": st2.index, "stage_m": st2.values, "E_med": np.nanmedian(Es, 0), "E_lo": np.nanquantile(Es, .05, 0), "E_hi": np.nanquantile(Es, .95, 0), "P_E_above_0.5": np.nanmean(Es > 0.5, 0)}); ser.to_csv(D / "nest_bayes_series.csv", index=False)
out["reach_mean_suitability"] = {str(yr): {"median": float(np.nanmedian(ser[ser.time_utc.dt.year == yr].E_med)), "share_of_hours_with_P(E>0.5)_over_0.9": float(np.mean(ser[ser.time_utc.dt.year == yr]["P_E_above_0.5"] > 0.9))} for yr in [2025, 2026]}
out["reach_mean_suitability"]["at_0.90_m"] = {"median": float(curve.E_med[np.argmin(np.abs(grid_st - 0.90))]), "lo": float(curve.E_lo[np.argmin(np.abs(grid_st - 0.90))]), "hi": float(curve.E_hi[np.argmin(np.abs(grid_st - 0.90))])}
json.dump(out, open(D / "nest_bayes.json", "w"), indent=1); print(json.dumps(out, indent=1))
# ---- figure
INK, AXIS, GRID, BLUE, GREEN, DARK, MUTED, ORANGE = "#0b0b0b", "#c3c2b7", "#e1e0d9", "#2a78d6", "#1a9e74", "#0d366b", "#898781", "#eb6834"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig = plt.figure(figsize=(16, 8.5)); gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.3)
ax = fig.add_subplot(gs[0, 0]); xx = np.linspace(-14, 4, 400); from scipy.stats import norm
for lab, (m, s), c in [("prior: Edwards et al. 1983 spawning curve, weak", out["velocity_coef_per_m_s"]["prior"], MUTED), ("data alone (59 nests, flat prior)", out["velocity_coef_per_m_s"]["data_alone"], GREEN), ("posterior", out["velocity_coef_per_m_s"]["posterior"], DARK)]: ax.plot(xx, norm.pdf(xx, m, s), color=c, lw=1.4, label=f"{lab}: {m:.1f} ± {s:.1f}")
ax.axvline(0, color=INK, lw=0.6); ax.set_xlabel("velocity coefficient (log-odds of use per m/s, at the same bank distance and depth)"); ax.set_ylabel("density"); ax.set_title("The Bayesian update on the velocity effect", fontsize=8.5); ax.legend(fontsize=7); ax.grid(True); ax.set_axisbelow(True)
ax = fig.add_subplot(gs[0, 1]); vv = np.linspace(0, 0.8, 80); Sv = np.exp(np.outer(bv, vv)); ax.fill_between(vv, np.quantile(Sv, .05, 0), np.quantile(Sv, .95, 0), color=DARK, alpha=0.15, lw=0); ax.plot(vv, np.median(Sv, 0), color=DARK, lw=1.5, label="posterior relative suitability S(v)"); ax.plot(vv, np.exp(pm[1] * vv), color=MUTED, lw=1, ls="--", label="prior curve (fit to the 1983 spawning SI curve)"); ax.plot(vv_lit, S_lit, color=ORANGE, lw=1, ls=":", label="Edwards et al. 1983 spawning SI curve as digitised")
ax.axhline(0.5, color=INK, lw=0.5, ls=":"); ax.axvline(0.2, color=ORANGE, lw=0.8); ax.set_xlabel("depth-averaged velocity (m/s)"); ax.set_ylabel("use relative to still water"); ax.set_title(f"The suitability curve: half use at v50 = {out['v50_m_s_posterior']['median']:.2f} m/s [{out['v50_m_s_posterior']['ci90'][0]:.2f}, {out['v50_m_s_posterior']['ci90'][1]:.2f}] (orange: the 0.2 m/s convention)", fontsize=8); ax.legend(fontsize=7); ax.grid(True); ax.set_axisbelow(True); ax.set_ylim(0, 1.05)
ax = fig.add_subplot(gs[0, 2]); ax.fill_between(curve.stage_m, curve.E_lo, curve.E_hi, color=BLUE, alpha=0.2, lw=0); ax.plot(curve.stage_m, curve.E_med, color=BLUE, lw=1.5, label="reach mean suitability E, 90% band"); ax2 = ax.twinx(); ax2.plot(curve.stage_m, curve["P_E_above_0.5"], color=ORANGE, lw=1, ls="--", label="P(E > 0.5)"); ax2.set_ylim(0, 1.05); ax2.set_ylabel("P(E > 0.5)", color=ORANGE)
ax.axvspan(0.74, 1.09, color=GRID, alpha=0.5, lw=0); ax.set_xlabel("stage (m)"); ax.set_ylabel("reach mean relative suitability E (0 to 1)"); ax.set_ylim(0, 1.05); ax.set_title("The threshold-free ecological number against stage\n(the curve's posterior and the velocity chain, jointly)", fontsize=8.5); ax.legend(loc="upper right", fontsize=7); ax.grid(True); ax.set_axisbelow(True)
for k, yr in enumerate([2025, 2026]):
    ax = fig.add_subplot(gs[1, k]); s = ser[ser.time_utc.dt.year == yr]; ax.fill_between(s.time_utc, s.E_lo, s.E_hi, color=BLUE, alpha=0.2, lw=0); ax.plot(s.time_utc, s.E_med, color=BLUE, lw=0.8); ax2 = ax.twinx(); ax2.plot(s.time_utc, s.stage_m, color=ORANGE, lw=0.5, alpha=0.8); ax2.set_ylim(0.5, 1.9); ax2.set_ylabel("stage (m)", color=ORANGE)
    ax.set_ylim(0, 1); ax.set_ylabel("reach mean suitability E"); ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b")); ax.set_title(f"{yr}: E every two hours (median over the record {out['reach_mean_suitability'][str(yr)]['median']:.2f})", fontsize=8.5); ax.grid(True); ax.set_axisbelow(True)
ax = fig.add_subplot(gs[1, 2]); ax.text(0.0, 0.95, "What this is\n\nS(v) is the odds that a cell holds a nest at velocity v relative to still water,\nat the same bank distance and depth (used against available, 59 nests, 20\navailable cells each on the nest's day; Gaussian priors; Laplace posterior).\n\nE is the wetted-width-weighted mean of S over the reach: 1 if the whole reach\nwere still water, 0 if all of it were too fast for a nest. No threshold enters.\n\nThe prior on the velocity effect is the USFWS spawning suitability curve\n(Edwards et al. 1983, judgement-based), fitted by an exponential and made weak:\nits sd puts the steep adult utilisation curve of the same report two sd away.\n\nThe data alone put the effect at " + f"{out['velocity_coef_per_m_s']['data_alone'][0]:.1f} ± {out['velocity_coef_per_m_s']['data_alone'][1]:.1f} per m/s; the posterior at\n{out['velocity_coef_per_m_s']['posterior'][0]:.1f} ± {out['velocity_coef_per_m_s']['posterior'][1]:.1f}. Bank distance stays the strongest selector\n({out['bank_distance_coef_per_m']['posterior'][0]:.2f} per m); depth enters positively ({out['depth_coef_per_m']['posterior'][0]:.1f} per m): nests\nsit in deep, slow water near the bank.", transform=ax.transAxes, va="top", fontsize=7.5); ax.axis("off")
fig.suptitle("The ecological number as a Bayesian suitability: a literature prior updated by the nests, carried through the velocity chain", fontsize=10, y=0.95); fig.savefig(F / "11_nest_bayes.png", bbox_inches="tight"); plt.close(fig); print("saved 11")
