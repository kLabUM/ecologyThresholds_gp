"""Step 0c (part 2). Clean the MaxBotix stage record.
Two parts. (1) The daily sensor artefact is measured on dry days (nights change by less than 2 cm) as the mean within-day
pattern by month, set to zero at night (23:00 to 05:00 EST) so that the nights are the anchor, and subtracted from every
reading. A fixed correction of at most a few cm cannot remove a storm peak. (2) The corrected series is smoothed with a GP,
stage(t) = seasonal(t) + hydrologic(t) + noise (RBF 10 to 200 days + Matern 3/2 0.02 to 5 days + white), hyperparameters
fitted by marginal likelihood on hourly medians of a 45-day window per season and applied in overlapping 20-day windows.
A marginal-likelihood decomposition with a periodic term was tried first and is kept in NOTES.md: it could not tell the
daily pattern from the storm term (fitted daily amplitude 0.3 cm against 1.5 cm measured), so the template is used.
Outputs: data_clean/stage_clean_10min.csv, stage_clean_hourly.csv, stage_clean_template.csv, stage_clean_kernels.json,
stage_clean_checks.json, figures/00f_stage_clean.png"""
import json, sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy import linalg
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel as C, RBF, Matern, WhiteKernel
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"; F = HERE / "figures"
t0 = time.time(); log = lambda *a: (print(f"[{time.time()-t0:5.0f}s]", *a), sys.stdout.flush())
TZ = "Etc/GMT+5"; FIT_WINDOWS = {2025: ("2025-06-01", "2025-07-16"), 2026: ("2026-05-20", "2026-07-05")}; DRY_NIGHT_CHANGE_M = 0.02
st = pd.read_csv(D / "stage_10min.csv", parse_dates=["time_utc"]).set_index("time_utc").stage_m.dropna().tz_convert(TZ)
df = pd.DataFrame({"raw_m": st.values}, index=st.index); df["hour"] = df.index.hour + df.index.minute / 60
df["night_date"] = (df.index + pd.Timedelta(hours=1)).normalize()      # the night 23:00 to 05:00 belongs to the date it ends on
df["is_night"] = (df.hour >= 23) | (df.hour < 5); df["month"] = df.index.strftime("%Y-%m")
# ---------------------------------------------------------------- (1) the daily artefact from dry days, zero at night
nightmed = df[df.is_night].groupby("night_date").raw_m.median(); d2d = nightmed.diff().abs(); dry_dates = d2d[d2d < DRY_NIGHT_CHANGE_M].index
df["dry"] = df.night_date.isin(dry_dates); df["anom"] = df.raw_m - df.raw_m.rolling("24h", center=True).median()
hb = (df.hour // 0.5 * 0.5)                                            # half-hour bins
tmpl_all = df[df.dry].groupby(hb[df.dry]).anom.mean(); tmpl_all -= tmpl_all[(tmpl_all.index >= 23) | (tmpl_all.index < 5)].mean()
tmpl = {}
for mth, g in df[df.dry].groupby("month"):
    n_dry = g.night_date.nunique()
    if n_dry >= 6: t = g.groupby(hb[g.index]).anom.mean().reindex(tmpl_all.index).interpolate(); t -= t[(t.index >= 23) | (t.index < 5)].mean(); tmpl[mth] = t
    else: tmpl[mth] = tmpl_all
template = pd.DataFrame(tmpl); template["all_months"] = tmpl_all; template.index.name = "hour_est"; template.to_csv(D / "stage_clean_template.csv")
def artefact_of(idx):
    m = idx.strftime("%Y-%m"); h = (idx.hour + idx.minute / 60); hbin = h // 0.5 * 0.5
    return np.array([tmpl.get(mm, tmpl_all).get(hh, tmpl_all.get(hh, 0.0)) for mm, hh in zip(m, hbin)])
df["artefact_m"] = artefact_of(df.index); df["corrected_m"] = df.raw_m - df.artefact_m
log(f"dry nights used: {len(dry_dates)}; template peak-to-trough by month (cm): " + ", ".join(f"{m} {(t.max()-t.min())*100:.1f}" for m, t in tmpl.items()))
# ---------------------------------------------------------------- (2) GP smoothing per season
def make_kernel(): return C(1.0, (1e-2, 1e2)) * RBF(30.0, (10.0, 200.0)) + C(0.3, (1e-3, 1e2)) * Matern(0.3, (0.02, 5.0), nu=1.5) + WhiteKernel(1e-3, (1e-6, 1e-1))
out = []; kern = {}
for yr, (f0, f1) in FIT_WINDOWS.items():
    s = df[df.index.year == yr].corrected_m; origin = s.index[0].normalize()
    hr = s.resample("1h").median().dropna(); th = ((hr.index - origin).total_seconds() / 86400).values; mu, sd = hr.mean(), hr.std(); yh = (hr.values - mu) / sd
    m = (hr.index >= pd.Timestamp(f0, tz=TZ)) & (hr.index < pd.Timestamp(f1, tz=TZ)); Xf, yf = th[m][:, None], yh[m]
    gp = GaussianProcessRegressor(make_kernel(), n_restarts_optimizer=3, random_state=0).fit(Xf, yf); res = yf - gp.predict(Xf); keep = np.abs(res) < 3 * res.std()
    gp = GaussianProcessRegressor(gp.kernel_, n_restarts_optimizer=2, random_state=0).fit(Xf[keep], yf[keep]); k_fit = gp.kernel_
    ls_season, ls_storm = k_fit.k1.k1.k2.length_scale, k_fit.k1.k2.k2.length_scale; noise_cm = float(np.sqrt(k_fit.k2.noise_level) * sd * 100)
    kern[str(yr)] = {"kernel": str(k_fit), "fit_window": [f0, f1], "n_fit": int(keep.sum()), "rejected": int((~keep).sum()), "seasonal_length_scale_days": float(ls_season), "storm_length_scale_days": float(ls_storm), "noise_sd_cm": noise_cm}
    log(f"{yr}: {k_fit}  -> storm scale {ls_storm*24:.1f} h, noise {noise_cm:.2f} cm")
    t10 = ((s.index - origin).total_seconds() / 86400).values; mean = np.full(len(s), np.nan); var = np.full(len(s), np.nan)
    ksig = k_fit.k1                                                    # signal part (seasonal + storm), white noise left out
    for c0 in np.arange(0, t10.max() + 1e-9, 10.0):
        c1 = c0 + 10.0; w = (th >= c0 - 5) & (th < c1 + 5); q = (t10 >= c0) & (t10 < c1)
        if w.sum() < 24 or q.sum() == 0: continue
        g = GaussianProcessRegressor(k_fit, optimizer=None, alpha=1e-10).fit(th[w][:, None], yh[w]); Xq = t10[q][:, None]
        Kq = ksig(Xq, th[w][:, None]); mean[q] = Kq @ g.alpha_; v = linalg.solve_triangular(g.L_, Kq.T, lower=True); var[q] = np.clip(ksig.diag(Xq) - (v ** 2).sum(axis=0), 0, None)
    o = df.loc[s.index, ["raw_m", "artefact_m", "corrected_m", "is_night", "dry"]].copy(); o["stage_clean_m"] = mean * sd + mu; o["stage_clean_sd_m"] = np.sqrt(var) * sd; o["residual_m"] = o.corrected_m - o.stage_clean_m; out.append(o)
all10 = pd.concat(out); all10.index.name = "time_est"; all10["time_utc"] = all10.index.tz_convert("UTC"); all10.to_csv(D / "stage_clean_10min.csv")
all10.resample("1h").mean(numeric_only=True).to_csv(D / "stage_clean_hourly.csv"); json.dump(kern, open(D / "stage_clean_kernels.json", "w"), indent=2)
# ---------------------------------------------------------------- checks
a = all10; hour = a.index.hour + a.index.minute / 60; dday = a.dry.values
def hourprof(x): b = x.rolling("24h", center=True).median(); an = (x - b)[dday]; return an.groupby((hour[dday] // 1)).mean() * 100
hp_raw, hp_clean = hourprof(a.raw_m), hourprof(a.stage_clean_m)
checks = {"dry_nights": int(len(dry_dates)), "dry_day_hour_profile_peak_to_trough_cm": {"raw": float(hp_raw.max() - hp_raw.min()), "cleaned": float(hp_clean.max() - hp_clean.min())},
          "artefact_cm": {"mean": float(a.artefact_m.mean() * 100), "min": float(a.artefact_m.min() * 100), "max": float(a.artefact_m.max() * 100)}, "residual_after_smoothing_sd_cm": float(a.residual_m.std() * 100)}
rise = a.raw_m.rolling("12h").max() - a.raw_m.rolling("12h").min(); ev = []
for t in rise.sort_values(ascending=False).index:
    if len(ev) >= 8: break
    if all(abs((t - e).total_seconds()) > 5 * 86400 for e in ev): ev.append(t)
pk = []
for t in ev:
    w = a[(a.index >= t - pd.Timedelta(hours=12)) & (a.index <= t + pd.Timedelta(hours=12))]; pk.append({"time": str(w.raw_m.idxmax()), "rise_m": float(rise[t]), "raw_max_m": float(w.raw_m.max()), "cleaned_max_m": float(w.stage_clean_m.max())})
checks["peaks"] = pk; checks["peak_change_cm"] = {"mean": float(np.mean([(p["cleaned_max_m"] - p["raw_max_m"]) * 100 for p in pk])), "worst": float(min((p["cleaned_max_m"] - p["raw_max_m"]) * 100 for p in pk))}
nm = a.raw_m[a.is_night].groupby(a.index[a.is_night].to_series().add(pd.Timedelta(hours=1)).dt.normalize().values).median()
cm = a.stage_clean_m[a.is_night].groupby(a.index[a.is_night].to_series().add(pd.Timedelta(hours=1)).dt.normalize().values).mean()
checks["cleaned_minus_raw_night_median_cm"] = {"mean": float((cm - nm).mean() * 100), "sd": float((cm - nm).std() * 100)}
json.dump(checks, open(D / "stage_clean_checks.json", "w"), indent=2); print(json.dumps({k: v for k, v in checks.items() if k != "peaks"}, indent=2)); print(pd.DataFrame(pk).round(3).to_string(index=False))
# ---------------------------------------------------------------- figure
INK, INK2, AXIS, GRID, BLUE, ORANGE, NIGHTC = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#dcdcf5"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig = plt.figure(figsize=(11, 7.6)); gs = fig.add_gridspec(3, 2, height_ratios=[1.3, 0.6, 1.1], hspace=0.55, wspace=0.2)
def shade(ax, t0_, t1_):
    for d in pd.date_range(t0_, t1_, freq="D", tz=TZ): ax.axvspan(d - pd.Timedelta(hours=1), d + pd.Timedelta(hours=5), color=NIGHTC, alpha=0.6, lw=0)
for j, (t0_, t1_, title) in enumerate([("2025-07-27", "2025-08-06", "Ten dry days, 2025"), ("2026-06-12", "2026-06-30", "June 2026, the two survey days")]):
    w = a[(a.index >= pd.Timestamp(t0_, tz=TZ)) & (a.index < pd.Timestamp(t1_, tz=TZ))]
    ax = fig.add_subplot(gs[0, j]); shade(ax, t0_, t1_); ax.plot(w.index, w.raw_m, color=INK, lw=0.5, alpha=0.7, label="raw, 10-min")
    ax.fill_between(w.index, w.stage_clean_m - 1.96 * w.stage_clean_sd_m, w.stage_clean_m + 1.96 * w.stage_clean_sd_m, color=BLUE, alpha=0.2, lw=0); ax.plot(w.index, w.stage_clean_m, color=BLUE, lw=1.0, label="cleaned (artefact removed, GP smoothed), 95% band")
    ax.set_title(title, fontsize=9); ax.set_ylabel("stage (m)"); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=6.5, loc="upper right" if j == 0 else "upper left")
    ax2 = fig.add_subplot(gs[1, j], sharex=ax); shade(ax2, t0_, t1_); ax2.plot(w.index, w.artefact_m * 100, color=ORANGE, lw=0.8); ax2.axhline(0, color=INK, lw=0.5); ax2.set_ylabel("artefact (cm)"); ax2.set_title("what was subtracted: the dry-day template, zero at night", fontsize=8); ax2.grid(True); ax2.set_axisbelow(True)
ax = fig.add_subplot(gs[2, 0]); ax.plot(hp_raw.index, hp_raw.values, color=INK, lw=1.2, label="raw"); ax.plot(hp_clean.index, hp_clean.values, color=BLUE, lw=1.2, label="cleaned"); ax.axvspan(0, 5, color=NIGHTC, alpha=0.6, lw=0); ax.axvspan(23, 24, color=NIGHTC, alpha=0.6, lw=0); ax.axhline(0, color=INK, lw=0.5)
ax.set_xlim(0, 24); ax.set_xticks(range(0, 25, 3)); ax.set_xlabel("hour of day (EST)"); ax.set_ylabel("stage minus 24 h median (cm)"); ax.set_title(f"Within-day pattern on dry days: {checks['dry_day_hour_profile_peak_to_trough_cm']['raw']:.1f} cm before, {checks['dry_day_hour_profile_peak_to_trough_cm']['cleaned']:.1f} cm after", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7)
ax = fig.add_subplot(gs[2, 1]); lo, hi = min(p["raw_max_m"] for p in pk) - 0.05, max(p["raw_max_m"] for p in pk) + 0.05; ax.plot([lo, hi], [lo, hi], color=AXIS, lw=0.8)
ax.scatter([p["raw_max_m"] for p in pk], [p["cleaned_max_m"] for p in pk], s=28, color=BLUE, edgecolor="white", lw=0.5, zorder=3); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
ax.set_xlabel("raw peak (m)"); ax.set_ylabel("cleaned peak (m)"); ax.set_title(f"The eight largest rises: peak change mean {checks['peak_change_cm']['mean']:+.1f} cm, worst {checks['peak_change_cm']['worst']:+.1f} cm", fontsize=9); ax.grid(True); ax.set_axisbelow(True)
fig.suptitle("Cleaning the stage record: dry-day template removed (nights as anchor), then GP smoothing", fontsize=9, y=0.96); fig.savefig(F / "00f_stage_clean.png", bbox_inches="tight"); plt.close(fig); log("saved 00f")
