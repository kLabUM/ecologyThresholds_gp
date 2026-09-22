"""Step 3b. Refine the day law and re-read the ADCP. (1) The reach-wide stage law a ∝ (stage - h0)^gamma over a grid of
h0 and gamma, judged by the day-out test, against the local-depth law. (2) The ADCP as a shape check: under each law the
ratio observed/predicted over the season, its trend with stage, and the exponent the ADCP itself shows under the point
and the column reading. Outputs: step3b_law_grid.csv, step3b_adcp_shape.json, figures/03b_law_refine.png"""
import json, sys, time, warnings
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import matplotlib.dates as mdates
from gpcore import *
warnings.filterwarnings("ignore"); t0 = time.time(); log = lambda *a: (print(f"[{time.time()-t0:5.0f}s]", *a), sys.stdout.flush())
df = load_points(); days = sorted(df.date.unique())
def dayout(law):
    d = with_target(df, law); MU, SD = np.full(len(d), np.nan), np.full(len(d), np.nan)
    for day in days:
        tr, te = np.where(d.date != day)[0], np.where(d.date == day)[0]; m = Model(n_restarts=2).fit(d.iloc[tr]); mu, sd = m.predict(d.iloc[te]); MU[te], SD[te] = mu, sd
    return scores(np.log(d.a.values), MU + np.log(d.g.values), np.sqrt(SD ** 2 + d.sd_y.values ** 2))
rows = []
for h0 in [0.35, 0.45, 0.55, 0.65]:
    for gamma in [0.5, 0.75, 1.0, 1.25]:
        sc = dayout(("stage", gamma, h0)); rows.append({"law": "stage", "h0": h0, "gamma": gamma, **sc}); log(f"stage law h0 {h0} gamma {gamma}: CRPS {sc['CRPS']:.3f} MAE {sc['MAE']:.3f} cov {sc['coverage95']:.2f}")
grid = pd.DataFrame(rows); grid.to_csv(D / "step3b_law_grid.csv", index=False); best = grid.sort_values("CRPS").iloc[0]; log(f"best stage law: h0 {best.h0}, gamma {best.gamma}, CRPS {best.CRPS:.3f}")
# ---- the ADCP as a shape check
ad = pd.read_csv(D / "adcp_amplitude_2h.csv", parse_dates=["time_utc"]); h0_adcp = json.load(open(D / "amplitude_inputs.json"))["h0_adcp_m_column_fit"]
ad["h"] = ad.stage_m - h0_adcp; ad = ad[ad.h > Z0 * 3].copy(); qd = pd.read_csv(D / "usgs_daily.csv").set_index("date").discharge_cfs; ad["Q"] = ad.time_utc.dt.strftime("%Y-%m-%d").map(qd)
res = {}
for name, g in [("reach stage law (stage-%.2f)^%.2f" % (best.h0, best.gamma), np.clip(ad.stage_m - best.h0, 0.02, None) ** best.gamma), ("local depth^0.5 at the ADCP", ad.h ** 0.5), ("discharge^0.75", (ad.Q / 350) ** 0.75)]:
    for reading, denom in [("column", np.log(ad.h / Z0) - 1), ("point at fixed height", np.ones(len(ad)))]:
        a_obs = ad.velocity_ms / denom; ratio = np.log(a_obs / g); slope = np.polyfit(np.log(ad.h), ratio, 1)[0]; res[f"{name} | ADCP read as {reading}"] = {"log_ratio_sd_over_season": float(ratio.std()), "slope_of_log_ratio_vs_log_depth": float(slope)}
for k, v in res.items(): log(f"{k}: after fixing the level, residual sd over the season {v['log_ratio_sd_over_season']:.2f} (log), remaining trend with depth {v['slope_of_log_ratio_vs_log_depth']:+.2f}")
adcp_exp = {"point_reading_exponent_vs_local_depth": float(np.polyfit(np.log(ad.h), np.log(ad.velocity_ms), 1)[0]), "column_reading_exponent_vs_local_depth": float(np.polyfit(np.log(ad.h), np.log(ad.velocity_ms / (np.log(ad.h / Z0) - 1)), 1)[0]), "transects_day_out_exponent_on_stage_minus_h0": float(best.gamma), "h0": float(best.h0)}
json.dump({"best_stage_law": {"h0": float(best.h0), "gamma": float(best.gamma), "CRPS": float(best.CRPS)}, "adcp_shape": res, "adcp_exponents": adcp_exp}, open(D / "step3b_adcp_shape.json", "w"), indent=2)
# ---- figure
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, DARK, MUTED = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#0d366b", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), gridspec_kw={"wspace": 0.3})
ax = axes[0]; piv = grid.pivot(index="gamma", columns="h0", values="CRPS"); im = ax.imshow(piv.values, cmap="Blues_r", origin="lower", aspect="auto"); ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns); ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index)
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]): ax.text(j, i, f"{piv.values[i, j]:.3f}", ha="center", va="center", fontsize=7, color=INK if piv.values[i, j] > piv.values.min() + 0.01 else "white")
ax.set_xlabel("h0 (m), the stage at which the reach's flow would stop"); ax.set_ylabel("exponent γ"); ax.set_title("Reach stage law a = A (stage − h0)^γ: day-out CRPS", fontsize=9)
ax = axes[1]; g_best = np.clip(ad.stage_m - best.h0, 0.02, None) ** best.gamma
for reading, denom, col in [("column average", np.log(ad.h / Z0) - 1, DARK), ("point at fixed height", np.ones(len(ad)), ORANGE)]:
    a_obs = ad.velocity_ms / denom; r = a_obs / g_best; r = r / np.median(r); ax.scatter(ad.h, r, s=3, color=col, alpha=0.4, lw=0, label=f"ADCP read as {reading}")
ax.axhline(1, color=INK, lw=0.6); ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("local depth at the ADCP (m)"); ax.set_ylabel("observed / law, level removed"); ax.set_title(f"Does the ADCP follow the reach law? flat = yes\n(ADCP exponent: point {adcp_exp['point_reading_exponent_vs_local_depth']:.2f}, column {adcp_exp['column_reading_exponent_vs_local_depth']:.2f}; transects {best.gamma:.2f})", fontsize=8.5); ax.grid(True, which="both"); ax.set_axisbelow(True); ax.legend(fontsize=7)
ax = axes[2]; a_obs = ad.velocity_ms; pred = g_best * np.median(a_obs / g_best); ax.plot(ad.time_utc, a_obs, color=INK, lw=0.6, alpha=0.7, label="ADCP record"); ax.plot(ad.time_utc, pred, color=BLUE, lw=1.0, label="reach stage law, level fitted")
ax.set_ylabel("velocity at the ADCP (m/s)"); ax.set_title(f"The season under the reach law with a free level: MAE {np.mean(np.abs(pred - a_obs))*1000:.0f} mm/s", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7); ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
fig.suptitle("Step 3b: the day law refined, and the ADCP as a shape check", fontsize=10, y=1.02); fig.savefig(F / "03b_law_refine.png", bbox_inches="tight"); plt.close(fig); log("saved 03b")
