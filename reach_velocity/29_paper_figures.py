"""The five paper figures, built from the saved tables (no refitting). Outputs paper/F1..F5 (.png at 300 dpi and .pdf).
F1 site, data, target. F2 the conveyance field and how well it is known. F3 one stage sensor carries the field through time.
F4 spatial rating curves for two seasons. F5 the ecological number and the nests."""
import json, sys, warnings
from pathlib import Path
from uncertainty import PRIMARY   # the primary stage law, read from the saved law grid (entry 30)
import numpy as np, pandas as pd, geopandas as gpd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import matplotlib.dates as mdates; from matplotlib.colors import ListedColormap
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]; HERE = Path(__file__).resolve().parent; D = HERE / "data_clean"; RAW = HERE / "data_raw"; P = HERE / "paper"; P.mkdir(exist_ok=True)
FT = 0.3048; ADCP_E, ADCP_N = 13244857.0, 319995.5; H0, GAM = PRIMARY; Z0 = json.load(open(D / "field_qc.json"))["z0_m"]; THR = 0.20; SURVEYED = (0.74, 1.09)
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, DARK, MUTED, BEIGE = "#0b0b0b", "#52514e", "#b9b8ae", "#e4e3dc", "#2a78d6", "#eb6834", "#1a9e74", "#0d366b", "#898781", "#e9e5dc"
C25, C26 = BLUE, ORANGE
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.edgecolor": AXIS, "axes.linewidth": 0.5, "grid.color": GRID, "grid.linewidth": 0.4, "legend.frameon": False, "axes.facecolor": "white", "axes.titlelocation": "left", "xtick.major.size": 2, "ytick.major.size": 2, "axes.titlepad": 3})
def letter(ax, s, x=-0.2, y=1.0): ax.text(x, y, s, transform=ax.transAxes, fontweight="bold", fontsize=9, va="bottom", ha="left")
def save(fig, name): fig.savefig(P / f"{name}.png", dpi=300, bbox_inches="tight"); fig.savefig(P / f"{name}.pdf", bbox_inches="tight"); plt.close(fig); print("saved", name)
def vbar_at(logA, bed, stage):
    h = np.clip(stage - bed, 0, None); a = np.exp(logA) * np.clip(stage - H0, 0.02, None) ** GAM; wet = h > np.e * Z0; return np.where(wet, a * (np.log(np.clip(h, np.e * Z0, None) / Z0) - 1), np.nan)
# ---------------- inputs
loc = pd.read_csv(D / "locations.csv"); loc = loc[loc.a > 0].copy(); loc["A"] = loc.a / np.clip(loc.stage_m - H0, 0.02, None) ** GAM
rd = pd.read_csv(D / "readings_transect.csv"); ser = pd.read_csv(D / "series_10min.csv", parse_dates=["time_utc"]); ser["time_utc"] = pd.to_datetime(ser.time_utc, utc=True)
grid = pd.read_csv(D / "grid.csv"); sec_along = loc.groupby("transect").along_m.median().sort_values(); order = sec_along.index.tolist()
bk = gpd.read_file(D / "banks_c" / "C_Middle.geojson"); banks = [(np.array(g.coords)[:, 0] - ADCP_E) * FT for g in bk.geometry], [(np.array(g.coords)[:, 1] - ADCP_N) * FT for g in bk.geometry]   # bank set C, fixed at the median survey-day level (19_, entry 31)
def draw_banks(ax, lw=0.5): [ax.plot(x, y, color="#8a877f", lw=lw) for x, y in zip(*banks)]
def draw_sections(ax, lw=0.5, color=INK, alpha=0.7):
    for t, g in loc.groupby("transect"): g = g.sort_values("across_m"); ax.plot(g.x_m, g.y_m, color=color, lw=lw, alpha=alpha)
survey_days = sorted(loc.date.unique()); adcp = pd.read_csv(D / "step3_adcp.csv", parse_dates=["time_utc"]); adcp["time_utc"] = pd.to_datetime(adcp.time_utc, utc=True)
# ================= F1: site, data, target
fig = plt.figure(figsize=(7.1, 6.2)); gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1], width_ratios=[1.35, 1, 1], hspace=0.5, wspace=0.42)
ax = fig.add_subplot(gs[0, 0]); draw_banks(ax)
for yr, c in [(2025, C25), (2026, C26)]:
    for t, g in loc[loc.year == yr].groupby("transect"): g = g.sort_values("across_m"); ax.plot(g.x_m, g.y_m, color=c, lw=1.0); ax.scatter(g.x_m, g.y_m, s=3, color=c, lw=0, zorder=3)
ax.scatter([0], [0], marker="*", s=70, color=INK, zorder=5); ax.annotate("stage sensor\nand ADCP", (0, 0), xytext=(28, -55), textcoords="offset points", fontsize=6.5, arrowprops=dict(arrowstyle="-", lw=0.5, color=INK))
ax.plot([], [], color=C25, lw=1, label="2025: 13 sections, 9 days"); ax.plot([], [], color=C26, lw=1, label="2026: 8 sections, 2 days"); ax.legend(loc="upper center", fontsize=6, bbox_to_anchor=(0.5, -0.12), ncol=1)
ax.annotate("", xy=(loc.x_m.min() - 25, loc.y_m.min() + 40), xytext=(loc.x_m.min() - 25, loc.y_m.min() + 110), arrowprops=dict(arrowstyle="->", lw=0.6, color=INK2)); ax.text(loc.x_m.min() - 20, loc.y_m.min() + 70, "flow", fontsize=6.5, color=INK2)
ax.set_aspect("equal"); ax.set_xlabel("east of the sensor (m)"); ax.set_ylabel("north of the sensor (m)"); ax.set_title("The reach and its 20 sections"); ax.grid(True); ax.set_axisbelow(True); letter(ax, "a", x=-0.42)
gsb = gs[0, 1:].subgridspec(2, 1, hspace=0.55)
for k, yr in enumerate([2025, 2026]):
    ax = fig.add_subplot(gsb[k]); s = ser[ser.time_utc.dt.year == yr]; ax.plot(s.time_utc, s.stage_m, color=INK, lw=0.5)
    if yr == 2025: ax.axvspan(adcp.time_utc.min(), adcp.time_utc.max(), color=GREEN, alpha=0.12, lw=0); ax.text(adcp.time_utc.min() + pd.Timedelta(days=2), 1.42, "ADCP recording", fontsize=6, color=GREEN)
    for d in [d for d in survey_days if d.startswith(str(yr))]: ax.axvline(pd.Timestamp(d, tz="UTC"), color=C25 if yr == 2025 else C26, lw=0.7, ymin=0.9)
    ax.axhspan(*SURVEYED, color=GRID, alpha=0.6, lw=0); ax.set_ylim(0.55, 1.75); ax.set_ylabel("stage (m)"); ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b")); ax.grid(True); ax.set_axisbelow(True)
    ax.set_title(f"{yr}: stage every 10 min (ticks: survey days; grey: the surveyed stage range)" if yr == 2025 else f"{yr}: the surveys caught 0.74 to 0.95 m of a season that reached 1.65 m")
    if k == 0: letter(ax, "b", x=-0.09)
ax = fig.add_subplot(gs[1, 0]); picks = [("2025-07-21|AT10|AT1014", C25, "slow"), ("2025-07-14|AT07|AT0710", DARK, "typical"), ("2026-06-25|0625-T03|0625-30", C26, "fast")]
for lid, c, lab in picks:
    p = rd[rd.location_id == lid].sort_values("height_m"); a = float(loc.set_index("location_id").a[lid]); z = np.linspace(Z0 * 1.05, p.height_m.max() * 1.1, 60); ax.plot(a * np.log(z / Z0), z, color=c, lw=1); ax.scatter(p.velocity_ms, p.height_m, s=14, color=c, edgecolor="white", lw=0.4, zorder=3, label=f"{lab}: a = {a:.3f} m/s")
ax.set_xlabel("velocity (m/s)"); ax.set_ylabel("height above the bed (m)"); ax.set_title(f"One number per location: the log-law\namplitude a, with z0 = {Z0*1000:.0f} mm for the reach"); ax.legend(loc="lower right"); ax.grid(True); ax.set_axisbelow(True); ax.set_xlim(0, None); ax.set_ylim(0, None); letter(ax, "c")
for k, (col, lab, ttl) in enumerate([("a", "log a, the amplitude on its day", "Before: the day's stage\nis in the target"), ("A", "log A = log a − log(stage − 0.45)", "After: the conveyance,\na property of the place")]):
    ax = fig.add_subplot(gs[1, 1 + k]); dm = loc.groupby("date").agg(st=("stage_m", "median"), med=(col, lambda v: np.log(v).median()), lo=(col, lambda v: np.log(v).quantile(.25)), hi=(col, lambda v: np.log(v).quantile(.75)), yr=("year", "first"))
    for d, r in dm.iterrows(): ax.errorbar(r.st, r.med, yerr=[[r.med - r.lo], [r.hi - r.med]], fmt="o", ms=3.5, color=C25 if r.yr == 2025 else C26, elinewidth=0.7, capsize=0)
    if col == "a": xx = np.linspace(0.72, 1.11, 30); ax.plot(xx, GAM * np.log(xx - H0) + dm.med.mean() - GAM * np.mean(np.log(dm.st - H0)), color=INK, lw=0.6, ls="--"); ax.text(0.86, -3.45, f"dashed: a = A (stage − {H0})^{GAM:g}", fontsize=6.2)
    ax.set_xlabel("stage on the survey day (m)"); ax.set_ylabel(lab); ax.set_title(ttl); ax.grid(True); ax.set_axisbelow(True); ax.set_ylim(-3.6, -1.2); letter(ax, "de"[k])
save(fig, "F1_site_data_target")
# ================= F2: what is learned (the field is at-a-station geometry, and its held-out error)
fig = plt.figure(figsize=(7.1, 5.8)); gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.32, width_ratios=[1, 1])
fe = pd.read_csv(D / "features.csv"); fe = fe[fe.A > 0]; fe["A_law"] = loc.set_index("location_id").A.reindex(fe.location_id).values
ax = fig.add_subplot(gs[0, 0]); sc = ax.scatter(fe.dist_bank_m, fe.A_law, c=fe.section_mean_depth_m, cmap="viridis", s=12, edgecolor=INK, lw=0.2); cb = fig.colorbar(sc, ax=ax, shrink=0.7, pad=0.02); cb.set_label("section mean depth at 0.90 m (m)", fontsize=6.5); cb.ax.tick_params(labelsize=6)
ax.set_yscale("log"); ax.set_xlabel("distance to the nearer bank (m)"); ax.set_ylabel(f"place factor A (m/s per (m of stage above {H0} m)^{GAM:g})"); ab = pd.read_csv(D / "step2_feature_ablation.csv").set_index("features"); ax.set_title("The place factor A against bank distance\n(alone it predicts a held-out section as well as six features)"); ax.grid(True, which="both"); ax.set_axisbelow(True); letter(ax, "a")
pr_all = pd.read_csv(D / "step2_predictions.csv"); pr_all = pr_all[pr_all.scheme == "transect_out"]; m2 = pd.read_csv(D / "step2_metrics.csv"); m_to = m2[m2.scheme == "transect_out"].set_index("model"); m_rf = m2[m2.scheme == "random_10fold"].set_index("model")
ax = fig.add_subplot(gs[0, 1]); tsec = "AT0401"; pf = pr_all[(pr_all.model == "full features + notes") & (pr_all.transect == tsec)].merge(loc[["location_id", "across_m"]], on="location_id").sort_values("across_m"); pm = pr_all[(pr_all.model == "reach mean") & (pr_all.transect == tsec)].merge(loc[["location_id", "across_m"]], on="location_id").sort_values("across_m")
ax.fill_between(pf.across_m, pf.mu - 1.96 * pf.sd, pf.mu + 1.96 * pf.sd, color=BLUE, alpha=0.2, lw=0, label="six features, 95% band"); ax.plot(pf.across_m, pf.mu, color=BLUE, lw=1.2); ax.plot(pm.across_m, pm.mu, color=MUTED, lw=1, ls="--", label="reach mean"); ax.scatter(pf.across_m, pf.y, s=16, color=INK, zorder=4, label="measured")
ax.set_xlabel("across the section (m)"); ax.set_ylabel("log A"); ax.set_title(f"Section {tsec} held out, predicted from the other 19\n(the field beats the reach mean on 16 of 20 sections)"); ax.legend(loc="lower left"); ax.grid(True); ax.set_axisbelow(True); letter(ax, "b")
ax = fig.add_subplot(gs[1, 0]); models = [("reach mean", "reach mean"), ("generic shape (bank distance)", "bank distance"), ("full features + notes", "six features + notes")]; xs = np.arange(3)
ax.bar(xs - 0.18, [m_to.CRPS_log[m] for m, _ in models], 0.36, color=BLUE, label="section held out"); ax.bar(xs + 0.18, [m_rf.CRPS_log[m] for m, _ in models], 0.36, color=MUTED, label="random folds (leak neighbours)")
for i, (m, _) in enumerate(models): ax.text(xs[i] - 0.18, m_to.CRPS_log[m] + 0.005, f"{m_to.CRPS_log[m]:.2f}", ha="center", fontsize=6); ax.text(xs[i] + 0.18, m_rf.CRPS_log[m] + 0.005, f"{m_rf.CRPS_log[m]:.2f}", ha="center", fontsize=6, color=INK2)
ax.set_xticks(xs); ax.set_xticklabels([n for _, n in models]); ax.set_ylabel("CRPS of log A (lower is better)"); ax.set_ylim(0, 0.5); ax.set_title("Held-out score by model\n(random folds leak neighbours and flatter every model)"); ax.legend(loc="upper right"); ax.grid(True, axis="y"); ax.set_axisbelow(True); letter(ax, "c")
ax = fig.add_subplot(gs[1, 1]); sc = ax.scatter(grid.x_m, grid.y_m, c=np.exp(grid.logA), cmap="Blues", s=1.2, lw=0, vmin=0, vmax=np.nanquantile(np.exp(grid.logA), .99)); draw_banks(ax); draw_sections(ax, lw=0.4, alpha=0.5)
cb = fig.colorbar(sc, ax=ax, shrink=0.75, pad=0.02); cb.set_label("place factor A", fontsize=6.5); cb.ax.tick_params(labelsize=6); ax.set_aspect("equal"); ax.set_title(f"The field on a 2 m grid, an illustration\n(sd of log A {np.median(grid.logA_sd[grid.dist_to_transect_m > 15]):.2f} between sections, {np.median(grid.logA_sd[grid.dist_to_transect_m < 3]):.2f} on them)"); ax.tick_params(labelsize=6); ax.set_xlabel("east (m)"); ax.set_ylabel("north (m)"); ax.grid(True); ax.set_axisbelow(True); letter(ax, "d", x=-0.45)
save(fig, "F2_conveyance_field")
# ================= F3: identification of the stage law, and closure
fig = plt.figure(figsize=(7.1, 5.8)); gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.32)
ax = fig.add_subplot(gs[0, 0]); off = pd.read_csv(D / "step3_day_offsets.csv"); o = off[off.law == "none"].copy(); xx = np.linspace(0.72, 1.11, 40); lawc = GAM * np.log(xx - H0) - GAM * np.mean(np.log(o.stage_m - H0)) + o.offset_log.mean(); after = o.offset_log - (np.log(o.stage_m - H0) - np.mean(np.log(o.stage_m - H0)) + o.offset_log.mean())
ax.plot(xx, lawc, color=INK, lw=0.8, ls="--", label=f"the law: a = A (stage − {H0} m)^{GAM:g}")
for _, r in o.iterrows(): ax.scatter(r.stage_m, r.offset_log, s=24, color=C25 if r.year == 2025 else C26, edgecolor=INK, lw=0.3, zorder=3)
ax.scatter(o.stage_m, after, s=20, facecolor="white", edgecolor=INK2, lw=0.6, zorder=3, label="after the law")
ax.set_xlabel("stage on the survey day (m)"); ax.set_ylabel("day held out: mean residual of log a"); ax.set_title(f"Each day held out: the day level before the law\n(filled, sd {o.offset_log.std(ddof=1):.2f}) and after it (open, sd {after.std(ddof=1):.2f})"); ax.legend(loc="upper left"); ax.grid(True); ax.set_axisbelow(True); letter(ax, "a")
ax = fig.add_subplot(gs[0, 1]); ab = pd.read_csv(D / "adcp_single_beam_series.csv", parse_dates=["time_utc"]); ab["time_utc"] = pd.to_datetime(ab.time_utc, utc=True); abj = json.load(open(D / "adcp_single_beam.json"))
ax.fill_between(ab.time_utc, ab.v_lo, ab.v_hi, color=BLUE, alpha=0.2, lw=0, label="the chain at the instrument's cell, 90% band"); ax.plot(ab.time_utc, ab.v_pred, color=BLUE, lw=0.7); ax.plot(ab.time_utc, ab.u_beam_ms, color=INK, lw=0.6, label="ADCP record read as one 25° beam")
ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y")); ax.set_ylabel("depth-averaged velocity at the ADCP (m/s)"); ax.set_title(f"The ADCP through 2025, never used in the fit: level ratio {abj['level_ratio_record_over_field_median']:.2f},\nMAE {abj['MAE_m_s']*1000:.0f} mm/s, the band covers {abj['coverage90']*100:.0f}% of {abj['n_2h']:,} two-hour means to 1.20 m"); ax.legend(loc="upper right"); ax.grid(True); ax.set_axisbelow(True); ax.set_ylim(0, 0.9); letter(ax, "b")
ax = fig.add_subplot(gs[1, 0]); rh = pd.read_csv(D / "step6_readings_heldout_adopted_law.csv"); ax.scatter(rh.velocity_ms, rh.v_pred, s=5, color=BLUE, alpha=0.45, lw=0); ax.plot([0, 0.8], [0, 0.8], color=AXIS, lw=0.7)
mae = np.mean(np.abs(rh.velocity_ms - rh.v_pred)); mad = np.mean(np.abs(rh.velocity_ms - rh.velocity_ms.median())); cov = np.mean((rh.velocity_ms >= rh.v_lo) & (rh.velocity_ms <= rh.v_hi)); below = rh.velocity_ms < THR; hit = np.mean((rh.v_pred < THR) == below); base = max(below.mean(), 1 - below.mean())
ax.set_xlim(0, 0.8); ax.set_ylim(0, 0.8); ax.set_xlabel("measured velocity (m/s), each reading"); ax.set_ylabel("predicted at its height, its day held out (m/s)"); ax.set_title(f"{len(rh):,} readings on held-out days: MAE {mae:.3f} m/s against\na deviation of {mad:.3f}; below/above {THR} m/s right {hit*100:.0f}% (base {base*100:.0f}%)"); ax.grid(True); ax.set_axisbelow(True); letter(ax, "c")
ax = fig.add_subplot(gs[1, 1]); cs = pd.read_csv(D / "step6_continuity_common_stage.csv"); cn = pd.read_csv(D / "step6_continuity.csv"); rr = json.load(open(D / "step6_reach_rating.json"))
for t, q in cs.groupby("section"): ax.plot(q.stage_m, q.Q_med, color=BLUE, lw=0.5, alpha=0.5)
xx = np.linspace(0.70, 1.30, 40); ax.plot(xx, rr["coef"] * (xx - rr["h0"]) ** rr["exponent"], color=INK, lw=1.0, ls="--", label=f"rating through the measured points: Q = {rr['coef']:.0f} (stage − 0.45)^{rr['exponent']:.2f}, {rr['scatter_pct']:.0f}% scatter")
ax.scatter(cn.stage_m, cn.Q_measured_m3s, s=16, marker="s", color=INK, zorder=4, label="measured from the profiles, 20 sections on 11 days"); ax.plot([], [], color=BLUE, lw=0.5, label="each section's modelled discharge at common stages")
sp = cs.groupby("stage_m").Q_med.agg(lambda x: np.std(np.log(x), ddof=1)); ax.set_xlabel("stage (m)"); ax.set_ylabel("discharge (m³/s)"); ax.set_title(f"Continuity at common stage: the 20 sections' discharges\nspread by {sp[0.90]*100:.0f}% at 0.90 m and {sp[1.30]*100:.0f}% at 1.30 m"); ax.legend(loc="upper left", fontsize=6); ax.grid(True); ax.set_axisbelow(True); ax.set_ylim(0, 30); letter(ax, "d")
save(fig, "F3_stage_law")
# ================= F4: spatial rating curves for two seasons
ss = pd.read_csv(D / "step6_section_series.csv", parse_dates=["time_utc"]); ss["time_utc"] = pd.to_datetime(ss.time_utc, utc=True)
fig = plt.figure(figsize=(7.1, 7.4)); gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.15], hspace=0.5, wspace=0.3); cax = fig.add_axes([0.955, 0.66, 0.012, 0.18])
for k, yr in enumerate([2025, 2026]):
    ax = fig.add_subplot(gs[k, :]); s = ss[ss.time_utc.dt.year == yr]; piv = s.pivot(index="section", columns="time_utc", values="mean_v_med").reindex(order); tt = piv.columns
    xe = mdates.date2num(tt.tz_convert(None).to_pydatetime()); xe = np.r_[xe, xe[-1] + (xe[-1] - xe[-2])]; im = ax.pcolormesh(xe, np.arange(len(order) + 1), piv.values, cmap="Blues", vmin=0, vmax=0.6, shading="flat"); ax.set_yticks(np.arange(len(order)) + 0.5); ax.set_yticklabels([f"{t} ({sec_along[t]:+.0f} m)" for t in order], fontsize=5.2); ax.invert_yaxis()
    st = s[s.section == order[0]]; out_ = (st.stage_m > SURVEYED[1]) | (st.stage_m < SURVEYED[0]); ax.scatter(mdates.date2num(st.time_utc[out_].dt.tz_convert(None).dt.to_pydatetime()), np.full(out_.sum(), len(order) + 0.4), s=3, marker="|", color=INK2, lw=0.5, clip_on=False)
    ax2 = ax.twinx(); ax2.plot(mdates.date2num(st.time_utc.dt.tz_convert(None).dt.to_pydatetime()), st.stage_m, color=ORANGE, lw=0.7); ax2.set_ylim(0.5, 1.9); ax2.set_ylabel("stage (m)", color=ORANGE, labelpad=1); ax2.tick_params(axis="y", colors=ORANGE, labelsize=6)
    ax.xaxis_date(); ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b")); ax.set_title(f"{yr}: mean velocity of every section every two hours from the stage record (rows: sections, upstream at the top; ticks below: stage outside the surveyed range)"); letter(ax, "ab"[k], x=-0.13)
    if k == 1: cb = fig.colorbar(im, cax=cax); cb.set_label("section mean velocity (m/s)", fontsize=6.5); cb.ax.tick_params(labelsize=6)
ax = fig.add_subplot(gs[2, 0])
for t, c in zip(["AT10", "0625-T02", "AT0401"], [BLUE, GREEN, DARK]):
    s = ss[ss.section == t].sort_values("stage_m"); b = pd.cut(s.stage_m, np.linspace(0.62, 1.7, 55)); q = s.groupby(b, observed=True).agg(st=("stage_m", "median"), med=("mean_v_med", "median"), lo=("mean_v_lo", "median"), hi=("mean_v_hi", "median"))
    ax.fill_between(q.st, q.lo, q.hi, color=c, alpha=0.18, lw=0); ax.plot(q.st, q.med, color=c, lw=1.1, label=f"{t} ({sec_along[t]:+.0f} m)"); g = loc[loc.transect == t]; ax.scatter(g.stage_m.median(), g.v_mean.mean(), s=28, marker="s", color=c, edgecolor="white", lw=0.5, zorder=4)
ax.axvspan(*SURVEYED, color=GRID, alpha=0.6, lw=0); ax.axvspan(SURVEYED[1], 1.7, facecolor="none", hatch="////", edgecolor=GRID, lw=0); ax.text(0.755, 0.02, "surveyed range", fontsize=6, color=INK2); ax.text(1.25, 0.02, "extrapolated: the law set diverges", fontsize=6, color=INK2)
ax.set_xlabel("stage (m)"); ax.set_ylabel("section mean velocity (m/s)"); ax.set_title("Three sections' rating curves, 90% bands\n(squares: measured on the survey day; hatched: extrapolated)"); ax.legend(loc="upper left"); ax.grid(True); ax.set_axisbelow(True); ax.set_xlim(0.65, 1.7); ax.set_ylim(0, 1.0); letter(ax, "c")
ax = fig.add_subplot(gs[2, 1]); at = pd.read_csv(D / "step4_attribution.csv"); srcs = ["A_field", "disc", "law", "stage", "bed", "z0"]; names = ["field", "discrep-\nancy", "law", "stage", "bed", "z0"]; W_ = 0.38
for k, (st, c) in enumerate([(0.90, BLUE), (1.30, DARK)]):
    a_ = at[(at.section == "AT10") & (at.stage_m == st)].set_index("source"); ax.bar(np.arange(6) + (k - 0.5) * W_, [a_.mean_v_iqr90[s_] * 100 for s_ in srcs], W_, color=c, label=f"at {st:.2f} m ({'inside' if st < 1.1 else 'above'} the surveyed range)")
ax.set_xticks(np.arange(6)); ax.set_xticklabels(names, fontsize=6.5); ax.set_ylabel("90% band on the section mean velocity\nfrom this source alone (cm/s)"); ax.set_title("What the band is made of (section AT10)\n(inside the range the field and the discrepancy; above it the law)"); ax.legend(loc="upper right"); ax.grid(True, axis="y"); ax.set_axisbelow(True); letter(ax, "d")
save(fig, "F4_spatial_rating_curves")
# ================= F5: the ecological number, learned from the nests
rs = pd.read_csv(D / "step7_reach_series.csv", parse_dates=["time_utc"]); rs["time_utc"] = pd.to_datetime(rs.time_utc, utc=True); nbs = pd.read_csv(D / "nest_bayes_series.csv", parse_dates=["time_utc"]); nbs["time_utc"] = pd.to_datetime(nbs.time_utc, utc=True); nb = json.load(open(D / "nest_bayes.json")); nh = json.load(open(D / "nest_habitat_model.json")); s7 = json.load(open(D / "step7_summary.json"))
fig = plt.figure(figsize=(7.1, 6.0)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=0.5, wspace=0.3)
for k, yr in enumerate([2025, 2026]):
    ax = fig.add_subplot(gs[0, k]); r = rs[rs.time_utc.dt.year == yr]; e = nbs[nbs.time_utc.dt.year == yr]
    ax.fill_between(e.time_utc, e.E_lo * 100, e.E_hi * 100, color=DARK, alpha=0.15, lw=0); ax.plot(e.time_utc, e.E_med * 100, color=DARK, lw=0.9, label="mean suitability E, learned from the nests")
    ax.fill_between(r.time_utc, r.share_lo * 100, r.share_hi * 100, color=BLUE, alpha=0.15, lw=0); ax.plot(r.time_utc, r.share_med * 100, color=BLUE, lw=0.7, ls="--", label=f"share below {THR} m/s (the convention)")
    oo = r[~r.in_surveyed_range]; ax.scatter(oo.time_utc, np.full(len(oo), 1.5), s=3, marker="|", color=INK2, lw=0.5); ax2 = ax.twinx(); ax2.plot(r.time_utc, r.stage_m, color=ORANGE, lw=0.5, alpha=0.9); ax2.set_ylim(0.5, 1.9); ax2.tick_params(axis="y", colors=ORANGE, labelsize=6)
    if k == 1: ax2.set_ylabel("stage (m)", color=ORANGE)
    for d in [d for d in survey_days if d.startswith(str(yr))]: ax.axvline(pd.Timestamp(d, tz="UTC"), color=INK, lw=0.5, ls=":", alpha=0.6)
    ax.set_ylim(0, 100); ax.xaxis.set_major_locator(mdates.MonthLocator()); ax.xaxis.set_major_formatter(mdates.DateFormatter("%b")); ax.grid(True); ax.set_axisbelow(True); ax.set_title(f"{yr}: the ecological number every two hours, 90% bands\n(E median {nb['reach_mean_suitability'][str(yr)]['median']*100:.0f}%; share below {THR} m/s median {s7['season'][str(yr)]['median_share']*100:.0f}%)")
    if k == 0: ax.set_ylabel("percent of the reach"); ax.legend(loc="upper left", fontsize=6)
    letter(ax, "ab"[k], x=-0.14)
sens = pd.read_csv(D / "step7_sensitivity.csv"); ec = pd.read_csv(D / "nest_bayes_curve.csv")
ax = fig.add_subplot(gs[1, 0]); cmap = plt.get_cmap("Blues")
for i, thr in enumerate([0.10, 0.20, 0.30, 0.40]):
    s_ = sens[(sens.reference == "depth-averaged") & (sens.threshold_m_s == thr)]; ax.plot(s_.stage_m, s_.share_med * 100, color=cmap(0.4 + 0.6 * i / 3), lw=1.0 if thr != THR else 1.5, ls="--", label=f"share below {thr:.2f} m/s")
ax.fill_between(ec.stage_m, ec.E_lo * 100, ec.E_hi * 100, color=DARK, alpha=0.18, lw=0); ax.plot(ec.stage_m, ec.E_med * 100, color=DARK, lw=2.0, label="mean nesting suitability E, 90% band")
ax.axvspan(*SURVEYED, color=GRID, alpha=0.6, lw=0); ax.axvspan(SURVEYED[1], 1.36, facecolor="none", hatch="////", edgecolor=GRID, lw=0); ax.text(0.755, 3, "surveyed range", fontsize=6, color=INK2); ax.text(1.12, 3, "extrapolated", fontsize=6, color=INK2)
ax.set_xlabel("stage at the sensor (m)"); ax.set_ylabel("percent of the reach's wetted width"); ax.set_title("The ecological rating curve: habitat against stage\n(about 8 points of slow water lost per 10 cm of stage)"); ax.legend(loc="upper right", fontsize=6); ax.grid(True); ax.set_axisbelow(True); ax.set_ylim(0, 100); ax.set_xlim(0.65, 1.36); letter(ax, "c")
ax = fig.add_subplot(gs[1, 1]); from scipy.stats import norm; xx = np.linspace(-13, 4, 400); vc = nb["velocity_coef_per_m_s"]
for lab, (m_, s_), c in [("prior: USFWS 1983 spawning curve, weak", vc["prior"], MUTED), (f"the {json.load(open(D / 'nest_bayes.json'))['n_used']} nests alone", vc["data_alone"], GREEN), ("posterior", vc["posterior"], DARK)]: ax.plot(xx, norm.pdf(xx, m_, s_), color=c, lw=1.4, label=f"{lab}: {m_:.1f} ± {s_:.1f}")
ax.axvline(0, color=INK, lw=0.5); ax.set_xlabel("velocity effect on nest use (log-odds per m/s, at the same bank distance and depth)"); ax.set_ylabel("density"); ax.set_title(f"The rule learned from the fish: a published curve\nupdated by the nests (half use at {nb['v50_m_s_posterior']['median']:.2f} m/s [{nb['v50_m_s_posterior']['ci90'][0]:.2f}, {nb['v50_m_s_posterior']['ci90'][1]:.2f}])"); ax.legend(fontsize=6); ax.grid(True); ax.set_axisbelow(True); letter(ax, "d")
save(fig, "F5_ecological_number")
