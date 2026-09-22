"""The model's variables in one figure, drawn from real data: a profile and its amplitude a, the stage law and the place factor A, a cross-section
with the GP's inputs, and the chain from the sensor to the ecological number. Output: figures/17_model_variables.png"""
import json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; from matplotlib.patches import FancyBboxPatch
from pathlib import Path
H = Path(__file__).resolve().parent; D = H / "data_clean"; F = H / "figures"
Z0 = json.load(open(D / "field_qc.json"))["z0_m"]; H0, GAM = 0.45, 1.0; LAWS = [(0.45, 1.0), (0.35, 1.25), (0.55, 0.75), (0.35, 1.0), (0.45, 0.75), (0.65, 0.5), (0.55, 0.5)]
loc = pd.read_csv(D / "locations.csv"); loc = loc[loc.a > 0]; rd = pd.read_csv(D / "readings_transect.csv"); fe = pd.read_csv(D / "features.csv")
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, DARK, MUTED, BEIGE = "#0b0b0b", "#52514e", "#b9b8ae", "#e4e3dc", "#2a78d6", "#eb6834", "#1a9e74", "#0d366b", "#898781", "#e9e5dc"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8.5, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.4, "legend.frameon": False, "savefig.dpi": 220, "axes.facecolor": "white"})
fig = plt.figure(figsize=(14, 9.5)); gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.28)
# (a) one profile, one number
lid = "2025-07-14|AT07|AT0710"; p = rd[rd.location_id == lid].sort_values("height_m"); a = float(loc.set_index("location_id").a[lid]); h = float(p.water_depth_m.iloc[0]); vbar = a * (np.log(h / Z0) - 1)
ax = fig.add_subplot(gs[0, 0]); z = np.linspace(Z0 * 1.001, h, 200); ax.plot(a * np.log(z / Z0), z, color=DARK, lw=2, label="the log law  v(z) = a · ln(z / z0)")
ax.scatter(p.velocity_ms, p.height_m, s=40, color=ORANGE, edgecolor=INK, lw=0.5, zorder=4, label=f"the {len(p)} readings at one location")
ax.axhline(h, color=BLUE, lw=1.2); ax.text(0.01, h + 0.015, f"water surface: depth h = {h:.2f} m", color=BLUE, fontsize=8.5); ax.axhline(0, color=INK2, lw=1.5); ax.text(0.01, -0.035, "bed", color=INK2, fontsize=8.5)
ax.axvline(vbar, color=GREEN, lw=1.2, ls="--"); ax.text(vbar + 0.015, 0.055, f"depth-averaged velocity\nv̄ = a (ln(h/z0) − 1) = {vbar:.2f} m/s", color=GREEN, fontsize=8)
ax.annotate(f"z0 = {Z0*1000:.0f} mm: the height where the\nlog law reaches zero (one value for the reach)", xy=(0.0, Z0), xytext=(0.34, 0.135), fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.7, color=INK2))
ax.annotate(f"a = {a:.3f} m/s: the amplitude, one number\nper location (a = u*/κ in the friction language)", xy=(a * np.log(0.2 / Z0), 0.2), xytext=(0.03, 0.44), fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.7, color=INK2))
ax.set_xlim(0, 0.5); ax.set_ylim(-0.05, h + 0.06); ax.set_xlabel("velocity v (m/s)"); ax.set_ylabel("height above the bed z (m)"); ax.set_title("(a) One profile becomes one number: the amplitude a", fontsize=10, loc="left"); ax.legend(loc="upper right", fontsize=8, bbox_to_anchor=(0.99, 0.9)); ax.grid(True); ax.set_axisbelow(True)
# (b) the stage law and A
ax = fig.add_subplot(gs[0, 1]); dm = loc.groupby("date").agg(st=("stage_m", "median"), a=("a", "median"), lo=("a", lambda v: v.quantile(.25)), hi=("a", lambda v: v.quantile(.75)))
for d_, r in dm.iterrows(): ax.errorbar(r.st, r.a, yerr=[[r.a - r.lo], [r.hi - r.a]], fmt="o", ms=5, color=ORANGE, mec=INK, mew=0.5, elinewidth=0.8, capsize=0)
xx = np.linspace(0.46, 1.65, 100); Aref = float(np.median(loc.a / np.clip(loc.stage_m - H0, 0.02, None)))
for (h0, g_) in LAWS[1:]: A_ = float(np.median(loc.a / np.clip(loc.stage_m - h0, 0.02, None) ** g_)); ax.plot(xx, A_ * np.clip(xx - h0, 0, None) ** g_, color=MUTED, lw=0.7, alpha=0.7)
ax.plot(xx, Aref * np.clip(xx - H0, 0, None) ** GAM, color=DARK, lw=2, label=f"the law  a = A · (stage − h0)^γ   with h0 = {H0} m, γ = {GAM:.0f}")
ax.plot([], [], color=MUTED, lw=0.7, label="six other (h0, γ) pairs that fit the survey days as well: the law set")
ax.axvspan(0.74, 1.09, color=GRID, alpha=0.6, lw=0); ax.text(0.75, 0.305, "surveyed stage range", fontsize=8, color=INK2); ax.axvline(H0, color=INK2, lw=0.8, ls=":"); ax.text(H0 + 0.012, 0.04, "h0: the stage at which\nthe amplitude vanishes", fontsize=8, color=INK2)
ax.annotate("A: the place factor, the amplitude per metre of stage\nabove h0. It does not depend on the day. The GP learns A.", xy=(1.0, Aref * 0.55), xytext=(0.62, 0.24), fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.7, color=INK2))
ax.set_xlim(0.4, 1.65); ax.set_ylim(0, 0.32); ax.set_xlabel("stage at the sensor (m)"); ax.set_ylabel("amplitude a (m/s); orange: the survey days' medians"); ax.set_title("(b) The day's stage is removed by a one-parameter law: the place factor A", fontsize=10, loc="left"); ax.legend(loc="lower right", fontsize=7.5); ax.grid(True); ax.set_axisbelow(True)
# (c) a cross-section with the GP's inputs
ax = fig.add_subplot(gs[1, 0]); t = "AT0401"; q = fe[fe.transect == t].merge(loc[["location_id", "bed_rtk_sensor_datum_m", "bed_elev_m", "stage_m"]], on="location_id").sort_values("across_m"); bed = q.bed_rtk_sensor_datum_m.fillna(q.bed_elev_m).values; y = q.across_m.values; ref = 0.90
ax.fill_between(y, bed, ref, color="#dbe9f8", lw=0); ax.plot(y, bed, color=INK2, lw=1.5); ax.axhline(ref, color=BLUE, lw=1.2); sc = ax.scatter(y, bed + 0.03, c=q.A, cmap="Blues", s=55, edgecolor=INK, lw=0.5, zorder=4, vmin=0)
cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02); cb.set_label("place factor A at each location", fontsize=8)
i = int(np.argmax(ref - bed)); ax.annotate("", xy=(y[i], bed[i]), xytext=(y[i], ref), arrowprops=dict(arrowstyle="<->", lw=0.9, color=DARK)); ax.text(y[i] + 0.4, (bed[i] + ref) / 2, "local depth h at the reference stage\n(relative depth = h / max h in the section)", fontsize=7.5, color=DARK)
j = 4; ax.annotate("", xy=(y[0], ref + 0.06), xytext=(y[j], ref + 0.06), arrowprops=dict(arrowstyle="<->", lw=0.9, color=GREEN)); ax.text(y[j] + 0.3, ref + 0.075, "distance to the nearer bank d", fontsize=7.5, color=GREEN)
ax.text(y.min(), ref + 0.16, f"the GP's inputs at each location: d, relative depth, section mean depth ({(ref - bed).mean():.2f} m here),\nlateral bed slope, position along the reach, relative position across (−1 to +1); notes on substrate and wake in the mean", fontsize=7.5, color=INK2)
ax.set_xlabel("across the section (m)"); ax.set_ylabel("elevation on the sensor datum (m)"); ax.set_ylim(bed.min() - 0.1, ref + 0.3); ax.set_title(f"(c) Where A lives: section {t} at the reference stage of {ref:.2f} m, with the GP's inputs", fontsize=10, loc="left"); ax.grid(True); ax.set_axisbelow(True)
# (d) the chain
ax = fig.add_subplot(gs[1, 1]); ax.axis("off"); ax.set_title("(d) From the sensor to the ecological number, every two hours", fontsize=10, loc="left")
def box(x, y, w, h_, t_, fc, fs=8.5): ax.add_patch(FancyBboxPatch((x, y), w, h_, boxstyle="round,pad=0.015", fc=fc, ec=INK, lw=0.7, transform=ax.transAxes)); ax.text(x + w / 2, y + h_ / 2, t_, ha="center", va="center", fontsize=fs, transform=ax.transAxes)
box(0.02, 0.84, 0.44, 0.13, "stage(t) from the sensor\n(every 10 min, cleaned)", "#fdebd0"); box(0.54, 0.84, 0.44, 0.13, "A(x) from the GP, with its posterior\n(one joint draw over the reach per sample)", "#dbe9f8")
box(0.10, 0.60, 0.80, 0.16, "a(x, t) = A(x) · (stage(t) − h0)^γ\nh(x, t) = stage(t) − bed(x)          v̄(x, t) = a (ln(h / z0) − 1)", "#e3f3ee")
box(0.10, 0.36, 0.80, 0.16, "for each section: the mean of v̄ over its wetted width (width-weighted)\nfor the reach: the length-weighted mean over the 20 sections", "#f1f0ea")
box(0.02, 0.10, 0.44, 0.18, "the share of the wetted width\nbelow a threshold (0.2 m/s: the convention)", "#f5e6e0"); box(0.54, 0.10, 0.44, 0.18, "E = mean over the reach of S(v̄),\nS(v) = exp(b_v v): the nesting suitability\nlearned from the nests (Bayesian, weak prior)", "#e6e6f5")
for (x0, y0, x1, y1) in [(0.24, 0.84, 0.35, 0.76), (0.76, 0.84, 0.65, 0.76), (0.5, 0.60, 0.5, 0.52), (0.35, 0.36, 0.24, 0.28), (0.65, 0.36, 0.76, 0.28)]: ax.annotate("", xy=(x1, y1), xytext=(x0, y0), xycoords="axes fraction", arrowprops=dict(arrowstyle="->", lw=0.9, color=INK))
ax.text(0.5, 0.02, "each sample also draws the bed (6 cm), the stage (1.5 cm), z0 (7 to 16 mm), the law from its set, and a discrepancy calibrated on held-out days: the 90% band", ha="center", fontsize=7.5, color=INK2, transform=ax.transAxes)
fig.suptitle("The model's variables, from real data: a and z0 from each profile, A from the stage law, the GP on the channel's geometry, the chain to the ecological number", fontsize=11, y=0.97); fig.savefig(F / "17_model_variables.png", bbox_inches="tight"); plt.close(fig); print("saved 17")
