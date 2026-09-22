"""Before and after: what the model's target variable was, and what it is now. Figure 00t.
Before: velocity at 1 ft above the bed, one value per reading, a function of (stage, along, across) learned by one GP.
After: the conveyance A = a / sqrt(h), one value per place, a function of (along, across) only; stage and depth are put
back by the law. A has a physical meaning: for uniform flow u* = sqrt(g h S), so a = u*/kappa and A = sqrt(g S)/kappa,
the local friction slope in disguise. The reach's measured water-surface slope gives the number to compare with."""
import json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"; F = HERE / "figures"
KAPPA, G, S_REACH = 0.41, 9.81, 5.0e-4
loc = pd.read_csv(D / "locations.csv"); loc = loc[loc.a > 0].copy(); loc["A"] = loc.a / np.sqrt(loc.water_depth_m)
A_from_slope = np.sqrt(G * S_REACH) / KAPPA; print(f"A implied by the reach slope sqrt(g S)/kappa = {A_from_slope:.3f}; measured A median {loc.A.median():.3f}, quartiles {loc.A.quantile(.25):.3f} to {loc.A.quantile(.75):.3f}; mid-channel (|across| < 8 m) median {loc[loc.across_m.abs() < 8].A.median():.3f}")
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, DARK, MUTED = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#0d366b", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig = plt.figure(figsize=(14, 10)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.35, wspace=0.25)
# ---- top left: BEFORE, the target against stage
ax = fig.add_subplot(gs[0, 0]); dd = loc.groupby("date").agg(stage=("stage_m", "median"), v=("v_1ft", "median"), lo=("v_1ft", lambda x: x.quantile(.25)), hi=("v_1ft", lambda x: x.quantile(.75))).sort_values("stage")
ax.scatter(loc.stage_m + np.random.default_rng(0).uniform(-0.004, 0.004, len(loc)), loc.v_1ft, s=6, color=MUTED, alpha=0.4, lw=0, label="every location")
ax.errorbar(dd.stage, dd.v, yerr=[dd.v - dd.lo, dd.hi - dd.v], fmt="o", color=INK, ms=5, capsize=2, label="day median and quartiles"); b = np.polyfit(np.log(dd.stage), np.log(dd.v), 1)[0]
ax.set_xlabel("stage on the survey day (m)"); ax.set_ylabel("velocity 1 ft above the bed (m/s)"); ax.set_title(f"BEFORE: the target moved with the day's stage (day medians: slope {b:.1f} in log terms)\nso the GP had to learn stage, depth and place at once, from 11 stages", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7, loc="upper left"); ax.set_ylim(0, 0.7)
# ---- top right: AFTER, the target against stage
ax = fig.add_subplot(gs[0, 1]); dd2 = loc.groupby("date").agg(stage=("stage_m", "median"), A=("A", "median"), lo=("A", lambda x: x.quantile(.25)), hi=("A", lambda x: x.quantile(.75))).sort_values("stage")
ax.scatter(loc.stage_m + np.random.default_rng(0).uniform(-0.004, 0.004, len(loc)), loc.A, s=6, color=BLUE, alpha=0.35, lw=0, label="every location")
ax.errorbar(dd2.stage, dd2.A, yerr=[dd2.A - dd2.lo, dd2.hi - dd2.A], fmt="s", color=DARK, ms=5, capsize=2, label="day median and quartiles"); b2 = np.polyfit(np.log(dd2.stage), np.log(dd2.A), 1)[0]
ax.axhline(A_from_slope, color=ORANGE, lw=1.3, ls="--"); ax.text(0.74, A_from_slope + 0.008, f"what the reach's measured water-surface slope implies: sqrt(g S)/κ = {A_from_slope:.2f}", fontsize=7, color=ORANGE)
ax.set_xlabel("stage on the survey day (m)"); ax.set_ylabel("conveyance A = a / sqrt(depth)  (m^0.5/s)"); ax.set_title(f"AFTER: the target should not move with stage (day medians: slope {b2:.1f}, left over to test)\nand it has a physical meaning: the local friction slope, sqrt(g S_local)/κ", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7, loc="upper left"); ax.set_ylim(0, 0.35)
# ---- bottom: the two chains as diagrams
def box(ax, x, y, w, h, text, fc, fs=8, ec=INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.02", fc=fc, ec=ec, lw=0.8)); ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=INK)
def arrow(ax, x0, y0, x1, y1, text=None):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=10, lw=0.8, color=INK2))
    if text: ax.text((x0 + x1) / 2, (y0 + y1) / 2 + 0.03, text, ha="center", fontsize=7, color=INK2)
ax = fig.add_subplot(gs[1, 0]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off"); ax.set_title("BEFORE: one black box learns everything", fontsize=9)
box(ax, 0.03, 0.55, 0.22, 0.3, "stage on the day\nalong\nacross", "#eeeeea"); box(ax, 0.36, 0.55, 0.28, 0.3, "Gaussian process\n3 inputs, 1,347 readings\n(rows of every kind, wrong heights)", "#f6e6df"); box(ax, 0.75, 0.55, 0.22, 0.3, "velocity at 1 ft\n(one per reading)", "#eeeeea")
arrow(ax, 0.25, 0.7, 0.36, 0.7); arrow(ax, 0.64, 0.7, 0.75, 0.7)
ax.text(0.5, 0.36, "The kernel had to discover, from eleven days never revisited:\n• how velocity changes with height (the profile)\n• how it changes with stage (the law)\n• where the bed is (the depth)\n• the pattern across and along the channel\nIt memorised days, and beyond them it fell back to the mean.", ha="center", va="top", fontsize=8, color=INK2)
ax = fig.add_subplot(gs[1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off"); ax.set_title("AFTER: measured and tested things do their jobs; one field is learned", fontsize=9)
box(ax, 0.02, 0.72, 0.20, 0.2, "stage record\n(the only live input)", "#eeeeea", 7); box(ax, 0.02, 0.44, 0.20, 0.2, "measured bed\n(RTK cross-sections)", "#eeeeea", 7); box(ax, 0.02, 0.16, 0.20, 0.2, "GP: conveyance A\n(along, across)\n260 points, known noise", "#dfe9f7", 7)
box(ax, 0.34, 0.58, 0.22, 0.2, "local depth\nh = stage − bed", "#eeeeea", 7); box(ax, 0.34, 0.24, 0.22, 0.2, "amplitude\na = A · h^β\n(β from the ADCP)", "#eeeeea", 7)
box(ax, 0.68, 0.58, 0.30, 0.2, "log law: v(z) = a ln(z/z0)\nany height, depth average", "#eeeeea", 7); box(ax, 0.68, 0.24, 0.30, 0.2, "habitat: share below a\nthreshold, area vs stage", "#eeeeea", 7)
arrow(ax, 0.22, 0.82, 0.34, 0.7); arrow(ax, 0.22, 0.54, 0.34, 0.66); arrow(ax, 0.45, 0.58, 0.45, 0.44); arrow(ax, 0.22, 0.26, 0.34, 0.32); arrow(ax, 0.56, 0.34, 0.68, 0.64); arrow(ax, 0.83, 0.58, 0.83, 0.44)
ax.text(0.5, 0.06, "The GP learns one static field with a physical meaning and error bars. Stage enters through the law, not the kernel.", ha="center", fontsize=8, color=INK2)
fig.suptitle("Before and after: what the model's target is", fontsize=10, y=0.96); fig.savefig(F / "00t_before_after_target.png", bbox_inches="tight"); plt.close(fig); print("saved 00t")
