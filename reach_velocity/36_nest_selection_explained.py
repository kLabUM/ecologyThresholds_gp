"""One figure to explain the nest selection model: what 'used against available' means, why velocity only shows once depth is in the
model, what the coefficients mean in plain units, and where the model sits in the chain (before: a 0.2 m/s cut; now: a learned suitability).
Output: figures/12_nest_selection_explained.png"""
import json, numpy as np, pandas as pd
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; from matplotlib.patches import FancyBboxPatch
from gpcore import *; from uncertainty import GridLevels, PRIMARY
rng = np.random.default_rng(0); Z0 = json.load(open(D / "field_qc.json"))["z0_m"]; H0, GAM = PRIMARY; nb = json.load(open(D / "nest_bayes.json"))
grid = pd.read_csv(D / "grid.csv"); nests = pd.read_csv(D / "step7_nests.csv").dropna(subset=["v_pred_ms"]).reset_index(drop=True); cell = grid.iloc[cKDTree(grid[["x_m", "y_m"]].values).query(nests[["x_m", "y_m"]].values)[1]].reset_index(drop=True); cix = cKDTree(grid[["x_m", "y_m"]].values).query(nests[["x_m", "y_m"]].values)[1]; GLV = GridLevels()   # the field and bank distance at each nest's own stage (bank set C, entry 30)
def vbar(logA, bed, st):
    h = np.clip(st - bed, 0, None); a = np.exp(logA) * np.clip(st - H0, 0.02, None) ** GAM; return np.where(h > np.e * Z0, a * (np.log(np.clip(h, np.e * Z0, None) / Z0) - 1), np.nan), h
U, A = [], []
for i, r in nests.iterrows():
    st = r.stage_use; gA, gD = GLV.at(st); vg, hg = vbar(gA, grid.bed_m.values, st); ok = ~np.isnan(vg); c = cell.iloc[i]; vc, hc = vbar(gA[cix[i]], c.bed_m, st); U.append([float(vc), gD[cix[i]], float(hc)])
    for j in rng.choice(np.where(ok)[0], 20, replace=False): A.append([vg[j], gD[j], hg[j]])
U, A = np.array(U), np.array(A); bv, bb, bh = nb["velocity_coef_per_m_s"]["posterior"][0], nb["bank_distance_coef_per_m"]["posterior"][0], nb["depth_coef_per_m"]["posterior"][0]
INK, INK2, AXIS, GRID, BLUE, GREEN, DARK, MUTED, ORANGE = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#1a9e74", "#0d366b", "#898781", "#eb6834"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8.5, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
fig = plt.figure(figsize=(17, 11)); gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 0.9], hspace=0.55, wspace=0.32)
# row 1: before and now
ax = fig.add_subplot(gs[0, 0]); vv = np.linspace(0, 0.8, 200); ax.step([0, 0.2, 0.2, 0.8], [1, 1, 0, 0], color=ORANGE, lw=2, where="post", label="before: a nest can be here if v < 0.2 m/s, else not")
ax.plot(vv, np.exp(bv * vv), color=DARK, lw=2, label=f"now: relative use exp({bv:.1f} v), learned from the nests"); ax.axhline(0.5, color=INK, lw=0.5, ls=":"); ax.text(nb["v50_m_s_posterior"]["median"] + 0.01, 0.52, f"half use at {nb['v50_m_s_posterior']['median']:.2f} m/s", fontsize=8, color=DARK)
ax.set_xlabel("depth-averaged velocity in the cell (m/s)"); ax.set_ylabel("suitability (use relative to still water)"); ax.set_title("Before and now: the velocity rule", fontsize=9.5); ax.legend(fontsize=7.5, loc="upper right"); ax.grid(True); ax.set_axisbelow(True); ax.set_ylim(0, 1.1)
ax = fig.add_subplot(gs[0, 1]); ax.axis("off"); ax.set_title("What 'used against available' means", fontsize=9.5)
ax.text(0.0, 0.98, "Each nest is one 'used' cell: the 2 m cell it sits in, with the\nvelocity, depth and bank distance the chain gives that cell\non the day the nest was found.\n\nFor each nest, 20 'available' cells are drawn at random from\nthe wet channel at the same stage. They are what the fish\ncould have chosen.\n\nA logistic model asks: what makes a cell more likely to be\na used one than an available one?\n\n   log-odds(used) = b0 + bv·velocity + bd·bank distance + bh·depth\n\nThe coefficients are the selection. A negative bv means nests\navoid faster cells, all else equal. 'All else equal' is the point:\neach coefficient is the effect of one variable when the other\ntwo are held fixed.", transform=ax.transAxes, va="top", fontsize=8.5, linespacing=1.35)
ax = fig.add_subplot(gs[0, 2]); ax.axis("off"); ax.set_title("The coefficients in plain units (posterior)", fontsize=9.5)
rows = [("velocity", f"{bv:.1f} per m/s", f"a cell 0.1 m/s faster is {np.exp(bv*0.1):.2f} times as likely to hold a nest"), ("bank distance", f"{bb:.2f} per m", f"a cell 1 m further from the bank: {np.exp(bb):.2f} times; 5 m further: {np.exp(5*bb):.2f} times"), ("depth", f"{bh:.1f} per m", f"a cell 0.1 m deeper is {np.exp(bh*0.1):.2f} times as likely (at the same velocity and bank distance)")]
for k, (n, c, t) in enumerate(rows): ax.text(0.0, 0.92 - k * 0.3, f"{n}: {c}", transform=ax.transAxes, fontsize=9, fontweight="bold", va="top"); ax.text(0.0, 0.83 - k * 0.3, t, transform=ax.transAxes, fontsize=8.5, va="top", wrap=True)
ax.text(0.0, 0.02, "Read: nests sit in deep, slow water near the bank.", transform=ax.transAxes, fontsize=9, color=DARK)
# row 2: used vs available, and the confound picture
for k, (col, lab, bins) in enumerate([(0, "velocity (m/s)", np.linspace(0, 0.8, 25)), (1, "distance to the nearer bank (m)", np.linspace(0, 27, 28)), (2, "depth (m)", np.linspace(0, 1.4, 29))]):
    ax = fig.add_subplot(gs[1, k]); ax.hist(A[:, col], bins=bins, density=True, color=MUTED, alpha=0.5, label="available (20 cells per nest, same stage)"); ax.hist(U[:, col], bins=bins, density=True, color=GREEN, alpha=0.75, label="used (59 nests)")
    ax.axvline(np.median(A[:, col]), color=MUTED, lw=1.2); ax.axvline(np.median(U[:, col]), color=GREEN, lw=1.2); ax.set_xlabel(lab); ax.set_ylabel("density"); ax.set_title(f"{lab.split(' (')[0]}: used median {np.median(U[:, col]):.2f}, available {np.median(A[:, col]):.2f}", fontsize=9); ax.grid(True); ax.set_axisbelow(True)
    if k == 0: ax.legend(fontsize=7.5)
# row 3: why velocity hides until depth is in the model, and where the model sits
ax = fig.add_subplot(gs[2, 0:2]); near = A[:, 1] < 8; nu = U[:, 1] < 8
ax.scatter(A[near, 2], A[near, 0], s=9, color=MUTED, alpha=0.5, lw=0, label="available cells within 8 m of a bank"); ax.scatter(U[nu, 2], U[nu, 0], s=26, color=GREEN, edgecolor=INK, lw=0.4, zorder=3, label="nests within 8 m of a bank")
hh = np.linspace(0.05, 1.3, 50); a_typ = np.median(np.exp(grid.logA[grid.dist_bank_m < 8])) * np.clip(np.median(nests.stage_use) - H0, 0.02, None) ** GAM; ax.plot(hh, a_typ * (np.log(np.clip(hh, np.e * Z0, None) / Z0) - 1), color=INK2, lw=1, ls="--", label="what the log law gives a typical near-bank cell: deeper is faster")
ax.annotate("at the same depth, nests are in the slower cells\n(the velocity effect, once depth is held fixed)", xy=(0.75, 0.13), xytext=(0.85, 0.52), fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.7, color=INK2)); ax.annotate("at the same bank distance, nests are in the deeper cells", xy=(0.55, 0.05), xytext=(0.05, 0.6), fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.7, color=INK2))
ax.set_xlabel("depth of the cell (m)"); ax.set_ylabel("depth-averaged velocity of the cell (m/s)"); ax.set_title("Why velocity only shows once depth is in the model: near the bank, deeper cells are faster, and nests take the deep-but-slow ones", fontsize=9.5); ax.legend(fontsize=7.5, loc="upper left"); ax.grid(True); ax.set_axisbelow(True); ax.set_xlim(0, 1.35); ax.set_ylim(0, 0.75)
ax = fig.add_subplot(gs[2, 2]); ax.axis("off"); ax.set_title("Where it sits in the chain", fontsize=9.5)
def box(x, y, w, h, t, fc): ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02", fc=fc, ec=INK, lw=0.6, transform=ax.transAxes)); ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=8, transform=ax.transAxes)
box(0.05, 0.78, 0.9, 0.16, "the velocity chain: v(x, t), depth h(x, t), bank distance d(x)\nfor every surveyed point, every two hours, with draws", "#e8f0fb")
box(0.05, 0.46, 0.42, 0.2, "before\n1[v < 0.2 m/s]\n(assumed)", "#fbe9e0"); box(0.53, 0.46, 0.42, 0.2, "now\nS = exp(bv v + bh h + bd d)\n(learned from the nests,\nprior from the literature)", "#e3f3ee")
box(0.05, 0.08, 0.9, 0.22, "the ecological number with its band\nbefore: the share of wetted width below 0.2 m/s\nnow: the reach's mean suitability E, or the share above any level of S,\nwith the curve's posterior and the chain's draws combined", "#f1f0ea")
for x0 in (0.26, 0.74): ax.annotate("", xy=(x0, 0.67), xytext=(x0, 0.78), xycoords="axes fraction", arrowprops=dict(arrowstyle="->", lw=0.8, color=INK)); ax.annotate("", xy=(x0, 0.31), xytext=(x0, 0.46), xycoords="axes fraction", arrowprops=dict(arrowstyle="->", lw=0.8, color=INK))
fig.suptitle("The nest selection model, explained: from an assumed 0.2 m/s cut to a suitability learned from where the fish nested", fontsize=11, y=0.97); fig.savefig(F / "12_nest_selection_explained.png", bbox_inches="tight"); plt.close(fig); print("saved 12")
