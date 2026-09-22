"""A theoretical picture of the conveyance: where A sits in the velocity profile, in the cross-section, and in the
formulas. Synthetic numbers, no data. Figure 00u_conveyance_schematic.png"""
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
F = Path(__file__).resolve().parent / "figures"
INK, INK2, AXIS, GRID, BLUE, ORANGE, GREEN, DARK, MUTED = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9", "#2a78d6", "#eb6834", "#1a9e74", "#0d366b", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
z0, A = 0.011, 0.15; rng = np.random.default_rng(2)
fig = plt.figure(figsize=(14, 9.5)); gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], hspace=0.4, wspace=0.32)
# ---- (1) the profile at one point, two stages
ax = fig.add_subplot(gs[0, 0])
for h, col, lab in [(0.4, GREEN, "low stage, depth 0.4 m"), (0.8, DARK, "high stage, depth 0.8 m")]:
    a = A * np.sqrt(h); z = np.linspace(z0, h, 80); ax.plot(a * np.log(z / z0), z, color=col, lw=1.6, label=f"{lab}: a = A√h = {a:.3f}"); ax.axhline(h, color=col, lw=0.7, ls="--")
    zr = np.array([0.05, 0.1, 0.2, 0.3, 0.5, 0.7]); zr = zr[zr < h]; ax.scatter(a * np.log(zr / z0) + rng.normal(0, 0.012, len(zr)), zr, s=18, color=col, edgecolor="white", lw=0.5, zorder=4)
ax.plot([0], [z0], "o", color=INK, ms=4); ax.annotate("z0: where the curve reaches zero,\nthe roughness height (11 mm)", (0, z0), xytext=(0.3, 0.04), fontsize=7, color=INK2, arrowprops=dict(arrowstyle="->", color=INK2, lw=0.6))
ax.text(0.02, 0.62, "dots: what the crew measured\n(velocity at several heights)", fontsize=7, color=INK2)
ax.set_xlim(0, 0.65); ax.set_ylim(0, 0.9); ax.set_xlabel("velocity (m/s)"); ax.set_ylabel("height above the bed, z (m)"); ax.set_title("1. One point, two stages: the log profile v(z) = a ln(z/z0)\nsteep near the bed, bending over upward", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7, loc="upper left")
# ---- (2) same, log axis: straight lines
ax = fig.add_subplot(gs[0, 1])
for h, col in [(0.4, GREEN), (0.8, DARK)]:
    a = A * np.sqrt(h); z = np.linspace(z0, h, 80); ax.plot(a * np.log(z / z0), np.log(z / z0), color=col, lw=1.6, label=f"slope = a = {a:.3f} m/s")
ax.annotate("the amplitude a is the slope of this line:\nhow fast velocity grows per unit of ln(z)", (0.25, 2.3), xytext=(0.02, 3.9), fontsize=7, color=INK2, arrowprops=dict(arrowstyle="->", color=INK2, lw=0.6))
ax.text(0.25, 0.25, "both lines start at the same z0", fontsize=7, color=INK2)
ax.set_xlim(0, 0.65); ax.set_ylim(0, 4.6); ax.set_xlabel("velocity (m/s)"); ax.set_ylabel("ln(z / z0)"); ax.set_title("2. The same two profiles with height on a log scale:\nstraight lines, and a is their slope", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7, loc="lower right")
# ---- (3) a against depth: where A comes from
ax = fig.add_subplot(gs[0, 2]); hh = np.linspace(0.05, 1.2, 100)
for Ax, col, lab in [(0.06, MUTED, "a bank point, A = 0.06"), (0.15, BLUE, "mid-channel, A = 0.15"), (0.22, ORANGE, "a riffle, A = 0.22")]: ax.plot(hh, Ax * np.sqrt(hh), color=col, lw=1.5, label=lab)
ax.plot([0.4, 0.8], [A * np.sqrt(0.4), A * np.sqrt(0.8)], "o", color=DARK, ms=6, zorder=5); ax.annotate("the two profiles of panel 1", (0.8, A * np.sqrt(0.8)), xytext=(0.55, 0.19), fontsize=7, color=INK2, arrowprops=dict(arrowstyle="->", color=INK2, lw=0.6))
ax.set_xlabel("local depth h (m), = stage − bed"); ax.set_ylabel("amplitude a (m/s)"); ax.set_title("3. At one place, a grows with depth as A√h (the ADCP's law).\nA is the curve the place sits on: its steepness at unit depth", fontsize=9); ax.grid(True); ax.set_axisbelow(True); ax.legend(fontsize=7, loc="upper left"); ax.set_xlim(0, 1.2); ax.set_ylim(0, 0.26)
# ---- (4) the cross-section with A varying across
ax = fig.add_subplot(gs[1, 0:2]); x = np.linspace(-20, 20, 200); bed = 0.9 - 0.85 * np.exp(-((x + 3) / 11) ** 2) - 0.1 * np.exp(-((x - 12) / 4) ** 2); stage = 0.95
ax.fill_between(x, -0.05, bed, color="#d9d3c7", lw=0); ax.plot(x, bed, color=INK, lw=1.2); ax.axhline(stage, color=BLUE, lw=1.0, ls="--"); ax.text(-19.5, stage + 0.02, "water surface at this stage", fontsize=7, color=BLUE)
Ax_prof = lambda xx: 0.16 * np.exp(-((xx + 3) / 12) ** 2) + 0.03
for xp in [-15, -9, -3, 3, 9, 14]:
    b = np.interp(xp, x, bed); h = stage - b
    if h <= 0.05: continue
    Ap = Ax_prof(xp); a = Ap * np.sqrt(h); z = np.linspace(z0, h, 40); v = a * np.log(z / z0); scale = 6.0
    ax.plot(xp + v * scale, b + z, color=DARK, lw=1.2); ax.plot([xp, xp], [b, b + h], color=AXIS, lw=0.5); ax.text(xp - 0.3, b + h + 0.03, f"A = {Ap:.2f}", fontsize=7, ha="center", color=DARK)
ax.text(12, 0.15, "at each point across the channel: one A,\nand the whole profile at any stage follows\nfrom A, the depth and the law", fontsize=7.5, color=INK2)
ax.set_xlim(-20, 20); ax.set_ylim(-0.05, 1.15); ax.set_xlabel("across the channel (m)"); ax.set_ylabel("elevation (m)"); ax.set_title("4. Across a surveyed cross-section: A is a number at each point (small at the banks, large mid-channel);\nthe little curves are the velocity profiles it generates at this stage, drawn sideways", fontsize=9); ax.grid(True); ax.set_axisbelow(True)
# ---- (5) the formulas
ax = fig.add_subplot(gs[1, 2]); ax.axis("off")
lines = [("Measured", INK), ("  readings v at heights z (crew)", INK2), ("  bed elevation (RTK), stage (sensor); h = stage − bed", INK2), ("Fitted once for the reach", INK), ("  z0 = 11 mm, from all profiles", INK2), ("  β = 1/2, from the ADCP's record", INK2), ("Learned per place (the GP)", INK), ("  A, from the profile: a = slope of v vs ln z,  A = a / h^β", INK2), ("Physics behind it", INK), ("  u* = √(g h S),  a = u*/κ,  so  A = √(g S)/κ", INK2), ("  so A is the local friction slope; the reach's", INK2), ("  measured slope, 5 cm/100 m, gives A = 0.17", INK2), ("Prediction at any stage", INK), ("  a = A √h,   v(z) = a ln(z/z0),   v̄ = a (ln(h/z0) − 1)", INK2)]
for i, (t, c) in enumerate(lines): ax.text(0.0, 0.97 - i * 0.068, t, fontsize=8.2 if c == INK else 7.6, color=c, va="top", fontweight="bold" if c == INK else "normal", family="sans-serif")
ax.set_title("5. What is measured, what is fitted once, what is learned", fontsize=9)
fig.suptitle("The conveyance, in theory", fontsize=10, y=0.97); fig.savefig(F / "00u_conveyance_schematic.png", bbox_inches="tight"); print("saved 00u")
