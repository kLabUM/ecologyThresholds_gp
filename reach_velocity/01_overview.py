"""Step 0b. Look at the data before modelling anything: where the transects are, when the surveys happened relative to
stage and discharge, what the ADCP recorded, and what a vertical profile looks like. Reads data_clean/, writes figures/00*.png."""
from pathlib import Path
import json, numpy as np, pandas as pd, geopandas as gpd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT = Path(__file__).resolve().parents[1]; HERE = ROOT / Path(__file__).resolve().parent.name; D = HERE / "data_clean"; F = HERE / "figures"; F.mkdir(exist_ok=True)
FT = 0.3048; ADCP_E, ADCP_N = 13244857.0, 319995.5
INK, INK2, AXIS, GRID = "#0b0b0b", "#52514e", "#c3c2b7", "#e1e0d9"; C25, C26, CADCP, CQ = "#2a78d6", "#1a9e74", "#eb6834", "#898781"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "font.size": 8, "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False, "savefig.dpi": 200, "axes.facecolor": "#fcfcfb"})
tr = pd.read_csv(D / "readings_transect.csv", parse_dates=["time_utc"]); nests = pd.read_csv(D / "nests.csv")
st = pd.read_csv(D / "stage_10min.csv", parse_dates=["time_utc"]).set_index("time_utc"); ad = pd.read_csv(D / "adcp.csv", parse_dates=["time_utc"]); q = pd.read_csv(D / "usgs_daily.csv", parse_dates=["date"])
loc = tr.groupby("location_id").agg(year=("year", "first"), date=("date", "first"), transect=("Transect ID", "first"), x=("x_m", "first"), y=("y_m", "first"), n=("velocity_ms", "size"), stage=("stage_m", "median")).reset_index()
days = tr.groupby("date").agg(year=("year", "first"), stage=("stage_m", "median"), n_loc=("location_id", "nunique"), n=("velocity_ms", "size")).reset_index()
days["q_cfs"] = days.date.map(q.set_index(q.date.dt.strftime("%Y-%m-%d")).discharge_cfs)
ycol = lambda y: C25 if y == 2025 else C26

# ---------------------------------------------------------------- 00a map
def to_local(gdf):
    if gdf.crs is None or gdf.crs.is_geographic: gdf = gdf.set_crs("EPSG:4269", allow_override=True).to_crs("ESRI:102690") if gdf.crs is None else gdf.to_crs("ESRI:102690")
    return gdf.translate(xoff=-ADCP_E, yoff=-ADCP_N).scale(FT, FT, origin=(0, 0))
RAW = HERE / "data_raw"
cl = to_local(gpd.read_file(D / "banklines_survey" / "centerline_union_Middle.geojson")); bk = to_local(gpd.read_file(D / "banks_c" / "C_Middle_reference.geojson"))   # 44_ and 19_ (bank set C at 0.90 m)
fig, ax = plt.subplots(figsize=(8.5, 6.5))
pad = 120; xlo, xhi, ylo, yhi = loc.x.min() - pad, loc.x.max() + pad, loc.y.min() - pad, loc.y.max() + pad
bk.plot(ax=ax, color="#b9b8b0", lw=0.8); cl.plot(ax=ax, color="#6f6e68", lw=1.0, linestyle="--")
for tid, g in loc.groupby("transect"):
    ax.scatter(g.x, g.y, s=12, color=ycol(g.year.iloc[0]), edgecolor="white", lw=0.4, zorder=4)
    ax.annotate(tid, (g.x.mean(), g.y.mean()), xytext=(6, 6), textcoords="offset points", fontsize=6, color=INK2)
ax.scatter([0], [0], marker="*", s=160, color=CADCP, edgecolor="white", lw=0.6, zorder=6); ax.annotate("ADCP (SonTek)", (0, 0), xytext=(8, -12), textcoords="offset points", fontsize=7, color=CADCP)
ax.axhline((321000 - ADCP_N) * FT, color=INK2, lw=0.6, ls=":"); ax.text(xlo + 5, (321000 - ADCP_N) * FT + 4, "Northing 321000 ft: cut used in the notebooks (transects above it were dropped)", fontsize=6.5, color=INK2)
ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, yhi); ax.set_aspect("equal"); ax.set_xlabel("east of the ADCP (m)"); ax.set_ylabel("north of the ADCP (m)"); ax.grid(True); ax.set_axisbelow(True)
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=C25, markersize=5, label=f"2025 transect locations ({(loc.year==2025).sum()}, 9 days)"), Line2D([], [], marker="o", ls="", color=C26, markersize=5, label=f"2026 transect locations ({(loc.year==2026).sum()}, 2 days)"),
                   Line2D([], [], color="#6f6e68", ls=(0, (5, 3)), label="centreline (44_banklines_survey)"), Line2D([], [], color="#b9b8b0", label="bank lines (set C at 0.90 m stage)")], loc="lower right", fontsize=7)
ax.set_title("Middle reach: where the hand-held readings were taken", fontsize=9); fig.savefig(F / "00a_map.png", bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- 00b stage and discharge with survey days
fig, axes = plt.subplots(2, 2, figsize=(10, 5), sharex="col", gridspec_kw={"hspace": 0.12, "wspace": 0.12})
for j, yr in enumerate([2025, 2026]):
    s = st[st.index.year == yr]; qq = q[q.date.dt.year == yr]
    axes[0, j].plot(s.index, s.stage_m, color=INK, lw=0.6); axes[1, j].plot(qq.date, qq.discharge_cfs, color=CQ, lw=1.0)
    if yr == 2025: axes[0, j].axvspan(ad.time_utc.min(), ad.time_utc.max(), color=CADCP, alpha=0.08, lw=0); axes[0, j].text(ad.time_utc.min(), 1.75, "ADCP recording", fontsize=6.5, color=CADCP)
    for _, d in days[days.year == yr].iterrows():
        for ax in axes[:, j]: ax.axvline(pd.Timestamp(d.date, tz="UTC"), color=ycol(yr), lw=0.9, alpha=0.9)
    axes[0, j].set_title(f"{yr}: {int((days.year==yr).sum())} survey days (vertical lines)", fontsize=9); axes[0, j].set_ylim(0.6, 1.9); axes[1, j].set_ylim(0, 2100)
    for ax in axes[:, j]: ax.grid(True); ax.set_axisbelow(True)
    axes[1, j].tick_params(axis="x", rotation=0)
axes[0, 0].set_ylabel("Stage at the Middle sensor (m)"); axes[1, 0].set_ylabel("USGS discharge at Ann Arbor (cfs)")
fig.suptitle("When the surveys happened, against stage and discharge", fontsize=9, y=0.95); fig.savefig(F / "00b_stage_discharge.png", bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- 00c ADCP
fig = plt.figure(figsize=(10, 4)); gs = fig.add_gridspec(2, 2, width_ratios=[1.6, 1], hspace=0.1, wspace=0.2)
a0 = fig.add_subplot(gs[0, 0]); a1 = fig.add_subplot(gs[1, 0], sharex=a0); a2 = fig.add_subplot(gs[:, 1])
s = st[(st.index >= ad.time_utc.min()) & (st.index <= ad.time_utc.max())]; a0.plot(s.index, s.stage_m, color=INK, lw=0.6); a0.set_ylabel("stage (m)"); a0.tick_params(labelbottom=False)
a1.plot(ad.time_utc, ad.velocity_ms, color=CADCP, lw=0.4, alpha=0.9); a1.set_ylabel("ADCP velocity (m/s)")
for _, d in days[days.year == 2025].iterrows():
    for ax in (a0, a1): ax.axvline(pd.Timestamp(d.date, tz="UTC"), color=C25, lw=0.8, alpha=0.8)
hr = ad.set_index("time_utc").resample("1h").mean(); mon = hr.index.month
sc = a2.scatter(hr.stage_m, hr.velocity_ms, c=mon, cmap="Blues", vmin=4, vmax=10, s=6, lw=0); cb = fig.colorbar(sc, ax=a2, pad=0.02); cb.set_label("month (2025)")
a2.set_xlabel("stage (m)"); a2.set_ylabel("ADCP velocity, hourly mean (m/s)"); a2.set_title("The one place where velocity is known at every stage", fontsize=9)
for ax in (a0, a1, a2): ax.grid(True); ax.set_axisbelow(True)
a0.set_title("ADCP record at one point, with the 2025 survey days (blue lines)", fontsize=9); fig.savefig(F / "00c_adcp.png", bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- 00d vertical profiles by survey day
dl = days.sort_values("date").reset_index(drop=True); n = len(dl); ncol = 4; nrow = int(np.ceil(n / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(10, 2.4 * nrow), sharex=True, sharey=True, gridspec_kw={"hspace": 0.45, "wspace": 0.08}); axes = axes.ravel()
for i, d in dl.iterrows():
    ax = axes[i]; g = tr[tr.date == d.date]
    for _, p in g.groupby("location_id"):
        p = p.sort_values("height_m"); ax.plot(p.velocity_ms, p.height_m, color=ycol(d.year), lw=0.5, alpha=0.45); ax.scatter(p.velocity_ms, p.height_m, s=4, color=ycol(d.year), alpha=0.7, lw=0)
    ax.axhline(FT, color=INK, lw=0.6, ls=(0, (4, 3))); ax.axvline(0.2, color=INK2, lw=0.5, ls=":")
    ax.set_title(f"{d.date}\nstage {d.stage:.2f} m, {d.q_cfs:.0f} cfs, {int(d.n_loc)} locations", fontsize=6.5, pad=3); ax.grid(True); ax.set_axisbelow(True)
for ax in axes[n:]: ax.axis("off")
for ax in axes[::ncol]: ax.set_ylabel("height above bed (m)")
for ax in axes[max(0, n - ncol):n]: ax.set_xlabel("velocity (m/s)")
axes[0].set_xlim(-0.02, 1.0); axes[0].set_ylim(0, 1.0); axes[0].text(0.97, FT + 0.02, "1 ft", fontsize=6.5, ha="right", color=INK)
fig.suptitle("Vertical profiles: every reading, one line per location, per survey day (dashed: 1 ft above the bed; dotted: 0.2 m/s)", fontsize=9, y=0.995)
fig.savefig(F / "00d_profiles.png", bbox_inches="tight"); plt.close(fig)
print(days.to_string(index=False)); print("saved 00a-00d")
