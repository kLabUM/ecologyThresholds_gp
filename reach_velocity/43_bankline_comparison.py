"""Step 10. Bank lines: three sets compared under the rule written before any result (logbook entry 27).

    python 43_bankline_comparison.py build    # the sets and the visual check (figures/18_bankline_sets.png)
    python 43_bankline_comparison.py score    # sections held out under A, B and C; the rule applied

A: the old lines (made by 44_banklines_survey.py, the original notebook code). B: surveyed edges (1 m out) with the 3DEP edges between them. C: stage-dependent edges at each survey
day's sensor water level, with 20 draws of their uncertainty (lidar error, survey-to-lidar gap, the sensor tie).
Inputs: data_raw/3dep_bank_transects_{r}.csv, bank_lines_{r}.geojson, centerline_union_{r}.geojson; data_clean point tables.
Outputs (data_clean/bankcmp/): sets_{r}.npz, banks_B_{r}.geojson, banks_C_{r}_{date}.geojson, anchors_{r}.csv, summary_build.json;
score: scores_points.csv, scores_sections.csv, paired.csv, decision.json, figures/18_bankline_scores.png"""
import json, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd, geopandas as gpd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import norm
import banklines_3dep as B
warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent; D = HERE / "data_clean"; F = HERE / "figures"; OUT = D / "bankcmp"; OUT.mkdir(exist_ok=True)
REACHES = ["Upper", "Middle", "Lower"]; FT = B.FT; SEED = 27


def build(reach):
    t = B.load_transects(reach); st = B.stations(reach, t); e = B.flat_edges(t); an = B.survey_anchors(reach, st); old = B.set_A(reach)
    an.to_csv(OUT / f"anchors_{reach}.csv", index=False)
    bneg, bpos = B.set_B(reach, st, e, an); banksB = B.offsets_to_lines(reach, st, bneg, bpos, old); banksB.to_file(OUT / f"banks_B_{reach}.geojson", driver="GeoJSON")
    C = B.StageBanks(reach, t, st, e, an); k = B.tie(reach); p = B.survey_points(reach)
    days = p.groupby("date").stage_m.median()
    best = {}
    for day, s in days.items():
        lev = B.water_level(reach, s, st.along_m.to_numpy()); ed = C.edges(lev)
        n_, p_ = B.push_out(st, ed["neg"][0].copy(), ed["pos"][0].copy(), reach, p[p.date == day])   # the day's own points inside
        best[day] = (n_, p_, ed["neg"][1], ed["pos"][1])
        B.offsets_to_lines(reach, st, n_, p_, old).to_file(OUT / f"banks_C_{reach}_{day}.geojson", driver="GeoJSON")
    rng = np.random.default_rng(SEED); draws = {}
    for d in range(B.N_DRAWS):
        zn, zp = B.correlated_normals(len(st), rng), B.correlated_normals(len(st), rng)
        for day, s in days.items():
            lev = B.water_level(reach, s, st.along_m.to_numpy()) + rng.normal(0, k["sd"])
            n_, p_ = C.draw(lev, zn, zp); draws[(d, day)] = B.push_out(st, n_, p_, reach, p[p.date == day])
    np.savez_compressed(OUT / f"sets_{reach}.npz", st=st.to_numpy(float), st_cols=np.array(st.columns), Bneg=bneg, Bpos=bpos,
                        days=np.array(list(days.index)), stage=days.to_numpy(),
                        Cbest=np.array([best[dy][:2] for dy in days.index]), Csd=np.array([best[dy][2:] for dy in days.index]),
                        Cdraws=np.array([[draws[(d, dy)] for dy in days.index] for d in range(B.N_DRAWS)]), valid=C.valid, W=C.W, hs=C.hs)
    return {"t": t, "st": st, "e": e, "an": an, "old": old, "banksB": banksB, "C": C, "days": days, "best": best, "draws": draws, "p": p}


def outside(pts, banks):
    poly = B.polygon(banks).buffer(0.5 / FT)
    from shapely.geometry import Point
    return int(sum(not poly.contains(Point(x, y)) for x, y in zip(pts.Easting, pts.Northing)))


def figure_sets(R):
    n = len(R); fig = plt.figure(figsize=(6.6 * n, 11.5)); gs = fig.add_gridspec(2, n, height_ratios=[3.3, 1], hspace=0.28, wspace=0.16)
    ORANGE, GREEN, INK, GREY = "#ea580c", "#15803d", "#1f2937", "#9ca3af"; BLUES = ["#93c5fd", "#3b82f6", "#1e3a8a"]
    summ = {}
    for j, (reach, r) in enumerate(R.items()):
        st, p, old, days = r["st"], r["p"], r["old"], r["days"]; x0, y0 = p.Easting.median(), p.Northing.median()
        X = lambda v: (np.asarray(v) - x0) * FT; Y = lambda v: (np.asarray(v) - y0) * FT
        def xy(off): return st.x.to_numpy() + off / FT * st.nx.to_numpy(), st.y.to_numpy() + off / FT * st.ny.to_numpy()
        ax = fig.add_subplot(gs[0, j]); ax.set_aspect("equal")
        order = days.sort_values().index; pick = [order[0], order[len(order) // 2], order[-1]]
        mid = pick[1]; D_ = np.array([r["draws"][(d, mid)] for d in range(B.N_DRAWS)])
        for k_ in (0, 1):
            lo, hi = np.quantile(D_[:, k_, :], 0.05, axis=0), np.quantile(D_[:, k_, :], 0.95, axis=0)
            xl, yl = xy(lo); xh, yh = xy(hi); ax.fill(np.r_[X(xl), X(xh)[::-1]], np.r_[Y(yl), Y(yh)[::-1]], color="#bfdbfe", lw=0, zorder=0)
        for g in old.geometry: gx, gy = g.xy; ax.plot(X(gx), Y(gy), color=ORANGE, lw=1.2, ls=(0, (4, 2)), zorder=2)
        for g in r["banksB"].geometry: gx, gy = g.xy; ax.plot(X(gx), Y(gy), color=GREEN, lw=1.4, zorder=3)
        for c, dy in zip(BLUES, pick):
            for k_ in (0, 1):
                gx, gy = xy(r["best"][dy][k_]); ax.plot(X(gx), Y(gy), color=c, lw=1.3, zorder=4)
        ax.scatter(X(p.Easting), Y(p.Northing), s=5, color=INK, lw=0, zorder=5)
        an = r["an"]
        for side in ("neg", "pos"):
            a = an.dropna(subset=[f"{side}_m"]); ii = a.station.to_numpy(); gx = st.x.to_numpy()[ii] + a[f"{side}_m"].to_numpy() / FT * st.nx.to_numpy()[ii]; gy = st.y.to_numpy()[ii] + a[f"{side}_m"].to_numpy() / FT * st.ny.to_numpy()[ii]
            ax.scatter(X(gx), Y(gy), s=20, facecolor="white", edgecolor=INK, lw=0.8, zorder=6)
        inval = ~r["C"].valid
        if inval.any(): ax.scatter(X(st.x[inval]), Y(st.y[inval]), s=10, marker="x", color="#b91c1c", lw=0.8, zorder=6)
        ax.set_xlim(X(p.Easting).min() - 40, X(p.Easting).max() + 40); ax.set_ylim(Y(p.Northing).min() - 40, Y(p.Northing).max() + 40)
        ax.set_xlabel("east (m)"); ax.set_ylabel("north (m)" if j == 0 else "")
        banksC_mid = gpd.read_file(OUT / f"banks_C_{reach}_{mid}.geojson")
        lv = {dy: float(np.median(B.water_level(reach, days[dy], st.along_m))) for dy in pick}
        o = {"A": outside(p, old), "B": outside(p, r["banksB"]), "C_own_day": int(sum(outside(p[p.date == dy], gpd.read_file(OUT / f"banks_C_{reach}_{dy}.geojson")) for dy in days.index))}
        summ[reach] = {"survey_days": int(len(days)), "stations": int(len(st)), "stations_surveyed_edges_used": int(inval.sum()), "lidar_water_m_median": float(np.nanmedian(r["C"].W)),
                       "survey_day_water_m": {dy: round(float(np.median(B.water_level(reach, days[dy], st.along_m))), 2) for dy in days.index}, "points_outside": o, "n_points": int(len(p))}
        ax.set_title(f"{reach}: C drawn at {pick[0]} ({lv[pick[0]]:.2f} m), {pick[1]} ({lv[pick[1]]:.2f} m), {pick[2]} ({lv[pick[2]]:.2f} m)\n"
                     f"lidar water {summ[reach]['lidar_water_m_median']:.2f} m; surveyed edges stand in at {inval.sum()} of {len(st)} stations\n"
                     f"points outside (0.5 m): A {o['A']}, B {o['B']}, C on each point's own day {o['C_own_day']} (of {len(p)})", fontsize=9, loc="left", color=INK)
        a2 = fig.add_subplot(gs[1, j]); s_ = st.station_m
        wA = B.bank_features(pd.DataFrame({"Easting": st.x, "Northing": st.y}), old).width_m
        a2.plot(s_, wA, color=ORANGE, lw=1.2, ls=(0, (4, 2)), label="A old")
        a2.plot(s_, r["banksB"].pipe(lambda b: B.bank_features(pd.DataFrame({"Easting": st.x, "Northing": st.y}), b).width_m), color=GREEN, lw=1.4, label="B hybrid")
        for c, dy in zip(BLUES, pick): a2.plot(s_, r["best"][dy][1] - r["best"][dy][0], color=c, lw=1.2, label=f"C at {lv[dy]:.2f} m")
        w_lid = r["e"].edge_pos_m_s - r["e"].edge_neg_m_s; a2.plot(s_, w_lid, color=GREY, lw=0.8, ls=":", label="3DEP alone")
        a2.set_ylim(bottom=0); a2.grid(alpha=0.3); a2.set_xlabel("distance along the 3DEP stations (m)"); a2.set_ylabel("width (m)" if j == 0 else ""); a2.legend(fontsize=7, frameon=False, ncol=2)
    fig.legend(handles=[Line2D([], [], color=ORANGE, lw=1.2, ls=(0, (4, 2)), label="A: old lines"), Line2D([], [], color=GREEN, lw=1.4, label="B: survey + 3DEP hybrid"),
                        Line2D([], [], color=BLUES[0], lw=1.3, label="C at the lowest survey-day level"), Line2D([], [], color=BLUES[1], lw=1.3, label="C at a middle level"),
                        Line2D([], [], color=BLUES[2], lw=1.3, label="C at the highest level"), plt.Rectangle((0, 0), 1, 1, color="#bfdbfe", label="C 5 to 95% (middle level, 20 draws)"),
                        Line2D([], [], marker="o", ls="", mfc="white", mec=INK, ms=6, label="surveyed edges (1 m out)"), Line2D([], [], marker="x", ls="", color="#b91c1c", ms=5, label="lidar not used (canopy rule)")],
               loc="upper center", bbox_to_anchor=(0.5, 0.965), ncol=4, frameon=False, fontsize=9)
    fig.suptitle("Bank-line sets for the comparison (logbook entry 27): water levels are the sensor's, carried to each station", y=0.995, fontsize=12)
    fig.savefig(F / "18_bankline_sets.png", dpi=160, bbox_inches="tight"); plt.close(fig)
    json.dump(summ, open(OUT / "summary_build.json", "w"), indent=1); return summ


# ---------------------------------------------------------------- scoring (after sign-off) ---
def model_points(reach):
    """The GP's point table with each reach's primary stage law applied (as in 20_/38_)."""
    from gpcore import load_points, with_target
    if reach == "Middle":
        import uncertainty as U
        df = load_points().merge(pd.read_csv(D / "locations.csv")[["location_id", "Easting", "Northing"]], on="location_id")
        h0, gam = U.LAW_RIDGE[0]
    else:
        df = pd.read_csv(D / f"reach_{reach}_chain_points.csv"); h0, gam = json.load(open(D / f"reach_{reach}_chain_summary.json"))["law"]["primary"]
    return with_target(df, ("stage", gam, h0)).reset_index(drop=True)


def lines_from(reach, z, st, neg, pos, old):
    return B.offsets_to_lines(reach, st, neg, pos, old)


def section_out(d):
    """Leave one section out; returns mu, sd (GP sd combined with the point's own sd, as in 38_)."""
    from gpcore import Model
    MU, SD = np.full(len(d), np.nan), np.full(len(d), np.nan)
    for t in sorted(d.transect.unique()):
        te = (d.transect == t).to_numpy(); m = Model(n_restarts=1).fit(d[~te]); mu, sd = m.predict(d[te])
        MU[te], SD[te] = mu, np.sqrt(sd ** 2 + d.sd_y.to_numpy()[te] ** 2)
    return MU, SD


def _A(m, s):
    return m * (2 * norm.cdf(m / s) - 1) + 2 * s * norm.pdf(m / s)


def mixture_scores(y, MU, SD):
    """CRPS and PIT of an equal-weight Gaussian mixture per point; MU, SD are (draws, points)."""
    k = MU.shape[0]; w = 1.0 / k
    e1 = w * _A(y[None, :] - MU, SD).sum(0)
    e2 = sum(w * w * _A(MU[i] - MU[j], np.sqrt(SD[i] ** 2 + SD[j] ** 2)) for i in range(k) for j in range(k))
    pit = (w * norm.cdf((y[None, :] - MU) / SD)).sum(0)
    return e1 - 0.5 * e2, pit


def score_reach(reach, log):
    d = model_points(reach); z = np.load(OUT / f"sets_{reach}.npz", allow_pickle=True)
    st = pd.DataFrame(z["st"], columns=list(z["st_cols"])); old = B.set_A(reach); days = list(z["days"])
    banksB = gpd.read_file(OUT / f"banks_B_{reach}.geojson")
    res = {}
    for arm, banks in (("A", old), ("B", banksB)):
        dd = d.copy(); dd[["dist_bank_m", "width_m", "rel_across"]] = B.bank_features(dd, banks)[["dist_bank_m", "width_m", "rel_across"]]
        mu, sd = section_out(dd); zz = (dd.y.to_numpy() - mu) / sd
        crps = sd * (zz * (2 * norm.cdf(zz) - 1) + 2 * norm.pdf(zz) - 1 / np.sqrt(np.pi)); res[arm] = (crps, norm.cdf(zz), mu, sd); log(f"  {reach} {arm}: CRPS {crps.mean():.3f}")
    MU, SD = [], []
    for k in range(B.N_DRAWS):
        dd = d.copy()
        for dy in sorted(dd.date.unique()):
            m_ = (dd.date == dy).to_numpy(); j = days.index(dy); neg, pos = z["Cdraws"][k, j, 0], z["Cdraws"][k, j, 1]
            dd.loc[m_, ["dist_bank_m", "width_m", "rel_across"]] = B.bank_features(dd[m_], B.offsets_to_lines(reach, st, neg, pos, old))[["dist_bank_m", "width_m", "rel_across"]].to_numpy()
        mu, sd = section_out(dd); MU.append(mu); SD.append(sd); log(f"  {reach} C draw {k + 1}/{B.N_DRAWS}")
    MU, SD = np.array(MU), np.array(SD); crps, pit = mixture_scores(d.y.to_numpy(), MU, SD); res["C"] = (crps, pit, MU.mean(0), np.sqrt(SD.var(0) + 0) )
    rows = []
    for arm, (crps, pit, *_ ) in res.items():
        rows.append(pd.DataFrame({"reach": reach, "set": arm, "location_id": d.location_id, "transect": d.transect, "crps": crps, "pit": pit,
                                  "in90": (pit >= 0.05) & (pit <= 0.95)}))
    return pd.concat(rows)


def decide(P):
    S = P.groupby(["reach", "set", "transect"]).agg(crps=("crps", "mean"), in90=("in90", "mean"), n=("crps", "size")).reset_index()
    S.to_csv(OUT / "scores_sections.csv", index=False)
    W = S.pivot_table(index=["reach", "transect"], columns="set", values="crps")
    rows = []
    for a, b in (("B", "A"), ("C", "A"), ("C", "B")):
        for scope, w in [("all three", W)] + [(r, W.loc[r]) for r in REACHES]:
            dif = (w[a] - w[b]).dropna(); rows.append({"pair": f"{a} - {b}", "scope": scope, "sections": int(len(dif)), "mean_diff": float(dif.mean()),
                                                     "se": float(dif.std(ddof=1) / np.sqrt(len(dif))), "better_in_sections": int((dif < 0).sum())})
    pr = pd.DataFrame(rows); pr["tie"] = pr.mean_diff.abs() < pr.se; pr.to_csv(OUT / "paired.csv", index=False)
    tot = S.groupby("set").agg(crps=("crps", "mean"), in90=("in90", "mean")); per = S.groupby(["reach", "set"]).agg(crps=("crps", "mean"), in90=("in90", "mean"))
    return S, pr, tot, per


def figure_scores(S, pr, tot, per):
    fig, axs = plt.subplots(1, 3, figsize=(17, 5.2)); COL = {"A": "#ea580c", "B": "#15803d", "C": "#2563eb"}
    ax = axs[0]; x = np.arange(len(REACHES)); w = 0.26
    for k, s in enumerate("ABC"):
        v = [per.loc[(r, s), "crps"] for r in REACHES]; ax.bar(x + (k - 1) * w, v, w, color=COL[s], label=s)
        for xi, vi in zip(x + (k - 1) * w, v): ax.text(xi, vi + 0.005, f"{vi:.3f}", ha="center", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(REACHES); ax.set_ylabel("CRPS of ln A, sections held out (lower is better)"); ax.legend(frameon=False, title="bank lines"); ax.set_title("Mean of section scores", loc="left")
    ax = axs[1]
    for k, s in enumerate("ABC"):
        v = [per.loc[(r, s), "in90"] for r in REACHES]; ax.bar(x + (k - 1) * w, v, w, color=COL[s])
    ax.axhline(0.9, color="#1f2937", lw=1, ls="--"); ax.set_ylim(0, 1); ax.set_xticks(x); ax.set_xticklabels(REACHES); ax.set_ylabel("share of points inside the 90% band"); ax.set_title("Coverage (target 0.90)", loc="left")
    ax = axs[2]; q = pr[pr.scope == "all three"].reset_index(drop=True)
    ax.errorbar(q.mean_diff, np.arange(len(q)), xerr=q.se, fmt="o", color="#1f2937", capsize=4)
    for i, r in q.iterrows(): ax.text(r.mean_diff, i + 0.18, f"{r.mean_diff:+.3f} +/- {r.se:.3f}" + (" (tie)" if r.tie else ""), ha="center", fontsize=9)
    ax.axvline(0, color="#9ca3af", lw=1); ax.set_yticks(range(len(q))); ax.set_yticklabels(q.pair); ax.set_ylim(-0.6, len(q) - 0.4)
    ax.set_xlabel("paired difference in section CRPS (negative: first set better)"); ax.set_title(f"All three reaches, {int(q.sections.iloc[0])} sections, mean +/- 1 se", loc="left")
    fig.suptitle("Bank-line sets scored under the rule in logbook entry 27 (each section held out in turn)", fontsize=12)
    fig.tight_layout(); fig.savefig(F / "18_bankline_scores.png", dpi=160); plt.close(fig)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "build"
    if mode == "build":
        R = {r: build(r) for r in REACHES}
        print(json.dumps(figure_sets(R), indent=1))
    elif mode == "score":
        import time
        T0 = time.time(); log = lambda m: print(f"[{time.time() - T0:6.0f}s] {m}", flush=True)
        parts = []
        for r in (sys.argv[2:] or REACHES):
            f = OUT / f"scores_points_{r}.csv"
            if not f.exists(): score_reach(r, log).to_csv(f, index=False)
            parts.append(pd.read_csv(f))
        P = pd.concat(parts); P.to_csv(OUT / "scores_points.csv", index=False)
        if set(P.reach) == set(REACHES):
            S, pr, tot, per = decide(P); figure_scores(S, pr, tot, per)
            print(per.round(3).to_string()); print(tot.round(3).to_string()); print(pr.round(4).to_string())
