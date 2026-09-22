"""Feature ablation for the conveyance GP (Middle, stage law), every section held out, with a paired test per section against the full set.
Outputs: step2_feature_ablation.csv (also the kernel family and restart check from step2_kernel_ablation.csv is left as is)."""
import numpy as np, pandas as pd, warnings, json; warnings.filterwarnings("ignore")
from scipy.stats import wilcoxon
from gpcore import *
df = load_points(); d = with_target(df, ("stage", 1.0, 0.45))
sets = {"full (6)": FULL, "no section depth": [f for f in FULL if f != "section_mean_depth_m"], "geometry only (rel_across, dist_bank, rel_depth)": ["rel_across", "dist_bank_m", "rel_depth"], "bank distance only": ["dist_bank_m"], "position only (along, across)": ["along_m", "across_m"], "position + rel_depth": ["along_m", "across_m", "rel_depth"]}
def crps_terms(y, mu, sd):
    from scipy.stats import norm; z = (y - mu) / sd; return sd * (z * (2 * norm.cdf(z) - 1) + 2 * norm.pdf(z) - 1 / np.sqrt(np.pi))
per = {}; rm = float(np.average(d.y, weights=1 / d.sd_y ** 2)); rsd = float(d.y.std())
for name, feats in sets.items():
    parts = []
    for t in sorted(d.transect.unique()):
        te = d[d.transect == t]; m = Model(feats=feats, n_restarts=1).fit(d[d.transect != t]); mu, sd = m.predict(te); parts.append(pd.Series(crps_terms(te.y.values, mu, np.sqrt(sd ** 2 + te.sd_y.values ** 2)), index=te.index))
    per[name] = pd.concat(parts)
per["reach mean"] = pd.Series(crps_terms(d.y.values, rm, np.sqrt(rsd ** 2 + d.sd_y.values ** 2)), index=d.index)
rows = []
for name, c in per.items():
    by_t = c.groupby(d.transect).mean(); full_t = per["full (6)"].groupby(d.transect).mean(); diff = by_t - full_t
    rows.append({"features": name, "CRPS": float(c.mean()), "CRPS_minus_full_per_section_mean": float(diff.mean()), "sections_better_than_full": int((diff < 0).sum()), "wilcoxon_p_vs_full": float(wilcoxon(diff).pvalue) if name != "full (6)" else np.nan})
R = pd.DataFrame(rows).sort_values("CRPS"); R.to_csv(D / "step2_feature_ablation.csv", index=False); print(R.round(4).to_string(index=False))
