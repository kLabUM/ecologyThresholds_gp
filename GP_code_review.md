# Code & Statistics Review — Huron River GP Models

**Reviewed:** `dataSorting.ipynb`, `gpModel_{Upper,Middle,Lower}.ipynb`, `modelVis_{Upper,Middle,Lower}.ipynb`, `main_csvInput.ipynb`, plus all files in `data/`, `sorted_data/`, `mapping/`, `model outputs/`.
**Date:** 2026-09-15
**Method:** every finding below was reproduced against your actual CSVs in a clean environment (scikit-learn 1.9.0, pandas 2.3.3, geopandas 1.1.4). Reported R² values were re-derived, not taken from the notebooks.

---

## Headline

Two things need to be settled before this becomes a paper:

1. **The reported cross-validation R² is not a generalization estimate.** Random `KFold` shuffles points that come from the same transect on the same survey day into both train and test. When folds are grouped properly, **Upper goes from R² = 0.36 to −0.22 and Lower from 0.49 to −0.45** — worse than predicting the mean. Middle survives best at +0.13.

2. **`upstream_depth` — the variable every management threshold depends on — takes only 12–16 distinct values per section** (one per field survey day). The effective sample size for learning the stage→velocity response is the number of survey days, not the number of points.

Everything else is fixable. These two are structural.

---

## Severity 1 — Fix before any further modeling

### 1.1 Upper and Lower depth series are built from the Middle sensor's file

`dataSorting.ipynb` cells 2b and 2c read only the 2026 file into `depth_data_mm_26`, then concatenate it with `depth_data_mm_25` — which is still bound to `25_depthSensorData_Middle.csv` from cell 2a:

```python
# cell 2b (Upper) — and identically in 2c (Lower)
depth_data_mm_26 = pd.read_csv('data/26_depthSensorData_Upper.csv', skiprows=3)
depth_data_mm = pd.concat([depth_data_mm_25, depth_data_mm_26], ignore_index=True)
#                          ^^^^^^^^^^^^^^^^ still the MIDDLE 2025 sensor
```

There is no `25_depthSensorData_Upper.csv` or `_Lower.csv` in `data/`, so the 2025 half of both series is Middle's gauge (ARB012). Consequences:

- **Upper** has zero 2025 observations, so no observation is matched to contaminated data — but the combined series becomes bimodal (~6450 mm vs ~1058 mm), which silently disables the 2σ outlier trim. My run confirms it: 1489 points in, 1489 out. Nothing was trimmed.
- **Lower** has 50 observations in 2025 that are currently matched to Middle's gauge, ~19 river miles upstream. With the correct file they have no depth data and should drop out (687 → 637 matched).

### 1.2 One site's stage datum is hard-coded for all three sites

Every depth CSV carries an `offset` column — the deployment datum in mm. The code ignores it and uses `7313.2` for all three sections. That value belongs to **ARB012, the Middle gauge**:

| File | site_id | `offset` values present |
|---|---|---|
| `25_depthSensorData_Middle.csv` | ARB012 | 6858.0 (May 30 – Jun 9 2025), 7324.34 (11 readings), **7313.2** (Jun 9 2025 →) |
| `26_depthSensorData_Middle.csv` | ARB012 | 7313.2 |
| `26_depthSensorData_Upper.csv` | ARB017 | **1587.5** |
| `26_depthSensorData_Lower.csv` | ARB112 | **5943.6** (Mar 1 – Jun 8 2026), **6705.6** (Jun 8 →) |

**Do not simply switch to applying the `offset` column per row.** I tested that first and it is wrong. At both sites where the offset changes, the `device_id` is *unchanged* across the change — so these are datum *revisions* entered retrospectively, not physical remountings. Applying them per row injects a step that isn't there:

| | step at the offset change | |
|---|---|---|
| Middle, 2025-06-09 (device 074A both sides) | per-row offset: **+1.27 ft** | any single constant: **−0.22 ft** |
| Lower, 2026-06-08 (device 067A both sides) | per-row offset: **+2.62 ft** | any single constant: **+0.12 ft** |

A correct datum removes steps; these create them. The constant-offset approach was right — the bug is that **the same constant was used for three different gauges.**

(Lower's actual sensor swap is later: device 067A → 147A between Jul 2 and Jul 6 2026, with the offset already at 6705.6 on both sides. That transition shows a +0.48 ft step regardless of treatment, which may well be real.)

### 1.3 Net effect on `upstream_depth`

Correcting 1.1 and 1.2 together — each section reads its own files, with one constant datum per site (ARB012 = 7313.2, ARB017 = 1587.5, ARB112 = 6705.6), plus error-code filtering:

| Section | As written (m) | Corrected (m) | sd | Obs matched | Verdict |
|---|---|---|---|---|---|
| Middle | 0.710 – 1.339 | **0.710 – 1.339** | unchanged | 1691 → 1691 | **No change.** 7313.2 is ARB012's correct datum |
| Upper | 6.024 – 6.426 | 0.299 – 0.701 | identical (0.0818) | 1168 → 1168 | **Constant 5.73 m shift**, correlation 1.0000 |
| Lower | 0.788 – 1.440 | 0.282 – 0.832 | 0.157 → 0.143 | 687 → **637** | **Real change** — 50 contaminated 2025 obs drop out |

So the three sections land very differently, and it's worth being precise about it:

- **Middle is fine.** Its depth pipeline is correct as written. Middle's published model does not need refitting on this account.
- **Upper's model is numerically unaffected.** `StandardScaler` centres the feature, so a constant shift cancels — the correlation between the as-written and corrected values is 1.0000 and the standard deviations match to four decimals. **But every reported number is wrong.** The paper would state a 6.0–6.4 m stage range for a reach whose true range is 0.30–0.70 m, and any threshold expressed as an absolute stage would be unusable by a water manager.
- **Lower does need refitting.** Not because the depth values shift much (correlation 0.9995), but because 50 observations currently carry a depth value sourced from a gauge ~19 river miles upstream, and should have none.

Please confirm ARB112's datum against your deployment notes — I picked 6705.6 as the later, presumably better-surveyed revision, but that choice moves Lower's absolute stage by 2.5 ft and only your field record can settle it. It does not affect the fitted model (constant shift), only the reported values.

### 1.4 One off-channel transect destroys Lower's `downstream_dist`

`gpModel_Lower.ipynb` never applies the 150 ft centerline filter that `dataSorting` assumes. 17 points survive — all of transect **LT0101, 2025-08-06**, roughly 200 ft off the main stem. They project to `dist_along_centerline` ≈ 16,350 ft against a reference of 99,494 ft, giving `downstream_dist` ≈ **−83,150 ft (−25,347 m)**.

`StandardScaler` then sees:

| | mean | sd |
|---|---|---|
| With the 17 outliers | −1,712 ft | **16,729 ft** |
| Without | 1,715 ft | 906 ft |

An **18.5× inflation** of the scale. The 404 genuine in-channel points, which span 0–3,499 ft, get compressed into **0.209 standardized units** — and the GP's RBF length-scale lower bound is 0.1. The entire study reach is about two length scales wide. The downstream dimension is effectively annihilated.

The `gpModel` notebooks even contain a check for this that fires and is ignored:

```python
negative_count = (valid_downstream < -0.01).sum()
if negative_count > 0:
    print(f"WARNING: {negative_count} points have negative distances (should be 0)")
```

**Fix:** apply `complete_data = complete_data[complete_data['distance_to_centerline'].abs() <= 150]` in all three `gpModel` notebooks, and make the warning raise instead of print.

---

## Severity 2 — The validation strategy

### 2.1 Random KFold does not measure what the paper needs it to measure

Your data is clustered three ways: by survey **date** (all points on a day share one `upstream_depth`), by **transect**, and by **location** (repeat visits — Upper has 741 points at 363 unique coordinates; Middle 1360 at 386; Lower 421 at 220). Random `KFold` puts siblings from every cluster in both train and test.

Re-run with identical kernel and settings, changing only the fold assignment:

| Section | Random KFold *(published)* | Group by **date** | Group by **location** | Group by **transect** |
|---|---|---|---|---|
| Upper | **+0.363** *(pub. 0.366)* | **−0.225** | +0.269 | — |
| Middle | **+0.622** *(pub. 0.664)* | **+0.131** | −2.62 | −2.83 |
| Lower | **+0.494** | +0.132 | +0.378 | +0.001 |
| Lower *(1.4 fixed)* | +0.512 | −0.451 | — | +0.049 |

My random-KFold numbers reproduce yours to within 0.04, so this is a like-for-like comparison.

**Group by date is the number that matters for your application.** Water managers will ask "at this stage, what will the velocity field look like?" — a stage the model has not seen. Under that test, Upper and Lower are at or below zero, and Middle is at +0.13.

Middle's catastrophic location/transect folds (one fold at **−15.0**) are the ADCP block; see 2.3.

### 2.2 `upstream_depth` has 12–16 unique values

| Section | n points | unique `upstream_depth` | unique locations | unique transects |
|---|---|---|---|---|
| Upper | 741 | **16** | 363 | 28 |
| Middle | 1360 | 110 *(475 of them ADCP-only)* | 386 | 21 |
| Lower | 404 | **12** | 220 | 20 |

The first GP input is functionally a **categorical survey-day index**. An ablation makes this concrete — `upstream_depth` *alone*, with no spatial information at all:

| Section | inputs | Random KFold | Group by date |
|---|---|---|---|
| Upper | all 3 | +0.363 | −0.217 |
| Upper | spatial only (2) | +0.223 | −0.106 |
| Upper | **`upstream_depth` only (1)** | **+0.221** | −0.120 |
| Lower | all 3 | +0.512 | −0.451 |
| Lower | spatial only (2) | +0.350 | −0.119 |
| Lower | **`upstream_depth` only (1)** | **+0.147** | −0.317 |

For Upper, a single 16-level variable scores **as well as both spatial coordinates combined** under random KFold, and collapses under grouped CV. That is the signature of day-mean memorization, not a learned hydraulic relationship.

Note also that adding inputs makes grouped-CV performance *worse* in both sections — the model is using the extra dimensions to overfit.

### 2.3 The ADCP block: 34 % of Middle at one point, measuring a different quantity

`dataSorting` section 3 appends 475 SonTek rows to Middle — **475 of 1378 points (34 %)**, all at the single hard-coded coordinate (319995.5, 13244857), all with identical `downstream_dist` and `distance_to_centerline`. Under GroupKFold-by-location they land in one fold, and that fold scores **R² = −15.0**.

Separately, these are not the same physical measurement as the rest of the target variable:

- `_measurement` is `sontek_Vb` at site **ARB111** — a side-looking index/beam velocity.
- The wading-rod data is a **point velocity**, and the `Depth <= 1` filter restricts it to the top foot of the water column.
- The code assigns the ADCP rows `Depth = 1` purely so they survive that filter.

Pooling an index velocity and near-surface point velocities into one `Velocity` response, with no measurement-type covariate, is hard to defend in review. At minimum add `Type` as a categorical input or a per-type noise term; the cleaner option is to exclude the ADCP block from the response and use it as validation of the stage–velocity relationship instead.

---

## Severity 3 — Model specification

### 3.1 The fitted kernels are degenerate — hyperparameters pinned at their bounds

```
Middle: 0.667² * RBF(length_scale=[0.71, 0.10, 0.10]) + RationalQuadratic(alpha=2.09, length_scale=0.848) + WhiteKernel(noise_level=0.195)
Upper:  0.782² * RBF(length_scale=[0.10, 0.10, 5.82]) + RationalQuadratic(alpha=0.10, length_scale=0.10) + WhiteKernel(noise_level=0.485)
```

Bold marks the lower bound of 0.1. **Middle has two length scales at the bound; Upper has four hyperparameters at the bound** (both RQ parameters plus two RBF dimensions). sklearn emits `ConvergenceWarning` for each, and the notebook output shows them being scrolled past.

A length scale at the lower bound means that kernel component has become a near-delta function — it contributes nothing but memorization. What's left doing the actual work is the `RationalQuadratic` term, which in sklearn is **isotropic**: a single length scale applied to all three standardized dimensions. So after the ARD RBF collapses, the model is isotropic over (stage, downstream, cross-stream) — treating a 718 m downstream axis and a 74 m cross-stream axis as interchangeable.

Upper is worse: its RQ is *also* at the bound in both parameters, and `WhiteKernel(noise_level=0.485)` means **48.5 % of the standardized variance is being called noise.**

Recommendations:

- Widen `length_scale_bounds` to something like `(1e-2, 1e3)` and refit. If they still hit the floor, the kernel is wrong for the data, not under-bounded.
- Replace `RBF + RationalQuadratic` with a single **ARD Matérn (ν = 3/2 or 5/2)**. It has one length scale per dimension, is far better behaved for rough spatial fields than RBF, and removes the RBF/RQ identifiability problem — the two terms are currently competing to explain the same structure, which is exactly why one collapses.
- Drop `n_restarts_optimizer=15` to 5 once the kernel is identifiable; the 30-minute fits are mostly the optimizer re-finding the same degenerate optimum.
- `alpha=1e-6` alongside a `WhiteKernel` is redundant but harmless.

### 3.2 Predictive uncertainty includes observation noise

`gp.predict(..., return_std=True)` returns a std that includes the `WhiteKernel` term. For a *threshold on the expected velocity field* you want the latent-function uncertainty, not the noise of a single wading-rod reading. With `noise_level` at 0.195 (Middle) and 0.485 (Upper), this is the dominant term — it's why `y_std_3v` is nearly flat across space (Upper: 0.078 to 0.091, a 17 % range over the whole reach). Those uncertainty maps in `modelVis` are showing you the noise floor, not spatial information content.

### 3.3 Negative predictions are clipped

`np.clip(velocity_grid, 0, None)` appears in every `modelVis` prediction path. Clipping means the GP is predicting physically impossible values and the fix is applied after the fact — which silently breaks the calibration of the uncertainty bands that the threshold analysis depends on. Model `sqrt(Velocity)` or `log(Velocity + ε)` and back-transform, or use a warped GP.

---

## Severity 4 — The threshold analysis specifically

### 4.1 "% area below threshold" is computed over a rectangle, not the wetted channel

`modelVis` cell 5 builds a 50 × 50 grid spanning `downstream_dist` × `distance_to_centerline` from the data's min/max, then reports the fraction of grid cells below 0.2 m/s. But the channel is not a rectangle in those coordinates — the half-width varies station to station, and the grid predicts everywhere out to the maximum half-width observed anywhere in the reach.

| Section | grid half-width | median sampled half-width | downstream bins (of 20) with **zero** data | grid cells >20 m from any training point |
|---|---|---|---|---|
| Upper | 23.4 m | 11.9 m | **15** | **86.9 %** |
| Middle | 37.2 m | 21.7 m | 0 | 7.3 % |
| Lower | 61.5 m | 16.0 m | **18** | **97.4 %** |

For Upper and Lower the reported area fraction is dominated by cells on dry land, in unsampled reaches, or both. **Middle is the only section where this statistic is defensible**, and even there roughly 40 % of cells are more than 10 m from any observation.

**Fix:** mask the grid to the wetted channel. The simplest defensible version is a per-station half-width envelope interpolated from the transect data, plus a convex hull or alpha-shape in (downstream, cross-stream) space, plus a hard mask where predictive variance exceeds some fraction of the prior variance. Then area-weight the cells rather than counting them.

### 4.2 The uncertainty band on the area curve is not a ±1σ band

```python
vel_upper = np.clip(vel_pred + vel_std, 0, None)
area_lower = ((vel_upper <= VELOCITY_THRESHOLD).sum() / len(vel_upper)) * 100
```

This shifts **every grid cell in the same direction simultaneously**, which assumes the GP's errors are perfectly correlated across the entire reach. It produces an extreme envelope, then labels it `±1σ uncertainty` on the plot in cell 6.

The correct quantity is the posterior distribution *of the area fraction*. Draw joint samples and compute the statistic per sample:

```python
# needs the joint posterior, not marginal stds
samples = gp_model_3v.sample_y(X_grid_scaled, n_samples=500, random_state=0)   # (n_grid, 500)
vel = scaler_y_3v.inverse_transform(samples.reshape(-1,1)).reshape(samples.shape)
area = (vel <= VELOCITY_THRESHOLD).mean(axis=0) * 100
mean, lo, hi = area.mean(), *np.percentile(area, [2.5, 97.5])
```

This will give a **much narrower and much more credible** band than the current one, because spatial errors partially cancel. It is also the number a management agency can actually act on.

### 4.3 Dropping negative velocities biases the threshold metric

`dataSorting` section 4 applies `Velocity >= 0`, removing **112 observations** across the three sections (28 Upper / 52 Middle / 32 Lower). Magnitudes run from −0.002 to −0.400 ft/s.

The small ones are instrument noise around zero and should be retained as ~0. The larger ones are real eddies and reverse flow — which is precisely the slow-water refugia habitat the ecological thresholds are meant to identify. Systematically deleting the low end of the velocity distribution biases the "% area below 0.2 m/s" statistic **downward**, in the non-conservative direction for habitat protection.

Suggested replacement: `Velocity = Velocity.clip(lower=0)` for |v| below the instrument noise floor, and retain larger reversals as |v| with a `reverse_flow` flag.

### 4.4 The threshold is a hard-coded magic number

`VELOCITY_THRESHOLD = 0.2` appears identically in all three `modelVis` notebooks with no citation and no sensitivity analysis. Given the project's framing — thresholds informing management decisions that prioritize ecological health — this needs a source (species, life stage, literature) and a sensitivity sweep across a plausible range. Since `Depth <= 1` makes this a **near-surface** velocity model, check that the threshold in your source is defined on the same quantity; most fish-habitat criteria are stated as depth-averaged or nose velocity.

---

## Severity 5 — Correctness issues worth fixing

### 5.1 Middle's `downstream_dist` no longer starts at zero, and the notebook notices

`gpModel_Middle.ipynb` applies `complete_data[complete_data["Northing"] < 321000]` — dropping exactly 15 rows, all transect `0729-BT01`. But `reference_distance` was chosen in `dataSorting` as the minimum along-line distance **including** those rows. After the filter, Middle's `downstream_dist` runs **1169.3 – 1887.1 m**, not 0 – 718 m.

`modelVis_Middle.ipynb` cell 3 contains a self-check for this:

```python
print(f"Expected: if downstream_dist was calculated correctly, min should be ~0")
print(f"Actual min downstream_dist value: {gp_data_3v['downstream_dist'].min():.6f}")
```

It prints 1169.3 and nothing acts on it. The back-projection still works (the reference is added back consistently), but any figure axis or table reporting "distance downstream" is offset by 1,169 m. Move the `Northing` filter into `dataSorting` before the reference is chosen.

### 5.2 `build_main_stem` is chaining far beyond the study reach

| Section | main stem length | data bbox diagonal | ratio | % of stem >500 ft from any data point |
|---|---|---|---|---|
| Upper | 44,906 ft | 16,014 ft | 2.8 | 84 % |
| Middle | 32,698 ft | 5,429 ft | 6.0 | 87 % |
| Lower | **153,501 ft** | 52,533 ft | 2.9 | **96 %** |

Lower's "main stem" is 29 miles for a reach whose data spans 3,500 ft. The function keeps extending through junctions as long as any branch connects, regardless of whether it goes anywhere near the study area. This is what put LT0101 83,000 ft upstream of the reference (1.4). Add a stopping criterion: stop extending once the chain is more than some buffer beyond the data's along-line extent.

### 5.3 Sensor error codes are not filtered

`_value` is exactly **500** in every depth file (the MaxBotix minimum-range code): 364 readings in `26_Upper` (1.4 %), 235 in `26_Lower` (1.2 %). Lower also has 100 readings above 9000 mm (0.5 %). These are treated as real stage. The 300-minute median resample hides most of them, but they should be dropped explicitly before aggregation.

### 5.4 The section-assignment gaps leave rows on their original CSV labels

```python
complete_data.loc[complete_data['Northing'] > 340000, 'Section'] = 'Upper'
complete_data.loc[(complete_data['Northing'] <= 323000) & (complete_data['Northing'] >= 280000), 'Section'] = 'Middle'
complete_data.loc[complete_data['Northing'] < 255000, 'Section'] = 'Lower'
```

Northings in 255,000–280,000 and 323,000–340,000 are never reassigned; they keep whatever the CSV said, and `fillna('Middle')` then sweeps up blanks. In the current data this affects 15 rows (all Northing 323,623–323,651, labelled Middle — the same `0729-BT01` transect that 5.1 later deletes) plus **9 rows with NaN Northing** that get silently assigned to Middle. It works today by coincidence. Make the boundaries exhaustive and assert that no row is left unassigned.

### 5.5 Pickle paths don't match what's on disk

`gpModel_*.ipynb` writes to `model outputs/gp_model_{section}.pkl`, and `modelVis_*` reads from there. On disk:

- `gp_model_Middle.pkl` — 15 MB, at the **project root**, not in `model outputs/`
- `model outputs/gp_model_Upper.pkl` — present
- `model outputs/gp_model_Lower.pkl` — **missing**; `gpModel_Lower.ipynb` has no saved output and its notebook shows no completed run

So `modelVis_Middle.ipynb` will raise `FileNotFoundError` as committed, and Lower has never been fit end to end. Also: `.gitignore` excludes two stale CSV names but not `*.pkl`, so a 15 MB binary is being committed on every model run.

### 5.6 Smaller items

- **`reference_distance` is labelled metres, stored as feet.** `modelVis` prints `"Reference distance along centerline: {reference_distance:.2f} meters"`. The arithmetic below it is correct (it converts to feet explicitly); only the label is wrong. Easy to trip over later.
- **`modelVis` cells 1 and 2 are near-duplicates** — cell 1 is an earlier copy without the `Type` colouring. Delete cell 1.
- **23 exact duplicate rows** in the concatenated ecohydro data.
- **`Actual depth` (2025) vs `Depth 2` (2026)** — the two years have different columns and `pd.concat` silently unions them into all-NaN-for-one-year columns. Harmless now since neither is used, but it will bite if either is ever added as an input.
- **`Velocity 2` is essentially empty** (1 / 12 / 0 non-null values) — not a lost-data concern, but worth deleting from the schema so it isn't mistaken for a usable second reading.
- **`Depth <= 1` removes 30–43 % of observations.** That's a legitimate modelling choice, but it makes this a near-surface velocity model and the paper should say so explicitly in the methods.

---

## Suggested order of work

1. **Rebuild the depth pipeline** — read each section's own files, use one constant datum *per gauge* (not per row, and not ARB012's value for everything), drop the 500/>9000 error codes. Confirm ARB112's datum against your field notes. Middle needs no refit from this; Lower does; Upper needs its reported values corrected.
2. **Apply the 150 ft centerline filter** in all three `gpModel` notebooks and make the negative-distance check raise.
3. **Move Middle's `Northing < 321000` filter into `dataSorting`** so the reference distance is computed on the final row set.
4. **Switch to `GroupKFold` on survey date** as the headline metric. Report random KFold alongside it if you like, but label it as in-sample.
5. **Refit with an ARD Matérn** and wider bounds; check no hyperparameter sits at a bound.
6. **Mask the prediction grid to the wetted channel** and recompute the area statistic with joint posterior samples.
7. Then — and only then — revisit whether 12–16 stage values can support the threshold claims.

### On point 7

This is the honest constraint. With 16 survey days, the stage→velocity relationship is estimated from 16 points, and Upper's observed stage range is only about 0.40 m wide. A GP will happily draw a smooth curve through the slider range and it will look convincing, but grouped CV says it does not predict an unseen day.

Two ways forward that don't require another field season:

- **Use the continuous depth-sensor record as the forcing.** You have 5-minute stage at all three sites across two seasons. If velocity can be related to stage through a physically-structured model (a rating relationship per location, with the GP modelling the spatial residual), you get the full stage range rather than 16 slices of it.
- **Report thresholds as interpolation within the observed stage envelope only**, with the envelope drawn on every figure, and state plainly that extrapolation beyond it is unsupported. This is defensible and publishable; the current framing is not, because the slider invites extrapolation the data can't back.

The ADCP record is a genuine asset here — a continuous velocity series at a known point. It would do more good as an independent check on the stage–velocity relationship than as 34 % of Middle's training set.
