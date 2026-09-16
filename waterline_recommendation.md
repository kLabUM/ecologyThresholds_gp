# Building an Explicit Channel Geometry for the Huron GP Models

**Question:** how to construct a stream waterline / edge geometry per site, given no edge GeoJSON, using transect edge points as a proxy.

**Short answer:** don't fit a curve to the transect endpoints — they are not on the waterline, and the resulting edge would be wrong by a variable 20–35%. You already hold something much better: **surveyed bed elevation plus a water-surface elevation for every transect**. The waterline is where those two intersect, and deriving it that way gets you a *stage-dependent* boundary, which is exactly what a model with stage as an input needs.

All numbers below were measured from `25_ecohydroData.csv` + `26_ecohydroData.csv`.

---

## 1. Why the endpoint-proxy plan won't work

The transects stop in open water, not at the bank:

| | median total water depth at first point | at last point | transects ending shallower than 0.5 ft at **both** ends |
|---|---|---|---|
| Upper | 1.05 ft | 1.40 ft | |
| Middle | 0.60 ft | 0.75 ft | |
| Lower | 0.58 ft | 2.10 ft | |
| **all** | | | **1 of 69** |

There is not a single zero in the `Total water depth` column across 2,946 transect readings; the minimum is 0.1 ft. Crews stopped where it got deep, unsafe, or simply where the last sensible reading was — which is a different distance from the bank on every transect and every visit.

I also tested the obvious repair, extrapolating each transect's near-bank bed slope down to zero depth:

| | median extra half-width implied | 25th–75th pct | max |
|---|---|---|---|
| Upper | 5.1 m | 2.5 – 25.7 m | 120 m |
| Middle | 14.2 m | 8.0 – 21.4 m | 53 m |
| Lower | 7.5 m | 3.3 – 15.4 m | 96 m |

As a fraction of surveyed width that's a median of 18–34%, but the interquartile range is enormous and the tail is absurd. The reason is in the geometry: these are wide, shallow channels (median surveyed width 19 m Upper / 41 m Middle / 33 m Lower, with under 1.5 m of bed relief), so near-bank slopes are only 0.03–0.10 m/m. A slope fitted from three points at 2.3–3.8 m spacing is too noisy to divide by. **Per-transect linear extrapolation is not usable.** A fitted cross-section *shape* is, which is what section 3 does.

---

## 2. What you actually have: bed topography and a water surface

`Elevation` is populated on 3,005 of 3,011 transect rows. Testing which combination is flat along a transect settles what it is:

| quantity | median within-transect sd (ft) — Upper / Middle / Lower |
|---|---|
| `Elevation` alone | 0.98 / 0.47 / 0.95 |
| **`Elevation + Total water depth`** | **0.34 / 0.17 / 0.19** |
| `Elevation − Total water depth` | 1.63 / 0.93 / 1.56 |

A water surface is flat across a channel; a bed is not. So:

- **`Elevation` is the surveyed streambed elevation.**
- **`Elevation + Total water depth` is the water-surface elevation (WSE)**, flat to ~0.2 ft within a transect, which is about RTK GNSS noise.

Corroborating: within-transect `corr(Elevation, Total water depth)` is −0.93 to −0.97, exactly what you get when both are measured from a common flat surface.

Plot three Middle cross-sections and it's unambiguous — they are textbook channel sections with a level water surface on top (`transect_bed_wse.png`).

**This is the whole ballgame.** You are not missing edge data. You have ~3,000 georeferenced bed points at 2.3–3.8 m cross-channel spacing, and a water-surface elevation per transect-day. The waterline at any stage is the contour where bed elevation equals water-surface elevation. Nothing needs to be fit to "edge points" at all.

---

## 3. Recommended construction

Work in the centerline frame you already have — `dist_along_centerline` (s) and `distance_to_centerline` (n) — so everything stays consistent with the GP inputs.

### Step 1 — bed surface `z_b(s, n)`

Interpolate the 3,005 bed points onto an (s, n) grid. Keep it simple and monotone-ish; this is a smooth, well-sampled surface within transects and a gap-filling problem between them.

- Across-channel (n): dense, 2.3–3.8 m spacing. A smoothing spline per transect is plenty.
- Along-channel (s): sparse. Transect spacing is ~34 m (Middle), ~36 m (Upper), ~53 m (Lower). **This is the hard limit on planform detail** — you cannot resolve a meander cutbank or a point bar shorter than that from this data. Say so in the methods.

A thin-plate spline or a 2-D GP with a long s length scale and a short n length scale is appropriate. Resist using the velocity GP's kernel here; bed elevation is a much smoother field.

### Step 2 — water-surface profile `WSE(s, h)`

You have one WSE per transect-day. Fit a longitudinal profile that shifts with gauge stage `h`:

```
WSE(s, h) ≈ a(s) + b · h
```

`a(s)` is the water-surface slope down the reach (fit as a low-order polynomial or monotone spline in s, constrained to be decreasing downstream). `b` converts gauge stage to local stage; it is near 1 for a reach that responds uniformly.

**Caveat, and it's a real one:** no transect in your dataset was visited on a second day. Every transect contributes exactly one (s, WSE, h) triple, so `b` is identified only by contrasts *between* transects surveyed on different days — it borrows strength across space to estimate a temporal effect. Fit it, report the uncertainty on `b`, and check the residuals against `a(s)`. If `b` comes out far from 1 or with a wide interval, fall back to `b = 1` (rigid translation of the water surface) and state that as an assumption.

### Step 3 — waterline and channel polygon

For a stage `h`, the wetted domain is `{(s, n) : z_b(s, n) ≤ WSE(s, h)}`. Walk each s-station outward from the thalweg to the first crossing on each side; that gives left and right offsets `n_L(s, h)` and `n_R(s, h)`. Smooth those two offset curves along s (a penalized spline, or a low-pass at the transect-spacing scale — do not smooth finer than ~35–55 m, you'd be inventing detail), then polygonize.

Export per site and per stage, e.g. `mapping/waterline_{Section}.geojson` with a `stage_m` attribute, so `dataSorting` and the visualization notebooks read one canonical file instead of each re-deriving a half-width.

### Step 4 — sanity-check against imagery

Digitize the waterline from one leaf-off aerial image and compare. This is a cross-check, not the primary source, because an image is a single unknown stage.

- **NAIP** (1 m or 60 cm, via [The National Map](https://apps.nationalmap.gov/lidar-explorer/)) — resolves a 20–40 m river easily.
- **USGS 3DEP lidar** (1 m DEM; [Michigan coverage](https://www.michigan.gov/dtmb/-/media/Project/Websites/dtmb/Services/GIS/3DEP.pdf), [OpenTopography](https://portal.opentopography.org/raster?opentopoID=OTNED.012021.4269.3)) — lidar does not penetrate water, so the flat no-return surface *is* the waterline at flight time. This is the best independent check, and it also gives you overbank topography your wading survey can't reach.
- **NHDPlus HR** ([USGS](https://www.usgs.gov/national-hydrography/nhdplus-high-resolution)) `NHDArea` polygons — worth a look, but rivers under roughly 30 m are often carried only as a flowline with no areal feature, and your median surveyed widths (19–41 m) straddle that. I could not confirm coverage for these reaches from the public service; check layer 8 of the [NHDPlus_HR MapServer](https://hydro.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/MapServer) against your bounding boxes before planning around it.

Also check [Washtenaw County open data](https://data-washtenaw.opendata.arcgis.com/) — counties often publish a hydrography layer digitized at a better scale than NHD.

---

## 4. The bigger payoff: this fixes the model, not just the mask

An explicit edge tidies up the map and the area denominator. But the same geometry gives you two changes to the **model inputs** that matter more, and that I think explain the pathologies in the current fits.

### 4a. Replace `distance_to_centerline` with normalized position

Right now cross-stream position is in absolute metres. But 15 m from the centerline is mid-channel at a 41 m-wide section and dry land at a 19 m-wide one. The model is being asked to learn one function of a coordinate that means different things at different places — so it can't, and the optimizer drives the cross-stream length scale to the floor.

Use instead:

```
η(s, n, h) = n / w_side(s, h)        # 0 at the centerline, ±1 at the waterline
```

Now every cross-section is on a common footing, and the physical boundary condition (velocity → 0 at the bank) sits at the same coordinate everywhere.

### 4b. Add local depth — probably the single biggest win

You can now compute, for any point and stage:

```
d(s, n, h) = WSE(s, h) − z_b(s, n)
```

Local water depth is a *far* better predictor of local velocity than distance from a centerline. For a wide channel, conveyance reasoning gives `v ∝ d^(2/3)`; the thalweg is fast because it is deep, not because it is in the middle. Two options, in increasing order of ambition:

1. Add `d` as a fourth GP input.
2. Model the ratio `v / d^(2/3)` (roughly a roughness/slope term) with the GP and multiply back. This puts the known physics in the mean function and leaves the GP to model only the departure from it — which is exactly the kind of structure a GP is good at, and it will let the length scales find sensible values.

**Look at what the current kernels are telling you.** Middle's fitted RBF sits on its lower bound in two dimensions: 15.0 m along a 718 m reach and **1.35 m across a 74 m channel**. Upper's is on the bound in the stage dimension at 8.6 mm. Those aren't length scales, they're the optimizer saying "this coordinate doesn't explain the data, let me memorize instead." Depth and normalized cross-stream position are the coordinates that should explain it.

### 4c. Optional: bank boundary anchors

With an explicit waterline you can add pseudo-observations of `v ≈ 0` along it, with a generous assigned variance (say 0.05 m/s). This stops the GP extrapolating a mid-channel velocity out to the bank. It is a physical prior, not data — label it as such in the methods, report how many anchors and what variance, and show the fit with and without.

Do this only *after* 4a and 4b. If the coordinates are right, the model may reach the bank correctly on its own, and an anchor you didn't need is an assumption you have to defend.

---

## 5. What this does and does not fix

**Fixes:** the rectangular prediction domain; the arbitrary half-width heuristic in the visualization cell; the area denominator; the cross-stream coordinate's inconsistency; and — via local depth — most likely the collapsed length scales.

**Does not fix:** the fundamental stage sampling. You still have 12–17 survey days per reach, and the waterline at an unobserved stage is an *extrapolation of the WSE model*, not an observation. A stage-dependent polygon will look authoritative in a figure; it carries the same uncertainty as `b` in step 2. Keep the shaded stage-gap band on the area plot.

**Also unchanged:** the ~34–53 m transect spacing caps planform resolution. If the paper shows a waterline with metre-scale wiggles, a reviewer will ask where they came from.

---

## 6. Suggested order

1. Confirm the vertical datum on `Elevation` (NAVD88?) and the GNSS method. The 0.17–0.34 ft within-transect WSE scatter suggests RTK, which is fine; autonomous GNSS would not be. Check the Lower section's WSE outlier too — one transect reads 609 ft against a section median of 578 ft, which is either the off-channel LT0101 transect or a bad fix.
2. Build `z_b(s, n)` and `WSE(s, h)`; write `mapping/waterline_{Section}.geojson` from `dataSorting`, not from the visualization notebooks.
3. Cross-check one waterline against lidar or NAIP.
4. Refit with `η` and local depth `d` as inputs, ARD Matérn kernel, wider `length_scale_bounds`. Compare under GroupKFold on survey date — that is the test that will show whether this helped.
5. Only then consider bank anchors.

Step 4 is where you'll find out if this was worth it. My expectation is that Middle's grouped-CV R² moves meaningfully off +0.13 once velocity is regressed on depth rather than on distance-from-a-line, and that the length scales come off their bounds. That is a testable prediction — worth running before committing to the rewrite.

---

**Sources for the external datasets mentioned:**

- [NHDPlus High Resolution — USGS](https://www.usgs.gov/national-hydrography/nhdplus-high-resolution)
- [Access National Hydrography Products — USGS](https://www.usgs.gov/national-hydrography/access-national-hydrography-products)
- [NHDPlus_HR MapServer (layer 8 = NHDArea, 9 = NHDWaterbody)](https://hydro.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/MapServer)
- [USGS Lidar Explorer](https://apps.nationalmap.gov/lidar-explorer/)
- [Michigan 3DEP fact sheet — Michigan DTMB](https://www.michigan.gov/dtmb/-/media/Project/Websites/dtmb/Services/GIS/3DEP.pdf)
- [OpenTopography — USGS 1 m DEM](https://portal.opentopography.org/raster?opentopoID=OTNED.012021.4269.3)
- [Washtenaw County Open Data](https://data-washtenaw.opendata.arcgis.com/)
