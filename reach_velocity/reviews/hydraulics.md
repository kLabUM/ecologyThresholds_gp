**Review of reach_velocity (Middle reach GP), hydraulics and measurements lens**

Numbers marked "recomputed" come from two read-only scripts that rebuild each section from `locations.csv`/`features.csv` via `gpcore` (no GP unless stated): `<reviewer scratch space, not kept>` and `review_checks2.py`.

## Findings, most severe first

**1. The ADCP is misread; its "column reading" rests on an impossible depth, and the field disagrees with it by a factor 2.3.**
(a) Entries 9, 12, 16: the instrument "measures the column above it", local depth = stage − 0.55 (fitted), "the level agrees at mid-stage (A 0.098 vs 0.094)", column exponent 0.46 "is u* = √(ghS)".
(b) Hand-held profile 0618-10 sits 1.7 m from the ADCP with bed −0.04 m on the sensor datum (`locations.csv`): the water over the instrument is ≈0.8 m deep at stage 0.75, not 0.20 m, so h0 = 0.55 (and 0.63, entry 8) contradicts the survey. That profile's column mean is 0.235 m/s; the ADCP at stage 0.73–0.77 reads 0.080 m/s (n 46), which is what the same log profile gives at z ≈ 3 cm. Under the adopted law the field's A at the ADCP is 0.258 against 0.113 from the record, obs/pred 0.44 (recomputed); under the Q law 0.58 (fig 03: bias +108 mm/s, coverage 0.25). The "0.098 vs 0.094 agreement" compares A defined under two different laws. The raw export's only field is `sontek_Vb`, 0.01 m/s steps, tag `stage: Unknown`.
(c) The β = ½ story (00n, 00o, 00t, 00u, entries 9–12) and the validation-table row "ADCP follows the reach law" depend on it. With a realistic depth the record scales as (stage − 0.45)^≈1.0–1.2 under either reading, i.e. it sits on the transects' ridge and identifies nothing about ½.
(d) Obtain the instrument model and configuration (IQ or Argonaut-SW; what Vb is; cell, blanking; its logged vertical-beam depth). Treat it as a near-bed cell or a wake point, redo 03b with the surveyed depth, and remove "β from the ADCP" everywhere.

**2. The stage exponent is a ridge, and the reach's own rating sits at its low edge, not at γ = 1.**
(a) Entries 16, 19: a ∝ (stage − 0.45)^1; Q = 28.1(stage − 0.45)^1.70 read as "consistent".
(b) Area grows almost linearly, so γ = 1 implies Q ∝ (stage − 0.45)^≈2. Carrying every section's measured a to common stages with the script's widths (recomputed): adopted law, median-Q exponent 1.95 over 0.70–1.10 m; local depth^(2/3), 1.49; the measured rating, 1.70, lies between. Model Q/rating: 0.91 at 0.70 m, 1.00 at 0.90, 1.26 at 1.30, 1.36 at 1.50 (entry 18's "≈1000 cfs at 1.3 m": the rating gives 753, 0.77 × gauge ≈ 830). `step3b_law_grid.csv` is flat from (0.35, 1.25) to (0.55, 0.75) at CRPS 0.309–0.310; depth laws 0.329/0.336; eleven days cannot separate them. The two revisits bracket the law: AT03/AT05 show no response where it predicts ×1.45, AT07 ×2.2 where it predicts ×1.39.
(c) Everything above 1.1 m inherits γ; the rating is the one GP-independent constraint and argues for γ ≈ 0.75.
(d) Fit γ jointly to amplitudes and the 20 measured discharges, or choose it by the between-section spread of Q at common stage (recomputed sd log Q at 0.90 m: 0.095 for the stage laws, 0.11–0.13 for the depth laws, 0.18 with none). Report the ridge, not a point.

**3. Above 1.1 m the band is too narrow by construction.**
(a) Entry 17: the law adds 0.06 m/s at 1.30 m; entry 20: May 2026 "5–20% slow water" at 1.4–1.6 m.
(b) `uncertainty.py` samples only the three ridge pairs (0.90→1.30 m multipliers 1.77–1.98); the depth laws give 1.29–1.41. Section share at 1.30 m: 0.14 adopted vs 0.18–0.24 depth laws; at 1.50 m 0.10 vs 0.15–0.18, at or beyond the chain's 0.12 [0.08, 0.18] (`step5_reach_totals.csv`). The discrepancy sd (0.19, iid) is calibrated at survey stages, yet the common-stage continuity spread grows from 0.095 (0.90 m) to 0.14–0.15 at 0.70 and 1.50 m. Flooded margins are absent above the surveyed points.
(c) The 2025-vs-2026 spawning contrast in fig 07 is the most quotable number and the least supported.
(d) Add the depth laws to the scenario set; let the discrepancy sd grow with |stage − 0.9|; hatch everything above 1.09 m in figs 06/07 and quote spring 2026 as a lower bound.

**4. "The slope does not change with stage" and assumption 5 are untested above 0.94 m.**
(a) Entries 10, 12.
(b) One offset per day absorbs any day-to-day sensor–reach mismatch (2025 total-station offsets span 0.41 m), so the 1.3 cm residual is within-day only. Within-day slopes (`rtk_sessions.csv`): 5.2 cm/100 m (06-18, 353 m, 0.74 m), 3.9 (06-25, 119 m, 0.94 m), 6.7, 3.6, 1.1 over 61, 42, 17 m. The two highest survey stages (1.09, 0.99 m) are single transects.
(c) Constant slope is the premise of no backwater and of any depth-based law.
(d) One high-stage RTK day at both reach ends; meanwhile quote 5 ± 2 cm/100 m for 0.74–0.94 m.

**5. The continuity panel (fig 06) is not the validation it is called.**
(a) Entry 19: "the modelled discharge reproduces the measured one (0.73 vs 0.77)".
(b) Q_model is evaluated at each section's own survey stage, where the GP was trained on those points; and `25_step6_rating_curves.py` uses mean_v × Σhw, which is 0.87 of Σvhw (median, recomputed) because v and h correlate; that explains most of the 0.73/0.77 gap.
(c) Readers will take it as independent.
(d) Show Q(stage) for all 20 sections under the law with the rating and 0.77 × gauge overlaid; the between-section spread is the real continuity test.

**6. The nest test validates bank proximity, not the velocity field.**
(a) Entry 20: AUC 0.72, "not an artefact".
(b) Distance-to-bank alone gives AUC 0.76 (KS p 3e-10); within ±1.5 m bank-distance strata the model's velocity gives AUC 0.52 (p 0.16, n 59). Nest cells sit 4.0 m from the bank, wet cells 11.4 m. The 30 m-along restriction cannot expose this because the generic cross-channel shape is the same everywhere.
(c) The ecological deliverable then has no independent check.
(d) Report the stratified AUC; use nests only as a plausibility check on the margin band.

**7. Terminology and physics narrative lag the adopted model.** Hydraulic conveyance is K = AR^(2/3)/n (m³/s); A = a/(stage − 0.45) is a stage-normalised log-law amplitude. Figs 00t/00u derive A = √(gS)/κ (0.17 vs 0.15) for β = ½, since rejected. The adopted a(x,t) = A(x)·g(stage) is a separable index-velocity model (USGS TM 3-A23); name it, redraw 00t/00u for it or drop them, and drop u*/τ_b columns.

**8. Profile reduction: state the limits.** 73/260 locations have one reading, 2026 profiles 3–4; z0 by substrate has sand (14 mm) above gravel (13), so "the ordering is physical" overstates; 70 readings in the bottom tenth run 18% low; the meter is unnamed. Implied Manning n at 0.90 m is 0.05–0.08, falling to 0.033–0.05 at 1.30 m under the chain (recomputed): the relative-submergence effect (h/D84 ≈ 5) is the hydraulic reason an exponent above ⅔ is plausible, and the paper should say so. Name the instrument, exclude z < 0.05 m from fits, report the substrate-z0 sensitivity.

**9. Lower reach: regulated and possibly backwatered.** It sits directly below Flat Rock dam, with two more dams upstream and Lake Erie ~10 km downstream; the 62 cm rise on 15 June at +22% gauge flow (`step8_bed_matching.csv`) may be a release or lake set-up. Check the dam schedule and Lake Erie levels before proposing a sensor there.

## Strong; do not touch
The RTK tie (per-day offsets, one slope, 1.3 cm; sensor zero 253.55 m NAVD88). The night-referenced diurnal stage correction. The transfer check (124 profiles, 30 mm/s). Day-out and transect-out with CRPS and coverage, random folds shown only as the warning. The calibrated band (0.946 before inflation) and readings on held-out days (0.083 m/s, coverage 0.91, n 1,192). Same-day section discharges within 3–12%, and the 0.77 gauge ratio (check it against the USGS drainage areas; it should land near 0.7–0.75). The two tributary-storm days as the case for a local sensor. The attribution design. The Upper-reach honesty. The near-bed reference at 0.10 m lying where the log law is within 3%. The separable law itself, which the common-stage check supports inside the surveyed range.

## The five figures
1. Premise and data: map with the 20 sections, both seasons' stage with survey days, the RTK tie (00q bottom). One campaign, one sensor.
2. Method: profile to amplitude (00m), the target before/after redrawn for the adopted law, the day-out law surface (03b left) with the reach rating overlaid. The spatially distributed at-a-station geometry in one figure.
3. Validation: a transect-out example (02 lower), readings on held-out days (03 lower right), continuity at common stage (new, finding 5), calibration text.
4. Deliverable 1: fig 06's two heatmaps and three rating curves with extrapolation hatched, plus the attribution bars (04).
5. Deliverable 2: fig 07's series, the ecological rating curve (two references, threshold sweep), nests stratified.
Map (05) and other reaches (08, 08a) to the supplement; the ADCP stays in text until finding 1 is resolved.