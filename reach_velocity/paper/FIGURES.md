# The five paper figures (revised after the reviewer pass, 22 September 2026)

Built by `29_paper_figures.py` from the saved tables in `data_clean/` (no refitting); PNG at 300 dpi and PDF in this folder. One message per figure. Captions are drafts in the paper's voice. Numbers changed between the first draft and this one where entry 22 of the logbook says so (the joint reach sampler, width-weighted shares, the widened law set, the re-read ADCP, the nest baselines).

## F1. The site, the data, and the target (`F1_site_data_target`)

**Message.** A hand-held campaign on 20 sections plus one stage record are the inputs; the target is one number per location that does not depend on the day.

**Caption.** The Middle reach of the Huron River (Michigan) and its inputs. (a) The 20 surveyed cross-sections (257 locations, 1,285 readings at three to five heights above the bed), the ultrasonic stage sensor and the bed-mounted acoustic Doppler instrument at the origin. (b) The stage record every 10 minutes in 2025 and 2026 with the survey days marked and the ADCP recording period shaded; the surveys covered 0.74 to 1.09 m of stage, the 2026 season reached 1.65 m. (c) Each location's readings reduce to the amplitude a of a logarithmic profile with one roughness scale for the reach (z0 = 11 mm); three locations are shown. (d) The amplitude on its survey day rises with stage (day medians with interquartile bars); the dashed line is the one-parameter law a = A (stage − 0.45 m). (e) Dividing by the law leaves the place factor A, a stage-normalised amplitude with no trend against stage, and it is A that the spatial model learns.

**Source.** `locations.csv`, `readings_transect.csv`, `series_10min.csv`, `step3_adcp.csv`, bank lines from `data_raw/`.

## F2. What is learned (`F2_conveyance_field`)

**Message.** The place factor is at-a-station geometry, mostly distance to the bank; a held-out section is predicted better than the reach mean, and random folds flatter every model.

**Caption.** The place factor and its held-out error. (a) A against distance to the nearer bank for every location, coloured by the section's mean depth; bank distance alone predicts a held-out section as well as six features (continuous ranked probability score 0.30 against 0.32, better on 15 of 20 sections, p = 0.14). (b) One section held out and predicted from the other nineteen by the six-feature Gaussian process, with its 95% band and the reach mean; the field beats the reach mean on 16 of 20 sections. (c) The held-out score by model, sections held out against random ten-fold cross-validation; random folds leak neighbouring points and flatter every model. (d) The field on a 2 m grid, an illustration: the posterior standard deviation of log A is 0.31 between sections against 0.23 on them, and every number in the paper comes from the sections.

**Source.** `features.csv`, `step2_predictions.csv`, `step2_metrics.csv`, `step2_feature_ablation.csv`, `grid.csv`.

## F3. The stage law: identification and closure (`F3_stage_law`)

**Message.** The day level is stage, the law's exponent is a ridge that three weak constraints bracket, readings on held-out days are predicted to 0.08 m/s, and the sections' discharges are consistent at common stage and lie on one measured rating.

**Caption.** The stage law and its tests. (a) Each survey day held out: the mean residual of log a from the spatial field alone against the day's stage (filled, sd 0.30 in log) and after the law a = A (stage − 0.45 m)^0.75 (open, sd 0.14). (b) The day-out score for every (h0, γ) pair on a grid: eleven days cannot separate the pairs within optimiser noise (grey band, the set the uncertainty chain draws from), the measured rating implies γ ≈ 0.75 for h0 = 0.45 m (orange), and the ADCP re-read at its surveyed depth gives 1.1 to 1.25 (green); the law is a ridge, not a point, and the band above the surveyed range says so. (c) The 1,208 readings on held-out days against their prediction at their own height, day held out: mean absolute error 0.082 m/s against a mean absolute deviation of 0.119 m/s, below or above 0.2 m/s called right for 74% against a base rate of 58%. (d) Continuity at common stage: every section's modelled discharge (Σ v h w per draw) over 0.70 to 1.30 m, the discharges measured from the profiles on the survey days (squares, 20 sections on 11 days) and the rating through them, Q = 28 (stage − 0.45)^1.70 m³/s with 6% scatter; the spread across sections at a common stage is 11% at 0.90 m and 14% at 0.70 and 1.30 m.

**Source.** `step3_day_offsets.csv`, `step3b_law_grid.csv`, `adcp_recheck.json`, `step6_readings_heldout_adopted_law.csv`, `step6_continuity.csv`, `step6_continuity_common_stage.csv`, `step6_reach_rating.json`.

## F4. Spatial rating curves for two seasons (`F4_spatial_rating_curves`)

**Message.** The deliverable: every section's mean velocity every two hours through two seasons from the stage record, with a band that is the field and the day-to-day discrepancy inside the surveyed range and the law set above it.

**Caption.** Spatial rating curves. (a, b) The width-weighted mean depth-averaged velocity of all 20 sections every two hours in 2025 and 2026 from the stage record (orange), sections ordered downstream; the ticks below each panel mark hours when the stage was outside the surveyed range (35% of the 2025 record, 50% of 2026). (c) Three sections' rating curves with 90% bands from the joint Monte Carlo chain and the section mean velocity measured on each one's survey day; hatched above the surveyed range, where the law set diverges. (d) The 90% band on one section's mean velocity when each source is switched on alone, at 0.90 m and at 1.30 m: inside the range the field and the day-to-day discrepancy dominate; above it the law set does.

**Source.** `step6_section_series.csv`, `step4_attribution.csv`, `locations.csv`.

## F5. The ecological number (`F5_ecological_number`)

**Message.** The share of the reach below an ecological velocity threshold every two hours for two seasons with its band, stated with its reference height and threshold, and the nests as a plausibility check, not a validation.

**Caption.** The ecological number. (a, b) The share of the reach's wetted width below 0.2 m/s every two hours in 2025 and 2026, width-weighted within sections and length-weighted over them, with 90% bands from the joint chain (median 43% in 2025 and 36% in 2026 over the records, 38% and 36% over the hours inside the surveyed range; band ±6.5 points inside the range); ticks below mark hours outside the surveyed range. (c) The ecological rating curve: the share against stage for four thresholds on the depth-averaged velocity and for 0.2 m/s at 0.10 m above the bed (dashed), with the band on 0.2 m/s and the extrapolated range hatched; the reference height moves the number by 8 to 10 points, as much as the band, and the margins beyond the outermost surveyed points are not counted, so the share is a lower bound at every stage. (d) The 57 smallmouth bass nests in wet cells ranked against every wet cell of the reach on their day by the model's velocity (area under the curve 0.71), by bank distance alone (0.75), by the place factor (0.73) and by depth (0.50), and by velocity within bank-distance strata of ±2 m (0.49): the nests sit in the slow margin, as bank distance alone predicts, and the 2 m field adds no discrimination within it.

**Source.** `step7_reach_series.csv`, `step7_sensitivity.csv`, `step7_summary.json`.

## What did not make the five, and where it goes

- The transfer to the Upper and Lower reaches (figure 08 in `figures/`): a table in the text from `step8_summary.csv` and one paragraph (logbook entry 21). The spatial field reproduces in the Lower (CRPS 0.37 against 0.54) and not in the slow Upper (0.43 against 0.38); without a local water level there is no season series, and the RTK water surface shows that no gauge supplies it.
- The ADCP through 2025 (figure 03, `adcp_recheck.json`): one paragraph. It is a near-bed value whose stage response sits on the transects' ridge; its level cannot be compared with the field until the instrument is known.
- The velocity maps at two stages and the sd map (figure 05): supplement.
- The data hygiene (the 'Other' rows, the height column, the transect relabels, the BT transects, the RTK sessions without a fix): methods and supplement.
- The stage record cleaning: supplement figure. The uncertainty attribution at four stages and four sections (figure 04): supplement. The kernel and feature ablations: one sentence each in results (`step2_kernel_ablation.csv`, `step2_feature_ablation.csv`).
