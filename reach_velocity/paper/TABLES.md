# Tables for the paper (Middle reach with the stage sensor and the ADCP; 22 September 2026)

Built by `32_paper_tables.py` from `data_clean/`. CSVs sit beside this file.

## Table 1. The surveyed sections

Twenty sections, eleven days. Velocity and share are width-weighted over the section's locations; the discharge is the sum of v h w over them; the last three columns are the chain's prediction for the section with its own day held out.

| section | survey day | along (m) | locations | width (m) | stage (m) | mean velocity measured (m/s) | share below 0.2 m/s measured | discharge measured (m3/s) | mean velocity, day held out (m/s) | 90% lo | 90% hi |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0625-T05 | 2026-06-25 | -215.830 | 6 | 37.430 | 0.938 | 0.319 | 0.206 | 7.893 | 0.345 | 0.158 | 0.732 |
| 0625-T04 | 2026-06-25 | -178.271 | 8 | 38.956 | 0.938 | 0.342 | 0.294 | 8.699 | 0.299 | 0.152 | 0.584 |
| 0625-T03 | 2026-06-25 | -146.835 | 8 | 45.781 | 0.939 | 0.312 | 0.110 | 8.370 | 0.407 | 0.234 | 0.768 |
| 0625-T02 | 2026-06-25 | -121.675 | 10 | 44.539 | 0.941 | 0.340 | 0.289 | 8.652 | 0.410 | 0.215 | 0.828 |
| 0625-T01 | 2026-06-25 | -97.055 | 12 | 48.842 | 0.952 | 0.304 | 0.164 | 8.003 | 0.411 | 0.212 | 0.774 |
| ACORNT01 | 2025-05-30 | -41.379 | 17 | 49.715 | 0.877 | 0.306 | 0.417 | 7.742 | 0.250 | 0.162 | 0.411 |
| 0618-T01 | 2026-06-18 | -1.915 | 12 | 44.152 | 0.748 | 0.145 | 0.684 | 3.519 | 0.213 | 0.111 | 0.384 |
| AT03 | 2025-06-03 | 24.645 | 18 | 47.344 | 0.894 | 0.301 | 0.254 | 6.600 | 0.396 | 0.239 | 0.640 |
| AT05 | 2025-06-10 | 28.342 | 18 | 54.181 | 1.092 | 0.360 | 0.278 | 13.873 | 0.369 | 0.271 | 0.508 |
| 0618-T02 | 2026-06-18 | 53.547 | 11 | 44.362 | 0.741 | 0.134 | 0.780 | 3.238 | 0.168 | 0.111 | 0.269 |
| AT0401 | 2025-06-06 | 81.209 | 18 | 44.949 | 0.987 | 0.294 | 0.277 | 9.601 | 0.261 | 0.190 | 0.372 |
| AT07 | 2025-07-09 | 130.726 | 14 | 36.921 | 0.925 | 0.313 | 0.323 | 7.696 | 0.310 | 0.226 | 0.427 |
| AT06 | 2025-07-14 | 157.119 | 9 | 50.963 | 0.792 | 0.182 | 0.613 | 4.808 | 0.220 | 0.147 | 0.342 |
| AT08 | 2025-07-14 | 174.534 | 12 | 47.679 | 0.795 | 0.222 | 0.427 | 4.346 | 0.205 | 0.128 | 0.316 |
| AT09 | 2025-07-16 | 231.472 | 17 | 53.383 | 0.853 | 0.264 | 0.325 | 5.590 | 0.374 | 0.251 | 0.632 |
| AT10 | 2025-07-21 | 275.193 | 13 | 49.062 | 0.827 | 0.234 | 0.528 | 5.939 | 0.263 | 0.127 | 0.444 |
| AT11 | 2025-07-21 | 317.413 | 15 | 47.247 | 0.830 | 0.233 | 0.438 | 5.324 | 0.272 | 0.154 | 0.494 |
| 0618-T03 | 2026-06-18 | 351.114 | 10 | 45.181 | 0.739 | 0.247 | 0.328 | 3.604 | 0.220 | 0.132 | 0.368 |
| AT12 | 2025-07-23 | 383.877 | 13 | 58.434 | 0.811 | 0.260 | 0.278 | 5.084 | 0.249 | 0.145 | 0.421 |
| AT13 | 2025-07-23 | 445.032 | 16 | 51.195 | 0.813 | 0.197 | 0.667 | 5.223 | 0.250 | 0.109 | 0.557 |

## Table 2. Parameters and settings of the chain

| block | item | value | how set |
|---|---|---|---|
| Readings to one number | vertical profile | v(z) = a ln(z / z0); readings at z >= 3 cm and v > 5 mm/s | fixed form |
|  | roughness scale z0 | 11.2 mm for the reach (flat between 7 and 16 mm), from 160 profiles with three or more readings | fitted, pooled grid search |
|  | amplitude a per location | 257 locations with a > 0; per-location sd from the spread of the readings' implied a (median 0.09 relative) | fitted per location |
| Stage record | sensor | ultrasonic range every 10 min; daytime sensor pattern removed against the night line, median filtered; 2 h means for the series | processing |
|  | stage error in the chain | 1.5 cm, one draw per time shared by every section | assumed |
| Stage law | form | a = A (stage - h0)^gamma, one law for the reach | fixed form |
|  | primary pair | h0 = 0.45 m, gamma = 1.0 (best day-out CRPS) | chosen on a grid |
|  | law set drawn in the chain | every grid pair within 0.01 CRPS of the best with h0 below 0.66 m: (0.45, 1.0), (0.35, 1.25), (0.55, 0.75), (0.35, 1.0), (0.45, 0.75), (0.65, 0.5), (0.55, 0.5) | uniform over the set |
| Spatial model | target | log A, A = a / (stage - h0)^gamma |  |
|  | inputs | along the reach, relative position across, distance to the nearer bank, relative depth in the section, section mean depth, lateral bed slope | measured |
|  | mean | weighted least squares on an intercept, six substrate classes from the notes and a wake flag | fitted |
|  | kernel | ARD Matern 3/2 plus a white nugget; per-point noise fixed from the readings' spread | fitted by marginal likelihood, 2 to 3 restarts |
|  | fitted scales | amplitude sd 0.84, nugget sd 0.12 (log); length scales: along 159 m, relative across 2.2, bank distance 16 m, relative depth 0.42, section mean depth 0.18 m, lateral slope 0.52 | fitted |
| Uncertainty chain | field | one joint Gaussian process draw over every point of every section | posterior |
|  | discrepancy | inflation 1.056 from held-out days (95% coverage before inflation 0.946, n 257); total 0.187 in log, split into a reach-level part 0.098 (one draw per reach) and a point-level part 0.159 | calibrated |
|  | bed | RTK cross-sections on the sensor datum through per-session offsets (one water-surface slope of 5 cm per 100 m, residual 1.3 cm, 11 sessions); 6 cm per point in the chain | measured; error assumed |
|  | roughness | z0 drawn between 7 and 16 mm, one per draw; a rescaled through the readings' typical height (0.25 m) | assumed range |
|  | draws | 100 per two-hourly stage for the series (3,106 stages), 300 for the rating curves | setting |
| Map (illustration only) | bed between sections | Gaussian process on the RTK sections and bank points (scales 117 m along, 29 m across, nugget 5 cm); 2 m grid inside the mapped bank lines | fitted |
| Ecological number | definition | share of the wetted width below 0.2 m/s on the depth-averaged velocity (alternative: 0.10 m above the bed); width-weighted within a section, length-weighted over the 20 sections; wetted width ends at the outermost surveyed points | choice, to be tied to the ecology |

## Table 3. Validation scores

| test | quantity | value |
|---|---|---|
| Section held out (20 folds, 257 points) | CRPS of log A: six features + notes / bank distance alone / reach mean | 0.307 / 0.320 / 0.414 |
|  | MAE of log A: six features + notes / reach mean | 0.417 / 0.567 |
|  | 95% band coverage, six features + notes | 0.93 |
|  | random ten-fold CRPS for the same three models (the cautionary comparison) | 0.280 / 0.315 / 0.412 |
| Day held out (11 folds, 257 points) | CRPS of log a: stage law (0.45, 1.0) / local depth^0.5 / no law | 0.309 / 0.329 / 0.394 |
|  | CRPS across the law set (seven pairs) | 0.309 to 0.318 |
|  | day level: sd of the held-out day means before / after the law (log) | 0.34 / 0.17 |
| Readings on held-out days (1,208 readings) | MAE / mean absolute deviation of the readings (m/s) | 0.082 / 0.119 |
|  | 95% band coverage | 0.88 |
|  | below or above 0.2 m/s called right / base rate | 76% / 58% |
| Sections on held-out days (20) | measured section mean velocity inside the 90% band | 100% (MAE 0.048 m/s, median relative error 18%) |
|  | measured share below 0.2 m/s inside the 90% band | 95% (MAE 12 points) |
| ADCP through 2025 (shape check, level fitted) | MAE of the law c (stage - 0.45) against the 10-min record | 8 mm/s |
|  | sd of the log ratio over the season | 0.110 |
|  | its own exponent on (stage - 0.45) at the surveyed depth: point / column reading | 1.25 / 1.12 |
| Continuity (measured) | 20 section discharges from 11 days on one rating | Q = 28 (stage - 0.45)^1.70 m3/s, 7% scatter |
|  | same-day sections, max/min of the measured discharge | 1.11, 1.12, 1.03, 1.11, 1.10 |
| Continuity (modelled, common stage) | sd of log Q across the 20 sections at 0.70 / 0.90 / 1.30 m | 0.135 / 0.092 / 0.140 |
| Reach share below 0.2 m/s at 0.90 m | sections (width- and length-weighted) / map with margins to the bank line | 0.342 / 0.378 |
|  | 90% band half-width inside the surveyed range (points) | 7.6 |
| Season medians | 2025 record / hours inside the surveyed range | 41% / 36% (65% of hours) |
|  | 2026 record / hours inside the surveyed range | 34% / 34% (50% of hours) |
|  | May 2026 (extrapolated: 24% of its hours above the range) | 24% [18, 32] |
| Nests (59 in wet cells) | AUC of nests against the wet channel: field velocity / bank distance / place factor / depth | 0.72 / 0.76 / 0.74 / 0.50 |
|  | field velocity within bank-distance strata (+-2 m) | 0.53 (p = 0.20) |

## Table 4. Ablations (section held out, stage law)

| kind | variant | CRPS | sections_better_than_full | wilcoxon_p_vs_full |
|---|---|---|---|---|
| features | bank distance only | 0.299 | 15.000 | 0.143 |
| features | geometry only (rel_across, dist_bank, rel_depth) | 0.310 | 12.000 | 0.701 |
| features | full (6) | 0.317 | 0.000 |  |
| features | position only (along, across) | 0.325 | 10.000 | 0.756 |
| features | no section depth | 0.354 | 8.000 | 0.058 |
| features | position + rel_depth | 0.357 | 6.000 | 0.033 |
| features | reach mean | 0.436 | 2.000 | 0.000 |
| kernel | Matern 3/2, 1 restarts | 0.317 |  |  |
| kernel | Matern 3/2, 3 restarts | 0.317 |  |  |
| kernel | Matern 3/2, 6 restarts | 0.317 |  |  |
| kernel | Matern 5/2, 3 restarts | 0.321 |  |  |
| kernel | squared exponential, 3 restarts | 0.341 |  |  |
