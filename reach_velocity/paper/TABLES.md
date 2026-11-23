# Tables for the paper (Middle reach with the stage sensor and the ADCP; 22 September 2026)

Built by `32_paper_tables.py` from `data_clean/`. CSVs sit beside this file.

## Table 1. The surveyed sections

Twenty sections, eleven days. Velocity and share are width-weighted over the section's locations; the discharge is the sum of v h w over them; the last three columns are the chain's prediction for the section with its own day held out.

| section | survey day | along (m) | locations | width (m) | stage (m) | mean velocity measured (m/s) | share below 0.2 m/s measured | discharge measured (m3/s) | mean velocity, day held out (m/s) | 90% lo | 90% hi |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0625-T05 | 2026-06-25 | -215.830 | 6 | 35.980 | 0.938 | 0.319 | 0.206 | 7.885 | 0.387 | 0.149 | 0.865 |
| 0625-T04 | 2026-06-25 | -178.271 | 8 | 36.241 | 0.938 | 0.342 | 0.294 | 8.625 | 0.346 | 0.159 | 0.719 |
| 0625-T03 | 2026-06-25 | -146.835 | 8 | 43.221 | 0.939 | 0.312 | 0.110 | 8.310 | 0.407 | 0.216 | 0.827 |
| 0625-T02 | 2026-06-25 | -121.675 | 10 | 45.009 | 0.941 | 0.340 | 0.289 | 8.648 | 0.422 | 0.229 | 0.771 |
| 0625-T01 | 2026-06-25 | -97.055 | 12 | 48.732 | 0.952 | 0.304 | 0.164 | 8.001 | 0.439 | 0.252 | 0.812 |
| ACORNT01 | 2025-05-30 | -41.379 | 17 | 45.102 | 0.877 | 0.306 | 0.417 | 7.601 | 0.269 | 0.161 | 0.459 |
| 0618-T01 | 2026-06-18 | -1.915 | 12 | 42.802 | 0.748 | 0.145 | 0.684 | 3.485 | 0.224 | 0.110 | 0.445 |
| AT03 | 2025-06-03 | 24.645 | 18 | 48.206 | 0.894 | 0.301 | 0.254 | 6.723 | 0.431 | 0.247 | 0.728 |
| AT05 | 2025-06-10 | 28.342 | 18 | 48.230 | 1.092 | 0.360 | 0.278 | 13.665 | 0.403 | 0.272 | 0.599 |
| 0618-T02 | 2026-06-18 | 53.547 | 11 | 46.261 | 0.741 | 0.134 | 0.780 | 3.240 | 0.187 | 0.124 | 0.322 |
| AT0401 | 2025-06-06 | 81.209 | 18 | 42.499 | 0.987 | 0.294 | 0.277 | 9.541 | 0.299 | 0.213 | 0.423 |
| AT07 | 2025-07-09 | 130.726 | 14 | 39.835 | 0.925 | 0.313 | 0.323 | 7.856 | 0.307 | 0.221 | 0.438 |
| AT06 | 2025-07-14 | 157.119 | 9 | 48.737 | 0.792 | 0.182 | 0.613 | 4.803 | 0.246 | 0.159 | 0.400 |
| AT08 | 2025-07-14 | 174.534 | 12 | 50.663 | 0.795 | 0.222 | 0.427 | 4.382 | 0.198 | 0.124 | 0.317 |
| AT09 | 2025-07-16 | 231.472 | 17 | 47.602 | 0.853 | 0.264 | 0.325 | 5.568 | 0.370 | 0.239 | 0.595 |
| AT10 | 2025-07-21 | 275.193 | 13 | 43.672 | 0.827 | 0.234 | 0.528 | 5.897 | 0.169 | 0.111 | 0.263 |
| AT11 | 2025-07-21 | 317.413 | 15 | 43.948 | 0.830 | 0.233 | 0.438 | 5.307 | 0.242 | 0.144 | 0.393 |
| 0618-T03 | 2026-06-18 | 351.114 | 10 | 48.192 | 0.739 | 0.247 | 0.328 | 3.611 | 0.224 | 0.129 | 0.393 |
| AT12 | 2025-07-23 | 383.877 | 13 | 53.464 | 0.811 | 0.260 | 0.278 | 5.063 | 0.268 | 0.143 | 0.512 |
| AT13 | 2025-07-23 | 445.032 | 16 | 51.416 | 0.813 | 0.197 | 0.667 | 5.179 | 0.260 | 0.112 | 0.565 |

## Table 2. Parameters and settings of the chain

| block | item | value | how set |
|---|---|---|---|
| Readings to one number | vertical profile | v(z) = a ln(z / z0); readings at z >= 3 cm and v > 5 mm/s | fixed form |
|  | roughness scale z0 | 11.2 mm for the reach (flat between 7 and 16 mm), from 160 profiles with three or more readings | fitted, pooled grid search |
|  | amplitude a per location | 257 locations with a > 0; per-location sd from the spread of the readings' implied a (median 0.09 relative) | fitted per location |
| Stage record | sensor | ultrasonic range every 10 min; daytime sensor pattern removed against the night line, median filtered; 2 h means for the series | processing |
|  | stage error in the chain | 1.5 cm, one draw per time shared by every section | assumed |
| Stage law | form | a = A (stage - h0)^gamma, one law for the reach | fixed form |
|  | primary pair | h0 = 0.35 m, gamma = 1.0 (best day-out CRPS) | chosen on a grid |
|  | law set drawn in the chain | every grid pair within 0.01 CRPS of the best with h0 below 0.66 m: (0.35, 1.0), (0.55, 0.75), (0.45, 0.75), (0.35, 1.25), (0.45, 1.0), (0.65, 0.5), (0.55, 0.5), (0.35, 0.75) | uniform over the set |
| Spatial model | target | log A, A = a / (stage - h0)^gamma |  |
|  | inputs | along the reach, relative position across, distance to the nearer bank, relative depth in the section, section mean depth, lateral bed slope | measured |
|  | mean | weighted least squares on an intercept, six substrate classes from the notes and a wake flag | fitted |
|  | kernel | ARD Matern 3/2 plus a white nugget; per-point noise fixed from the readings' spread | fitted by marginal likelihood, 2 to 3 restarts |
|  | fitted scales | amplitude sd 0.87, nugget sd 0.16 (log); length scales: along 198 m, relative across 1.9, bank distance 19 m, relative depth 0.42, section mean depth 0.18 m, lateral slope 1.54 | fitted |
| Uncertainty chain | field | one joint Gaussian process draw over every point of every section | posterior |
|  | discrepancy | inflation 1.241 from held-out days (95% coverage before inflation 0.899, n 257); total 0.390 in log, split into a reach-level part 0.057 (one draw per reach) and a point-level part 0.386 | calibrated |
|  | bed | RTK cross-sections on the sensor datum through per-session offsets (one water-surface slope of 5 cm per 100 m, residual 1.3 cm, 11 sessions); 6 cm per point in the chain | measured; error assumed |
|  | roughness | z0 drawn between 7 and 16 mm, one per draw; a rescaled through the readings' typical height (0.25 m) | assumed range |
|  | draws | 100 per two-hourly stage for the series (3,106 stages), 300 for the rating curves | setting |
| Map (illustration only) | bed between sections | Gaussian process on the RTK sections and bank points (scales 117 m along, 29 m across, nugget 5 cm); 2 m grid inside the mapped bank lines | fitted |
| Ecological number | definition | share of the wetted width below 0.2 m/s on the depth-averaged velocity (alternative: 0.10 m above the bed); width-weighted within a section, length-weighted over the 20 sections; wetted width ends at the outermost surveyed points | choice, to be tied to the ecology |

## Table 3. Validation scores

| test | quantity | value |
|---|---|---|
| Section held out (20 folds, 257 points) | CRPS of log A: six features + notes / bank distance alone / reach mean | 0.319 / 0.321 / 0.414 |
|  | MAE of log A: six features + notes / reach mean | 0.437 / 0.567 |
|  | 95% band coverage, six features + notes | 0.93 |
|  | random ten-fold CRPS for the same three models (the cautionary comparison) | 0.285 / 0.318 / 0.412 |
| Day held out (11 folds, 257 points) | CRPS of log a: stage law (0.35, 1.0) / local depth^0.5 / no law | 0.325 / 0.337 / 0.383 |
|  | CRPS across the law set (8 pairs) | 0.325 to 0.332 |
|  | day level: sd of the held-out day means before / after the law (log) | 0.34 / 0.17 |
| Readings on held-out days (1,208 readings) | MAE / mean absolute deviation of the readings (m/s) | 0.087 / 0.119 |
|  | 95% band coverage | 0.84 |
|  | below or above 0.2 m/s called right / base rate | 73% / 58% |
| Sections on held-out days (20) | measured section mean velocity inside the 90% band | 100% (MAE 0.055 m/s, median relative error 23%) |
|  | measured share below 0.2 m/s inside the 90% band | 95% (MAE 11 points) |
| ADCP through 2025 (shape check, level fitted) | MAE of the law c (stage - 0.35)^1.0 against the 10-min record | 12 mm/s |
|  | sd of the log ratio over the season | 0.159 |
|  | its own exponent on (stage - 0.45) at the surveyed depth: point / column reading | 1.25 / 1.12 |
| Continuity (measured) | 20 section discharges from 11 days on one rating | Q = 25 (stage - 0.45)^2.10 m3/s, 6% scatter |
|  | same-day sections, max/min of the measured discharge | 1.10, 1.11, 1.02, 1.11, 1.10 |
| Continuity (modelled, common stage) | sd of log Q across the 20 sections at 0.70 / 0.90 / 1.30 m | 0.133 / 0.090 / 0.143 |
| Reach share below 0.2 m/s at 0.90 m | sections (width- and length-weighted) / map with margins to the bank line | 0.361 / 0.419 |
|  | 90% band half-width inside the surveyed range (points) | 6.6 |
| Season medians | 2025 record / hours inside the surveyed range | 43% / 38% (65% of hours) |
|  | 2026 record / hours inside the surveyed range | 36% / 36% (50% of hours) |
|  | May 2026 (extrapolated: 24% of its hours above the range) | 27% [21, 34] |
| Nests (59 in wet cells) | AUC of nests against the wet channel: field velocity / bank distance / place factor / depth | 0.71 / 0.75 / 0.73 / 0.50 |
|  | field velocity within bank-distance strata (+-2 m) | 0.49 (p = 0.10) |

## Table 4. Ablations (section held out, stage law)

| kind | variant | CRPS | sections_better_than_full | wilcoxon_p_vs_full |
|---|---|---|---|---|
| features | bank distance only | 0.291 | 12.000 | 0.053 |
| features | geometry only (rel_across, dist_bank, rel_depth) | 0.295 | 13.000 | 0.165 |
| features | full (6) | 0.324 | 0.000 |  |
| features | position only (along, across) | 0.325 | 12.000 | 0.985 |
| features | position + rel_depth | 0.357 | 8.000 | 0.189 |
| features | no section depth | 0.365 | 3.000 | 0.011 |
| features | reach mean | 0.436 | 2.000 | 0.000 |
| kernel | Matern 3/2, 1 restarts | 0.320 |  |  |
| kernel | Matern 3/2, 3 restarts | 0.320 |  |  |
| kernel | Matern 3/2, 6 restarts | 0.320 |  |  |
| kernel | Matern 5/2, 3 restarts | 0.332 |  |  |
| kernel | squared exponential, 3 restarts | 0.348 |  |  |
