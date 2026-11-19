# Tables for the paper (Middle reach with the stage sensor and the ADCP; 22 September 2026)

Built by `32_paper_tables.py` from `data_clean/`. CSVs sit beside this file.

## Table 1. The surveyed sections

Twenty sections, eleven days. Velocity and share are width-weighted over the section's locations; the discharge is the sum of v h w over them; the last three columns are the chain's prediction for the section with its own day held out.

| section | survey day | along (m) | locations | width (m) | stage (m) | mean velocity measured (m/s) | share below 0.2 m/s measured | discharge measured (m3/s) | mean velocity, day held out (m/s) | 90% lo | 90% hi |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0625-T05 | 2026-06-25 | -215.830 | 6 | 36.837 | 0.938 | 0.319 | 0.206 | 7.896 | 0.389 | 0.159 | 0.827 |
| 0625-T04 | 2026-06-25 | -178.271 | 8 | 36.343 | 0.938 | 0.342 | 0.294 | 8.632 | 0.344 | 0.165 | 0.741 |
| 0625-T03 | 2026-06-25 | -146.835 | 8 | 43.403 | 0.939 | 0.312 | 0.110 | 8.316 | 0.410 | 0.220 | 0.816 |
| 0625-T02 | 2026-06-25 | -121.675 | 10 | 45.904 | 0.941 | 0.340 | 0.289 | 8.678 | 0.428 | 0.227 | 0.741 |
| 0625-T01 | 2026-06-25 | -97.055 | 12 | 48.732 | 0.952 | 0.304 | 0.164 | 8.001 | 0.447 | 0.255 | 0.868 |
| ACORNT01 | 2025-05-30 | -41.379 | 17 | 45.102 | 0.877 | 0.306 | 0.417 | 7.601 | 0.252 | 0.153 | 0.451 |
| 0618-T01 | 2026-06-18 | -1.915 | 12 | 40.841 | 0.748 | 0.145 | 0.684 | 3.484 | 0.206 | 0.101 | 0.400 |
| AT03 | 2025-06-03 | 24.645 | 18 | 47.390 | 0.894 | 0.301 | 0.254 | 6.693 | 0.420 | 0.237 | 0.721 |
| AT05 | 2025-06-10 | 28.342 | 18 | 49.768 | 1.092 | 0.360 | 0.278 | 13.758 | 0.411 | 0.266 | 0.590 |
| 0618-T02 | 2026-06-18 | 53.547 | 11 | 45.648 | 0.741 | 0.134 | 0.780 | 3.239 | 0.181 | 0.109 | 0.297 |
| AT0401 | 2025-06-06 | 81.209 | 18 | 43.782 | 0.987 | 0.294 | 0.277 | 9.567 | 0.307 | 0.208 | 0.441 |
| AT07 | 2025-07-09 | 130.726 | 14 | 40.642 | 0.925 | 0.313 | 0.323 | 7.901 | 0.307 | 0.227 | 0.437 |
| AT06 | 2025-07-14 | 157.119 | 9 | 47.600 | 0.792 | 0.182 | 0.613 | 4.798 | 0.245 | 0.154 | 0.427 |
| AT08 | 2025-07-14 | 174.534 | 12 | 48.220 | 0.795 | 0.222 | 0.427 | 4.350 | 0.199 | 0.121 | 0.336 |
| AT09 | 2025-07-16 | 231.472 | 17 | 47.419 | 0.853 | 0.264 | 0.325 | 5.567 | 0.354 | 0.230 | 0.549 |
| AT10 | 2025-07-21 | 275.193 | 13 | 42.961 | 0.827 | 0.234 | 0.528 | 5.886 | 0.160 | 0.098 | 0.265 |
| AT11 | 2025-07-21 | 317.413 | 15 | 42.796 | 0.830 | 0.233 | 0.438 | 5.305 | 0.225 | 0.119 | 0.380 |
| 0618-T03 | 2026-06-18 | 351.114 | 10 | 42.133 | 0.739 | 0.247 | 0.328 | 3.604 | 0.210 | 0.122 | 0.360 |
| AT12 | 2025-07-23 | 383.877 | 13 | 53.464 | 0.811 | 0.260 | 0.278 | 5.063 | 0.260 | 0.137 | 0.516 |
| AT13 | 2025-07-23 | 445.032 | 16 | 49.116 | 0.813 | 0.197 | 0.667 | 5.137 | 0.257 | 0.112 | 0.562 |

## Table 2. Parameters and settings of the chain

| block | item | value | how set |
|---|---|---|---|
| Readings to one number | vertical profile | v(z) = a ln(z / z0); readings at z >= 3 cm and v > 5 mm/s | fixed form |
|  | roughness scale z0 | 11.2 mm for the reach (flat between 7 and 16 mm), from 160 profiles with three or more readings | fitted, pooled grid search |
|  | amplitude a per location | 257 locations with a > 0; per-location sd from the spread of the readings' implied a (median 0.09 relative) | fitted per location |
| Stage record | sensor | ultrasonic range every 10 min; daytime sensor pattern removed against the night line, median filtered; 2 h means for the series | processing |
|  | stage error in the chain | 1.5 cm, one draw per time shared by every section | assumed |
| Stage law | form | a = A (stage - h0)^gamma, one law for the reach | fixed form |
|  | primary pair | h0 = 0.45 m, gamma = 0.75 (best day-out CRPS) | chosen on a grid |
|  | law set drawn in the chain | every grid pair within 0.01 CRPS of the best with h0 below 0.66 m: (0.45, 0.75), (0.35, 1.0), (0.55, 0.5), (0.35, 0.75), (0.55, 0.75), (0.45, 1.0), (0.65, 0.5), (0.35, 1.25) | uniform over the set |
| Spatial model | target | log A, A = a / (stage - h0)^gamma |  |
|  | inputs | along the reach, relative position across, distance to the nearer bank, relative depth in the section, section mean depth, lateral bed slope | measured |
|  | mean | weighted least squares on an intercept, six substrate classes from the notes and a wake flag | fitted |
|  | kernel | ARD Matern 3/2 plus a white nugget; per-point noise fixed from the readings' spread | fitted by marginal likelihood, 2 to 3 restarts |
|  | fitted scales | amplitude sd 0.84, nugget sd 0.14 (log); length scales: along 197 m, relative across 2.1, bank distance 17 m, relative depth 0.40, section mean depth 0.17 m, lateral slope 1.54 | fitted |
| Uncertainty chain | field | one joint Gaussian process draw over every point of every section | posterior |
|  | discrepancy | inflation 1.159 from held-out days (95% coverage before inflation 0.911, n 257); total 0.312 in log, split into a reach-level part 0.070 (one draw per reach) and a point-level part 0.304 | calibrated |
|  | bed | RTK cross-sections on the sensor datum through per-session offsets (one water-surface slope of 5 cm per 100 m, residual 1.3 cm, 11 sessions); 6 cm per point in the chain | measured; error assumed |
|  | roughness | z0 drawn between 7 and 16 mm, one per draw; a rescaled through the readings' typical height (0.25 m) | assumed range |
|  | draws | 100 per two-hourly stage for the series (3,106 stages), 300 for the rating curves | setting |
| Map (illustration only) | bed between sections | Gaussian process on the RTK sections and bank points (scales 117 m along, 29 m across, nugget 5 cm); 2 m grid inside the mapped bank lines | fitted |
| Ecological number | definition | share of the wetted width below 0.2 m/s on the depth-averaged velocity (alternative: 0.10 m above the bed); width-weighted within a section, length-weighted over the 20 sections; wetted width ends at the outermost surveyed points | choice, to be tied to the ecology |

## Table 3. Validation scores

| test | quantity | value |
|---|---|---|
| Section held out (20 folds, 257 points) | CRPS of log A: six features + notes / bank distance alone / reach mean | 0.311 / 0.316 / 0.414 |
|  | MAE of log A: six features + notes / reach mean | 0.425 / 0.567 |
|  | 95% band coverage, six features + notes | 0.93 |
|  | random ten-fold CRPS for the same three models (the cautionary comparison) | 0.283 / 0.314 / 0.412 |
| Day held out (11 folds, 257 points) | CRPS of log a: stage law (0.45, 0.75) / local depth^0.5 / no law | 0.321 / 0.328 / 0.373 |
|  | CRPS across the law set (8 pairs) | 0.321 to 0.325 |
|  | day level: sd of the held-out day means before / after the law (log) | 0.34 / 0.17 |
| Readings on held-out days (1,208 readings) | MAE / mean absolute deviation of the readings (m/s) | 0.084 / 0.119 |
|  | 95% band coverage | 0.85 |
|  | below or above 0.2 m/s called right / base rate | 74% / 58% |
| Sections on held-out days (20) | measured section mean velocity inside the 90% band | 100% (MAE 0.055 m/s, median relative error 24%) |
|  | measured share below 0.2 m/s inside the 90% band | 95% (MAE 12 points) |
| ADCP through 2025 (shape check, level fitted) | MAE of the law c (stage - 0.45)^0.75 against the 10-min record | 14 mm/s |
|  | sd of the log ratio over the season | 0.169 |
|  | its own exponent on (stage - 0.45) at the surveyed depth: point / column reading | 1.25 / 1.12 |
| Continuity (measured) | 20 section discharges from 11 days on one rating | Q = 28 (stage - 0.45)^1.70 m3/s, 6% scatter |
|  | same-day sections, max/min of the measured discharge | 1.10, 1.11, 1.01, 1.11, 1.10 |
| Continuity (modelled, common stage) | sd of log Q across the 20 sections at 0.70 / 0.90 / 1.30 m | 0.143 / 0.113 / 0.137 |
| Reach share below 0.2 m/s at 0.90 m | sections (width- and length-weighted) / map with margins to the bank line | 0.362 / 0.454 |
|  | 90% band half-width inside the surveyed range (points) | 7.7 |
| Season medians | 2025 record / hours inside the surveyed range | 45% / 39% (65% of hours) |
|  | 2026 record / hours inside the surveyed range | 35% / 35% (50% of hours) |
|  | May 2026 (extrapolated: 24% of its hours above the range) | 23% [18, 30] |
| Nests (59 in wet cells) | AUC of nests against the wet channel: field velocity / bank distance / place factor / depth | 0.69 / 0.72 / 0.71 / 0.47 |
|  | field velocity within bank-distance strata (+-2 m) | 0.46 (p = 0.00) |

## Table 4. Ablations (section held out, stage law)

| kind | variant | CRPS | sections_better_than_full | wilcoxon_p_vs_full |
|---|---|---|---|---|
| features | bank distance only | 0.291 | 14.000 | 0.027 |
| features | geometry only (rel_across, dist_bank, rel_depth) | 0.294 | 14.000 | 0.070 |
| features | full (6) | 0.320 | 0.000 |  |
| features | position only (along, across) | 0.325 | 10.000 | 1.000 |
| features | no section depth | 0.351 | 6.000 | 0.097 |
| features | position + rel_depth | 0.357 | 7.000 | 0.090 |
| features | reach mean | 0.436 | 2.000 | 0.000 |
| kernel | Matern 3/2, 1 restarts | 0.320 |  |  |
| kernel | Matern 3/2, 3 restarts | 0.320 |  |  |
| kernel | Matern 3/2, 6 restarts | 0.320 |  |  |
| kernel | Matern 5/2, 3 restarts | 0.332 |  |  |
| kernel | squared exponential, 3 restarts | 0.348 |  |  |
