# Bank-line methods explored (logbook entries 27 to 31)

Scores: each surveyed section held out in turn (the rule of entry 27; lower CRPS is better). Ecological numbers: the Middle reach at the sensor stage 0.90 m, 90% bands in brackets. Blank: not run for that method. Figure: figures/20_bank_methods_comparison.png. Written by 46_bank_methods_comparison.py.

| | A: old lines | B: fixed hybrid | C moving with stage, outline at record high | C moving with stage | C fixed at 0.90 m | C fixed at the median survey level (in use) |
|---|---|---|---|---|---|---|
| What it is | the original notebook's lines, one set for all stages | surveyed edges, 3DEP edges between transects, one set | edges at the sensor's water level each time step; map outline at the record's highest stage | edges at the sensor's water level each time step; map outline at the highest survey-day stage | C's edges at one level, 0.90 m on the sensor, all points inside | C's edges at one level per reach, the median of the survey-day levels, all points inside |
| Level of the fixed lines (Upper / Middle / Lower) |  |  |  |  | Middle only: 0.90 m sensor stage | 262.62 m NAVD88 / 0.877 m sensor / 176.06 m NAVD88 |
| Surveyed points outside the lines (Upper / Middle / Lower) | 6 / 16 / 4 | 0 / 0 / 0 | 0 / 0 / 0 (each point on its own day) | 0 / 0 / 0 (each point on its own day) | 0 (Middle) | 0 / 0 / 0 |
| Upper: section held out, CRPS of ln A (share of points in the 90% band) | 0.450 (0.79) | 0.421 (0.83) |  | 0.405 (0.86) |  | 0.407 (0.84) |
| Middle: section held out, CRPS of ln A (share of points in the 90% band) | 0.314 (0.88) | 0.331 (0.85) |  | 0.310 (0.88) | 0.311 (0.89) | 0.312 (0.89) |
| Lower: section held out, CRPS of ln A (share of points in the 90% band) | 0.357 (0.83) | 0.365 (0.79) |  | 0.335 (0.84) |  | 0.339 (0.85) |
| All sections: CRPS minus A's, mean +/- se (sections better than A) |  | -0.003 +/- 0.011 (34 of 56) |  | -0.024 +/- 0.008 (41 of 56) | -0.003 +/- 0.014 (10 of 20) | -0.021 +/- 0.009 (39 of 56) |
| Run through the whole model? | yes | no (scored under the entry-27 rule only) | yes (Middle) | yes | yes | yes |
| Middle: readings on held-out days, MAE / spread (m/s) | 0.082 / 0.119 |  |  | 0.084 / 0.119 | 0.086 / 0.119 | 0.087 / 0.119 |
| Middle: share of the width below 0.2 m/s at 0.90 m | 0.342 [0.260, 0.423] |  | 0.362 [0.290, 0.439] | 0.362 [0.290, 0.439] | 0.359 [0.294, 0.420] | 0.361 [0.296, 0.425] |
| Middle: mean suitability E at 0.90 m | 0.331 [0.205, 0.567] |  | 0.524 [0.303, 1.299] | 0.431 [0.240, 0.857] | 0.374 [0.233, 0.681] | 0.398 [0.236, 0.738] |
| Middle: E median over 2025 / 2026 | 0.397 / 0.343 |  | 0.517 / 0.446 | 0.545 / 0.483 | 0.471 / 0.416 | 0.494 / 0.442 |
| Middle: velocity effect on nest use (per m/s) | -4.59 +/- 1.49 |  | -2.50 +/- 1.76 | -3.20 +/- 1.84 | -3.81 +/- 1.68 | -3.55 +/- 1.72 |
| Middle: share of draws with velocity helping nest use | below 0.1% |  | 7.7% | 4.0% | 1.1% | 1.7% |
| Upper: E median per season | 2026 0.495 |  |  | 2026 0.586 |  | 2026 0.560 |
| Lower: E median per season | 2025 0.289, 2026 0.232 |  |  | 2025 0.384, 2026 0.319 |  | 2025 0.358, 2026 0.302 |
