# Working in this folder (instructions for a coding agent)

This folder is a complete, reproducible analysis: a season-long velocity record for a river reach from one survey campaign, one stage sensor and one ADCP, and an ecological number learned from nest sites. A person new to it will ask you to explain the method, rerun it, check a number, change a setting, or extend it to another reach. Read in this order before doing anything: `README.md`, `paper/METHODS_AND_RESULTS.md`, `figures/17_model_variables.png` (the variables), then the parts of `LOGBOOK.md` the question touches (the entries are dated and numbered; the status section at the end lists the numbers to quote).

## Running and verifying

- Environment: Python 3.11 with `requirements.txt` (pinned; `requirements_frozen.txt` is the full environment). The scripts expect the interpreter at `../.venv/bin/python` by default; `PY=/path/to/python zsh run_all.sh` overrides it.
- `zsh run_all.sh` reruns everything from `data_raw/` in about ten minutes and rewrites `data_clean/`, `figures/` and `paper/`. Seeds are fixed; a rerun reproduces the numbers in `paper/TABLES.md` to the printed precision.
- To verify a claim: find its number in `paper/TABLES.md` or the logbook, open the `data_clean/` table named there, and read the script that wrote it (each script's header docstring names its inputs and outputs). `RERUN.md` lists the headline numbers and a verification procedure.
- Every script is standalone once its inputs exist. To rerun from step 2 after changing model code: `20_` to `26_`, then `30_` to `36_`, then `41_`, `29_`, `32_`, `40_`.

## The pipeline in one look

| scripts | what they do | outputs |
|---|---|---|
| `00_load.py` (runs twice) | raw sheets and sensor exports to clean tables; units converted here and nowhere else; named corrections listed in the code | `readings_transect.csv`, `nests.csv`, `stage_10min.csv`, `adcp.csv`, `usgs_*daily.csv` |
| `02_` to `07_` | the stage record: daytime artefact diagnosed, nights from the sun's elevation, the correction (night line, monthly template, median), the 10-min grid shared with the ADCP | `series_10min.csv` (stage_m, adcp_ms) |
| `08_`, `12_`, `13_` | log profiles: z0 for the reach, one amplitude a per location with its error, checks by substrate, year and crew | `locations.csv`, `field_qc.json`, `amplitude_*.json` |
| `14_`, `15_`, `19_` | transect labels checked by geometry; RTK sections tied to the sensor datum; bank set C and the candidate features (bank distance and relative position from C at each point's survey-day level, entry 30) | `features.csv`, `rtk_session_offsets.csv`, `banks.csv` |
| `gpcore.py` | `load_points()`, `destage()`, `with_target()` (A = a / g, y = log A, per-point sd), `Model` (WLS mean on notes + ARD Matérn GP), `scores()` (MAE, RMSE, CRPS, coverage) | |
| `uncertainty.py` | `LAW_RIDGE` (the law set), `Chain` (GPs for every law and every C draw, calibration on held-out days, the discrepancy split), `Field` (the field on the 5 cm stage grid, read at each draw's stage), `sample`, `sample_series`, `sample_reach` (one joint draw of the whole reach per sample), `GridLevels` | `uncertainty_calibration.json` |
| `20_` to `26_` | steps 2 to 7: the GP on A with sections held out; the day effect and the laws; the (h0, γ) grid; the chain's calibration and attribution; the bed GP and 2 m grid; the section series, continuity and section-level coverage; the ecological share | `step2_*` to `step7_*`, `grid.csv` |
| `30_`, `31_`, `34_`, `33_`, `35_`, `36_` | feature and kernel ablations; the ADCP re-read and single-beam test; nest-calibrated thresholds and the selection model; the Bayesian suitability; the explanatory figure | `step2_*_ablation.csv`, `adcp_*.json`, `nest_*.json/csv` |
| `41_`, `29_`, `32_`, `40_` | the variables figure, the paper's F1 to F5, the tables, the methods document (docx and md), all from saved outputs | `paper/` |
| `27_`, `39_`, `37_`, `28_`, `38_` | the other two reaches, in this order (paused; not in `run_all.sh`); `37_` builds their bank set C from the sensor tie | `reach_*`, `step8_*`, figures 08 to 16 |
| `44_banklines_survey.py` | the centrelines and the old bank lines (set A): the main branch's dataSorting.ipynb (commit 2346281) carried over; fingerprint check against the copies used until November 2026 | `data_clean/banklines_survey/` |
| `banklines_3dep.py` | bank lines from USGS 3DEP plus the survey: the 3DEP transects (`transect_points`, `fetch_transects`, cache checked against the geometry), sets B and C, and `StageFeatures`, set C as the model uses it (edges, distance to the bank and relative position at any stage, best estimate or one of 20 draws) | `data_raw/3dep_bank_transects_*.csv`, `data_clean/banks_c/` (via `19_`, `37_`) |
| `43_bankline_comparison.py` | the comparison of A, B and C under the rule in logbook entry 27 (the record of entries 27 and 28) | `data_clean/bankcmp/`, figures 18 |
| `provenance.py` | every line or polygon file and the script that writes it; `00_load.py` stops if one has none | |

## Rules that keep the analysis honest

1. Numbers come from the surveyed sections; the 2 m map is illustration. Do not quote a map value as a result.
2. Claims are scored with sections held out, days held out, readings on held-out days, or an independent record; never with random folds (they leak neighbours and are shown once as the cautionary comparison).
3. The per-point noise is the readings' own spread and stays fixed; the nugget is learned on top of it.
4. The stage law's exponent is a ridge. If you change `LAW_RIDGE` or the primary law, rerun steps 4 to 7 and the paper scripts, and say so in the logbook.
5. Units are converted once, in `00_load.py`. Nothing downstream sees feet.
6. No personal names anywhere: crews are coded in the raw sheets and in every output; people are "the PI", "the field lead", "the crew".
7. Every change of method gets a numbered logbook entry with a figure and the numbers before and after; every figure has a caption that says what it shows; plain language, numbers with units, no em-dashes.
8. Paper figures and tables are built only from saved tables (`29_`, `32_`, `40_`), so a rerun of the chain and a rebuild of the figures are separate, checkable steps.
9. Every bank line, channel polygon and centreline in the folder is written by a script here (`provenance.py` lists them; `00_load.py` stops otherwise). The model uses bank set C (entry 30); A and B stay only as the record of the comparison.

## Common requests, and where to go

- "Explain A / the place factor": `figures/17_model_variables.png`, logbook entries 11 to 13, methods document Section 3.
- "Why is the band wide above 1.1 m?": entries 16, 17 and 22 (item 3); `data_clean/step3b_law_grid.csv`; Figure 4d.
- "How was the stage record cleaned?": entries 3 to 6; `05_stage_simple.py`; figures 00e to 00k.
- "What does the ADCP show?": entry 23; `34_adcp_single_beam.py`; `data_clean/adcp_single_beam.json`; the open question about what 'Vb' is in `OPEN_QUESTIONS.md`.
- "Change the ecological threshold or the prior": the threshold is `THR` in `26_step7_ecology.py` and `29_paper_figures.py`; the prior is built in `35_nest_bayes.py` from `data_raw/literature_priors.json` (keep it weak and keep its provenance).
- "Run another reach": the Middle's chain needs, in this order, a stage sensor inside the surveyed stretch (or a tie to one), the RTK sections on one datum, transect profiles with three to five readings per location, and nests held aside. `27_`/`37_` prepare a reach the Middle's way; `39_` diagnoses its sensor record (clean it the way `05_stage_simple.py` cleans the Middle's before use); `28_`/`38_` run the tests and the chain. Read entries 21 and 24 first: they say why the two other reaches are paused.
- "Reproduce the paper's numbers": `RERUN.md`.

## What is pending

`OPEN_QUESTIONS.md`: the ADCP's instrument model and the meaning of 'Vb' (the paper's Figure 3b assumes a single 25° beam), the hand-held meter's make, a source for the 0.2 m/s convention, and four items for the other reaches.
