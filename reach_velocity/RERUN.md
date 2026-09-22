# Rerunning and verifying this analysis

## Environment

Python 3.11 with the packages in `requirements.txt` (pinned to the versions used: scikit-learn 1.9.0, numpy 2.4.6, scipy 1.17.1, pandas 2.3.3, geopandas, shapely, pyproj, matplotlib, python-docx). `requirements_frozen.txt` is the complete environment for an exact match. The scripts expect the interpreter at `../.venv/bin/python`, that is, a venv at the repository root:

```bash
uv venv --python 3.11 .venv
uv pip install --python .venv/bin/python -r reach_velocity/requirements.txt
```

Any other interpreter works with `PY=/path/to/python zsh reach_velocity/run_all.sh`.

## Run everything

```bash
zsh reach_velocity/run_all.sh
```

About ten minutes on a laptop. The script moves into its own folder, runs the scripts in the logbook's order, and rewrites `data_clean/`, `figures/` and `paper/`. Random seeds are fixed, so a rerun reproduces every number in `paper/TABLES.md` to the printed precision; with other package versions the third decimal of a GP score can move, which is not a concern.

## What a fresh clone reproduces

Checked on 22 September 2026: the branch was cloned into an empty directory, a new environment was built from `requirements.txt`, `run_all.sh` was run, and the five other-reach scripts were run by hand. Of the 254 committed files in this folder, 247 came back byte-identical (every table, every JSON, every PNG figure, the markdown and the tables); the other seven differ only in timestamps (the PDFs and the docx container), with their content identical. Nothing outside the folder was read.

## Verify

After a rerun, compare against the committed outputs:

1. `paper/TABLES.md`: Table 3 (scores) and Table 4 (ablations) match to the printed precision.
2. The headline numbers: section held out CRPS 0.307 against 0.414 for the reach mean (`data_clean/step2_metrics.csv`); readings on held-out days MAE 0.082 m/s against a deviation of 0.119 (`step6_validation_table.csv`); section means on held-out days inside the 90% band 20 of 20 (`step6_section_coverage.csv`); reach share below 0.2 m/s at 0.90 m 0.342 (`step5_reach_totals.csv`); mean suitability E at 0.90 m 0.33 with the velocity effect −4.6 ± 1.5 per m/s (`nest_bayes.json`); ADCP level ratio 0.91 (`adcp_single_beam.json`).
3. `paper/F1_*.png` to `F5_*.png` regenerate identically.

A prompt for a coding-agent session that verifies the work independently:

> Read reach_velocity/README.md, CLAUDE.md and LOGBOOK.md. Set up the environment as RERUN.md says, run run_all.sh, and compare every number in paper/TABLES.md with what the rerun wrote. Then pick three claims from paper/STORY.md and check each against the data_clean tables that support it, reading the script that produced them. Report any number that does not reproduce and any claim you cannot trace to a table.

## Rerun a single step

Every script is standalone once its inputs exist; the header docstring names its inputs and outputs. After changing the model code, rerun `20_` to `26_`, then `30_` to `36_`, then `41_`, `29_`, `32_` and `40_`.

## The paused scripts

`27_step8_prep.py`, `28_step8_chain.py`, `39_other_sensors_diagnose.py`, `37_other_reaches_prep.py` and `38_other_reaches_chain.py` (the other two reaches) are not in `run_all.sh`. Their inputs are in `data_raw/`; run them by hand from this folder once the sensor records are cleaned and the sensor placement is settled (logbook entries 21 and 24).
