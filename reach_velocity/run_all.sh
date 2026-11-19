#!/bin/zsh
# Reproduces the Middle-reach analysis end to end from the raw copies in data_raw/. Run from anywhere: the script moves into its own folder.
# Python: ../.venv/bin/python by default (a venv at the repository root, see RERUN.md); override with PY=/path/to/python zsh run_all.sh.
# About ten minutes on a laptop. Every script is standalone once its inputs exist; the header docstring of each names its inputs and outputs.
set -e
cd "$(dirname "$0")"
PY="${PY:-../.venv/bin/python}"
# --- the inputs and their cleaning (logbook entries 1 to 13)
$PY 00_load.py                     # raw sheets and sensor files to clean tables (first pass: raw stage); checks every line file has a script (provenance.py)
$PY 44_banklines_survey.py         # the centrelines and the old bank lines (set A) from the survey (entry 29)
$PY 02_stage_diagnose.py           # the sensor's daytime pattern, diagnosed (figure 00e)
$PY 04_stage_night.py              # nights from the sun's elevation (figure 00g)
$PY 05_stage_simple.py             # the adopted stage correction: night line, monthly template, median (figures 00h, 00i)
$PY 06_align_10min.py              # stage and ADCP on one 10-min grid (figure 00j)
$PY 00_load.py                     # second pass: readings get the corrected stage at their own time
$PY 07_full_series.py              # both seasons in one view (figure 00k)
$PY 08_field_prep.py               # log profiles: z0 for the reach, one amplitude per location (figures 00l, 00m)
$PY 10_step1_adcp_gp.py            # step 1 (superseded, kept for the record): the ADCP alone as a GP of stage, cross-validated (figure 01)
$PY 11_step1b_law.py               # step 1b (superseded): laws for the ADCP's stage response (figure 01b); 12_ records its h0
$PY 12_amplitude_inputs.py         # amplitude, its noise, bed elevation on the sensor datum (figure 00n)
$PY 13_amplitude_checks.py         # profile checks: z0 by substrate, year, crew (figure 00o, figures/profiles)
$PY 14_transect_id_check.py        # transect labels checked by geometry
$PY 15_rtk_elevations.py           # RTK cross-sections tied to the sensor datum (figure 00q)
$PY 16_what_we_learn.py            # the target explained (figures 00r, 00s)
$PY 17_before_after_target.py      # before/after the stage law (figure 00t)
$PY 18_conveyance_schematic.py     # the place factor, schematic (figure 00u)
$PY 19_features.py                 # bank set C and the candidate input features (figure 00v; entry 30)
$PY 01_overview.py                 # the data at a glance (figures 00a to 00d; after 19_, which writes the bank lines it draws)
# --- the chain (entries 15 to 20, 22, 23)
$PY 20_step2_conveyance_gp.py      # step 2: the GP on the place factor, sections held out (figure 02)
$PY 21_step3_stage_transfer.py     # step 3: the day effect, laws, days held out, readings on held-out days (figure 03)
$PY 22_step3b_law_refine.py        # step 3b: the (h0, gamma) grid and the ADCP shape (figure 03b)
$PY 23_step4_uncertainty.py        # step 4: the uncertainty chain, calibration, attribution (figure 04)
$PY 24_step5_map.py                # step 5: the bed GP and the 2 m grid, map against sections (figure 05)
$PY 25_step6_rating_curves.py      # step 6: spatial rating curves, continuity, section-level coverage (figure 06)
$PY 26_step7_ecology.py            # step 7: the share below the threshold, nests (figure 07)
$PY 30_feature_ablation.py         # feature ablation with a paired test
$PY 30b_kernel_ablation.py         # kernel family and restart check (Table 4)
$PY 31_adcp_recheck.py             # the ADCP re-read with the surveyed depth
$PY 34_adcp_single_beam.py         # the ADCP as a single 25-degree beam against the chain (figure 09)
$PY 33_nest_habitat_model.py       # nest-calibrated thresholds and the selection model (figure 10)
$PY 35_nest_bayes.py               # the Bayesian suitability with the USFWS 1983 curve as a weak prior (figure 11)
$PY 36_nest_selection_explained.py # the selection model explained (figure 12)
# --- the paper's figures, tables and the methods document
$PY 41_model_variables_figure.py   # Figure M: the model's variables (figure 17)
$PY 29_paper_figures.py            # F1 to F5 in paper/
$PY 32_paper_tables.py             # Tables 1 to 4 in paper/
$PY 40_methods_docx.py             # paper/METHODS_AND_RESULTS.docx and .md
# --- paused: the other two reaches (entries 21 and 24). Run by hand once their sensor records are cleaned and the sensor placement is settled:
#   $PY 27_step8_prep.py; $PY 39_other_sensors_diagnose.py; $PY 37_other_reaches_prep.py; $PY 28_step8_chain.py; $PY 38_other_reaches_chain.py
#   (37_ builds bank set C for the Upper and Lower, so 28_ now runs after it; entry 30)
# --- the bank-line comparison as recorded (entries 27 and 28): $PY 43_bankline_comparison.py build; $PY 43_bankline_comparison.py score
