# A season-long velocity record for a river reach from one survey campaign, one stage sensor and one ADCP

This folder is the complete analysis behind the paper: the raw inputs, every script in run order, every table and figure they produce, the experiment logbook, the reviewer reports, and the paper's figures, tables and methods document. It is self-contained and reruns from the raw files in about ten minutes with one command (`RERUN.md`). The study site is the Middle reach of the Huron River, Michigan; the other two reaches are set up but paused (see the end).

## Start here

1. `paper/METHODS_AND_RESULTS.md` (also as `.docx`): the approach, where it comes from, the model with its variables, the data and its cleaning, the results with their takeaways, the five figures and the tables, and directions for the paper. Read this first.
2. `figures/17_model_variables.png`: every variable of the model drawn from real data (a, z0, h, A, h0, γ, the GP's inputs, the chain). Keep it open while reading anything else.
3. `LOGBOOK.md`: the experiment notebook, entries 0 to 24 in the order the work was done, with a figure per step and every decision and number. Entry 22 answers three reviewer reports (`reviews/`). The table "Where we stand" compares this approach with the original notebooks.
4. `paper/FIGURES.md`, `paper/TABLES.md`, `paper/STORY.md`: the captions, the tables, the story and an abstract draft.
5. `OPEN_QUESTIONS.md`: what the analysis still needs from the field team.

## The method on one page

The question: how far can one stage record and one survey campaign be pushed toward a season-long, uncertainty-bounded velocity record for a reach, and an ecological number from it?

1. **One profile becomes one number.** Each location's three to five readings at heights z above the bed follow the log law v(z) = a ln(z/z0), with one roughness scale z0 for the reach (11 mm). The amplitude a is one number per location, with a measurement error from the readings' own spread. The depth-averaged velocity is a (ln(h/z0) − 1).
2. **A one-parameter law separates the day from the place.** On a survey day the amplitude rises with stage; a = A (stage − h0)^γ removes it. A, the place factor, does not depend on the day. Eleven survey days cannot pin (h0, γ): several pairs fit equally well, so the whole set is carried into the uncertainty and the band widens above the surveyed stage range.
3. **A Gaussian process learns A from channel geometry.** Six measured inputs (position along and across, distance to the bank, relative depth, section mean depth, lateral bed slope), substrate notes in the mean, the readings' spread as fixed per-point noise, an ARD Matérn 3/2 kernel. A held-out section is predicted better than the reach mean (CRPS 0.31 against 0.41), and bank distance alone does nearly as well: the field is at-a-station hydraulic geometry.
4. **A joint uncertainty chain carries everything through the stage record.** One draw of the whole reach at a time: the GP posterior, the bed, the stage, the roughness scale, the law from its set, and a discrepancy term calibrated on held-out days and split into a reach-level and a point-level part. The output is the depth-averaged velocity of every section every two hours over two seasons, with 90% bands.
5. **The ecological rule is learned from the fish.** The 59 nests are asked what they chose relative to what was available (a used-against-available logistic model on velocity, bank distance and depth), with the published USFWS spawning curve as a weak prior. The reach's mean nesting suitability E follows through the same chain, with its band; the share below 0.2 m/s is kept as the plain-language companion.

Validation: sections held out one at a time, days held out one at a time, the readings on held-out days at their own heights, the section means on held-out days against the band, the ADCP through a season at its own cell, and continuity of discharge between sections. Scores are MAE, CRPS and coverage in physical units; random folds are shown once, as the cautionary comparison.

| symbol | meaning | where it comes from |
|---|---|---|
| z, v(z) | height above the bed and the velocity read there | measured, 3 to 5 per location |
| z0 | the height at which the log law reaches zero | 11 mm, fitted for the reach (7 to 16 mm carried) |
| a | the log-law amplitude at a location on its day | fitted per location, with its own error |
| h, v̄ | local depth, and the depth-averaged velocity a (ln(h/z0) − 1) | stage minus bed; derived |
| stage | the sensor's water level every 10 minutes | measured, cleaned (night line, monthly template, median) |
| h0, γ | the stage law's offset and exponent | (0.45 m, 1.0) adopted; seven equivalent pairs carried |
| A | the place factor a / (stage − h0)^γ | the GP's target: one static field |
| S(v), E | relative nesting suitability exp(b_v v) and its reach mean | b_v learned from the nests with a weak prior |

## What is in the folder

| path | what it holds |
|---|---|
| `data_raw/` | the raw inputs as used: survey sheets (2025, 2026; the crew column is coded), stage-sensor and ADCP exports, centreline and bank lines for all three reaches, USGS daily discharge (reference only), the published suitability curve with its provenance (`literature_priors.json`) |
| `00_` to `19_` | the inputs, cleaned and checked step by step (logbook entries 1 to 13) |
| `20_` to `26_` | the chain: the GP on A, the stage law, the uncertainty chain, the map, the spatial rating curves, the ecological number (entries 15 to 20) |
| `30_` to `36_` | the reviewer-driven additions: ablations, the ADCP re-read and its single-beam test, the nest selection model and the Bayesian suitability (entries 22, 23) |
| `29_`, `32_`, `40_`, `41_` | the paper's figures and tables, the methods document, the variables figure; all from saved outputs, no refitting |
| `gpcore.py`, `uncertainty.py` | shared code: the GP, the target, the scores; the Monte Carlo chain and the joint reach sampler |
| `data_clean/` | every table the scripts write; the logbook's provenance section maps each to its script |
| `figures/` | every diagnostic figure (00a to 17) and the per-location profile pages in `figures/profiles/` |
| `paper/` | `F1` to `F5` (PNG and PDF), `FIGURES.md`, `TABLES.md` with `table1..4_*.csv`, `STORY.md`, `METHODS_AND_RESULTS.md` and `.docx` |
| `reviews/` | the brief and the three reviewer reports (hydraulics, statistics, story); entry 22 answers each finding |
| `27_`, `28_`, `37_`, `38_`, `39_` | the other two reaches: paused (below) |
| `run_all.sh`, `requirements.txt`, `requirements_frozen.txt`, `RERUN.md` | reproduction |
| `CLAUDE.md` | instructions for a coding agent working in this folder |

## Questions people ask

**Why not learn velocity directly from stage and position, as the notebooks did?** Because velocity at a point depends on the day (through stage) and on the place, and eleven days of surveys at different places cannot separate the two by brute force. The log law turns a profile into one number, the stage law removes the day, and what is left is a property of the place that 257 locations can pin. Entry 11 and 12 of the logbook, and Section 2 of the methods document.

**Where does z0 = 11 mm come from?** A pooled grid search over every profile with three or more readings; the fit is flat between 7 and 16 mm, and the chain draws from that range. Entry 7; `08_field_prep.py`; figure 00m.

**Why hold out whole sections and whole days rather than random points?** Random folds leave a point's neighbours in the training set, so the score reflects interpolation over a few metres, not prediction of a new place or a new day. On this data random folds flatter every model (Figure 2c). Entry 15.

**What is the 'law set' and why does the band widen above 1.1 m?** The day-out score is flat across (h0, γ) pairs from (0.35, 1.25) to (0.65, 0.5). Inside the surveyed range they agree; above it they diverge. The chain draws from the set, so the band says what the data cannot decide. Entries 16, 17, 22; figure 03b; Figure 4d.

**What does the ADCP add?** It was never used in any fit. Read as a single 25° beam, its record through 2025 matches the chain at its own cell (ratio 0.87, 37 mm/s, the band covers every two-hour mean to 1.20 m). That reading is to be confirmed with the instrument's configuration (`OPEN_QUESTIONS.md`). Entry 23; figure 09; Figure 3b.

**Where does the ecological rule come from?** From the nests: a used-against-available model, with the USFWS 1983 spawning curve as a weak prior. The nests sit in deep, slow water near the bank; the velocity effect is −3.2 ± 1.8 log-odds per m/s, half use at 0.21 m/s (weaker with bank set C, entry 30). The 0.2 m/s convention sits at the 60th percentile of the nests. Entry 23; figures 10 to 12; Figure 5.

**What are the numbers to quote?** `paper/TABLES.md`, Table 3, and the "Status" section at the end of the logbook.

**What is the RTK tie?** The RTK water surface on each section minus the sensor stage at the time gives one offset per session and one slope for the reach (5 cm per 100 m, residual 1.3 cm). It puts every bed on the sensor's datum and shows that one sensor describes the whole reach. Entry 10; figure 00q.

## Conventions

Units are converted once at load; nothing downstream sees feet. Random seeds are fixed. Numbers come from the surveyed sections; the 2 m map is illustration. Crew identities in the survey sheets are coded (`crew01`, ...); no personal names appear anywhere in the folder. Every step has a script, a figure and a logbook entry.

## The other two reaches

Their survey data, RTK, sensor records and geometry are in `data_raw/` and the scripts exist (`27_step8_prep.py`, `28_step8_chain.py`, `39_other_sensors_diagnose.py`, `37_other_reaches_prep.py`, `38_other_reaches_chain.py`; figures 08, 08a, 13 to 16). The work is paused: the Lower's sensor record needs the Middle's cleaning (a daytime artefact of 2 to 9 cm and spikes), the Upper's sensor sits on another branch of the network far from its sections, and both sensors sit outside the surveyed stretches. Entries 21 and 24 of the logbook and `OPEN_QUESTIONS.md` say what to do next. Nothing in the Middle-reach results depends on them.
