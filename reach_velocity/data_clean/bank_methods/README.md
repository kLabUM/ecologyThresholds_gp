# The bank-line methods explored: the tables kept for the comparison

Each folder holds the few tables that one earlier run of the pipeline wrote, kept so that `46_bank_methods_comparison.py` can put every method
side by side (figure 20, `comparison_table.md`). They are records, not inputs to the model; the version in use (C fixed at each reach's median
survey-day level, logbook entry 31) is the live `data_clean/`.

| folder | method | from |
|---|---|---|
| `set_A/` | A, the old lines | commit 4298715 (entry 29) |
| `C_stage_record_outline/` | C moving with stage, the map's outline at the record's highest stage (Middle ecology only) | commit 87bcb46 (entry 30) |
| `C_stage/` | C moving with stage, the map's outline at the highest survey-day stage | commit 5f876e4 (entry 30) |

Set B was scored under the entry-27 rule only (`data_clean/bankcmp/`, entry 28), not run through the whole model. `run_logs/` holds the console logs of the fixed-line run and its scoring.
