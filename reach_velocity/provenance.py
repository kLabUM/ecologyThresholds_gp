"""Every line or polygon file in this folder and the script that writes it (logbook entries 30 and 31). 00_load.py runs check() first and stops if a
line or polygon file has no script here: no bank line, channel outline or centreline is used unless the folder can make it again.
    python provenance.py      # lists the files and their scripts"""
import fnmatch
from pathlib import Path
HERE = Path(__file__).resolve().parent
WRITERS = {   # pattern (relative to this folder) -> the script that writes it
    "data_clean/banklines_survey/bank_lines_*.geojson": "44_banklines_survey.py (set A: dataSorting.ipynb at 2346281, carried over)",
    "data_clean/banklines_survey/channel_extent_*.geojson": "44_banklines_survey.py (set A's channel polygon)",
    "data_clean/banklines_survey/centerline_union_*.geojson": "44_banklines_survey.py (the main stem from the river network)",
    "data_clean/bankcmp/banks_B_*.geojson": "43_bankline_comparison.py build (set B)",
    "data_clean/bankcmp/banks_C_*.geojson": "43_bankline_comparison.py build (set C as scored in entry 28)",
    "data_clean/banks_c/C_Middle.geojson": "19_features.py via banklines_3dep.StageFeatures.write (set C in the model, fixed at the median survey-day level, entry 31)",
    "data_clean/banks_c/C_outline_Middle.geojson": "19_features.py via banklines_3dep.StageFeatures.write",
    "data_clean/banks_c/C_Upper.geojson": "37_other_reaches_prep.py via banklines_3dep.StageFeatures.write (entry 31)",
    "data_clean/banks_c/C_Lower.geojson": "37_other_reaches_prep.py via banklines_3dep.StageFeatures.write (entry 31)",
    "data_clean/banks_c/C_outline_*.geojson": "37_other_reaches_prep.py via banklines_3dep.StageFeatures.write",
}
EXTERNAL = {  # inputs that are not lines the analysis draws: where they come from
    "data_raw/Huron_RiverSystem.geojson": "the Huron River network (stream centrelines) from the main branch; 44_ takes the main stem from it",
}
LINE_SUFFIXES = (".geojson", ".shp", ".gpkg", ".kml")


def files():
    return sorted(p.relative_to(HERE).as_posix() for p in HERE.rglob("*") if p.suffix.lower() in LINE_SUFFIXES and "__pycache__" not in p.parts)


def writer(rel):
    if rel in EXTERNAL: return "external: " + EXTERNAL[rel]
    for pat, w in WRITERS.items():
        if fnmatch.fnmatch(rel, pat): return w
    return None


def check():
    orphans = [f for f in files() if writer(f) is None]
    if orphans:
        raise SystemExit("line or polygon files with no script that writes them (see provenance.py):\n  " + "\n  ".join(orphans))
    return True


if __name__ == "__main__":
    for f in files(): print(f"{f}\n    <- {writer(f) or 'NO SCRIPT'}")
