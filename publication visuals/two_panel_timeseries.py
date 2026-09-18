"""
Four-panel time series figure: 2026 daily depth at each site (Upper, Middle,
Lower) plus 2025 daily ADCP velocity, one panel per row, top to bottom in
that order. Tall-and-narrow layout and enlarged text, since this figure is
meant to be viewed small (e.g. embedded at reduced size in a document).

Depth panels (Upper, Middle, Lower)
------------------------------------
Relative water depth (m) for 2026, from
data/26_depthSensorData_{Site}.csv. Raw `_value` is a distance in mm from
the sensor down to the water surface; converting to depth needs a per-site
datum offset:

    depth_m = (site_offset_mm - _value_mm) / 1000

Offsets are hard-coded per site rather than read from the file's own
`offset` column, per claude/GP_code_review_2026-09-15.md secs 1.2-1.3
(Lower's file offset column changes mid-year; the corrected constant is
used for the whole series):
    Upper  (ARB017): 1587.5 mm
    Middle (ARB012): 7313.2 mm
    Lower  (ARB112): 6705.6 mm

Each panel's series is aggregated to one value per day using the same
method the modeling pipeline (dataSorting.ipynb, cells 2a-2c) applies
before fitting: a 300-minute (5-hour) median resample, a global 2-standard-
deviation outlier trim on that resampled series, then a per-day mean.

Depth smoothing (plot only)
----------------------------
After the daily aggregation above, each site's depth series gets a
6-point centered moving median (`smooth_daily`) -- a further despike on
top of what the modeling pipeline does, applied only for this figure. It
is not part of dataSorting.ipynb and does not touch the ADCP series.

ADCP panel
----------
Sontek-brand ADCP velocity, 2025, from data/25_ADCPData.csv, filtered to
`site_id == 'ARB111'` (its instrument tag in that file -- physically in
the Middle reach). ARB111 is not shown on the plot itself: it's an
internal site code, not a name a reader outside this dataset would
recognize, so it stays in the filter/comments but out of the title and
legend. `_value` is already m/s; the x1/0.3048 conversion in
dataSorting.ipynb is only there to match the ft/s model pipeline and is
not applied here. Same 300-min-median / 2-std-trim / daily-mean
aggregation as the depth panels. Restricted to the calendar months
spanned by the depth panels' 2026 daily series (month-of-year match,
since the two years' field seasons don't run the same calendar dates).
Labeled "ADCP" rather than "Sontek" (Sontek is the instrument brand;
ADCP -- acoustic Doppler current profiler -- is what it is) on the plot
itself; internal names still say "sontek" since that is what the source
files and pipeline call it.

Survey day markers
-------------------
Each panel gets a dashed vertical line for the calendar dates that
section was actually surveyed, taken from sorted_data/complete_data_
{Section}.csv (the dataSorting.ipynb output, which carries a corrected
`Section` label -- unlike the raw ecohydro CSVs, where a large share of
rows have no `Section` at all; see GP_code_review sec 5.4). Sontek/ADCP
rows (`Type == 'Sontek'`) are excluded since those are continuous
instrument readings, not field survey visits. Upper and Lower only have
2026 dates; Middle's sorted file carries both years, so it's filtered to
the year that panel (or the ADCP panel, at 2025) actually shows.
"""

from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "data"
SORTED_DATA_DIR = SCRIPT_DIR.parent / "sorted_data"

SITE_OFFSET_MM = {
    "Upper": 1587.5,
    "Middle": 7313.2,
    "Lower": 6705.6,
}
DEPTH_FILES = {
    "Upper": "26_depthSensorData_Upper.csv",
    "Middle": "26_depthSensorData_Middle.csv",
    "Lower": "26_depthSensorData_Lower.csv",
}
SONTEK_FILE = "25_ADCPData.csv"
SONTEK_SITE_ID = "ARB111"

SORTED_DATA_FILES = {
    "Upper": "complete_data_Upper.csv",
    "Middle": "complete_data_Middle.csv",
    "Lower": "complete_data_Lower.csv",
}

RESAMPLE_BIN = "300min"   # matches dataSorting.ipynb depth_interval_minutes / sontek_interval_minutes
OUTLIER_STD = 2.0         # matches dataSorting.ipynb outlier_std_threshold
DEPTH_SMOOTH_WINDOW = 6   # centered moving-median window, plot only

# dataviz skill categorical palette: red (slot 8) for all depth/stage
# curves, green (slot 6) for ADCP -- one color per PANEL now that each
# site has its own subplot, rather than one color per site on a shared axis
DEPTH_COLOR = "#e34948"
ADCP_COLOR = "#008300"

INK = "#000000"          # pure black text, for contrast against pure white
SECONDARY_INK = "#000000"
MUTED = "#000000"        # axis ticks/labels also go black per the contrast ask
GRID = "#d9d9d9"
SURFACE = "#ffffff"      # pure white background
BOX_COLOR = "#000000"    # outline box around each panel's plot area

SURVEY_LINE_COLOR = MUTED
SURVEY_LINE_STYLE = "--"
SURVEY_LINE_WIDTH = 1.0

# Text sized up relative to line weight: this figure is meant to be viewed
# small (e.g. shrunk into a document), so labels need to survive that
# downscaling while the data lines themselves stay thin and fine-grained.
TITLE_FONTSIZE = 15
LABEL_FONTSIZE = 13
TICK_FONTSIZE = 11
LEGEND_FONTSIZE = 11
LINE_WIDTH = 1.3
MARKER_SIZE = 2.5  # unused now that panels are line-only; kept for easy reversion


def _resample_and_trim(raw_value: pd.Series) -> pd.Series:
    """300-min median resample, then a global 2-std outlier trim. Matches
    dataSorting.ipynb's depth_agg / sontek_agg steps exactly."""
    agg = raw_value.resample(RESAMPLE_BIN).median().dropna()
    mean, std = agg.mean(), agg.std()
    return agg[(agg - mean).abs() <= OUTLIER_STD * std]


def load_depth_daily(site: str) -> pd.Series:
    df = pd.read_csv(DATA_DIR / DEPTH_FILES[site], skiprows=3)
    df["_time"] = pd.to_datetime(df["_time"])
    df = df.set_index("_time")

    trimmed_mm = _resample_and_trim(df["_value"])
    depth_m = (SITE_OFFSET_MM[site] - trimmed_mm) / 1000.0

    daily = depth_m.groupby(depth_m.index.date).mean().sort_index()
    daily.index = pd.to_datetime(daily.index)
    print(f"  {site}: {len(daily)} daily values "
          f"({daily.index.min().date()} to {daily.index.max().date()})")
    return daily


def smooth_daily(daily: pd.Series, window: int = DEPTH_SMOOTH_WINDOW) -> pd.Series:
    """Centered `window`-point moving median, applied to the plotted depth
    series only -- a further despike on top of the modeling pipeline's own
    daily aggregation, not a replacement for it."""
    return daily.rolling(window=window, center=True, min_periods=1).median()


def load_sontek_daily(months: range) -> pd.Series:
    df = pd.read_csv(DATA_DIR / SONTEK_FILE, skiprows=3)
    df = df[df["site_id"] == SONTEK_SITE_ID]
    df["_time"] = pd.to_datetime(df["_time"])
    df = df.set_index("_time")

    trimmed = _resample_and_trim(df["_value"])

    daily = trimmed.groupby(trimmed.index.date).mean().sort_index()
    daily.index = pd.to_datetime(daily.index)
    daily = daily[daily.index.month.isin(months)]
    print(f"  ADCP: {len(daily)} daily values for months {sorted(months)}")
    return daily


def load_section_survey_dates(section: str, year: int) -> list:
    """Unique field-survey calendar dates for one section/year, from the
    dataSorting.ipynb output in sorted_data/ (which has a corrected,
    complete `Section` label, unlike the raw ecohydro CSVs). Excludes
    Sontek/ADCP rows -- those are continuous instrument readings, not
    field survey visits."""
    df = pd.read_csv(SORTED_DATA_DIR / SORTED_DATA_FILES[section],
                     usecols=["Collection Date", "Type"])
    dt = pd.to_datetime(df["Collection Date"], format="mixed", utc=True,
                        errors="coerce").dt.tz_localize(None)
    mask = (df["Type"] != "Sontek") & (dt.dt.year == year)
    dates = sorted(pd.to_datetime(dt[mask].dt.normalize().dropna().unique()))
    print(f"  {section} {year}: {len(dates)} unique survey dates")
    return dates


def add_survey_lines(ax, dates):
    for d in dates:
        ax.axvline(d, color=SURVEY_LINE_COLOR, linestyle=SURVEY_LINE_STYLE,
                   linewidth=SURVEY_LINE_WIDTH, alpha=0.7, zorder=1)


def style_axis(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    # full outline box around the plot area, all four spines
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color(BOX_COLOR)
        ax.spines[side].set_linewidth(1.1)
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax.tick_params(colors=MUTED, labelsize=TICK_FONTSIZE)
    for label in ax.get_xticklabels():
        label.set_rotation(30)
        label.set_ha("right")
    ax.yaxis.label.set_color(SECONDARY_INK)


def add_legend(ax, survey_proxy):
    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(handles=handles + [survey_proxy], labels=labels + ["Survey day"],
                       frameon=False, loc="upper right", fontsize=LEGEND_FONTSIZE)
    for text in legend.get_texts():
        text.set_color(INK)


def main():
    print("Loading and daily-aggregating 2026 depth data...")
    depth_daily = {site: smooth_daily(load_depth_daily(site)) for site in DEPTH_FILES}
    print(f"  applied a {DEPTH_SMOOTH_WINDOW}-point centered moving median to each site's daily series")

    months_present = sorted(set(pd.concat([pd.Series(d.index) for d in depth_daily.values()]).dt.month))
    month_range = range(min(months_present), max(months_present) + 1)

    print("Loading and daily-aggregating 2025 ADCP data...")
    sontek_daily = load_sontek_daily(month_range)

    print("Loading per-section survey days from sorted_data...")
    survey_dates = {site: load_section_survey_dates(site, 2026) for site in DEPTH_FILES}
    survey_dates_adcp = load_section_survey_dates("Middle", 2025)

    fig, (ax_upper, ax_middle, ax_lower, ax_adcp) = plt.subplots(
        4, 1, figsize=(6.5, 14), facecolor=SURFACE, constrained_layout=True
    )
    survey_proxy = Line2D([0], [0], color=SURVEY_LINE_COLOR, linestyle=SURVEY_LINE_STYLE,
                          linewidth=SURVEY_LINE_WIDTH, label="Survey day")

    depth_axes = {"Upper": ax_upper, "Middle": ax_middle, "Lower": ax_lower}
    for site, ax in depth_axes.items():
        d = depth_daily[site]
        ax.plot(d.index, d.values, color=DEPTH_COLOR, linewidth=LINE_WIDTH,
                label=site, zorder=2)
        add_survey_lines(ax, survey_dates[site])
        ax.set_title(f"Water Stage - {site} Reach", color=INK,
                     fontsize=TITLE_FONTSIZE, loc="left")
        ax.set_ylabel("Depth (m)", fontsize=LABEL_FONTSIZE)
        add_legend(ax, survey_proxy)

    ax_adcp.plot(sontek_daily.index, sontek_daily.values, color=ADCP_COLOR,
                linewidth=LINE_WIDTH, label="ADCP velocity", zorder=2)
    add_survey_lines(ax_adcp, survey_dates_adcp)
    lo_name = pd.Timestamp(2025, month_range.start, 1).strftime("%b")
    hi_name = pd.Timestamp(2025, month_range.stop - 1, 1).strftime("%b")
    ax_adcp.set_title(
        f"ADCP Velocity - Middle Reach", color=INK,
        fontsize=TITLE_FONTSIZE, loc="left"
    )
    ax_adcp.set_ylabel("Velocity (m/s)", fontsize=LABEL_FONTSIZE)
    add_legend(ax_adcp, survey_proxy)

    for ax in (ax_upper, ax_middle, ax_lower, ax_adcp):
        style_axis(ax)

    out_path = SCRIPT_DIR / "two_panel_timeseries.png"
    fig.savefig(out_path, dpi=200, facecolor=SURFACE)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
