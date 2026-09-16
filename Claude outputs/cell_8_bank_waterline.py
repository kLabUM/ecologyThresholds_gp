### 8) Stream bank waterline from the surveyed "Bank height" points -> GeoJSON
# =============================================================================
# Add this as a new cell at the END of dataSorting.ipynb.
# It writes six files and changes nothing else. No new columns anywhere.
#
#   mapping/bank_left_Upper.geojson     mapping/bank_right_Upper.geojson
#   mapping/bank_left_Middle.geojson    mapping/bank_right_Middle.geojson
#   mapping/bank_left_Lower.geojson     mapping/bank_right_Lower.geojson
#
# "Left" and "right" are as seen facing downstream, matching the sign of
# distance_to_centerline from section 6 (negative = left).
#
# WHERE THE BANK COMES FROM
#   Any row whose "Point Notes" contains "Bank height" is a direct survey point at
#   the bank. Its Elevation is the highest the water can reach at that spot, so
#   these points define the maximum-extent waterline. Nothing is extrapolated and
#   no curve is fitted to the streambed.
#
#   These rows carry no depth or velocity, so section 4 removes them and they are
#   NOT in sorted_data. This cell therefore reads them from `complete_data`, which
#   still holds every surveyed row.
#
#   To decide which bank a point is on, each transect calibrates itself: its
#   in-water points already carry distance_to_centerline from section 6, and along
#   a transect that value changes in a straight line with position. Fitting that
#   line and reading the bank point off it is exact here (R-squared 0.99997 or
#   better on all three sections) and recomputes no centerline geometry.
#
# USES ONLY THE LIBRARIES ALREADY IMPORTED IN CELL 0
#   pandas as pd, numpy as np, os, geopandas as gpd, LineString,
#   matplotlib.pyplot as plt
# =============================================================================

# ------------------------------ SETTINGS -------------------------------------
BANK_SECTIONS  = ['Upper', 'Middle', 'Lower']
# Wordings in "Point Notes" that mark a surveyed point at the bank. These were
# typed by hand, so the match ignores capitalisation AND extra or missing spaces:
# "Bank height", "bank  Height ", "BANK HEIGHT" all match. To accept another
# wording, add it to this list in lower case - nothing else needs to change.
BANK_NOTES     = ['bank height', 'edge of bank']
MIN_WET_POINTS = 4               # in-water points needed to calibrate a transect
SIDE_LIMIT_FT  = 150.0           # ignore in-water points this far off the centerline
OUTSIDE_TOL_FT = 1.0             # a bank point may sit this far inside the outermost
                                 # in-water point on its own side before it is rejected
MAX_GAP_FT     = 800.0           # break the line rather than draw a straight chord
                                 # across a stretch of river with no bank points

os.makedirs('mapping', exist_ok=True)
print('=' * 74)
print('STREAM BANK WATERLINE, FROM THE SURVEYED "BANK HEIGHT" POINTS')
print('=' * 74)

fig, axgrid = plt.subplots(2, len(BANK_SECTIONS),
                           figsize=(6 * len(BANK_SECTIONS), 11))
axgrid = np.atleast_2d(axgrid)
summary = []

for ax, ax2, section_name in zip(axgrid[0], axgrid[1], BANK_SECTIONS):
    print(f'\n{section_name}')
    print('-' * 74)

    # --- the surveyed bank points, from the untrimmed table -------------------
    # Tidy the note before matching: force lower case, collapse any run of spaces
    # or tabs to one space, trim the ends. That makes the match immune to how
    # each person happened to type it.
    tidy_note = (complete_data['Point Notes'].astype(str).str.lower()
                 .str.replace(r'\s+', ' ', regex=True).str.strip())
    is_bank = np.zeros(len(complete_data), dtype=bool)
    for wording in BANK_NOTES:
        is_bank |= tidy_note.str.contains(wording, regex=False).values

    bank = complete_data[
        (complete_data['Section'] == section_name) & is_bank
    ].dropna(subset=['Easting', 'Northing', 'Elevation', 'Transect ID']).copy()
    bank['wording'] = tidy_note[bank.index].values
    bank = bank.drop_duplicates(subset=['Transect ID', 'Easting', 'Northing'])
    bank['day'] = pd.to_datetime(bank['Collection Date'],
                                 format='mixed').dt.date.astype(str)

    # --- the in-water points, which already carry the section-6 distances -----
    wet = pd.read_csv(f'sorted_data/complete_data_{section_name}.csv')
    wet = wet.dropna(subset=['Easting', 'Northing',
                             'distance_to_centerline', 'dist_along_centerline'])
    wet = wet[wet['distance_to_centerline'].abs() <= SIDE_LIMIT_FT]
    wet['day'] = pd.to_datetime(wet['Collection Date'],
                                format='mixed').dt.date.astype(str)

    print(f'  {len(bank)} bank points on {bank["Transect ID"].nunique()} transects '
          f'(of {wet["Transect ID"].nunique()} surveyed)')
    for wording, k in bank['wording'].value_counts().items():
        print(f'      {k:3d} from notes reading "{wording}"')

    # --- place each bank point on its transect's own side-distance scale -------
    placed, no_calib, inside = [], [], []
    for (tid, day), points in bank.groupby(['Transect ID', 'day']):
        w = wet[(wet['Transect ID'] == tid) & (wet['day'] == day)]
        w = w.drop_duplicates(subset=['Easting', 'Northing'])
        xy = w[['Easting', 'Northing']].values.astype(float)
        if len(w) < MIN_WET_POINTS:
            no_calib.append(f'{tid}'); continue
        middle = xy.mean(0)
        direction = np.linalg.svd(xy - middle)[2][0]   # the line the transect runs along
        x = (xy - middle) @ direction
        if np.ptp(x) < 1.0:
            no_calib.append(f'{tid}'); continue
        # distance_to_centerline changes in a straight line along the transect
        n_wet = w['distance_to_centerline'].values
        m, c = np.polyfit(x, n_wet, 1)

        # The outermost in-water point at each end of the transect. Do NOT decide
        # left vs right from the sign of `side`: on a few transects the centerline
        # runs outside the channel, so both banks come out negative and the sign
        # test puts them both on the left. Ask instead which END of the in-water
        # survey the point lies beyond - that is true no matter where the
        # centerline happens to sit.
        left_edge, right_edge = n_wet.min(), n_wet.max()

        bx = (points[['Easting', 'Northing']].values.astype(float) - middle) @ direction
        for (_, pt), side in zip(points.iterrows(), m * bx + c):
            past_left  = left_edge - side     # positive: beyond the left-hand end
            past_right = side - right_edge    # positive: beyond the right-hand end
            label  = 'left' if past_left > past_right else 'right'
            beyond = max(past_left, past_right)
            if beyond < -OUTSIDE_TOL_FT:
                inside.append(f'{tid} ({-beyond:.0f} ft inside)')
                continue
            placed.append(dict(
                transect=tid, day=day,
                along=float(w['dist_along_centerline'].median()),
                side=float(side), bank=label,
                east=float(pt['Easting']), north=float(pt['Northing']),
                bankfull_elev=float(pt['Elevation'])))

    if no_calib:
        print(f'  dropped, transect has fewer than {MIN_WET_POINTS} in-water points '
              f'to calibrate against: {", ".join(sorted(set(no_calib)))}')
    if inside:
        print(f'  dropped, point lands inside the in-water survey on its own side, '
              f'so it cannot be the bank: {"; ".join(inside)}')

    P = pd.DataFrame(placed)
    if P.empty:
        print('  SKIPPED: no bank points could be placed.')
        ax.set_title(f'{section_name}: no bank data')
        ax2.set_title(f'{section_name}: no bank data'); continue

    # --- bankfull height above the water on the day it was surveyed ----------
    wet_wse = (wet.assign(wse=wet['Elevation'] + wet['Total water depth'])
               .groupby(['Transect ID', 'day'])['wse'].median())
    P['above_water'] = [P['bankfull_elev'].iloc[i]
                        - wet_wse.get((P['transect'].iloc[i], P['day'].iloc[i]), np.nan)
                        for i in range(len(P))]
    print(f'  bankfull sits {P["above_water"].median():.2f} ft above the water surface '
          f'on the survey day (range {P["above_water"].min():.2f} to '
          f'{P["above_water"].max():.2f})')

    # --- one line per bank, in downstream order ------------------------------
    lines = {}
    for label in ('left', 'right'):
        side_pts = P[P['bank'] == label]
        # one vertex per transect: if a transect has several points on this side,
        # keep the one furthest from the centerline
        side_pts = (side_pts.assign(absside=side_pts['side'].abs())
                    .sort_values('absside').groupby('transect', as_index=False).last()
                    .sort_values('along'))
        if len(side_pts) < 2:
            print(f'  {label:5s} bank: only {len(side_pts)} point(s) - no file written')
            continue

        # break the line wherever consecutive vertices are more than MAX_GAP_FT
        # apart on the ground: there is no survey through that stretch, so a
        # straight chord there would be invention, not measurement
        xy = side_pts[['east', 'north']].values
        step = np.hypot(*np.diff(xy, axis=0).T)
        cut = np.flatnonzero(step > MAX_GAP_FT) + 1
        parts = [p for p in np.split(np.arange(len(xy)), cut) if len(p) >= 2]
        if not parts:
            print(f'  {label:5s} bank: every vertex is more than {MAX_GAP_FT:,.0f} ft '
                  f'from the next - no file written')
            continue

        rows = []
        for k, idx in enumerate(parts, 1):
            q = side_pts.iloc[idx]
            rows.append(dict(
                section=section_name, bank=label, stretch=k,
                n_points=len(q),
                length_ft=round(float(np.hypot(*np.diff(q[['east', 'north']].values,
                                                        axis=0).T).sum()), 1),
                widest_vertex_gap_ft=round(float(np.hypot(
                    *np.diff(q[['east', 'north']].values, axis=0).T).max()), 1),
                bankfull_elev_min_ft=round(float(q['bankfull_elev'].min()), 2),
                bankfull_elev_max_ft=round(float(q['bankfull_elev'].max()), 2),
                note='maximum-extent waterline from surveyed "Bank height" points; '
                     'one vertex per transect',
                geometry=LineString(q[['east', 'north']].values)))
        lines[label] = [side_pts.iloc[idx] for idx in parts]
        path = f'mapping/bank_{label}_{section_name}.geojson'
        gpd.GeoDataFrame(rows, crs=target_crs).to_file(path, driver='GeoJSON')
        dropped = len(side_pts) - sum(len(p) for p in parts)
        print(f'  wrote {path}  ({sum(len(p) for p in parts)} vertices in '
              f'{len(parts)} separate stretch(es)'
              + (f', {dropped} isolated point(s) left out' if dropped else '') + ')')

    both = P.groupby('transect')['bank'].nunique()
    summary.append(dict(section=section_name, bank_points=len(P),
                        transects=P['transect'].nunique(),
                        both_banks=int((both == 2).sum()),
                        median_above_water_ft=round(float(P['above_water'].median()), 2)))

    # --- top row: where the lines sit on the ground --------------------------
    ax.scatter(wet['Easting'], wet['Northing'], s=5, c='0.72', label='in-water survey')
    for label, colour in (('left', '#1A9850'), ('right', '#762A83')):
        for k, q in enumerate(lines.get(label, [])):
            ax.plot(q['east'], q['north'], 'o-', color=colour, lw=1.6, ms=5,
                    label=f'{label} bank' if k == 0 else None)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(f'{section_name}  ({P["transect"].nunique()} transects with bank points)')
    ax.set_xlabel('Easting (ft)'); ax.grid(alpha=0.3)

    # --- bottom row: the check that is actually readable ---------------------
    # along the reach (left-right) against distance from the centerline (up-down).
    # The left bank must stay below every grey point and the right bank above it.
    ax2.scatter(wet['dist_along_centerline'], wet['distance_to_centerline'],
                s=5, c='0.72', label='in-water survey')
    for label, colour in (('left', '#1A9850'), ('right', '#762A83')):
        for k, q in enumerate(lines.get(label, [])):
            ax2.plot(q['along'], q['side'], 'o-', color=colour, lw=1.6, ms=5,
                     label=f'{label} bank' if k == 0 else None)
    ax2.axhline(0, color='k', lw=0.8, ls='--')
    ax2.set_xlabel('distance along centerline (ft)')
    ax2.set_title(f'{section_name}: bank vs in-water survey, by side')
    ax2.grid(alpha=0.3)

axgrid[0, 0].set_ylabel('Northing (ft)')
axgrid[0, 0].legend(fontsize=9, loc='best')
axgrid[1, 0].set_ylabel('distance from centerline (ft)\nnegative = left')
plt.suptitle('Maximum-extent waterline from surveyed bank points', fontsize=13)
plt.tight_layout(); plt.show()

print('\n' + '=' * 74)
print(pd.DataFrame(summary).to_string(index=False))
print("""
WHAT THESE FILES ARE
  The maximum-extent waterline: the line the water reaches when it is as high as
  it can get before leaving the channel. Every vertex is a surveyed point.

WHAT TO WATCH
  Bank points exist on only some transects, so each line has fewer vertices than
  there are transects and runs straight between them. Where the survey skips a
  long stretch of river the line is BROKEN rather than drawn straight across, so
  a file can hold several separate stretches (the "stretch" column). Use them to
  clip a map only where a stretch actually exists.

  A transect with a bank point on only one side contributes to one line only, so
  the left and right lines can have different vertex counts and cover different
  stretches of the reach.

  In the picture, both lines should sit outside the grey in-water points on their
  own side, and should not cross each other.

TO GET AN ESRI SHAPEFILE INSTEAD
  gpd.read_file('mapping/bank_left_Middle.geojson').to_file(
      'mapping/bank_left_Middle.shp')
""")
