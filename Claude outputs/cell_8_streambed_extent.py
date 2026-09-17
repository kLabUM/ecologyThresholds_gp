### 8) Streambed extent: surveyed bank points -> continuous waterline -> polygon
# =============================================================================
# Add this as a new cell at the END of dataSorting.ipynb. It REPLACES the two
# earlier cells (the bank-waterline cell and the complete-bank-lines cell).
# It adds no columns to any table.
#
# WRITES, per section:
#   mapping/channel_extent_{Section}.geojson   the polygon - this is the output
#   mapping/bank_lines_{Section}.geojson       the two edges, for checking
#
# NEW IMPORT (everything else is already in cell 0)
#   from shapely.geometry import Polygon
#   from reachframe import to_map_frame            # project root
#
# ---------------------------------------------------------------------------
# HOW IT WORKS, IN FOUR STEPS
#
#   1. FIND THE BANK. Rows whose "Point Notes" name the bank are direct survey
#      points at the water's maximum extent. They carry no depth or velocity, so
#      section 4 drops them - this cell reads them from `complete_data`.
#
#   2. PUT EACH ONE ON THE REACH'S OWN SCALE. Along a transect,
#      distance_to_centerline changes in a straight line with position. Fitting
#      that line from the transect's in-water points and reading the bank point
#      off it gives its offset exactly (R-squared 0.99997 or better), using the
#      columns section 6 already wrote and recomputing no geometry.
#
#   3. TURN EACH BANK INTO A CURVE ALONG THE REACH. A bank described as "offset
#      from the centerline at each distance downstream" is a 1-D series, so the
#      gaps between transects are filled by interpolation rather than geometry.
#      Where no bank point was surveyed, the curve follows the outermost survey
#      point pushed outward by a margin - and that margin is measured, as the
#      median gap between the two where both exist.
#
#   4. CLOSE THE TWO CURVES INTO A POLYGON and check that every survey point
#      falls inside it. Any that do not push the nearby edge outward, and the
#      polygon is rebuilt until none are left.
#
# WHY IT IS ONE CELL
#   Steps 2 and 3 both work in (distance downstream, offset from centerline).
#   Writing the bank out as map coordinates in between and reading it back would
#   mean recovering those numbers through the curvilinear frame, which loses a
#   few feet on Upper and Lower and about fifteen on Middle. Keeping it in one
#   cell keeps the exact values.
# =============================================================================

from shapely.geometry import Polygon
from reachframe import to_map_frame

# ------------------------------ SETTINGS -------------------------------------
SECTIONS       = ['Upper', 'Middle', 'Lower']

# Wordings in "Point Notes" that mark a surveyed point at the bank. Typed by
# hand, so the match ignores capitalisation and any extra or missing spaces.
# Add a wording here in lower case; nothing else needs to change.
BANK_NOTES     = ['bank height', 'edge of bank']

MIN_WET_POINTS = 4       # in-water points a transect needs before it can place a bank point
CALIB_LIMIT_FT = 150.0   # in-water points further off the centerline than this are
                         # not used to calibrate a transect (they distort the fit)
OUTSIDE_TOL_FT = 1.0     # a bank point may sit this far inside the in-water survey
                         # on its own side before it is rejected as a bad note

GRID_FT        = 50.0    # spacing of the stations the curve is defined at
BIN_FT         = 100.0   # length of reach each survey-edge sample summarises
WINDOW_FT      = 200.0   # a station counts as surveyed if data lies this close
BANK_REACH_FT  = 400.0   # how far along the reach a surveyed bank point is trusted
SMOOTH_FT      = 300.0   # moving-average length applied to the curve
SAFETY_FT      = 2.0     # keep the edge at least this far outside the survey
MIN_WIDTH_FT   = 10.0    # never let the two edges meet
MAX_GAP_FT     = 800.0   # break the polygon across a stretch with no data at all
MAX_PASSES     = 12      # attempts to pull the edge outside every point on the ground

os.makedirs('mapping', exist_ok=True)
print('=' * 74)
print('STREAMBED EXTENT AT MAXIMUM STAGE')
print('=' * 74)

fig, axgrid = plt.subplots(2, len(SECTIONS), figsize=(6 * len(SECTIONS), 10))
axgrid = np.atleast_2d(axgrid)
extent_summary = []

for ax, ax2, section_name in zip(axgrid[0], axgrid[1], SECTIONS):
    print(f'\n{section_name}')
    print('-' * 74)

    centerline = gpd.read_file(
        f'mapping/centerline_union_{section_name}.geojson').geometry.iloc[0]

    # =========================================================================
    # 1. the data, and the surveyed bank points
    # =========================================================================
    wet = pd.read_csv(f'sorted_data/complete_data_{section_name}.csv')
    wet = wet.dropna(subset=['Easting', 'Northing',
                             'dist_along_centerline', 'distance_to_centerline'])
    wet['day'] = pd.to_datetime(wet['Collection Date'], format='mixed').dt.date.astype(str)
    s_water = wet['dist_along_centerline'].values.astype(float)
    n_water = wet['distance_to_centerline'].values.astype(float)

    tidy_note = (complete_data['Point Notes'].astype(str).str.lower()
                 .str.replace(r'\s+', ' ', regex=True).str.strip())
    is_bank = np.zeros(len(complete_data), dtype=bool)
    for wording in BANK_NOTES:
        is_bank |= tidy_note.str.contains(wording, regex=False).values

    bank = complete_data[
        (complete_data['Section'] == section_name) & is_bank
    ].dropna(subset=['Easting', 'Northing', 'Elevation', 'Transect ID']).copy()
    bank = bank.drop_duplicates(subset=['Transect ID', 'Easting', 'Northing'])
    bank['day'] = pd.to_datetime(bank['Collection Date'], format='mixed').dt.date.astype(str)
    print(f'  {len(bank)} surveyed bank points on {bank["Transect ID"].nunique()} of '
          f'{wet["Transect ID"].nunique()} transects')

    # =========================================================================
    # 2. place each bank point on the reach's scale
    # =========================================================================
    placed, no_calib, rejected = [], [], []
    calib_pool = wet[wet['distance_to_centerline'].abs() <= CALIB_LIMIT_FT]

    for (tid, day), points in bank.groupby(['Transect ID', 'day']):
        t = calib_pool[(calib_pool['Transect ID'] == tid) & (calib_pool['day'] == day)]
        t = t.drop_duplicates(subset=['Easting', 'Northing'])
        if len(t) < MIN_WET_POINTS:
            no_calib.append(str(tid)); continue
        xy = t[['Easting', 'Northing']].values.astype(float)
        middle = xy.mean(0)
        direction = np.linalg.svd(xy - middle)[2][0]     # the line the transect runs along
        along_transect = (xy - middle) @ direction
        if np.ptp(along_transect) < 1.0:
            no_calib.append(str(tid)); continue
        n_t = t['distance_to_centerline'].values
        slope, intercept = np.polyfit(along_transect, n_t, 1)

        # The outermost in-water point at each END of the transect. Do not use the
        # sign of the offset to tell left from right: on a few transects the
        # centerline runs outside the channel and both banks come out negative.
        left_end, right_end = n_t.min(), n_t.max()

        proj = (points[['Easting', 'Northing']].values.astype(float) - middle) @ direction
        for (_, pt), offset in zip(points.iterrows(), slope * proj + intercept):
            past_left, past_right = left_end - offset, offset - right_end
            side = 'left' if past_left > past_right else 'right'
            if max(past_left, past_right) < -OUTSIDE_TOL_FT:
                rejected.append(f'{tid} ({-max(past_left, past_right):.0f} ft inside)')
                continue
            placed.append(dict(side=side,
                               station=float(t['dist_along_centerline'].median()),
                               offset=float(offset),
                               elev=float(pt['Elevation'])))

    if no_calib:
        print(f'  dropped, transect has fewer than {MIN_WET_POINTS} in-water points to '
              f'calibrate against: {", ".join(sorted(set(no_calib)))}')
    if rejected:
        print(f'  dropped, point lands inside the in-water survey on its own side: '
              f'{"; ".join(rejected)}')

    P = pd.DataFrame(placed)
    if P.empty:
        print('  SKIPPED: no bank points could be placed.')
        ax.set_title(f'{section_name}: no bank data')
        ax2.set_title(f'{section_name}: no bank data'); continue
    print(f'  {len(P)} placed; bank elevation {P["elev"].min():.1f} to '
          f'{P["elev"].max():.1f} ft')

    # =========================================================================
    # 3. turn each bank into a curve along the reach
    # =========================================================================
    grid = np.arange(s_water.min(), s_water.max() + GRID_FT, GRID_FT)

    # The edge of the survey. Collect the data into bins along the reach, take the
    # outermost point on each side in each bin, and draw a straight line between
    # bins. A rolling maximum would pick up and drop whole transects as it slid
    # along and come out as a staircase.
    cuts = np.arange(s_water.min(), s_water.max() + BIN_FT, BIN_FT)
    which = np.clip(np.digitize(s_water, cuts) - 1, 0, len(cuts) - 1)
    b_s, b_lo, b_hi = [], [], []
    for b in np.unique(which):
        m = which == b
        b_s.append(s_water[m].mean()); b_lo.append(n_water[m].min()); b_hi.append(n_water[m].max())
    survey_edge = {'left': np.interp(grid, b_s, b_lo),
                   'right': np.interp(grid, b_s, b_hi)}
    has_data = np.min(np.abs(grid[:, None] - s_water[None, :]), axis=1) <= WINDOW_FT

    curves = {}
    for side, outward in (('left', -1.0), ('right', +1.0)):
        pts_side = P[P['side'] == side].sort_values('station')
        edge = survey_edge[side]
        if len(pts_side) >= 2:
            bs = pts_side['station'].values
            bn = pts_side['offset'].values
            near = np.min(np.abs(grid[:, None] - bs[None, :]), axis=1) <= BANK_REACH_FT
            measured = np.where(near, np.interp(grid, bs, bn), np.nan)
            overlap = np.isfinite(measured)
            margin = max(float(np.median(outward * (measured[overlap] - edge[overlap]))), 0.0)
            print(f'  {side:5s}: {len(pts_side)} bank points; beyond them the edge is the '
                  f'survey plus {margin:.1f} ft (measured)')
        else:
            measured = np.full(grid.shape, np.nan)
            margin = 10.0
            print(f'  {side:5s}: too few bank points to measure a margin, using '
                  f'{margin:.1f} ft')

        curve = np.where(np.isfinite(measured), measured, edge + outward * margin)
        curve = np.convolve(curve, np.ones(int(round(SMOOTH_FT / GRID_FT))), 'same') / \
                np.convolve(np.ones_like(curve), np.ones(int(round(SMOOTH_FT / GRID_FT))), 'same')
        limit = edge + outward * SAFETY_FT                    # stay outside the survey
        curves[side] = np.where(outward * (limit - curve) > 0, limit, curve)

    narrow = (curves['right'] - curves['left']) < MIN_WIDTH_FT
    mid = 0.5 * (curves['left'] + curves['right'])
    curves['left'] = np.where(narrow, mid - MIN_WIDTH_FT / 2, curves['left'])
    curves['right'] = np.where(narrow, mid + MIN_WIDTH_FT / 2, curves['right'])

    # =========================================================================
    # 4. close into a polygon, and make it hold every survey point
    # =========================================================================
    gaps = np.flatnonzero(np.diff(grid[has_data]) > MAX_GAP_FT) + 1
    blocks = [b for b in np.split(np.flatnonzero(has_data), gaps) if len(b) >= 2]
    print(f'  {len(blocks)} continuous stretch(es) of river with data')

    survey_pts = gpd.GeoSeries(gpd.points_from_xy(wet['Easting'], wet['Northing']),
                               crs=target_crs)

    def build_shapes():
        edges, shapes = {}, []
        for side in ('left', 'right'):
            edges[side] = [LineString(np.column_stack(
                to_map_frame(centerline, grid[idx], curves[side][idx]))) for idx in blocks]
        for k in range(len(blocks)):
            ring = np.vstack([np.array(edges['left'][k].coords),
                              np.array(edges['right'][k].coords)[::-1]])
            poly = Polygon(ring)
            shapes.append(poly if poly.is_valid else poly.buffer(0))
        return edges, shapes

    # Containment was enforced on the curve, in reach coordinates. That is not the
    # same as containment on the ground, because the reach frame distorts where the
    # centerline bends and the polygon edge is a straight chord between stations.
    # So build it, find what fell outside, push those stations out, and repeat.
    for attempt in range(MAX_PASSES):
        bank_edges, polys = build_shapes()
        extent = gpd.GeoSeries(polys, crs=target_crs).union_all()
        outside = ~survey_pts.within(extent.buffer(0.01))
        if not outside.any():
            break
        push = survey_pts[outside].distance(extent) + 0.5
        # take the LARGEST push each station needs, once - adding them up would let
        # neighbouring points compound and throw the edge far past anything measured
        need = {'left': np.zeros(grid.shape), 'right': np.zeros(grid.shape)}
        for ps, pn, dist in zip(s_water[outside.values], n_water[outside.values], push.values):
            i = int(np.argmin(np.abs(grid - ps)))
            side = 'left' if abs(pn - curves['left'][i]) <= abs(pn - curves['right'][i]) else 'right'
            near = np.abs(grid - ps) <= WINDOW_FT
            need[side][near] = np.maximum(need[side][near], dist)
        curves['left'] -= need['left']
        curves['right'] += need['right']
    if attempt:
        print(f'  {attempt} extra pass(es) to pull the edge outside every point on the ground')

    kept = np.concatenate(blocks)          # only the stations that make it into a polygon
    width = (curves['right'] - curves['left'])[kept]
    acres = sum(p.area for p in polys) / 43560.0
    print(f'  width {width.min():.0f} to {width.max():.0f} ft (median {np.median(width):.0f})')
    print(f'  CHECK: {int(outside.sum())} of {len(survey_pts)} survey points outside the '
          f'polygon' + ('' if not outside.any() else '   <-- investigate'))

    gpd.GeoDataFrame(
        [dict(section=section_name, stretch=k + 1,
              area_acres=round(p.area / 43560.0, 2),
              bank_elev_min_ft=round(float(P['elev'].min()), 2),
              bank_elev_max_ft=round(float(P['elev'].max()), 2),
              note='streambed at maximum stage; surveyed bank where available, '
                   'otherwise the survey edge plus the measured margin',
              geometry=p) for k, p in enumerate(polys)],
        crs=target_crs).to_file(f'mapping/channel_extent_{section_name}.geojson',
                                driver='GeoJSON')
    gpd.GeoDataFrame(
        [dict(section=section_name, bank=side, stretch=k + 1, geometry=g)
         for side in ('left', 'right') for k, g in enumerate(bank_edges[side])],
        crs=target_crs).to_file(f'mapping/bank_lines_{section_name}.geojson',
                                driver='GeoJSON')
    print(f'  wrote mapping/channel_extent_{section_name}.geojson  ({acres:.1f} acres)')
    print(f'  wrote mapping/bank_lines_{section_name}.geojson')

    extent_summary.append(dict(section=section_name, bank_points=len(P),
                               stretches=len(blocks),
                               median_width_ft=round(float(np.median(width)), 1),
                               acres=round(acres, 1),
                               points_outside=int(outside.sum())))

    # --- pictures ------------------------------------------------------------
    ax.scatter(s_water, n_water, s=4, c='0.75', label='survey')
    for side, colour in (('left', '#1A9850'), ('right', '#762A83')):
        q = P[P['side'] == side]
        ax.plot(q['station'], q['offset'], 'o', color=colour, ms=5,
                label=f'{side} surveyed')
        for k, idx in enumerate(blocks):
            ax.plot(grid[idx], curves[side][idx], '-', color=colour, lw=2,
                    label=f'{side} edge' if k == 0 else None)
    ax.set_title(f'{section_name}: the bank as a curve along the reach')
    ax.set_xlabel('distance along centerline (ft)'); ax.grid(alpha=0.3)

    for p in polys:
        xs, ys = p.exterior.xy
        ax2.fill(xs, ys, color='#9ECAE1', alpha=0.6, zorder=1)
    ax2.scatter(wet['Easting'], wet['Northing'], s=4, c='0.35', zorder=2)
    ax2.set_aspect('equal', adjustable='box')
    ax2.set_title(f'{section_name}: streambed extent ({acres:.1f} acres)')
    ax2.set_xlabel('Easting (ft)'); ax2.grid(alpha=0.3)

axgrid[0, 0].set_ylabel('distance from centerline (ft)')
axgrid[0, 0].legend(fontsize=8, loc='best')
axgrid[1, 0].set_ylabel('Northing (ft)')
plt.tight_layout(); plt.show()

print('\n' + '=' * 74)
print(pd.DataFrame(extent_summary).to_string(index=False))
print("""
USING THE POLYGON TO BOUND A MODEL
  extent = gpd.read_file('mapping/channel_extent_Middle.geojson').union_all()
  inside = gpd.GeoSeries(gpd.points_from_xy(grid_E, grid_N),
                         crs=target_crs).within(extent)
  Predict on the full grid, then drop the cells where `inside` is False.

WHAT THIS IS AND IS NOT
  It IS one stage: the streambed area available when the water is as high as it
  gets before leaving the channel, taken from bank points surveyed at that height.

  It is NOT a floodplain boundary, and it is NOT a measured edge everywhere. Where
  no bank point was surveyed, the edge is the outermost survey point plus a margin
  measured elsewhere in the same reach - a lower bound on the true width there.
  The printed margin is the size of that assumption.

  Where the survey skips a long stretch the polygon BREAKS rather than spanning
  it. The reach between stretches is not described, and predictions there are
  unbounded by this file.

  The curvilinear frame limits precision: a few feet on Upper and Lower, about
  fifteen on Middle, whose centerline bends more tightly than the channel is wide.
""")
