### 9) Complete bank lines: extend the surveyed bank into continuous waterlines
# =============================================================================
# Add this as a new cell AFTER the bank-waterline cell in dataSorting.ipynb.
# It reads the six measured bank files and writes nine new ones. No new columns.
#
#   mapping/waterline_left_{Section}.geojson
#   mapping/waterline_right_{Section}.geojson
#   mapping/channel_extent_{Section}.geojson      <- the polygon between them
#
# NEW IMPORTS (everything else is already in cell 0)
#   from shapely.geometry import Polygon
#   from reachframe import to_reach_frame, to_map_frame      # project root
#
# THE IDEA, IN ONE SENTENCE
#   Describe each bank as a curve n(s) - how far it sits from the centerline at
#   each distance downstream - because in that form a bank is just two numbers
#   per station and can be interpolated, smoothed and extended like any other
#   1-D series. Then convert back to Easting/Northing.
#
# WHAT FILLS THE GAPS
#   Two kinds of evidence, in order of preference:
#     1. A measured bank point. Used directly wherever one is within BANK_REACH_FT.
#     2. The outermost in-water point, pushed outward by MARGIN. The margin is not
#        assumed - it is measured, as the median distance between the bank and the
#        outermost in-water point on the stations where both exist.
#   Where neither exists the curve is interpolated between the nearest stations
#   that do, and held constant beyond the last one.
#
# THE GUARANTEE
#   After smoothing, the curve is forced back outside every in-water point. No
#   training point can end up outside the waterline. The count of stations that
#   needed that correction is printed - if it is large, the smoothing is fighting
#   the data and SMOOTH_FT should come down.
# =============================================================================

from shapely.geometry import Polygon
from reachframe import to_reach_frame, to_map_frame

# ------------------------------ SETTINGS -------------------------------------
LINE_SECTIONS = ['Upper', 'Middle', 'Lower']
GRID_FT       = 50.0    # spacing of the stations along the reach
BIN_FT        = 100.0   # length of reach each survey-edge sample summarises
WINDOW_FT     = 200.0   # a station counts as surveyed if data lies this close
BANK_REACH_FT = 400.0   # how far a measured bank point is trusted along the reach
SMOOTH_FT     = 300.0   # moving-average length applied to the curve
SAFETY_FT     = 2.0     # keep the waterline at least this far outside the data
MIN_HALF_FT   = 5.0     # never let the two banks meet
MAX_GAP_FT    = 800.0   # break the line across a stretch with no data at all

os.makedirs('mapping', exist_ok=True)
print('=' * 74)
print('COMPLETE BANK LINES (maximum water extent)')
print('=' * 74)

fig, axgrid = plt.subplots(2, len(LINE_SECTIONS),
                           figsize=(6 * len(LINE_SECTIONS), 10))
axgrid = np.atleast_2d(axgrid)
line_summary = []


def _moving_average(y, n_nodes):
    """Moving average that ignores the NaNs and keeps the array length."""
    if n_nodes < 2:
        return y
    k = np.ones(int(n_nodes))
    ok = np.isfinite(y).astype(float)
    filled = np.where(np.isfinite(y), y, 0.0)
    num = np.convolve(filled, k, mode='same')
    den = np.convolve(ok, k, mode='same')
    out = np.where(den > 0, num / np.maximum(den, 1e-9), np.nan)
    return out


def _fill_ends(x, y):
    """Interpolate the NaNs, and hold the first/last real value beyond the ends."""
    good = np.isfinite(y)
    if not good.any():
        return y
    return np.interp(x, x[good], y[good])


for ax, ax2, section_name in zip(axgrid[0], axgrid[1], LINE_SECTIONS):
    print(f'\n{section_name}')
    print('-' * 74)

    centerline = gpd.read_file(
        f'mapping/centerline_union_{section_name}.geojson').geometry.iloc[0]

    # --- the data the model sees ---------------------------------------------
    w = pd.read_csv(f'sorted_data/complete_data_{section_name}.csv')
    w = w.dropna(subset=['dist_along_centerline', 'distance_to_centerline'])
    s_water = w['dist_along_centerline'].values.astype(float)
    n_water = w['distance_to_centerline'].values.astype(float)

    # --- the measured bank, turned into (station, offset) --------------------
    # The left bank is the low side of the channel and the right bank the high
    # side. Do not use the sign of the offset to tell them apart: on a few Upper
    # transects the centerline runs outside the channel and both banks come out
    # negative. Membership of the left/right file is the reliable label.
    measured = {}
    for label in ('left', 'right'):
        path = f'mapping/bank_{label}_{section_name}.geojson'
        if not os.path.exists(path):
            measured[label] = (np.array([]), np.array([])); continue
        parts = gpd.read_file(path).geometry
        xy = np.vstack([np.array(p.coords) for p in parts])
        bs, bn = to_reach_frame(centerline, xy[:, 0], xy[:, 1])
        order = np.argsort(bs)
        measured[label] = (bs[order], bn[order])
        print(f'  {label:5s} measured bank: {len(bs)} points, '
              f'station {bs.min():,.0f} to {bs.max():,.0f} ft')

    # --- stations along the reach, over the span the data covers -------------
    grid = np.arange(s_water.min(), s_water.max() + GRID_FT, GRID_FT)

    # The edge of the survey, as a curve. Collect the data into bins of BIN_FT
    # along the reach, take the outermost point on each side within each bin, and
    # draw a straight line between bins. Do NOT take a rolling maximum over a
    # window: transects are 50-400 ft apart, so a window picks up and drops whole
    # transects as it slides and the result comes out as a staircase.
    edges = np.arange(s_water.min(), s_water.max() + BIN_FT, BIN_FT)
    which = np.clip(np.digitize(s_water, edges) - 1, 0, len(edges) - 1)
    b_s, b_lo, b_hi = [], [], []
    for b in np.unique(which):
        m = which == b
        b_s.append(s_water[m].mean())
        b_lo.append(n_water[m].min())
        b_hi.append(n_water[m].max())
    b_s, b_lo, b_hi = np.array(b_s), np.array(b_lo), np.array(b_hi)
    lo = np.interp(grid, b_s, b_lo)
    hi = np.interp(grid, b_s, b_hi)

    # a station counts as covered if real data sits within WINDOW_FT of it
    has_data = np.min(np.abs(grid[:, None] - s_water[None, :]), axis=1) <= WINDOW_FT

    curves, needed_fix = {}, {}
    for label, edge, outward in (('left', lo, -1.0), ('right', hi, +1.0)):
        bs, bn = measured[label]

        # 1. the measured bank, but only near a real measurement
        if len(bs):
            near = np.min(np.abs(grid[:, None] - bs[None, :]), axis=1) <= BANK_REACH_FT
            from_bank = np.where(near, np.interp(grid, bs, bn), np.nan)
        else:
            near = np.zeros(grid.shape, dtype=bool)
            from_bank = np.full(grid.shape, np.nan)

        # 2. the margin, MEASURED where the bank and the data overlap
        both = np.isfinite(from_bank) & np.isfinite(edge)
        if both.sum() >= 3:
            margin = float(np.median(outward * (from_bank[both] - edge[both])))
            margin = max(margin, 0.0)
            print(f'  {label:5s} margin beyond the outermost data: {margin:.1f} ft '
                  f'(median of {both.sum()} stations)')
        else:
            margin = 10.0
            print(f'  {label:5s} margin: too little overlap to measure, using '
                  f'{margin:.1f} ft')

        # 3. prefer the measurement, fall back to the data plus the margin
        curve = np.where(np.isfinite(from_bank), from_bank, edge + outward * margin)
        curve = _fill_ends(grid, curve)

        # 4. smooth
        curve = _moving_average(curve, max(1, round(SMOOTH_FT / GRID_FT)))
        curve = _fill_ends(grid, curve)

        # 5. force it back outside every data point
        limit = edge + outward * SAFETY_FT
        inside = np.isfinite(limit) & (outward * (limit - curve) > 0)
        curve = np.where(inside, limit, curve)
        needed_fix[label] = int(inside.sum())
        curves[label] = curve

    # never let the banks touch or cross
    middle = 0.5 * (curves['left'] + curves['right'])
    too_narrow = (curves['right'] - curves['left']) < 2 * MIN_HALF_FT
    curves['left'] = np.where(too_narrow, middle - MIN_HALF_FT, curves['left'])
    curves['right'] = np.where(too_narrow, middle + MIN_HALF_FT, curves['right'])

    for label in ('left', 'right'):
        pct = 100.0 * needed_fix[label] / len(grid)
        print(f'  {label:5s}: {needed_fix[label]} of {len(grid)} stations ({pct:.0f}%) '
              f'had to be pushed back out after smoothing')

    # --- keep only the stretches that have data ------------------------------
    cut = np.flatnonzero(np.diff(grid[has_data]) > MAX_GAP_FT) + 1
    blocks = [b for b in np.split(np.flatnonzero(has_data), cut) if len(b) >= 2]
    print(f'  {len(blocks)} continuous stretch(es) of river with data')

    # --- build the lines and the polygon, then check IN MAP SPACE ------------
    # Containment was enforced on the curve, in reach coordinates. That is not
    # quite the same as containment on the ground: the reach frame is not a
    # perfect map where the centerline bends tightly, and the polygon edge is a
    # straight chord between stations. So build it, look for survey points that
    # ended up outside, nudge the nearby stations out, and repeat. This usually
    # finishes in one or two passes and costs a foot or two of width where it
    # acts, instead of widening the whole reach by the frame's error.
    pts = gpd.GeoSeries(gpd.points_from_xy(w['Easting'], w['Northing']),
                        crs=target_crs)

    def build():
        lines = {}
        for label in ('left', 'right'):
            lines[label] = [LineString(np.column_stack(
                to_map_frame(centerline, grid[idx], curves[label][idx])))
                for idx in blocks]
        shapes = []
        for k in range(len(blocks)):
            ring = np.vstack([np.array(lines['left'][k].coords),
                              np.array(lines['right'][k].coords)[::-1]])
            poly = Polygon(ring)
            shapes.append(poly if poly.is_valid else poly.buffer(0))
        return lines, shapes

    MAX_PASSES = 12
    for attempt in range(MAX_PASSES):
        seg_lines, polys_geom = build()
        extent = gpd.GeoSeries(polys_geom, crs=target_crs).union_all()
        outside = ~pts.within(extent.buffer(0.01))
        if not outside.any():
            break
        push = pts[outside].distance(extent) + 0.5
        # Collect the biggest push each station needs and apply it once. Adding
        # them up instead would let neighbouring points compound and throw the
        # edge tens of feet past anything that was measured.
        need = {'left': np.zeros(grid.shape), 'right': np.zeros(grid.shape)}
        for point_s, point_n, dist in zip(s_water[outside.values],
                                          n_water[outside.values], push.values):
            i = int(np.argmin(np.abs(grid - point_s)))
            side = ('left' if abs(point_n - curves['left'][i])
                    <= abs(point_n - curves['right'][i]) else 'right')
            near = np.abs(grid - point_s) <= WINDOW_FT
            need[side][near] = np.maximum(need[side][near], dist)
        curves['left'] -= need['left']
        curves['right'] += need['right']
    n_outside = int(outside.sum())
    if attempt:
        print(f'  {attempt} extra pass(es) needed to pull the edge outside every '
              f'survey point on the ground')

    written = {}
    for label in ('left', 'right'):
        rows = []
        for k, idx in enumerate(blocks, 1):
            rows.append(dict(section=section_name, bank=label, stretch=k,
                             station_from_ft=round(float(grid[idx].min()), 1),
                             station_to_ft=round(float(grid[idx].max()), 1),
                             n_stations=len(idx),
                             note='maximum water extent; measured bank where '
                                  'available, otherwise the outermost survey '
                                  'point plus the measured margin',
                             geometry=seg_lines[label][k - 1]))
        path = f'mapping/waterline_{label}_{section_name}.geojson'
        gpd.GeoDataFrame(rows, crs=target_crs).to_file(path, driver='GeoJSON')
        written[label] = rows
        print(f'  wrote {path}')

    polys = [dict(section=section_name, stretch=k + 1,
                  area_acres=round(g.area / 43560.0, 2),
                  note='maximum water extent, one stage', geometry=g)
             for k, g in enumerate(polys_geom)]
    path = f'mapping/channel_extent_{section_name}.geojson'
    gpd.GeoDataFrame(polys, crs=target_crs).to_file(path, driver='GeoJSON')
    total_acres = sum(p['area_acres'] for p in polys)
    print(f'  wrote {path}  ({total_acres:.1f} acres over {len(polys)} stretch(es))')
    print(f'  CHECK: {n_outside} of {len(pts)} survey points fall outside the polygon'
          + ('' if n_outside == 0 else '   <-- investigate'))
    width = curves['right'] - curves['left']
    print(f'  width across the reach: {width.min():.0f} to {width.max():.0f} ft '
          f'(median {np.median(width):.0f})')

    line_summary.append(dict(section=section_name, stretches=len(blocks),
                             median_width_ft=round(float(np.median(width)), 1),
                             acres=round(total_acres, 1), points_outside=outside))

    # --- pictures ------------------------------------------------------------
    ax.scatter(s_water, n_water, s=4, c='0.75', label='survey')
    for label, colour in (('left', '#1A9850'), ('right', '#762A83')):
        bs, bn = measured[label]
        ax.plot(bs, bn, 'o', color=colour, ms=5, label=f'{label} measured')
        for k, idx in enumerate(blocks):
            ax.plot(grid[idx], curves[label][idx], '-', color=colour, lw=2,
                    label=f'{label} waterline' if k == 0 else None)
    ax.set_title(f'{section_name}: bank as a curve along the reach')
    ax.set_xlabel('distance along centerline (ft)'); ax.grid(alpha=0.3)

    for p in polys:
        xs, ys = p['geometry'].exterior.xy
        ax2.fill(xs, ys, color='#9ECAE1', alpha=0.55, zorder=1)
    ax2.scatter(w['Easting'], w['Northing'], s=4, c='0.35', zorder=2, label='survey')
    ax2.set_aspect('equal', adjustable='box')
    ax2.set_title(f'{section_name}: maximum water extent')
    ax2.set_xlabel('Easting (ft)'); ax2.grid(alpha=0.3)

axgrid[0, 0].set_ylabel('distance from centerline (ft)')
axgrid[0, 0].legend(fontsize=8, loc='best')
axgrid[1, 0].set_ylabel('Northing (ft)')
plt.tight_layout(); plt.show()

print('\n' + '=' * 74)
print(pd.DataFrame(line_summary).to_string(index=False))
print("""
HOW TO USE THE POLYGON
  extent = gpd.read_file('mapping/channel_extent_Middle.geojson').union_all()
  keep   = gpd.points_from_xy(grid_E, grid_N)
  mask   = gpd.GeoSeries(keep, crs=target_crs).within(extent)
  Predict everywhere, then drop the cells where `mask` is False.

WHAT THIS IS AND IS NOT
  It IS one stage: the widest the water gets before it leaves the channel, taken
  from the bank points that were surveyed at that height.

  It is NOT a floodplain boundary, and it is NOT a measured edge everywhere. Where
  no bank point was surveyed the edge is the outermost survey point plus a margin
  measured elsewhere in the same reach, so it is a lower bound on the true width
  in those stretches. The printed margin tells you how big that assumption is.

  The reach frame itself limits precision: converting back to Easting/Northing
  carries a few feet of error on Upper and Lower and about fifteen feet on Middle,
  whose centerline bends more tightly than the channel is wide.
""")
