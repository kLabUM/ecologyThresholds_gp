### 4) One velocity per survey point, transferred to a common height above the bed
# =============================================================================
#
# WHY THIS REPLACES A FILTER
#   Each survey point was measured at several heights above the bed, and the 2-D
#   model needs one velocity per point. Keeping only the readings taken at exactly
#   one height throws away every point that was not sampled there. Instead the
#   readings are transferred to a common height.
#
# THE TRANSFER: a one-stage logarithmic profile
#
#       v(h) = v(z) * ln(h / z_p) / ln(z / z_p)
#
#   z_p is the PROFILE SCALE PARAMETER: one fitted constant per section, the height
#   at which the fitted logarithmic curve extrapolates to zero velocity.
#
#   *** z_p IS NOT A BED ROUGHNESS LENGTH. *** It is named z_p, not z0, to keep that
#   clear. Three things rule out the roughness reading, all measured from this data:
#     - Fitted per substrate class it runs BACKWARDS: 27.7 mm on silt and muck,
#       14.5 mm on sand, 8.9 mm on gravel. Roughness must rise with grain size.
#       It is tracking hydraulic setting, since fine material settles in slow pools
#       and gravel is exposed in fast riffles.
#     - Where the notes state a grain size, the grain-based value (k_s/30 with
#       k_s ~ 2.5d) is 0.6 to 2.9 mm. The fit gives 3 to 14 times more, consistent
#       with form drag from bedforms, vegetation and woody debris.
#     - It lands far below the shallowest measurement, so nothing in the survey
#       constrains it directly.
#   Shear velocity is never needed: it cancels in the ratio, and each point is
#   anchored on its own measured reading.
#
# READINGS BELOW FIT_MIN_HEIGHT_FT ARE NOT USED
#   The near-bed layer is where a profile law is not expected to hold. Readings
#   below the cutoff are excluded from the fit AND from use as a transfer source,
#   so the law is never applied outside the range it was fitted on. In open-channel
#   flow the logarithmic region reaches roughly half the flow depth, and this cutoff
#   sits at about a fifth to a quarter of it, so the retained data lie inside that
#   region.
#
# UNITS
#   "Depth" is the height above the bed in FEET, so the settings below are in feet.
#   Output is in the same units as "Velocity" (ft/s) - convert to SI later, in the
#   model notebook, exactly as "Velocity" already is.
#
# WHAT IT ADDS  (names come from VEL_COLUMN)
#   vel_1     the velocity at TARGET_HEIGHT_FT above the bed
#   vel_1_sd  how uncertain that number is: 0 where it was measured directly, the
#              section's measured transfer error where it was converted. Pass it to
#              the GP as a per-point `alpha` so a converted value is not trusted
#              like a measured one.
#
#   One row survives per survey point - the reading taken closest to the target
#   height - so every other column stays a real measurement from that point.
#
# USES ONLY THE LIBRARIES ALREADY IMPORTED IN CELL 0 (pandas as pd, numpy as np)
# =============================================================================

# ------------------------------ SETTINGS -------------------------------------
TARGET_HEIGHT_FT  = 1.0      # height above the bed that every reading is moved to
FIT_MIN_HEIGHT_FT = 0.5      # ignore readings below this, for fitting and as sources
VEL_COLUMN        = 'vel_1'  # output column name, matching TARGET_HEIGHT_FT above -
                             # rename both together if the target height changes
MEASURED_ONLY     = []       # sections to keep ONLY direct 1 ft measurements for,
                             # discarding converted points. Tested result: converted
                             # points did not improve spatially blocked prediction on
                             # any reach, and on Upper the 10 of them sit within 21 ft
                             # of a point already measured. Middle is the exception
                             # worth keeping - its 94 converted points are the only
                             # data in 3,800 ft of reach. Set to ['Upper','Lower'] to
                             # act on that; left empty here so nothing changes silently.
Z_P_BOUNDS_FT     = (1e-4, 0.30)   # search range for the profile scale parameter
Z_P_STEPS         = 60       # points in the log-spaced search
MIN_LEVELS        = 2        # heights a point needs before it can inform the shape
POINT_KEY         = ['Collection Date', 'Measurement Point ID']   # one survey point

print('=' * 74)
print(f'VELOCITY TRANSFERRED TO {TARGET_HEIGHT_FT} ft ABOVE THE BED  ->  "{VEL_COLUMN}"')
print(f'one-stage log profile, fitted on readings at or above {FIT_MIN_HEIGHT_FT} ft')
print('=' * 74)


def _shape(z, z_p):
    """log of the log-profile shape. The per-point scale is an additive constant
    here, which is why it drops out of both the fit and the transfer."""
    return np.log(np.log(z / z_p))


def _fit_profile_scale(df):
    """Search z_p for the best within-point fit of the logarithmic profile.

    Each point's mean log-velocity and mean shape are removed before scoring, so
    z_p is chosen on the SHAPE of the profiles only and is not pulled around by how
    fast or slow individual locations are.

    Returns z_p, the number of points used, the fraction of within-point variation
    explained, and the range of z_p values that score within 2% of the best.
    """
    q = df[(df['Velocity'] > 0) & (df['Depth'] >= FIT_MIN_HEIGHT_FT)].copy()
    q = q[q.groupby(q['_pt'])['Depth'].transform('nunique') >= MIN_LEVELS]
    if q['_pt'].nunique() < 10:
        return np.nan, 0, np.nan, (np.nan, np.nan)

    lv = np.log(q['Velocity'].values)
    grp = q['_pt'].values
    lv_c = lv - pd.Series(lv).groupby(grp).transform('mean').values
    z = q['Depth'].values

    lo, hi = Z_P_BOUNDS_FT
    hi = min(hi, 0.9 * z.min())                 # z_p must stay below every reading
    if not (hi > lo):
        return np.nan, 0, np.nan, (np.nan, np.nan)

    scores = []
    for z_p in np.logspace(np.log10(lo), np.log10(hi), Z_P_STEPS):
        F = _shape(z, z_p)
        F_c = F - pd.Series(F).groupby(grp).transform('mean').values
        # The shape enters with a coefficient of exactly ONE: v = A * ln(z/z_p)
        # becomes log v = log A + log(ln(z/z_p)), and log A is the per-point offset
        # already removed by centring. Letting a coefficient float here would trade
        # off against z_p and leave both unidentified.
        scores.append((np.sum((lv_c - F_c) ** 2), z_p))
    if not scores:
        return np.nan, 0, np.nan, (np.nan, np.nan)

    sse, z_p = min(scores, key=lambda t: t[0])
    total = np.sum(lv_c ** 2)
    r2 = 1.0 - sse / max(total, 1e-12)
    band = [zp for sc, zp in scores if sc <= sse * 1.02]
    at_bound = (z_p <= lo * 1.01) or (z_p >= hi * 0.99)
    if at_bound:
        print(f'  WARNING: z_p ran to the edge of Z_P_BOUNDS_FT - widen the range '
              f'or treat the value as unidentified')
    return z_p, q['_pt'].nunique(), r2, (min(band), max(band))


def _transfer(heights, speeds, z_p, target):
    """Move a point's reading to `target` height using the fitted log profile.

    Only the reading closest to the target is used, and only readings at or above
    FIT_MIN_HEIGHT_FT qualify. Averaging all of a point's readings was tested and is
    clearly worse - a reading far away says much less about the target than a near
    one, and averaging lets the distant ones drag the answer. Ties are averaged.

    A reading of exactly zero cannot be logged, so zeros are set aside; if every
    qualifying reading is zero the answer is zero, which is right - still water at
    one height is still water at another.
    """
    heights = np.asarray(heights, dtype=float)
    speeds = np.asarray(speeds, dtype=float)
    ok = (np.isfinite(heights) & np.isfinite(speeds) & (speeds >= 0)
          & (heights >= FIT_MIN_HEIGHT_FT))
    heights, speeds = heights[ok], speeds[ok]
    if len(heights) == 0:
        return np.nan
    moving = speeds > 0
    if not moving.any():
        return 0.0
    heights, speeds = heights[moving], speeds[moving]
    gap = np.abs(np.log(heights) - np.log(target))
    closest = gap == gap.min()
    moved = (np.log(speeds[closest])
             + _shape(np.full(closest.sum(), target), z_p) - _shape(heights[closest], z_p))
    return float(np.exp(moved.mean()))


complete_data_by_section = {'Middle': complete_data_Middle,
                            'Upper':  complete_data_Upper,
                            'Lower':  complete_data_Lower}
trimmed_by_section = {}

for section_name, raw in complete_data_by_section.items():
    print(f'\n{section_name}')
    print('-' * 74)
    print(f'  starting rows: {len(raw)}')

    # --- rows that describe a velocity measured at a height ------------------
    work = raw.dropna(subset=['Velocity', 'Depth', 'Measurement Point ID']).copy()
    work = work[work['Depth'] > 0]
    print(f'  rows with a velocity at a known height: {len(work)}')

    # --- the velocity >= 0 filter, kept from the previous version -------------
    before = len(work)
    work = work[work['Velocity'] >= 0]
    print(f'  after velocity >= 0 filter: {len(work)} (removed {before - len(work)})')

    work['_pt'] = work[POINT_KEY].astype(str).agg('|'.join, axis=1)
    usable = work[work['Depth'] >= FIT_MIN_HEIGHT_FT]
    print(f'  survey points: {work["_pt"].nunique()}  '
          f'({usable["_pt"].nunique()} have a reading at or above {FIT_MIN_HEIGHT_FT} ft; '
          f'{work["_pt"].nunique() - usable["_pt"].nunique()} do not and are dropped)')

    # --- 1. the profile scale parameter, fitted on this section's own profiles
    z_p, n_shape_pts, shape_r2, band = _fit_profile_scale(work)
    if not np.isfinite(z_p):
        print('  NOT ENOUGH multi-height points above the cutoff to fit the profile '
              '- section skipped'); continue
    FT_TO_MM = 304.8
    print(f'  profile scale z_p = {z_p:.5f} ft ({z_p*FT_TO_MM:.1f} mm), fitted on '
          f'{n_shape_pts} points measured at {MIN_LEVELS}+ heights above the cutoff')
    print(f'      explains {100*shape_r2:.0f}% of the within-point variation; values from '
          f'{band[0]*FT_TO_MM:.1f} to {band[1]*FT_TO_MM:.1f} mm fit within 2% as well')
    print(f'      report this as a fitted profile constant, NOT as a bed roughness')

    # --- 2. how wrong is a transferred value? measure it, do not assume -------
    # For every point measured AT the target height and also elsewhere above the
    # cutoff: hide the target reading, transfer the others to it, compare.
    errors, truths = [], []
    for _, g in work.groupby('_pt'):
        at_target = g[np.isclose(g['Depth'], TARGET_HEIGHT_FT)]
        others = g[~np.isclose(g['Depth'], TARGET_HEIGHT_FT)]
        if len(at_target) and len(others):
            guess = _transfer(others['Depth'], others['Velocity'], z_p, TARGET_HEIGHT_FT)
            if np.isfinite(guess):
                truth = at_target['Velocity'].mean()
                errors.append(guess - truth); truths.append(truth)
    transfer_sd = float(np.sqrt(np.mean(np.square(errors)))) if len(errors) >= 5 else np.nan
    if np.isfinite(transfer_sd):
        spread = float(np.std(truths)) if len(truths) > 1 else np.nan
        print(f'  transfer error: {transfer_sd:.3f} ft/s, measured on {len(errors)} '
              f'points sampled at the target height and elsewhere')
        print(f'      that is {transfer_sd/spread:.2f} of the spread of the thing being '
              f'predicted ({spread:.3f} ft/s); above ~0.8 a converted value would '
              f'carry little information')
    else:
        print(f'  transfer error: too few points to measure it - {VEL_COLUMN}_sd '
              f'left blank on converted points')

    # --- 3. one value per point ----------------------------------------------
    rows, n_measured, n_converted, n_too_shallow, n_no_source, n_skipped_conv = [], 0, 0, 0, 0, 0
    for _pt, g in work.groupby('_pt'):
        water = g['Total water depth'].median()
        if np.isfinite(water) and water < TARGET_HEIGHT_FT:
            n_too_shallow += 1          # no water at that height - nothing to predict
            continue

        at_target = g[np.isclose(g['Depth'], TARGET_HEIGHT_FT)]
        if len(at_target):
            value, sd = float(at_target['Velocity'].mean()), 0.0
            n_measured += 1
        elif section_name in MEASURED_ONLY:
            n_skipped_conv += 1         # this section keeps measured values only
            continue
        else:
            value = _transfer(g['Depth'], g['Velocity'], z_p, TARGET_HEIGHT_FT)
            if not np.isfinite(value):
                n_no_source += 1        # every reading here is below the cutoff
                continue
            sd = transfer_sd
            n_converted += 1

        # keep the reading taken closest to the target height as this point's row,
        # so every other column is still a real measurement from this point
        keep = g.iloc[(g['Depth'] - TARGET_HEIGHT_FT).abs().values.argmin()].copy()
        keep[VEL_COLUMN] = value
        keep[VEL_COLUMN + '_sd'] = sd
        rows.append(keep)

    out = pd.DataFrame(rows).drop(columns=['_pt']).reset_index(drop=True)
    print(f'  {VEL_COLUMN} set on {len(out)} points: {n_measured} measured directly, '
          f'{n_converted} converted')
    print(f'  dropped, water shallower than {TARGET_HEIGHT_FT} ft: {n_too_shallow} points')
    print(f'  dropped, no reading at or above {FIT_MIN_HEIGHT_FT} ft to transfer from: '
          f'{n_no_source} points')
    if n_skipped_conv:
        print(f'  dropped, MEASURED_ONLY is set for this section: {n_skipped_conv} '
              f'points that would have been converted')
    print(f'  {VEL_COLUMN} range {out[VEL_COLUMN].min():.3f} to {out[VEL_COLUMN].max():.3f} '
          f'ft/s (median {out[VEL_COLUMN].median():.3f})')
    trimmed_by_section[section_name] = out

# --- hand back the three names the rest of the notebook expects --------------
complete_data_Middle = trimmed_by_section['Middle']
complete_data_Upper = trimmed_by_section['Upper']
complete_data_Lower = trimmed_by_section['Lower']

complete_data_Middle.to_csv('sorted_data/complete_data_Middle.csv', index=False)
complete_data_Upper.to_csv('sorted_data/complete_data_Upper.csv', index=False)
complete_data_Lower.to_csv('sorted_data/complete_data_Lower.csv', index=False)

print('\n' + '=' * 74)
print(pd.DataFrame([{'section': s, 'points': len(v),
                     'measured': int((v[VEL_COLUMN + '_sd'] == 0).sum()),
                     'converted': int((v[VEL_COLUMN + '_sd'] != 0).sum()),
                     'median_ft_s': round(float(v[VEL_COLUMN].median()), 3)}
                    for s, v in trimmed_by_section.items()]).to_string(index=False))
print(f"""
FEEDING THE UNCERTAINTY TO THE MODEL
  Train on {VEL_COLUMN} instead of Velocity, and pass {VEL_COLUMN}_sd as the GP's
  per-point noise so a converted value is not trusted like a measured one:

      alpha = (gp_data['{VEL_COLUMN}_sd'] * FT_TO_M) ** 2 + 1e-6
      GaussianProcessRegressor(kernel=..., alpha=alpha)

  alpha is a variance, so it is squared, and it must be converted to SI along with
  the velocity itself. The 1e-6 keeps the directly measured points numerically
  stable rather than exactly noiseless.

WHAT THIS IS
  Velocity {TARGET_HEIGHT_FT} ft above the bed, at survey points where the water is at least that
  deep, from a one-stage logarithmic profile fitted to readings at or above
  {FIT_MIN_HEIGHT_FT} ft. It is not depth-averaged velocity.

  Points in shallower water were dropped - they still exist in the raw data, and
  they are slow, so account for them separately rather than letting the model
  extrapolate into them.

  z_p is a fitted profile constant. It is not a bed roughness length and should not
  be reported as one.
""")
