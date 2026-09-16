"""reachframe.py - the single source of truth for the reach coordinate frame.

WHY THIS FILE EXISTS
--------------------
The forward transform (Easting, Northing) -> (s, n) lived in dataSorting.ipynb and
the inverse (s, n) -> (Easting, Northing) lived, separately and in triplicate, in
modelVis_{Upper,Middle,Lower}.ipynb. The two used opposite signs for n, so every
spatial map came out MIRRORED left-for-right. Nothing caught it, because nothing
ever compared the two.

The fix is structural, not a sign flip:

  1. ONE function decides the sign - `tangent_normal`. Nothing else in the project
     is allowed to build a perpendicular vector.
  2. The forward and inverse transforms are ALGEBRAIC INVERSES by construction:
         forward:  n = (P - p(s)) . n_hat
         inverse:  P = p(s) + n * n_hat
     Both call `tangent_normal` for n_hat. There is no second place to get it
     wrong. The old code derived n from a cross product and the inverse from a
     hand-written perpendicular - two independent sign decisions, which is exactly
     one more than the problem has.
  3. `self_test()` runs a KNOWN-ANSWER check and a ROUND-TRIP check, and is called
     automatically on import.

Note carefully: a round-trip test ALONE is not enough. If both directions were
flipped together the round trip would still close perfectly and the map would
still be mirrored. The known-answer test is what pins the absolute convention to
something a human can verify; the round-trip test is what pins the two directions
to each other. You need both.

THE CONVENTION, stated once
---------------------------
    s  = distance along the centerline, increasing DOWNSTREAM (feet)
    n  = signed perpendicular offset from the centerline (feet)
         n < 0  ->  LEFT  of the downstream direction
         n > 0  ->  RIGHT of the downstream direction

"Left" and "right" are as seen by an observer facing downstream, which is the
standard river convention (river-left / river-right).

USAGE
-----
    from reachframe import to_reach_frame, to_map_frame, CONVENTION

    s, n = to_reach_frame(centerline, df['Easting'], df['Northing'])
    E, N = to_map_frame(centerline, s, n)

Everything is in the units of the centerline geometry - FEET in ESRI:102690 for
this project. Convert to SI at the model boundary, not here.
"""

import numpy as np
from shapely.geometry import LineString, Point

__all__ = ['CONVENTION', 'tangent_normal', 'to_reach_frame', 'to_map_frame',
           'round_trip_error', 'curvature_check', 'self_test']

CONVENTION = 'n_positive_right_of_downstream'
_EPS = 1.0          # feet, for the finite-difference tangent


def tangent_normal(line, s):
    """Point, unit tangent and unit normal at along-line distance(s) `s`.

    THIS IS THE ONLY PLACE IN THE PROJECT THAT DECIDES WHICH WAY +n POINTS.
    Do not write `(-ty, tx)` or `(ty, -tx)` anywhere else.

    The normal is the tangent rotated -90 degrees, i.e. RIGHT of the downstream
    direction. Check it by hand: a line heading east has tangent (1, 0); rotating
    -90 degrees gives (0, -1), which points south; south is on your right when you
    face east. That is the whole derivation, and `self_test` asserts it.

    Returns px, py, tx, ty, nx, ny  (arrays, same shape as `s`).
    """
    s = np.atleast_1d(np.asarray(s, dtype=float))
    L = line.length
    sc = np.clip(s, 0.0, L)
    px = np.empty(sc.shape); py = np.empty(sc.shape)
    tx = np.empty(sc.shape); ty = np.empty(sc.shape)
    for i, sv in enumerate(sc.ravel()):
        p = line.interpolate(float(sv))
        p1 = line.interpolate(max(float(sv) - _EPS, 0.0))
        p2 = line.interpolate(min(float(sv) + _EPS, L))
        dx, dy = p2.x - p1.x, p2.y - p1.y
        d = np.hypot(dx, dy)
        if d < 1e-12:
            dx, dy, d = 1.0, 0.0, 1.0
        idx = np.unravel_index(i, sc.shape)
        px[idx], py[idx] = p.x, p.y
        tx[idx], ty[idx] = dx / d, dy / d
    # rotate the tangent by -90 degrees -> right of travel
    nx, ny = ty, -tx
    return px, py, tx, ty, nx, ny


def to_reach_frame(line, easting, northing):
    """(Easting, Northing) -> (s, n). Inverse of `to_map_frame`.

    n is the projection of the offset vector onto the unit normal, so it is the
    exact algebraic inverse of `to_map_frame`, not a separately-signed distance.
    Points with non-finite coordinates return NaN.
    """
    E = np.atleast_1d(np.asarray(easting, dtype=float))
    N = np.atleast_1d(np.asarray(northing, dtype=float))
    shape = E.shape
    flat_E, flat_N = E.ravel(), N.ravel()

    s_flat = np.full(flat_E.shape, np.nan)
    n_flat = np.full(flat_E.shape, np.nan)

    finite = np.isfinite(flat_E) & np.isfinite(flat_N)
    for i in np.flatnonzero(finite):
        s_flat[i] = line.project(Point(float(flat_E[i]), float(flat_N[i])))

    ok = np.isfinite(s_flat)
    if ok.any():
        px, py, _, _, nx, ny = tangent_normal(line, s_flat[ok])
        # n is the offset vector projected onto the SAME unit normal that
        # to_map_frame adds back. That is what makes them exact inverses.
        n_flat[ok] = (flat_E[ok] - px) * nx + (flat_N[ok] - py) * ny

    return s_flat.reshape(shape), n_flat.reshape(shape)


def to_map_frame(line, s, n):
    """(s, n) -> (Easting, Northing). Inverse of `to_reach_frame`."""
    s = np.asarray(s, dtype=float)
    n = np.asarray(n, dtype=float)
    px, py, _, _, nx, ny = tangent_normal(line, s)
    return px + nx * n, py + ny * n


def curvature_check(line, max_abs_n, n_samples=400):
    """Smallest radius of curvature on `line`, and whether `max_abs_n` is safe.

    A curvilinear (s, n) frame stops being one-to-one once |n| reaches the local
    radius of curvature R: inside the centre of curvature, several (s, n) pairs
    map to the same point and `project` cannot invert. This is a property of the
    coordinate system, not a coding error, so it gets reported rather than fixed.

    Returns (R_min, ratio, ok) with ratio = max_abs_n / R_min.
    """
    L = line.length
    ss = np.linspace(0.0, L, n_samples)
    px, py, tx, ty, _, _ = tangent_normal(line, ss)
    ds = np.gradient(ss)
    dtx, dty = np.gradient(tx) / ds, np.gradient(ty) / ds
    kappa = np.hypot(dtx, dty)                    # |dT/ds|
    kappa = kappa[2:-2]                           # trim finite-difference ends
    k_max = float(np.nanmax(kappa)) if kappa.size else 0.0
    R_min = np.inf if k_max <= 0 else 1.0 / k_max
    ratio = float(max_abs_n / R_min) if np.isfinite(R_min) else 0.0
    return R_min, ratio, ratio < 0.5


def round_trip_error(line, easting, northing):
    """Max distance between each input point and its (s,n)->map reconstruction.

    A left/right sign error shows up as almost exactly 2*|n| - i.e. the reported
    error is the channel width, not a rounding wobble. It is impossible to miss.
    """
    s, n = to_reach_frame(line, easting, northing)
    E2, N2 = to_map_frame(line, s, n)
    d = np.hypot(np.asarray(easting, float) - E2, np.asarray(northing, float) - N2)
    return float(np.nanmax(d))


def self_test(verbose=False):
    """Known-answer + round-trip checks. Raises AssertionError on failure."""
    # --- known answer: pins the ABSOLUTE convention to a hand-checkable fact ---
    east = LineString([(0.0, 0.0), (100.0, 0.0)])          # heading east
    s, n = to_reach_frame(east, [50.0], [-10.0])           # a point to the SOUTH
    assert abs(s[0] - 50.0) < 1e-6, f"s wrong on the straight-line case: {s[0]}"
    assert n[0] > 0, (
        "CONVENTION BROKEN: a point south of an east-flowing centerline is on "
        f"RIVER-RIGHT and must have n > 0, got n = {n[0]:.3f}. Someone has "
        "flipped the normal in tangent_normal().")
    assert abs(n[0] - 10.0) < 1e-6, f"|n| wrong on the straight-line case: {n[0]}"

    s2, n2 = to_reach_frame(east, [50.0], [+10.0])         # a point to the NORTH
    assert n2[0] < 0, f"north of an east-flowing line must be river-LEFT, got {n2[0]:.3f}"

    # --- round trip on a curved line: pins the two directions to each other ---
    # Meander geometry chosen to resemble a real reach: amplitude 150 ft on a
    # 2000 ft wavelength gives R_min ~ 680 ft, and offsets of +/-60 ft keep
    # |n|/R below 0.1. Bending the test curve harder than the offsets allow makes
    # the frame non-invertible and the test fails for a reason that has nothing
    # to do with the sign - see curvature_check().
    x = np.linspace(0.0, 4000.0, 800)
    curve = LineString(np.column_stack([x, 150.0 * np.sin(2 * np.pi * x / 2000.0)]))
    rng = np.random.default_rng(0)
    ss = rng.uniform(100.0, curve.length - 100.0, 400)
    nn = rng.uniform(-60.0, 60.0, 400)
    Ex, Ny = to_map_frame(curve, ss, nn)
    err = round_trip_error(curve, Ex, Ny)
    R_min, ratio, ok = curvature_check(curve, 60.0)
    assert ok, f"test curve is too tightly bent for these offsets (|n|/R = {ratio:.2f})"
    assert err < 1.0, (
        f"ROUND TRIP BROKEN: max reconstruction error {err:.2f} units against a "
        f"channel half-width of 60. An error near 2*|n| means the forward and "
        f"inverse transforms disagree about the sign of n.")
    if verbose:
        print(f"reachframe self-test OK | convention={CONVENTION} | "
              f"round-trip max error {err:.4f} | test curve R_min {R_min:,.0f}")
    return True


self_test()     # runs on import - a broken convention fails loudly and immediately
