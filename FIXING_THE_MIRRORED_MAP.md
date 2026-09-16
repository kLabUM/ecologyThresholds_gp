# Fixing the Mirrored River Map at the Root

**Bug:** every spatial map produced by `modelVis_{Upper,Middle,Lower}.ipynb` has river-left and river-right swapped.

**Evidence:** take the Middle survey points, convert to `(downstream_dist, distance_to_centerline)` with the columns `dataSorting` already wrote, push them back through `model_to_spatial()`, and compare to the original Easting/Northing.

| | median round-trip error |
|---|---|
| `model_to_spatial` as written | **67.02 ft** |
| same, with the sign of `n` flipped | **0.00 ft** |

The median `|distance_to_centerline|` in that sample is 33.5 ft. The error is exactly twice the offset, which is the signature of a reflection rather than a drift.

**Where it is:** `dataSorting.ipynb` section 6 defines `n > 0` as *right* of the downstream direction. `model_to_spatial` builds `perp = (-ty, tx)`, which points *left*, and adds `perp * cl_dist`. Its own comment says "This points to the left of the centerline direction" — the comment is an accurate description of the code and still wrong, because it documents the implementation instead of asserting the contract.

---

## Why flipping the sign is not the fix

Change `perp` in three notebooks and the maps come out right today. Nothing stops it coming back, because the conditions that produced it are all still there:

1. **Two independent implementations of one transform.** The forward `(E,N) → (s,n)` is in `dataSorting`; the inverse `(s,n) → (E,N)` is in `modelVis`. Each makes its own sign decision. A problem with one degree of freedom is being expressed with two.
2. **Three copies of the inverse**, one per section notebook. A fix applied to one does not reach the others, and they drift.
3. **The convention lives in comments.** Comments do not fail.
4. **Nothing ever compared the two directions.** The check that finds this bug takes four lines and did not exist.

So the fix has to remove the second degree of freedom, remove the copies, and make the convention executable.

---

## The fix, in five parts

### Part 1 — one module, one sign decision

`reachframe.py` (attached) goes in the project root. It contains the forward transform, the inverse transform, and the tests. Exactly one function, `tangent_normal`, decides which way `+n` points:

```python
    # rotate the tangent by -90 degrees -> right of travel
    nx, ny = ty, -tx
```

**Nothing else in the project may construct a perpendicular vector.** If you find `(-ty, tx)` or `(ty, -tx)` anywhere outside this function, that is the bug returning.

The derivation is written so it can be checked without trusting the code: a line heading east has tangent `(1, 0)`; rotating −90° gives `(0, −1)`, which points south; south is on your right when you face east.

### Part 2 — make the two directions inverses by construction

This is the part that actually kills the bug class. Instead of deriving `n` from a cross product and the inverse from a hand-written perpendicular, both use the *same* normal vector:

```python
forward:  n = (P − p(s)) · n̂
inverse:  P = p(s) + n · n̂
```

These are algebraic inverses. There is no second sign to get wrong — not "two places that agree", but one place, used twice.

Verified against your data: `reachframe.to_reach_frame` reproduces the existing `distance_to_centerline` column with **median difference 0.0000 ft and 100% sign agreement** across all three sections. `dataSorting` was right; `modelVis` was wrong.

### Part 3 — two tests, and you need both

`reachframe.self_test()` runs on import, so a broken convention fails at the first cell rather than in a figure three notebooks later.

**Known-answer test** — pins the *absolute* convention:

```python
east = LineString([(0, 0), (100, 0)])          # heading east
s, n = to_reach_frame(east, [50.0], [-10.0])   # a point to the SOUTH
assert n[0] > 0    # south of east-flowing = river-RIGHT
```

**Round-trip test** — pins the two directions *to each other*:

```python
E2, N2 = to_map_frame(curve, s, n)
assert max_distance(P, (E2, N2)) < tol
```

> **This is the part people get wrong.** A round-trip test alone would have passed on the buggy code if both directions had been flipped together — the map would still be mirrored and the test would still be green. The round trip proves consistency; only the known-answer test proves correctness. Keep both.

### Part 4 — delete the duplicates, import instead

In **cell 0** of each `modelVis_{S}.ipynb`, add:

```python
import sys, os
sys.path.insert(0, os.path.abspath('.'))       # project root, where reachframe.py lives
from reachframe import to_map_frame, to_reach_frame, CONVENTION, round_trip_error
```

In **cell 2**, delete the whole body of `model_to_spatial` and replace it with a shim, so existing call sites keep working and become correct at the same time:

```python
def model_to_spatial(downstream_dist, distance_to_centerline):
    """(s, n) -> (Easting, Northing). Delegates to reachframe - do not reimplement.

    gp_data_3v carries downstream_dist and distance_to_centerline in METRES while
    Easting/Northing stay in FEET, so convert on the way in.
    """
    M_TO_FT = 1.0 / FT_TO_M
    s_ft = reference_distance + np.atleast_1d(downstream_dist) * M_TO_FT
    n_ft = np.atleast_1d(distance_to_centerline) * M_TO_FT
    return to_map_frame(centerline_geom, s_ft, n_ft)
```

This also fixes the cell 14 threshold figure, which calls `model_to_spatial` and therefore inherits the mirror.

Do the same in `dataSorting.ipynb`: replace the body of `distances_to_line` with a call to `to_reach_frame`. One transform, one definition, three notebooks importing it.

**Guard against a stale local copy.** Put this at the end of cell 2, where it runs every session:

```python
_E, _N = model_to_spatial(gp_data_3v['downstream_dist'], gp_data_3v['distance_to_centerline'])
_err = np.nanmedian(np.hypot(_E - gp_data_3v['Easting'], _N - gp_data_3v['Northing']))
assert _err < 5.0, (
    f"REACH FRAME BROKEN: {_err:.1f} ft median round-trip error. If this is close to "
    f"the channel width, the map is mirrored - something has reimplemented the "
    f"perpendicular instead of calling reachframe.")
print(f"reach frame OK ({CONVENTION}), round-trip {_err:.2f} ft")
```

That assertion is the thing that makes this permanent. It would have caught the original bug the first time a map was drawn.

### Part 5 — stamp the convention into the data

Any file carrying `distance_to_centerline` should say what the sign means, so a future reader verifies instead of guessing. When writing the GeoJSONs in `dataSorting`:

```python
gdf.attrs['frame_convention'] = CONVENTION      # 'n_positive_right_of_downstream'
```

and add a `frame_convention` column to the section CSVs, or a one-line sidecar. Costs nothing; removes the ambiguity that created this.

---

## Verification checklist

Run these in order. Each has a definite expected answer.

1. `python -c "import reachframe; reachframe.self_test(verbose=True)"` → prints OK.
2. Re-run `modelVis_Middle` cell 2 → the new assertion prints a round-trip under 5 ft.
3. Deliberately break it: flip `nx, ny = ty, -tx` to `-ty, tx` in `reachframe.py`. The known-answer test must fail **on import**, before any figure is drawn, with:

   ```
   AssertionError: CONVENTION BROKEN: a point south of an east-flowing centerline
   is on RIVER-RIGHT and must have n > 0, got n = -10.000.
   ```

   Then put it back. **Restart the kernel (or `rm -rf __pycache__`) afterwards** — a stale `.pyc` will keep serving the broken version and make you think the restore failed. This caught me while testing.
4. Overlay the training points on a spatial map. They must land on the coloured field, not its mirror.
5. `grep -rn "ty, tx\|tangent_y, tangent_x" --include=*.ipynb --include=*.py .` → only `reachframe.py` should match.

---

## What has to be regenerated

Every spatial figure ever produced by the modelVis notebooks is mirrored. That includes the cell 14 threshold maps I gave you. **No model needs refitting** — the GP was fitted in `(s, n)` space and never saw the map coordinates, so `gp_data_3v`, the kernels and all the CV numbers are unaffected. Only the projection back to Easting/Northing was wrong. Re-run the visualization notebooks and discard the old PNGs.

---

## A second problem this exposed

Building the round-trip check surfaced something separate: **the curvilinear frame is close to degenerate on these reaches.**

A `(s, n)` frame stops being one-to-one once `|n|` approaches the local radius of curvature `R` — inside the centre of curvature, several `(s, n)` pairs map to the same ground point and no inverse exists.

| section | R_min (ft) | max &#124;n&#124; / R | round-trip error |
|---|---|---|---|
| Upper | 117 | 0.66 | 2.4 ft |
| Middle | 151 | **0.80** | **15.3 ft** |
| Lower | 385 | 0.52 | 4.3 ft |

Middle's centerline bends to a 151 ft radius while the channel half-width reaches ~120 ft. That 15 ft round-trip error is *not* the sign bug — it is the frame itself, and it is a floor on the spatial precision of every map and every masked area statistic.

Resampling and smoothing the centerline helps Middle substantially and the others not at all:

| section | as built | resampled 25 ft + 9-point smooth |
|---|---|---|
| Upper | 2.4 ft | 3.3 ft |
| Middle | 15.3 ft | **2.2 ft** |
| Lower | 4.3 ft | 7.7 ft |

So Middle's error was largely NHD vertex noise; Upper's and Lower's is real meander curvature. `reachframe.curvature_check()` reports this, and the medial centerline from the channel-geometry cell should improve it further, since a medial line is smoother than a digitized flowline.

Worth a sentence in the methods: cross-stream position is accurate to a few feet, not to inches, and the curvilinear frame is the reason.

---

## Summary

| | before | after |
|---|---|---|
| implementations of the transform | 2 (forward, inverse) | 1 (inverse by construction) |
| copies of the inverse | 3 notebooks | 1 module |
| sign decisions in the codebase | 4+ | 1 |
| convention recorded as | a comment | an assertion that runs on import |
| test comparing the two directions | none | known-answer + round-trip, every session |
