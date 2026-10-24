# GED-V1-R2 — Development audit (2026-10-24)

## Design

- New experiment, independent of the earlier arbitrary V1-R1 `alpha` toy coupling.
- Physics baseline: two uniform-density concentric layers with identical constituent inventories; analytic Newtonian self-potential energy and E/c² accounting for energy crossing the body's boundary.
- Decimal precision fixed at 60 significant digits in model calculation; source has no external dependencies.
- Tests include independent numerical Simpson quadrature on each density profile and uniform-density analytic `U = -3GM²/(5R)`.

## Preserved failed initial attempt (development, not science)

Initial `python -m unittest -q test_ged_v1_r2.py` returned **4 failures / 26 tests**, due to test-side operations using the default Decimal precision (28 digits), while the model intentionally uses local Decimal precision 60. Failed tests:

- `test_absorb_vs_emit`: a 60-digit quantity compared against a 28-digit negation.
- `test_delta_g_inverse_square`: a last-place default-context rounding difference.
- `test_mass_parts_sum`: adding 60-digit component masses at default 28-digit context.
- `test_uniform_sphere_limit`: reference formula evaluated at default 28 digits.

**Correction:** Added `setUp()` in `test_ged_v1_r2.py` to use `getcontext().prec=60` for each test. No physical formula changed. Rerun: **26/26 pass**. Original failure output is summarized above; the original full failing terminal transcript was not retained as a separate file. Do not upgrade this to a clean-room scientific replication claim.

## Gates and limitations

Final local execution and input/output hash verification remain acceptance gates. No experiment provides empirical evidence for a new material-information gravity channel. No new physics should be inferred from Decimal precision or from a numerical code agreeing with its stated model equations. Retain failures and the accepted version with hashes under Manifest of Destiny.
