# GED-V1-R4 — Development and failed-test audit

## Source inheritance

- `ged_v1_r4.py` is the sole R4 implementation.
- The only local source dependency is pinned `ged_v1_r2.py`, SHA-256 `c1c8bd268fcfad2f8dea23bba703ef04284ced2590680584ccfa263cd9010813`. Its existing spherical volume, decimal parsing and E/c² conversion functions are reused.
- No new arbitrary entropy-to-gravity term was inserted.
- Each run captures exact source and environment receipts. This does not substitute for the authoritative standalone Turnip Run wrapper, which has not been installed or executed here.

## Initial development failure (preserved)

The first internal suite executed **55 tests, 15 failed**. Exact output is retained as `GED_V1_R4_PRECISION_FAILURE_A.log`.

Root causes were entirely in **test comparisons**, not changes in science formulas:

1. A few expected values were computed in the ambient 28-digit `Decimal` context rather than the model's explicit 60-digit precision.
2. Exact zero was represented by `Decimal` as strings such as `0E-45`, which compare mathematically equal to zero but are not the literal string `'0'`.

## Corrective action

- Fixed every arithmetic comparison to use the same explicit high-precision context, or used a physically appropriate relative tolerance.
- Compared numeric zero as a `Decimal` quantity rather than a text representation.
- Retained strict zero energy ledgers and exterior-null tests. No scientific model parameter or equation was changed to make the suite pass.
- Added a third reverse-gradient control to distinguish equal entropy from internal mass-energy arrangement.

## Independent checks

The final test suite includes independent sphere mass, heat capacity, analytical Boltzmann entropy and interior-gravity difference calculations. They do not simply compare the script against its own output. Remaining audit limits: the tests verify a mathematical model, not a measurement or a full relativistic result.
