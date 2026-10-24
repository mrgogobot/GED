# GED-V1-R5 Development and Acceptance Audit

- **Purpose:** Elastic material response and conventional mass-energy equivalent, without introducing hypothetical gravity laws.
- **Source:** `ged_v1_r5.py`; **tests:** `test_ged_v1_r5.py`; **input:** `GED_V1_R5_config.json`.
- **Dependency:** `ged_v1_r2.py` SHA-256 `c1c8bd268fcfad2f8dea23bba703ef04284ced2590680584ccfa263cd9010813` (already preserved in GED-V1-R4 run receipts).
- **Provenance:** Each run captures its own source/config/tests/dependency, environment, Git identity (if present), results and SHA-256 manifest. Existing Git records and earlier frozen versions are not modified.
- **Scientific boundary:** Linear uniaxial elasticity and far-field Newtonian monopole bookkeeping are not a complete stress-energy tensor or general relativity model.
- **Status:** Development acceptance pending clean-extraction tests, numerical cross-checks and user-side execution.
- **Failures:** Record any development failures below; do not delete or conceal them.

## Development test failure preserved

- **R1 initial developer run:** 47 tests executed, 4 assertions failed. Original stderr/console output preserved unchanged in `GED_V1_R5_PRECISION_FAILURE_A.log`.
- **Cause:** several independent verification calculations used Decimal's default 28-digit context while the model used 60 digits; one inverse-square assertion used exact equality after a rounded multiplication.
- **Correction:** test calculations now use explicit 60-digit local Decimal contexts; inverse-square comparison uses a tight relative tolerance. The model's physics formulas were **not changed** in response to this failure.
- **R2 developer run:** 47/47 tests PASS. Original successful developer log `test_dev_second.log` retained during development, not required in the distributed experiment inputs because official runs generate `tests.txt`.
- **Acceptance:** development PASS; user-side and independent audit still pending.
