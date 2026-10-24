# GED-V1-R1 Development Audit — preserved machinery failure

An initial 20-test local run had **19 PASS / 1 FAIL**:
`test_radial_falloff` required a ratio computed from two differences of very
close floating-point numbers to equal exactly 4 to eight decimal places.
Observed `4.000000014122946`, theoretical `4`. This is finite-precision
subtraction error, not a physical discrepancy. Fixed the *test tolerance*
only, to `math.isclose(actual,4,rel_tol=1e-8)`; the scientific formula and
model values were unchanged. Preserved the initial failure in this note as
a development/machinery outcome rather than silently rewriting history.

All official output runs must pass the corrected test suite, which is
preserved byte-for-byte in their input snapshots.
