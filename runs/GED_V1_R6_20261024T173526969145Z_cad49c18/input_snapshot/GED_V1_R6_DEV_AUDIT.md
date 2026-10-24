# GED-V1-R6 — Development and Scientific Audit

**Scope:** analytic Newtonian radial mass-energy proxy for two constructed loading profiles of the same iron-core/ice-shell sphere. Neither mechanical equilibrium nor full GR is solved.

## Preserved failures

- Initial run: the end-to-end clean-directory test recursively invoked itself from its own snapshot. Parent execution timed out. Captured incomplete original output is in `GED_V1_R6_FAILURE_A.log`.
- Second run: recursion guard succeeded, but two test assertions evaluated at Python's default 28-digit Decimal precision; no physics equation changed. A further assertion that attempted equality of independently rounded geometrical volumes was made relative-tolerance-based within 70-digit arithmetic. See the original second-run failure transcript in the same log.
- Additionally, the initial analytical test expectation `(75³-50³)/(100³-50³) = 297/875` was **arithmetically incorrect**; corrected to `19/56`. This was a test-reference error, not a model correction. Important: do not reclassify it as a mere precision issue.

## Final corrected tests

- Fixed the test-only recursive invocation by marking the inner snapshot execution and skipping the outer clean-room test only in that inner environment.
- Isolated independent Decimal arithmetic in 70-digit local contexts; used a physically meaningful relative tolerance where a rational volume relation is rounded.
- Corrected the independently derived analytic radial-shell fraction.
- Developer final test command: `python -m unittest -v test_ged_v1_r6.py`.
- **51 tests PASS**, including standalone clean-directory output/provenance verification; snapshot rerun deliberately reports one skipped meta-test to prevent infinite nesting.
- Every calculation is a model calculation, not an observation. Pinned R2 and R5 dependency hashes are required.

## Acceptance status

Developer computation/clean-room PASS; user-local reproduction and independent audit pending. No scientific evidence for extra gravitational information or entropy coupling. **Local only until GED-V1-R10 or a separately confirmed definitive finding; no GitHub push or Zenodo submission.**
