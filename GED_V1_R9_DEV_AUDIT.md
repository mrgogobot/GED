# GED-V1-R9 — Development and Numerical Audit

**Date:** 2026-10-24 | **Status:** developer computational verification; independent user-machine run pending.

## Design and claims

- Objective: infer six unknown nonnegative energy fractions at six assumed-known sites from external vector-gravity measurements. This is a restricted finite-dimensional Newtonian inverse problem.
- Source code: self-contained Python standard library. Physics: Newtonian point-equivalent mass-energy parcels and the ideal spherical shell theorem.
- Numerical scaling: `g_ref = G E/(c² r_ref²)`; solve dimensionless linear algebra rather than comparing tiny SI accelerations directly.
- Publication hold: no GitHub commits/pushes or Zenodo submission until GED-V1-R10 unless there is a separately vetted definite finding and explicit approval.

## Analytical and numerical controls

- **One exterior vector probe + exact total energy:** augmented rank 4/6, nullity 2. A second valid six-source allocation differs by L1=0.250044562244 while yielding the same three measured components within `2.8e-17 g_ref` (floating-point residual).
- **Six x-axis probes, scalar x-component only + exact total energy:** augmented rank 3/6, nullity 3. Identical y/z source columns provide a direct symmetry explanation.
- **Eight distributed 3-component vector probes + total:** augmented rank 6/6, with normalized matrix singular values independently evaluated by NumPy 2.3.5 as approximately `[4.97264, 0.83723, 0.65145, 0.61320, 0.23090, 0.15447]`. This indicates full rank for the **specified** basis and geometry, not universally recoverable structure.
- **Noiseless reconstruction:** active-set finite-dimensional inverse solver returns source shares with L1≈`1.11e-15`. Independent NumPy `lstsq` reconstruction returns the same shares to `5.6e-16` max absolute discrepancy.
- **Exact uniform concentric thin-shell degeneracy:** different shell radii and allocations but identical total mass have analytically identical Newtonian exterior fields. At radius 35m they differ, because the probe lies between certain shells.
- **Noise experiment:** fixed RNG seed `20261024`, 100 independent synthetic trials per level and injected Gaussian noise relative to `g_ref`; errors grow with noise but are contingent on the source basis, probe geometry and constraints. Simulated sensitivities are not instrument specifications.

## Verification

- Python tested in container: **3.13.5**. Intended user local interpreter: **3.13.17**.
- `python -m unittest -v test_ged_v1_r9.py`: **88 tests, 88 passed**; no failure generated in the development run.
- `python ged_v1_r9.py`: executes, writes separate results, config/source snapshots, test output, provenance and SHA-256 manifest.
- Independent external-to-main-solver cross-check: NumPy 2.3.5 SVD ranks match 4/6, 3/6, 6/6, and ordinary least squares recovers the same noiseless fractions. NumPy is **not a dependency** of GED-V1-R9.
- Fresh ZIP extraction and checksum checks: see final packaging verification; file checks are recorded by the run itself.

## Limitations / not earned

- No measurement, derivation of additional gravity, information-field channel, atomic-binding inference, continuum stress tensor, or full Einstein equation solution.
- Point-equivalent energy parcels are not a mechanically equilibrated physical asteroid. The full body's background gravity is removed perfectly by fiat. Six source positions and total energy are stipulated known.
- The inverse uses normal equations and active-set enumeration suitable for a tiny controlled problem; a production inverse would require numerical singular-value analysis, regularisation and uncertainty handling.
- R9 results are a *numerical illustration of established gravity inverse problems*, not new scientific observations or discoveries.
