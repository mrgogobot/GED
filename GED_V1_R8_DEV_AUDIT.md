# GED-V1-R8 development and verification audit

Date: 2026-10-24 (UTC). Scientific scope: synthetic conventional Newtonian gravimetry, exact shell-theorem non-uniqueness, seeded restricted three-class inversion and noise sensitivity.

## Design decisions

- Three equal-energy, zero-dipole arrangements from GED-V1-R7 (pair X, pair Y, cube) used as **preselected source templates**; no source-parameter fitting or arbitrary internal reconstruction claimed.
- Eight cube points share 100,000,000 J exactly as eight 12,500,000 J parcels (floating-point representation of the number is exact); pair sources each have 50,000,000 J. All share the same total notional mass `E/c²`.
- Measurements are 24 exterior positions; the 72 Cartesian acceleration components are tested as a single feature vector. No attempt to claim current observational feasibility.
- Two ideal uniform spherical shells at 30 m vs 70 m with identical mass demonstrate **exact**, not merely noise-limited, exterior non-uniqueness. An interior probe at 50 m serves as an internal positive control.
- Pseudorandom `random.Random(seed).gauss(...)` provides repeatable numerical Gaussian noise. The seed is frozen in JSON and results contain confusion matrices.
- Source standalone under `ged_v1_r8.py` rather than rewriting or importing GED-V1-R7. Only built-in Python libraries required. Previous accepted frozen experiments remain untouched.

## Checks and limitations

- Independent one-axis analytic formulas for Pair X, Pair Y and Cube used in separate verification tests; checks also cover orientation, symmetries, the far-field monopole limit and exact spherical shell theorem.
- All control inputs validated, including non-finite or boolean numerics and unsorted/inside probe radii.
- Development tests initially passed without a scientific-equation correction. Package metadata is finalized after clean-extraction verification; no failures are concealed or assumed.
- A successful source classification is only conditional on a known three-class model family. Full gravity inversions are non-unique; observable gravity generally does not encode complete material/atomic information.
- Third-party dependencies: none. Environment: Python standard library. Git provenance captured if available; no Git changes/publication performed.

## Publication rule

LOCAL HOLD through GED-V1-R10 unless the user chooses an earlier independently vetted definite finding. This scout's mathematical non-uniqueness is established physics, not a new discovery.
