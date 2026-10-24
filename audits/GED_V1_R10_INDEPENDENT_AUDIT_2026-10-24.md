# GED-V1-R10 — Independent Audit Receipt

Date: 2026-10-24 (UTC). Scope: uploaded Mac run archive `GED_V1_R10_audit.zip`, run ID `GED_V1_R10_20261024T200445115754Z_375dfe2a`.

## Integrity

- Uploaded archive SHA-256: `3d54e6b50558c771857813dca65b08749d1ede848015b3de82638d98b08b1049`
- ZIP CRC/integrity: PASS (22 entries, including directory entries).
- Internal run `SHA256SUMS`: 12/12 entries verified against extracted bytes.
- Seven working source/configuration/document files: identical to their frozen input snapshots.
- Saved Mac environment: Python 3.13.17, macOS arm64; local test exit code 0.
- Independent rerun environment: Python 3.13.5, Linux; 69/69 tests PASS; experiment program exit code 0.
- Independent clean rerun: `results.json` and `radial_samples.csv` byte-identical to saved Mac output.

## Independent scientific calculation

For spherical shells with fixed boundaries 0, 20, 50 and 100 metres, Newton's shell theorem gives
`g(r) = G * M_enclosed(r) / r**2`, with a uniform-density layer contributing
`M_shell * (r**3 - lower**3) / (upper**3 - lower**3)` inside its radial interval.

- Original shell masses A: 2e9, 3e9, 5e9 kg. Total = 1e10 kg.
- At r = 10 m: A encloses 2.5e8 kg; B encloses 3.75e8 kg. Therefore
  `g_A - g_B = -0.00008342875 m/s^2`, independently calculated.
- At r = 150 m: each has `g = 0.0000296635555555... m/s^2` under the stipulated model.
- All six sampled radii checked using independently calculated exact-rational enclosed masses; saved accelerations agree to the reported numerical precision.
- Exterior inverse matrix rows are each `[1, 1, 1]`: rank 1 of 3.
- With interior probe radii 10, 35, 75 m, the known-boundary response matrix is triangular with diagonal `1/8`, `31/104`, `19/56`; determinant `589/46592` is nonzero: rank 3 of 3.
- A and C have explicitly identical shell source masses with different *synthetic* binding/rest-energy ledgers, so their Newtonian fields agree by construction.

## Accepted scope

**Accepted:** Integrity-verified, independently reproduced and independently equation-checked *numerical Newtonian spherical-layer/inverse-problem scout*.

**Not established:** Full atomic CCMEI recovery, a new gravitational information/material-instruction field, physical realizability of matched A/C inventories, or full GR stress-energy and pressure treatment. The prior-evidence ledger is a retrospective status summary, **not** an independent audit of R5–R9 on-Mac run directories.

## Provenance and release gate

The Mac run record reports `git.commit = null`, so no commit maps exactly to this snapshot. The input hashes and this independent audit supply byte-level provenance; add a Git release tag only after the whole checkpoint is verified. R5–R9 archives remain to be independently inspected prior to a publication-wide provenance PASS. This audit neither commits nor publishes anything.
