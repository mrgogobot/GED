# Turnip Research Laboratory — GED-V1-R5 through GED-V1-R9 Independent Audit

**Audit date:** 2026-10-24 (UTC)  
**Audit scope:** Uploaded local execution archive and its preserved source, inputs, outputs, tests, and integrity receipts.  
**Archive:** `GED_V1_R5_to_R9_audit.zip`  
**Archive SHA-256:** `c611f9f32d5b6fa01895e7f861cc076740eb78f5e22ce4ca86a660c97255916d`  
**Overall determination:** ACCEPTED as reproducible, qualified **computational scouts**, **not** evidence of a new gravitational information channel or observational detection. Local publication readiness is conditional on versioning and review of the complete R10 release candidate.

## 1. Integrity and provenance

- ZIP central-directory/CRC integrity: **PASS**.
- Five expected distinct local run directories: **PRESENT** (one each, R5–R9).
- Every run's `SHA256SUMS` entry successfully verifies against its preserved bytes: **58/58 PASS** (R5=12, R6=13, R7=12, R8=10, R9=11).
- Working files present in the archive are byte-identical to their corresponding frozen `input_snapshot/` copies, including tests and configurations: **PASS**.
- Additional hashes in provenance metadata (where explicitly supplied) agree with the referenced frozen inputs/outputs: **44/44 PASS** (R5=11, R6=12, R7=11, R8=5, R9=5).
- R5 and R6 each snapshot exact source dependency `ged_v1_r2.py`: SHA-256 `c1c8bd268fcfad2f8dea23bba703ef04284ced2590680584ccfa263cd9010813`; R6 also pins R5 source SHA-256 `608b67e91756e26667818743cdb906d487eeae65ca3d3f1de350bd87af50749d`. These match preserved copies.
- Earlier test failures and development corrections are retained as logs and notes for R5–R7. No failed-run log is asserted for R8–R9.
- **Limitation:** R5–R9 do not have per-experiment authoritative Git commits/tags in the local provenance records. The frozen byte-level source hashes do identify the exact executed code. The publication release must record the Git commit/tag and explicitly link it to these earlier local frozen snapshots without retroactively pretending the earlier experiments were committed when run.

## 2. Independent execution

The audit extracted an isolated working copy rather than modifying the user's source folder. Since R5/R6 depend on R2, the pinned `ged_v1_r2.py` was copied **unchanged** from the R5 frozen snapshot into the isolated working directory. The audit interpreter was **Python 3.13.5**, versus **Python 3.13.17** recorded for the user's Mac. The same standard-library-only source was used; an independent NumPy-based SVD calculation was used solely for an additional audit, not for any GED run.

| Experiment | User tests | Independently rerun | Full experiment | Scientific scope |
|---|---:|---:|---|---|
| R5 — Material stress and elasticity | 47 PASS | 47 PASS | PASS | Hookean elastic energy, `E/c²` monopole proxy |
| R6 — Radial stress-energy distribution | 51 PASS | 51 PASS | PASS | Newtonian spherical enclosed-energy proxy |
| R7 — Asymmetric energy multipoles | 62 PASS | 62 PASS | PASS | Newtonian point-source multipoles |
| R8 — Known-template information recovery | 71 PASS | 71 PASS | PASS | Synthetic 3-template classification and shell-theorem ambiguity |
| R9 — Inverse reconstruction limits | 88 PASS | 88 PASS | PASS | Fixed-site inverse problem and noise sensitivity |
| **Total** | **319 PASS** | **319 PASS** | **5/5 PASS** | **No observations or new interaction law** |

Reproduced `results.json` and numerical CSV files from all five runs are **byte-for-byte identical** to the saved outputs (11/11 numerical data files). Text reports are additionally identical for R5, R6, R8 and R9; R7's `report.txt` differs **only in the automatically printed local run-directory path**, which changes for each execution. R7's numerical content matches.

## 3. Independent physics and numerical spot-checks

No GED calculation module was imported for these spot-checks. The check equations were separately implemented from the configuration assumptions.

- **R5:** For a uniaxial linear-elastic sample of volume 1 m³, independent evaluations of `U = 0.5 × stress² × volume / Young_modulus` and `U = 0.5 × Young_modulus × strain² × volume` match the stored iron/ice energy values; `Δg = G U/(c² r²)` matches the far-field proxy. This is a *chosen-system* increment, not a complete model including its loading apparatus.
- **R6:** At `r=50 m`, independently computed difference is `G × (5.0×10⁷ J)/(c² × (50 m)²) = 1.4852320538237332×10⁻²³ m/s²`, matching the recorded result. At exterior `r=200 m` the stipulated equal total energies give zero difference.
- **R7:** Independent exact Newtonian expressions for an opposed pair and eight cube-corner point-energy parcels give pair-minus-cube radial field `−6.14160004271731×10⁻²⁵ m/s²` at +x=200 m, agreeing with the archived calculation. The quadrupole truncation error is approximately 18.38% at that distance.
- **R8:** The fixed-seed known-template sweep reproduces the saved classification records exactly: 600/600 noiseless classifications, 293/600 at noise sigma `10⁻²⁴ m/s²`, and 189/600 at sigma `10⁻²³ m/s²`. This is a restricted, simulated classifier, not general interior-structure recovery.
- **R9:** An independent design matrix with SVD gives rank **4/6** (one 3-component vector probe with total-energy normalization), **3/6** (same-axis scalar measurements plus normalization), and **6/6** (distributed vector probes plus normalization). Separate least-squares fitting recovers the specified noiseless six fractions with L1 error approximately `1.07×10⁻¹⁵`.

## 4. Earned and unearned findings

**Earned:** Each source reproduces its stipulated mathematics and expected outcomes, including material-dependent elastic responses, internal/external spherical field distinctions, asymmetric multipoles, exact exterior shell non-uniqueness, and inverse-reconstruction limits under known source geometries.

**Not earned:** Measurement of a real gravitational signal, a new gravity law, proof of an independent material-binding information channel, complete recovery of atomic composition, or a full stress-energy/Einstein-equation solution. The models intentionally omit many real-world boundary, stress, equilibrium, thermal, and/or detector constraints as documented in the configurations.

**Scientific verdict:** ACCEPTED QUALIFIED NUMERICAL SCOUTS (R5–R9). These are not confirmations or falsifications of GED's additional gravitational-information hypothesis. Exact degeneracy under an assumed observational model does not rule out unrelated, independently detectable channels.

## 5. Release checkpoint instructions

1. Keep previously published `ged-v1-r4` immutable. Organise and include accepted snapshots, results, checksum manifests, preserved failures, this external audit and the separately audited R10 evidence in the **next** release candidate.
2. Before publication, perform a whole-release manifest check and rerun the test suites from the exact release-tree snapshot. Record the final authoritative Git commit/tag, dependency/toolchain versions, and any differences in original/packaged paths. Do not claim earlier local experiments had Git commits if they did not.
3. The R10 local numerical audit was handled separately. **This report does not itself re-audit R10.** Refer to its separately issued report and hashes.
4. Publish GitHub and Zenodo only after the user's review and explicit approval of the release contents and scientific claims. Neither service was modified as part of this audit.

**Archive disposition:** Preserve the uploaded ZIP and its hash as the audited local evidence source. Do not silently alter previously frozen experiments.
