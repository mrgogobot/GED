# GED V1-R10 publication candidate

Proposed tag: `ged-v1-r10`. Status: LOCAL CANDIDATE; publication requires explicit user approval. This supersedes no frozen experiment. Previous tag `ged-v1-r4` remains immutable.

## Contents and provenance

This checkpoint adds V1-R5 through R10 source, tests, configurations, original run inputs/outputs, SHA-256 receipts, failure logs, both original audit ZIPs, and both independent audit reports. The existing repository history and all earlier tracked evidence are retained. The authoritative publication codebase is https://github.com/mrgogobot/GED; this local checkout is a release preparation copy. The Desktop GED directory remains the original local evidence source and was not changed.

Historical runs recorded no exact Git commit for these experiments. Their original provenance remains unchanged; the candidate commit packages those exact bytes retrospectively. The final external candidate receipt records the commit and proposed tag without a self-referential commit hash inside the commit itself. No public tag has been created.

## Verification on 24 October 2026

Python 3.13.17, Clang 23.1.3, macOS 27 arm64; standard library only, no third-party dependencies. Complete runtime details and per-experiment results appear in `verification/r10-checkpoint/verification.json`.

| Experiment | Fresh tests | Original manifest entries | Matching numerical output files |
|---|---:|---:|---:|
| R5 | 47 | 12 | 2 |
| R6 | 51 | 13 | 2 |
| R7 | 62 | 12 | 2 |
| R8 | 71 | 10 | 2 |
| R9 | 88 | 11 | 3 |
| R10 | 69 | 12 | 2 |

All tests passed, all 70 manifest entries verified, and all 13 numerical output files matched byte for byte. Each working input matches its original snapshot; both audit ZIP hashes match their respective reports. The fresh executions used disposable snapshot copies. These checks confirm the reported counts; they do not repeat every additional independent equation check described in the preserved audit reports. Fresh test/execution logs are separate from original run receipts. Execution log paths refer to temporary verification locations, not new authoritative source locations.

## Accepted findings and limits

Accepted qualified numerical scouts reproducing conventional physics: elastic energy bookkeeping; spherical enclosed-energy fields; asymmetric Newtonian multipoles; restricted known-template classification; geometry-dependent inverse reconstruction and noise sensitivity; spherical-layer degeneracy and known-boundary interior inversion.

No observational detection, new gravitational information mechanism, atomic composition recovery, or full GR solution is established. R10 A/C binding ledgers are synthetic bookkeeping controls, not demonstrated realizable materials with identical complete stress-energy. Pressure, apparatus, dynamics, geometry assumptions and other limitations remain as documented in each original experiment. Numerical precision is not detector sensitivity. Development failures and their corrections remain archived; passing final tests does not erase those failures.

## Reproduction

After checking out the exact candidate commit (recorded in the external receipt), use Python 3.13.17. No package installation is required. Verify `RELEASE_SHA256SUMS` from the repository root with `shasum -a 256 -c RELEASE_SHA256SUMS`. Existing per-run manifests must be checked from their own run directories.

Copy the checkout to a disposable directory, excluding `.git`, before executing runners, which generate new `runs/` directories. From its root:

```sh
python3 -B -m unittest -v test_ged_v1_r5.py test_ged_v1_r6.py test_ged_v1_r7.py test_ged_v1_r8.py test_ged_v1_r9.py test_ged_v1_r10.py
python3 -B ged_v1_r5.py
python3 -B ged_v1_r6.py
python3 -B ged_v1_r7.py
python3 -B ged_v1_r8.py
python3 -B ged_v1_r9.py
python3 -B ged_v1_r10.py
```

R5/R6 require the unchanged R2 dependency already in the repository. Compare newly generated `results.json` and CSV files to original run paths in `verification.json`. Reports and provenance may contain new paths, times and Git state; preserve the originals.

## Publication plan

After approval, tag the approved commit `ged-v1-r10`, publish its GitHub release, and create a new version of the existing Zenodo record 23622641. Preserve the old DOI and use the new version DOI for this checkpoint. Verify uploaded archive hashes, release target, DOI resolution and reciprocal links before reporting completion. Existing `PUBLICATION.json` describes the earlier R4 publication only. Do not replace it with an unissued DOI or claim R10 is already public.

Creator and existing license scopes remain as recorded in `LICENSE.md`: Malin Hess; MIT for code, CC BY-NC 4.0 for research materials.
