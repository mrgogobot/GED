# GED — Gravitational Energy Distribution

Exploratory computational research by Malin Hess. This repository preserves numerical experiments, configurations, tests, run receipts and development failure logs. It does **not establish a new gravitational mechanism**, detect an information field, or provide observational evidence for modified gravity.

## Scientific scope

The experiments investigate which proposed structure, energy and entropy signatures are already explained by conventional mass-energy accounting. V1-R2 computes Newtonian spherical self-binding and weak-field energy bookkeeping. V1-R3 extends this to prescribed rearrangement and idealised thermal entropy ledgers. V1-R4 compares state-counting entropy and radial thermal energy redistribution through three numerical controls: different entropy with identical assumed gravitational sources; thermal relaxation at fixed total modeled energy; and reversed temperature gradients with equal entropy but different interior fields.

These are mathematical controls under explicit assumptions. Geometry and material properties are idealised. Pressure, mechanical dynamics, thermal expansion, full stress-energy and a self-consistent general-relativistic solution are outside the models. No independent entropy-gravity coupling is derived or measured. Numerical precision is not measurement sensitivity.

## Contents and evidence policy

- `ged_*.py`, `test_*.py`, and versioned JSON configurations: historical experiment implementations and tests.
- Versioned READMEs and `DEV_AUDIT` documents: preserved historical development descriptions; their audit wording is not a certification of independent review.
- `runs/`: original saved receipts and exact input snapshots, with their original relative paths.
- `*PRECISION_FAILURE_A.log`: deliberately retained development failures.
- Original audit ZIPs: preserved byte for byte.
- `verification/`: new integrity checks and fresh snapshot test logs, separate from the historical evidence.

Original evidence is preserved without rewriting embedded paths or historical statements. Receipts include the local username and absolute execution paths. Virtual environments, IDE settings, caches, macOS metadata, the empty Goblin placeholder and the unrelated launcher are excluded. Directory layout is retained so historical imports, relative paths and manifests remain usable.

## Verification performed on 24 October 2026

All 44 entries in the V1-R1 through V1-R4 run manifests match their SHA-256 digests. Recorded input/output hashes also match; top-level copies match the corresponding saved snapshot files. Fresh execution of the saved snapshot suites passed: V1-R1 20 tests, V1-R2 26, V1-R3 46, V1-R4 55.

V1-R2 source SHA-256 is `c1c8bd268fcfad2f8dea23bba703ef04284ced2590680584ccfa263cd9010813`, matching the dependency pinned by V1-R3 and V1-R4. The V1-R4 audit ZIP matches the digest recorded in the earlier conversation: `d165ce79eeefe4481f029ee15b5cce7588189307f95421c3d3ef1ae7a8e4183d`.

This is a local integrity and test-reproducibility check, **not independent peer review or an independent scientific audit**. Hashes do not authenticate original authorship, timestamps or scientific correctness. Historical Git provenance fields are null and remain null. V0-R1 lacks a comparable run SHA256SUMS manifest; its preserved record has a weaker evidence trail. Fresh tests do not constitute a fresh reproduction of every result file.

## Run locally

Python standard library only. Historical runs record Python 3.13.17; the current verification runtime is recorded in `verification/integrity.json`.

```sh
python3 -B -m unittest -v test_ged_v1_r2.py test_ged_v1_r3.py test_ged_v1_r4.py
python3 -B ged_v1_r4.py
```

Run from a disposable copy if you want to leave the release tree unchanged: runners create new directories under `runs/`. Keep the exact pinned V1-R2 module alongside V1-R3 and V1-R4. Consult each versioned README for its assumptions and configuration.

## Publication status

Prepared release candidate: `ged-v1-r4`. Public release and Zenodo deposit require author approval. No DOI is claimed until actually issued. See `PUBLICATION.json` for publication tracking. Licensing is pending author choice; no new reuse license is granted by this preparation.
