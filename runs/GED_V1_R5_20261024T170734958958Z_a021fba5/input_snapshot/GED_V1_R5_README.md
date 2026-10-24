# GED-V1-R5 — Stress-Energy and Material Response

**Turnip Research Laboratory · numerical scout · 2026-10-24**

This is a **conventional physics baseline** for the Gravitational Energy Distribution (GED) programme. **It is not experimental evidence of a gravitational information field, and it does not derive a new gravitational law.** It preserves the user's foundational suggestion that the collective properties of matter, including binding relationships, might help explain interactions, while testing conventional mechanisms first.

## Central question

Do materials with different elastic stiffness respond differently to identical loading conditions? Yes, according to Hooke's law. Does this difference alone demand an additional gravitational information channel? No: elastic energy and its transfer can already be accounted for using conventional energy conservation and mass-energy equivalence.

## Models and exact units

Uniaxial infinitesimal elastic deformation:

- **Stress** `sigma = E * epsilon` (Pa = N/m²).
- **Strain** `epsilon` (dimensionless).
- **Energy density** `u = (1/2) * sigma * epsilon = sigma²/(2E) = E*epsilon²/2` (J/m³ = Pa).
- **Stored energy** `U = u * V` (J).
- **Mass-equivalent for energy kept by the sample** `delta_m = U/c²` (kg).
- **Far-field weak-source monopole bookkeeping** `delta_g_proxy = G * U/(c²*r²)` (m/s²).

The last formula is deliberately labelled a **proxy**. A loaded uniaxial body is not spherically symmetric. Its supports, geometry changes, anisotropic stresses and full general-relativistic gravitational field are **not** solved. In a closed system, work transferred to a sample is lost from the reservoir that did the loading; total energy is conserved.

**Three controls:** (A) equal imposed stress, (B) equal imposed strain, and (C) equal stored elastic energy. Control C is especially helpful: iron and ice require different stresses and strains to store the same energy, but the *incremental* monopole mass-energy proxy is the same. Their baseline masses and baseline gravitational fields are **not** equal because the densities differ.

## Source values and assumptions

- Polycrystalline iron: **208.2 GPa** illustrative Young's modulus, from Takaki, Masumura & Tsuchiyama (2020), DOI [10.2355/tetsutohagane.TETSU-2019-129](https://doi.org/10.2355/tetsutohagane.TETSU-2019-129).
- Polycrystalline ice: **9.17 GPa** approximate dynamic Young's modulus near −5°C, compiled in [McCarthy, Elastic Properties, Caltech KISS workshop](https://www.kiss.caltech.edu/workshops/oceanworlds/presentations/McCarthy.pdf), citing Ewing et al. (1934).
- Mass density values **7,870** and **917 kg/m³** are nominal illustrative inputs rather than measured samples.
- CODATA 2022: `c = 299792458 m/s` (exact); `G = 6.67430e-11 m³/(kg s²)` (measured, uncertain); see [NIST constants](https://physics.nist.gov/cuu/pdf/all.pdf).
- Nominal reference sample volume **1 m³**; equal load **100,000 Pa**; equal strain **1e-5**; target elastic energy **1 J**; probe distance **100 m**.
- Purely elastic, isotropic, uniaxial, infinitesimal, reversible response: **no fracture, failure, creep, plasticity, porosity, thermal expansion, melting, or chemical/nuclear binding energy.** The two modulus references are from different sample conditions and serve as illustrative comparisons, not temperature-matched measurements.
- Densities already reflect real material mass-energy; do **not** add microscopic binding energy a second time.
- The figures are not predicted measurable gravity signals. Real apparatus stress-energy and mass distributions require explicit full boundary conditions.
- No entropy number is given for dissipation because its value requires actual temperatures and pathway.

## Files and authoritative dependency

- `ged_v1_r5.py` — main model, CLI and automatic provenance collector; WHAT/WHY code comments.
- `test_ged_v1_r5.py` — independent formula and conservation tests.
- `GED_V1_R5_config.json` — fixed input plan and limitations.
- `GED_V1_R5_DEV_AUDIT.md` and `GED_V1_R5_PRECISION_FAILURE_A.log` — development record and original failed test output.
- `ged_v1_r2.py` — **existing authoritative source** for gravitational mass equivalent and delta-g, **not bundled a second time**. Required SHA-256: `c1c8bd268fcfad2f8dea23bba703ef04284ced2590680584ccfa263cd9010813`.

## Run on the GED Python environment

```bash
cd ~/Desktop/Cosmology/GED
source .venv/bin/activate
unzip -n ~/Downloads/GED_V1_R5_bundle.zip -d .
python -m unittest -v test_ged_v1_r5.py
python ged_v1_r5.py
```

No third-party dependencies. Recommended Python 3.13, using existing project venv.

## Receipts and acceptance

Each execution produces `runs/GED_V1_R5_<UTC>_<ID>/` with `input_snapshot/`, `results.json`, `elastic_cases.csv`, `report.txt`, `tests.txt`, `provenance.json`, and `SHA256SUMS`. Frozen snapshots include the exact pinned R2 source; files have SHA-256 hashes. The run records the Python runtime and local Git commit/tag/dirty status when available.

**Acceptance ladder:** developer clean run → user execution → upload audit ZIP → independent audit → scientific numerical scout acceptance → optional separately approved public release. **Do not publish GED-V1-R5 on GitHub or Zenodo until review and explicit approval.**
