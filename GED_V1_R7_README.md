# GED-V1-R7 — Asymmetric Energy Distribution / Gravitational Multipoles

Turnip Research Laboratory · local computational scout · 24 October 2026

**Scientific classification:** conventional Newtonian multipole calculation; **not** evidence for an additional gravitational-information channel or a novel interaction law.

## Research question

GED-V0 through GED-V1-R6 used spherically symmetric controls. Those have identical external Newtonian fields when total enclosed masses are equal. Does ordinary gravity distinguish different arrangements **outside** the object once spherical symmetry is removed? **Yes:** higher multipoles record geometric mass-energy distribution even when total energy (monopole) and energy centroid (dipole) match.

## Numerical construction

- Common background: 100-m-radius sphere, iron core radius 50 m, iron density 7,870 kg/m³, ice shell density 917 kg/m³. This is the same *material geometry* used in GED-V1-R6. Its external monopole is shared among all cases.
- Equal additional retained energy: `100000000 J` for each constructed configuration. The energy equivalent is `E/c²`, with `c=299792458 m/s`.
- `pair_x`: two equal parcels at `(±60, 0, 0) m`.
- `pair_y`: the identical pair rotated by 90 degrees onto the y axis.
- `cube_eight`: eight equal parcels at `(±60/√3, ±60/√3, ±60/√3) m`. Each lies approximately 60 m from the centre, inside the ice shell. The eight shares are exactly `12,500,000 J` each; dipole and traceless quadrupole cancel by cubic symmetry. The cube has higher multipoles, so its full field is not spherical.
- External sample positions: along x, y, z, at radii 200, 400, 800, 1600 and 3200 m.
- Models are idealized *fixed source arrangements*, not automatically achievable stable elastic load cases. Point parcels are mathematical proxies, not literal infinite energy density or a stress-free construction.

## Equations and units

With energy parcels `E_i` joules, position vectors `x_i` metres, `G` in m³/(kg s²):

- Equivalent source mass: `m_i = E_i / c²`, kg.
- Exact Newtonian added acceleration: `g(r) = -G Σ_i m_i (r-x_i)/|r-x_i|³`, m/s².
- Mass-energy monopole: `M = Σ_i m_i`, kg.
- Dipole: `D_i = Σ_a m_a x_(a,i)`, kg·m; all constructed cases have zero dipole.
- Trace-free Cartesian quadrupole: `Q_ij = Σ_a m_a [3 x_(a,i)x_(a,j) - |x_a|² δ_ij]`, kg·m².
- Quadrupole potential: `Φ_Q = -G (rᵀ Q r)/(2 |r|^5)`, m²/s².
- Quadrupole acceleration correction: `g_Q = G[(Q r)/r^5 - (5/2)(rᵀ Q r)r/r^7]`, m/s².

The code compares the *exact parcel sum* with monopole-plus-quadrupole and measures truncation error versus distance. The 70-digit Decimal arithmetic protects tiny differences from cancellation against the large common asteroid field.

## Scientific limitations

This is a **Newtonian weak-field bookkeeping proxy** in which `E/c²` acts as added source mass. It does **not** solve the full Einstein field equations, elastic compatibility, stress or pressure equilibrium, source motion, binding-energy changes, apparatus, or conservation laws during hypothetical rearrangements. External symmetry findings are exact for the stipulated point-mass Newtonian model only. The full gravity of a real loaded body would include all stress-energy and supporting structures.

An exact source difference is expected from standard gravitational multipole theory and is **not** an information field beyond standard gravity. The cube's zero quadrupole does not imply absence of all higher-order moments. No measurement, claimed detectability, observational discovery, or extra coupling is offered.

## Run in the user's existing GED Python 3.13 virtual environment

```bash
cd ~/Desktop/Cosmology/GED
source .venv/bin/activate
unzip -n ~/Downloads/GED_V1_R7_bundle.zip -d .
python -m unittest -v test_ged_v1_r7.py
python ged_v1_r7.py
```

No added pip dependencies. Each execution generates `runs/GED_V1_R7_<UTC>_<ID>/` with a snapshot of the seven input files, `tests.txt`, `field_samples.csv`, `results.json`, `report.txt`, `provenance.json`, and `SHA256SUMS`. The runner also executes 62 verification tests **from the frozen snapshot**, independently of the user's initial command.

The intended release milestone is **GED-V1-R10**, or an earlier definite finding reviewed and explicitly approved by the user. **Do not push/commit or publish GED-V1-R7** as part of this run; preserve local evidence instead. Older GED experiments are not modified.

## Development failures preserved

`GED_V1_R7_PRECISION_FAILURE_A.log` preserves four first-run Decimal exact-zero assertions, corrected with explicit physical-scale tolerances. `GED_V1_R7_REPORT_FORMAT_FAILURE_A.log` records a 67-times repeated header discovered in an otherwise passing run; corrected through explicit string concatenation. Neither involved modifying the physics equations.

## Preserved development failures

`GED_V1_R7_PRECISION_FAILURE_A.log` contains the original four exact-zero test failures; `GED_V1_R7_REPORT_FORMAT_FAILURE_A.log` preserves a malformed development report header. Both were corrected without changing physics equations.
