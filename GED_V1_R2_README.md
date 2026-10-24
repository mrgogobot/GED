# GED-V1-R2 — Gravitational Binding and Mass-Energy Accounting

**Classification:** Established-physics, weak-field numerical scout. **No observation, no detection of a new field, no full GR calculation.**

## Research question

When the constituents of a physical object rearrange and the binding energy changes, must the object's external gravitational influence change? The crucial control is **whether energy crosses the boundary of the system**.

Continuing V1-R1, compare two spherical bodies of outer radius 100 m, the same quantities of iron and ice, and equal constituent mass inventories:

- **A:** iron core, ice shell.
- **B:** ice core, iron shell.

Core and shell each occupy half the volume. Rearranging them changes the **Newtonian gravitational self-binding potential**, even when their zeroth-order Newtonian external accelerations agree. All matter quantities are idealized mass inventories, *not* measured relativistic gravitating masses. GR would also account for stresses, other energies, and energy transfer mechanisms; this is only a weak-field energy-equivalent illustration.

## Mathematics and units

- `M0 = rho_iron V_iron + rho_ice V_ice` (kg)
- `U = -integral_0^R [G * M_enclosed(r) / r] dM(r)` (J)
- `delta_M = delta_E / c^2` (kg)
- `delta_g(r) = G * delta_M / r^2` (m/s²; fixed-radius, weak-field approximation)

`U` here is specifically **gravitational self-binding**. It is *not* nuclear binding, chemical bonding, fracture strength, melting threshold, or a complete material-response descriptor.

Energy transfer cases:

1. **Absorbed heating retained:** total body energy increases by input Q; mass equivalent +Q/c².
2. **Equivalent radiation emitted:** body energy decreases; mass equivalent −Q/c².
3. **B → A rearrangement with release escaping:** gravitational potential energy decreases; body mass-equivalent decreases by `(U_A-U_B)/c²`.
4. **B → A rearrangement with release retained as heat:** negative potential-energy change is compensated by positive internal thermal energy; net mass-energy change zero.

A closed system including radiation and external reservoirs conserves energy. Real rearrangement mechanisms can require external work; those are *not modeled*. Symmetric/quasistatic energy transfer is assumed, so recoil is neglected. Fixed-radius fields describe model snapshots, not dynamical field propagation.

**High precision warning:** mass changes of order 10⁻¹⁰ kg added to masses of order 10¹⁰ kg are lost by ordinary float arithmetic. The code uses Python `decimal` arithmetic (60-digit context). This prevents a numerical artefact; it does not imply the change is measurable.

## Run on Mac (Python 3.13.17, stdlib only)

From `~/Desktop/Cosmology/GED` with `.venv` active:

```bash
unzip -n ~/Downloads/GED_V1_R2_bundle.zip -d .
python -m unittest -v test_ged_v1_r2.py
python ged_v1_r2.py
```

The source code has extensive `# WHAT / WHY / LIMIT` commentary. The new `runs/GED_V1_R2_<UTC>_<id>` folder contains:

- frozen input source, configuration, README and development audit (`input_snapshot/`)
- `tests.txt` — independent unit and numerical quadrature test results
- `results.json` — full-precision quantities in decimal **strings**
- `energy_field_samples.csv` — per-case radial field differences
- `report.txt` — concise readable summary
- `provenance.json` — commands, interpreter, optional Git identity and SHA-256s
- `SHA256SUMS` — complete manifest of all run files

To verify the receipt hashes: `cd runs/<run-id> && shasum -a 256 -c SHA256SUMS`.

No prior experiment or `main.py` is modified. A failed regression suite prevents normal results from being emitted. The separate Turnip Run v1 execution wrapper is **not** used; the standalone runner records equivalent essential execution receipts. Git commit/tag is recorded only if actually present. Versions without Git can be identified using the individual file hashes.

## Interpretation and limitations

Conventional physics already predicts that aggregate internal energy contributes to gravitating mass-energy. Rearrangement affects the distribution and binding energy; whether it changes the **body's net total energy** depends on the full closed/open energy budget. Gravity does not thereby acquire a demonstrated additional instruction channel. **GED's new material-information mechanism remains untested.**

This scout uses ideal spheres, fixed radius, constant density, Newtonian potential energy and a first-order mass-equivalent correction. There is no GR solution, elastic stress calculation, chemistry, material rupture, measurement error budget, or real detection.

## Source and version links

- Precursor: `GED_V1_R1_20261024T142603147532Z_e1deca20`.
- User-supplied V1-R1 audit ZIP SHA-256: `8c85d0a3c2496ee43e62d6ac3b19ca825409e16cf41fe79786325796cd4ee091`.
- NIST CODATA 2022 values for c and G: https://physics.nist.gov/constants
