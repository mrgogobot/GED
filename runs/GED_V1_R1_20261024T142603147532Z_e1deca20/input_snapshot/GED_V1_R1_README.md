# GED-V1-R1 — Gravitational Energy Distribution

**Classification:** conventional Newtonian null + *deliberately hypothetical* information-coupling toy model. This is a numerical scout, **not** a detection of a new force or validation of information-carrying gravity.

## Experiment question

Do two idealized spherical objects with equal total mass and outer radius, but with their iron/ice layers swapped, have different **external** Newtonian gravitational fields? **No.** Their **internal** gravitational fields and **material responses** do differ.

Two bodies each have **R = 100 m**, half their volume ice and half iron:

- **A:** iron core, ice shell.
- **B:** ice core, iron shell.

Each has the same total mass. Because of spherical symmetry, external Newtonian gravity is `g = GM/r²`, regardless of layer order. The exterior spacetime of static spherical bodies in general relativity likewise depends on their matched gravitational mass rather than material labels; a fully relativistic mass match must account for stresses and all energies.

## Material-response proxy

Given a **uniformly heated** ideal material initially below its melting temperature, the energy per kilogram to bring it to the melting point and melt it fully is:

`Q_melt = cp * (Tm - T0) + L_f` in **J/kg**.

This describes an idealized thermal threshold. It is **not total atomic binding energy, nuclear binding energy or fracture strength**. The calculation also assumes a pressure regime permitting melting: ice in vacuum may **sublimate** instead. We do not attempt to predict actual impacts, object destruction, or thermodynamic gravitational mass changes.

## Explicitly hypothetical external coupling

We *invent* a test formula to give the GED idea a distinguishable prediction:

`g_TOY(r) = (G*M/r²) * (1 + alpha*chi_shell*(R/r)^p)` for `r >= R`.

`alpha` is an adjustable **dimensionless** hypothetical coupling. `p` is the chosen radial exponent. `chi_shell` is a dimensionless proxy obtained by normalizing each shell's heat-to-full-melt threshold; here **iron=+1, ice=−1**. The choice to use this proxy and radial law is purely a modelling decision, **not derived from the proposed binding-energy mechanism or from physics**.

Setting `alpha=0` restores the conventional null exactly. `alpha=1e-6`, `p=2`, `r=2R` gives an invented between-body field difference of **0.5 parts per million** of the conventional field. No device measured this; the program simply evaluates our stipulated equation.

**Unsolved:** an internally consistent field theory, conservation laws, relativity/covariance, actual binding-structure information channel, existing constraints, and any observed measurable effect. A null observation would constrain only the *chosen* toy equation, not all interpretations of gravity as information.

## Run (Python 3.13 standard library only)

No `pip install` required. Use your already-activated GED virtual environment:

```bash
cd ~/Desktop/Cosmology/GED
source .venv/bin/activate
unzip -n ~/Downloads/GED_V1_R1_bundle.zip -d .
python -m unittest -v test_ged_v1_r1.py
python ged_v1_r1.py
```

`ged_v1_r1.py` runs the tests **again** against a newly copied input snapshot before emitting science results. Do **not** modify the test suite to obtain passing output; failures are important receipts.

## Official per-run record

Each execution creates `runs/GED_V1_R1_<UTC>_<suffix>/` holding:

- `input_snapshot/` — exact script, test suite, configuration, README and preserved development audit
- `tests.txt` — verification output, including failures
- `results.json` — machine-readable scientific outputs and verdict
- `exterior_samples.csv` — sample field values at multiple distances
- `report.txt` — succinct human-readable report
- `provenance.json` — source/output SHA-256, toolchain, platform, Git information if available
- `SHA256SUMS` — checksum manifest of all files in the run directory, except itself

From that run directory, check all receipts with:

```bash
shasum -a 256 -c SHA256SUMS
```

Use the existing GED authoritative working codebase. This package adds **version-specific** filenames, preserving GED-V0-R1 and its results. You may additionally use your own Turnip Run execution wrapper; none is assumed or silently invoked. The per-run input snapshot records exact versions even when no Git repository is present.

## Acceptance boundary

**Accepted if local verification passes:** consistent Newtonian exterior null, differing internal acceleration and material heating threshold. **Not accepted:** new gravitational information-carrying physics. Preserve source, configuration, test receipts, outputs and scientific failures under the Manifest of Destiny.
