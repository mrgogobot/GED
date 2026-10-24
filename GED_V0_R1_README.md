# GED-V0-R1 — Gravitational Energy Distribution: controlled baseline

**Status: SCOUT ONLY — not a new-physics confirmation, not an accepted research result.**

## Original hypothesis to preserve

> Matter possesses physical properties arising from its constituent parts and binding energies. Gravitational fields might encode information about those collective properties, allowing physical interactions to follow from that information.

This is the research proposal, **not** what the baseline equations establish. The experiment first separates gravitational observables from material responses, then identifies what an additional information channel would have to predict.

## Files

- `ged_v0_r1.py` — heavily commented Python 3.13 standard-library experiment
- `GED_V0_R1_config.json` — all numerical input values and units
- `test_ged_v0_r1.py` — core analytic and numerical checks
- `GED_V0_R1_README.md` — purpose, assumptions, and execution protocol

Nothing overwrites an existing `main.py`. No external packages required.

## Install in your existing project

Extract this package into `~/Desktop/Cosmology/GED/` alongside your existing `main.py`, **not over it**. Activate your already-created Python 3.13 virtual environment:

```bash
cd ~/Desktop/Cosmology/GED
source .venv/bin/activate
python --version     # Expected: Python 3.13.17 on your system
python -m unittest -v test_ged_v0_r1.py
python ged_v0_r1.py
```

Each run creates a separate directory under `runs/GED_V0_R1_*` and writes:

- `results.json`: SI-unit-labelled numerical results
- `report.txt`: human-readable summary
- `provenance.json`: exact Python interpreter, invocation, UTC start timestamp, source/config SHA-256 and output SHA-256, and available Git commit

**For an official experiment:** use your existing **Turnip Run v1** wrapper and established Manifest of Destiny. This script's receipt is only a supplement; it cannot replace your official run receipt, source commit/tag, full environment and accepted/failed evidence policy. We did not guess your Turnip Run JSON plan schema. If Git hasn't been initialized or committed, the script reports missing commit information instead of inventing one.

## Three controlled comparisons

**Control A — Equal radius, different materials**

Iron: density 7,870 kg/m³. Ice: density 917 kg/m³. Both have radius 100 m. The model computes mass via `M = ρ(4πR³/3)` and surface gravity via `g = GM/R²`. Different **mass at the same radius** gives different surface gravity. Atomic weight alone does **not** determine a gravitational field.

**Control B — Equal radius and total mass, different arrangements**

Construct two exact spherical layered objects from identical volumes of ice and iron. In one, iron is inside; in the other, ice is inside. The two have equal total mass and radius; therefore their **external** Newtonian gravitational fields are exactly equal for the idealized static, spherically symmetric situation. Interior gravitational accelerations differ. This is a meaningful information-limit control, not proof that an alternative channel exists.

**Control C — Fragmentation and material response**

A pure-ice sphere is rearranged into 1,000 equal spherical fragments on a regular, spaced 10×10×10 grid. This is a *postulated configuration*, not a collision simulation. Total mass and material density are held constant. Each fragment has a gravitational field; near the group, spatial distribution changes the acceleration, while far away the group's field approaches that of its total mass concentrated at its centre.

Separately, uniform heating deposits 300,000 J/kg into initially 250 K iron or ice, at pressure permitting melting, with no heat loss. The calculation uses `Q/m = c_p(T_m-T_0) + fL_f` to estimate melt fraction `f`. It is **not** a strength, fracture, nuclear binding, or vacuum sublimation model. It is not a complete description of how asteroids behave under stress.

## Binding energy is not just a number describing stress

Binding energy is the **energy difference** between a bound configuration and appropriately separated constituents. Material strength, thermal phase transitions, mechanical fracture and nuclear binding are distinct physical quantities. They can be included in a fuller state model but must not be conflated. Gravity couples to mass-energy and stresses; it is not established to convey separate instructions for chemical or mechanical behaviour.

A body's Newtonian gravitational field has no strict outer boundary. If retained thermal energy is added, its mass-energy contribution is `Δm = Q/c²` (usually negligible here); ignoring it in Control C's simple mass conservation is an intentional approximation. Energy carried off by escaping fragments, radiation or vapor must be tracked in a more realistic model.

## Scientific decision gate

The script checks whether the **conventional model** is internally consistent. It does not compare a new information-field equation with data, because no additional equation or observable has yet been specified. A future GED-V1 proposal must predict a departure from conventional gravitational physics **without** merely relabelling known density, mass, binding energy, and electromagnetic/nuclear response.

## Provenance and version discipline

- Experiment ID: `GED-V0-R1`
- Code version: exactly the four files listed above as identified by the manifest's SHA-256 hashes
- Runtime: Python 3.13.17 intended; the actual runtime is recorded by each execution
- Third-party dependencies: none
- Units: SI, explicitly labelled in config and output
- Results: local generated files, each hashed, never overwritten
- Interpretation: **baseline/scout**, not hypothesis validated
- Prior output directories: never delete or silently replace; preserve failures

If you send back `results.json`, `report.txt`, `provenance.json` and test output, we can audit the run before proposing GED-V1.
