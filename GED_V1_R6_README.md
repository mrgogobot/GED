# GED-V1-R6 — Radial Stress-Energy Distribution Scout

Turnip Research Laboratory · numerical physics control · 24 October 2026

**Status:** exploratory conventional-physics model. **No gravitational information channel is derived, observed or confirmed.** R6 is local only; next planned GitHub/Zenodo publication milestone is GED-V1-R10, or an earlier definite finding after separate review and approval.

## Research question

If two spherically symmetric objects have *identical material inventory and total retained elastic energy* but put that energy at different radii, does a conventional Newtonian field calculation already predict a difference? **Yes, inside the idealized object. Outside, the two fields agree.** That is the control needed before suggesting any new information-field mechanism.

## Model

An illustrative 100-m-radius sphere has an iron core of radius 50 m, surrounded by an ice shell. All reference properties and quantities use SI units. Two independent **constructed profiles** each retain `1e8 J` of elastic energy: `core_weighted` puts 75% in the iron core; `shell_weighted` puts 25% there. The rest is in the shell. The underlying core/shell composition, size and rest masses are identical.

- `E = 208.2 GPa` iron, `E = 9.17 GPa` ice: illustrative moduli inherited from R5; not temperature- or provenance-matched asteroid measurements.
- For a uniform elastic energy density zone of volume `V`, `u=U/V` in J/m³, which is equivalent to Pa.
- Ideal Hookean strain from assigned energy: `epsilon=sqrt(2U/(E V))` (dimensionless), using the SHA-pinned R5 helper.
- For each spherical zone `a < r < b`, enclosed fraction is `(r^3-a^3)/(b^3-a^3)`; outside it is one and inside its inner boundary zero.
- Energy-equivalent enclosed mass: `delta M(<r)=U(<r)/c²` kg.
- Newtonian radial acceleration correction: `delta g_proxy(r)=G U(<r)/(c² r²)` m/s² for `r>0`; zero at the centre.

**Physical caveat:** These assigned elastic energy profiles are a *bookkeeping construction*. They do not solve elastic displacement compatibility, tensor stress equilibrium, boundary forces, pressure, evolving geometry, gravitational self-energy or Einstein's field equations. Uniform scalar energy density is not a complete stress-energy tensor, and a real body's actual external field cannot be predicted from these terms alone. Even the associated supports and work reservoirs have their own stress-energy.

## Contents

- `ged_v1_r6.py`: heavily commented scientific computation and receipt writer.
- `test_ged_v1_r6.py`: boundary-value, analytic and provenance tests.
- `GED_V1_R6_config.json`: machine-readable frozen plan, assumptions and limits.
- `GED_V1_R6_DEV_AUDIT.md` and `GED_V1_R6_FAILURE_A.log`: developer verification notes and preserved first-run failures.
- `ged_v1_r2.py` (**already in GED, not duplicated in ZIP**): SHA-256 `c1c8bd268fcfad2f8dea23bba703ef04284ced2590680584ccfa263cd9010813`.
- `ged_v1_r5.py` (**already in GED, not duplicated in ZIP**): SHA-256 `608b67e91756e26667818743cdb906d487eeae65ca3d3f1de350bd87af50749d`.

R6 checks both exact source hashes at run time. If the working R2 or R5 code has changed, it fails rather than silently accepting a new code version.

## Execute using the existing project Python 3.13 venv

```bash
cd ~/Desktop/Cosmology/GED
source .venv/bin/activate
unzip -n ~/Downloads/GED_V1_R6_bundle.zip -d .
python -m unittest -v test_ged_v1_r6.py
python ged_v1_r6.py
```

No new pip packages required. Each execution writes `runs/GED_V1_R6_<UTC>_<uniqueID>/` containing snapshots of **both** pinned dependencies, `tests.txt`, `results.json`, `radial_energy_gravity.csv`, `report.txt`, `provenance.json` and `SHA256SUMS`. The snapshot test run is executed separately from the parent working files.

Local test success is not external scientific validation. Preserve numerical failures, code snapshots, per-run checksums and the exact Git state. Please do not commit/push/release R6 to GitHub or Zenodo before the agreed R10 milestone or a separately approved definitive result.
