# GED-V1-R3 — Entropic Exhaust / Structural Persistence

**Status:** Numerical thermodynamic control scout; **not** an observational result, an entropy-driven gravitational law, or a detection of an extra information channel.

## Scientific question

Does a gravitationally bound system necessarily spend energy to remain organised, and does increased entropy require energy to escape the body?

**This model answers:** No. A stable idealised bound configuration requires no perpetual exhaust. During a prescribed reconfiguration, thermal entropy may increase even when all released gravitational potential energy stays inside the object. Conversely, when energy crosses the body boundary, the body's mass-energy drops by `E_export/c²`. Energy and entropy are separate accounting quantities.

## Calculation

**Same constituent inventory as GED-V1-R2:** ideal uniform-density iron and ice in two spherical half-volume layers, 100 m outer radius. B is ice core/iron shell; A is iron core/ice shell. We import V1-R2's already-audited formulas for Newtonian potential binding energies, and compute the available decrease `Q=U_B−U_A > 0`.

We assume a *prescribed* B→A rearrangement, with no dynamic simulation and no guarantee that it occurs spontaneously. The model divides released potential energy between retained internal heat `Q_body=(1−f)Q` and energy deposited into an infinite thermal bath `Q_bath=fQ`, where `0≤f≤1`.

**Energy ledger** (joules, J):

- `ΔU = −Q` — gravitational potential energy decreases.
- `ΔE_internal = (1−f)Q` — internal thermal energy rises.
- `ΔE_outside = fQ` — surroundings receive energy.
- `ΔU + ΔE_internal + ΔE_outside = 0` — closed total ledger.

**Thermal entropy ledger** (joules per kelvin, J/K), under constant positive heat capacity `C` in J/K and an ideal equilibrium bath at `T_bath` in K:

- `T_final = T_initial + Q_body/C` in K.
- `ΔS_body = C ln(T_final/T_initial)` in J/K.
- `ΔS_bath = Q_bath / T_bath` in J/K (AFTER full thermalisation in the bath).
- Only *modeled thermal* entropy is represented; structural/configurational entropy of rearrangement is **not** calculated.

**Gravitational mass change** (kg) of the body *only* when energy is exported:

- `ΔM_body = −Q_bath / c²`
- `Δg = G ΔM_body / r²` at fixed external probe radius `r`.
- Radiation and the environment still carry the exported energy and therefore continue to gravitate as part of the larger system. No energy or gravity "vanishes."

The weak-field mass correction is too tiny for ordinary asteroid measurements and is **not a GR field-equation solution**. The input constants are simplified bulk approximations. The tested binding is gravitational self-binding, **not chemical, nuclear, or fracture binding**.

## Environment

- Python **3.13.17** on your GED Mac `.venv` is appropriate; compatible Python 3.13 used for development verification.
- No third-party pip installations.
- **Requires the existing `ged_v1_r2.py` alongside this R3 source**, pinned to SHA-256 `c1c8bd268fcfad2f8dea23bba703ef04284ced2590680584ccfa263cd9010813`. The R3 runner refuses a mismatched predecessor; this is intentional to prevent version drift.
- R3 files do not overwrite the existing `main.py` or earlier experiments.

## Local instructions

```bash
cd ~/Desktop/Cosmology/GED
source .venv/bin/activate
unzip -n ~/Downloads/GED_V1_R3_bundle.zip -d .
python -m unittest -v test_ged_v1_r3.py
python ged_v1_r3.py
```

The runner reruns its verification tests using a frozen copy of the source and records:

```text
runs/GED_V1_R3_<UTC>_<ID>/
  input_snapshot/                     source, config, README, development audit, failed test log and pinned R2
  tests.txt                            local verification results
  results.json                         precise model numbers as strings
  entropy_energy_ledger.csv           stable and energy-partition cases
  report.txt                           human-readable summary
  provenance.json                      Python, Git if available, commands and SHA-256 receipts
  SHA256SUMS                           hashes for all included files except itself
```

This is a project-specific receipt generator, **not a claim that `turnip_run.py` was invoked**. If the user wants to integrate this with Turnip Run later, first audit the existing authoritative wrapper and JSON plan structure; do not create parallel provenance frameworks silently.

## Acceptance boundaries

- **Can pass:** code tests; conservation of energy; stable baseline; temperature/entropy response under toy thermal assumptions; tiny weak-field `Δg` from energy export.
- **Cannot establish:** whether actual asteroids spontaneously rearrange, whether they shed entropy continuously, whether gravity transmits mechanical/material-response instructions, or whether any extra gravitational information channel exists.
- This should be preserved as a **conventional physics control**, not a new-science result.

## Provenance

R2 baseline values supplied by the preceding GED-V1-R2 run, 2026-10-24: gravitational binding release `3.971322106352E+7 J`. Direct R3 computation reproduces this to the expected precision. All earlier errors and corrections are retained in `GED_V1_R3_DEV_AUDIT.md` with the original first test failure log.
