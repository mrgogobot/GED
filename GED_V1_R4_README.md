# GED-V1-R4 — Equal Energy, Different Entropy

**Status:** Verified *numerical controls* under the stated assumptions. This is **not** observational evidence, an independent entropy-based gravity law, an information field detection, or a full general-relativistic stress-energy calculation.

## Hypothesis under examination

The GED idea suggests gravitational fields might encode information about constituents' collective structure and responses to physical stresses. We must first determine which apparent information signatures can already be explained by mass-energy distribution, and which are only entropy descriptions of material states.

## What is tested

There are **three controls**, all using SI units:

**Control A: Different macrostate entropy, identical model stress-energy.** An ideal collection of `n=256` independent exactly degenerate binary labels has macrostates with `k=0` or `k=128` occupied labels. Their Boltzmann macrostate entropies are `k_B ln(comb(n,k))`, so the counted entropies differ. *By construction*, their energy and local stresses are identical, meaning there is no gravitational difference in the conventional model. This illustrates that a statistical state count is **not automatically an extra gravitational source**. The system is deliberately idealised; ordinary magnetic spins can produce magnetic fields and energy differences, and there is no simulation of dynamics from one macrostate to the other.

**Control B: Different entropy, same total energy, changed radial distribution.** A homogeneous `R=100 m` sphere of reference density `3000 kg/m³` is split into equal-volume inner/outer zones, with identical constant specific heat `1000 J/(kg K)` and fixed geometry. At the same total modeled thermal energy, temperatures `300 K / 200 K` relax to `250 K / 250 K`. The entropy production is `ΔS=C ln[(250²)/(300×200)] > 0`, where `C` is the heat capacity **of one zone**. No energy crosses the boundary. `g(r)` changes inside because radial mass-energy redistributes, but the exterior Newtonian field remains identical.

**Control C: Equal entropy, different radial distribution.** Swapping the hot and cold zones from `300/200 K` to `200/300 K` leaves the total thermal energy and thermal entropy identical (both zones have equal heat capacities), but reverses the interior gravitational difference. The exterior Newtonian monopole remains the same.

These controls prevent an accidental claim that a gravitational difference tracks entropy independently of energy distribution.

### Equations (what and why)

- `M_ref = ρ (4πR³/3)` (kg): identical base material mass by construction.
- `C = c_p M_ref/2` (J/K): fixed equal heat capacity in the two half-volume zones.
- `ΔE_i = C (T_i−T_ref)` (J): thermal energy relative to the common reference state.
- `ΔM_i = ΔE_i/c²` (kg): mass-equivalent of retained thermal energy, in weak-field bookkeeping.
- `ΔS_thermal = C[ln(T_core/T_ref)+ln(T_shell/T_ref)]` (J/K): idealised equilibrium entropy differences, even though the two zones may not be in mutual equilibrium.
- `M_enc(r)`: integrated radial mass of core and shell; `g(r)=G M_enc(r)/r²` (m/s²) in the Newtonian spherical model.
- `S_B−S_A = k_B ln[binom(n,k_Bmacro)/binom(n,k_Amacro)]` (J/K) for Control A. Here `k_B` the Boltzmann constant should **not** be confused with the macrostate count variable.

## Scientific limitations (read before interpreting)

The 100-m object is only an illustrative numerical control, not a realistic equilibrium asteroid. The model fixes the volume/densities even when temperature differs and omits thermal expansion, pressure/stress gradients, binding/strain/work changes, gravitating field self-energy, transport rates, and exact Einstein-field-equation solutions. A hot-core/cold-shell state is not stationary thermodynamic equilibrium without constraints; heat flow is prescribed, and the ideal heat-conduction path is not simulated. More complete general relativity requires the full stress-energy tensor and consistent boundary conditions. The reported exterior equality applies to **this Newtonian spherical monopole with fixed total modeled mass-energy**, not an unrestricted claim about arbitrary real bodies. **Entropy is not independent gravitational source data in the conventional model.** Whether GED can derive a new lawful observable remains open; no such equation is introduced here.

## Requirements

- Python **3.13.17** recommended (standard library only; no `pip install`).
- Existing exact `ged_v1_r2.py` in `~/Desktop/Cosmology/GED`, SHA-256:
  `c1c8bd268fcfad2f8dea23bba703ef04284ced2590680584ccfa263cd9010813`
- Existing `.venv` activated.
- **Do not replace/overwrite** the pinned V1-R2 source or previous experiments.

## Local Mac execution

```bash
cd ~/Desktop/Cosmology/GED
source .venv/bin/activate
unzip -n ~/Downloads/GED_V1_R4_bundle.zip -d .
python -m unittest -v test_ged_v1_r4.py
python ged_v1_r4.py
```

The runner independently reruns all tests *inside its input snapshot* before declaring PASS. Expect `runs/GED_V1_R4_<UTC>_<ID>/` containing:

```text
input_snapshot/             exact source, config, docs and pinned dependency
  ged_v1_r2.py
  ged_v1_r4.py
  test_ged_v1_r4.py
  GED_V1_R4_config.json
  GED_V1_R4_README.md
  GED_V1_R4_DEV_AUDIT.md
  GED_V1_R4_PRECISION_FAILURE_A.log
results.json                all model parameters, controls and limitations
radial_gravity.csv          selected inner, boundary and exterior probes
tests.txt                   snapshot test output (PASS or FAIL)
report.txt                  human-readable summary
provenance.json             runtime, Git if present, commands and hashes
SHA256SUMS                  all captured evidence files except itself
```

The failed development test output is retained intentionally; the error was in the original **test arithmetic/context and zero string representation**, not a demonstrated scientific failure. This package is a versioned scout; the ledger is generated locally by this script, **not** by invoking the separate Turnip Run wrapper.

### Acceptance boundaries

**May be accepted as:** mathematical demonstration, under explicitly stated conventional physics approximations, that state-counting entropy and gravitational mass-energy distribution are not the same variable; local numerical code verification and reproducibility.

**Must not be accepted as:** a GED entropy-gravity coupling, observational detection, or a proof that entropy can never affect gravity indirectly through state-dependent stress-energy.
