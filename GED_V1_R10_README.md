# GED-V1-R10 — CCMEI, radial layers and gravitational information retention

**Turnip Research Laboratory | Gravitational Energy Distribution | R10 milestone candidate | local numerical scout**

## The thought being preserved

> If gravity is a collective product of constituent matter, how much information about the constituent relationships survives in the resulting field?

**CCMEI** = Complete Constituent Mass-Energy Inventory. The idea under investigation is that a body's collective gravitational field originates from its constituents and their energy relationships, and **might** carry additional recoverable information about how its material responds. The italicized portion is a hypothesis, not an established consequence of the models.

**Do not mistake the term “layered” for a signal physically traveling from the core through successive shells.** In the model, every spherical mass-energy layer contributes according to Newton's shell theorem; time-dependent propagation is not modeled.

## What R10 actually tests

Three idealized spherical shells have edges **0, 20, 50, 100 metres**. They contain total mass-energy equivalents, respectively:

| Case | Inner | Middle | Outer | Total |
|---|---:|---:|---:|---:|
| A | 2 × 10⁹ kg | 3 × 10⁹ kg | 5 × 10⁹ kg | 1 × 10¹⁰ kg |
| B | 3 × 10⁹ kg | 2 × 10⁹ kg | 5 × 10⁹ kg | 1 × 10¹⁰ kg |
| C | 2 × 10⁹ kg | 3 × 10⁹ kg | 5 × 10⁹ kg | 1 × 10¹⁰ kg |

- **A versus B:** Same total mass-energy, distributed differently among shells. Their *interior* Newtonian fields differ, while their *exterior* fields are identical.
- **A versus C:** Same radial shell mass-energy; different **synthetic binding vs constituent rest-equivalent bookkeeping**, each arranged to reproduce the same total energy in every shell. Their Newtonian fields agree **at all radii**. **This is not a physical construction of actual iron or ice with equal stress-energy**—it is a mathematical demonstration that one scalar gravitational-source profile does not identify its components.
- **Inverse geometry:** Using *only* three exterior acceleration magnitudes, we can identify just **one** of the three independent shell masses (rank 1 of 3). We must not infer the other two from those data. At three specified interior radii (10, 35, 75 m), the model is full rank (3 of 3), so we can exactly reconstruct all three unknown shell masses **if shell boundaries and uniform densities are already known**. This says **nothing** about reconstructing quantum microstates or chemical identities.

## Mathematical model (with units)

For shell `i` with inner and outer radii `a_i, b_i` (metres) and total mass-energy equivalent `M_i` (kilograms), the enclosed fraction at radius `r` is:

`f_i(r) = 0 (r<=a_i); (r³-a_i³)/(b_i³-a_i³) (a_i<r<b_i); 1 (r>=b_i)`.

`M_enclosed(r) = Σ_i M_i f_i(r)` (kg).

`g(r) = G M_enclosed(r)/r²` (m/s²), for `r>0`, with `g(0)=0`.

A shell's **synthetic** aggregate source-energy ledger is

`M_i = M_rest_equiv,i + E_binding,i/c²`.

`E_binding,i` is a signed *contribution* to energy in joules; this is not a claim that total nuclear/chemical binding can be inferred from gravity. The constituent-rest-equivalent term is chosen separately in A and C to conserve shell total. It avoids double-counting binding energy. `G` is in m³ kg⁻¹ s⁻², `c` in m s⁻¹.

For inverse calculations rescale an acceleration to the **enclosed mass** `g(r) r²/G` (kg). The corresponding three columns of the linear measurement map are shell volume fractions. A **rank** indicates how many independent combinations can be identified from the selected positions, not the total information in a full dynamical gravitational field. The inverse test uses exact rational arithmetic, and high-precision Decimal is used for field strengths.

## Limitations and what this does NOT show

- Weak-field, static, spherical **Newtonian** gravity only; no pressure, stress, elastic equilibrium, full general relativity, time-dependent wave propagation, rotational or electromagnetic effects.
- The binding ledger is deliberately constructed; it does not model actual atomic species, microscopic bonding, material strength, or melting.
- It is incorrect to conclude that internal microscopic information is absent in nature merely because this restricted model does not recover it.
- R10 derives **no additional gravitational channel** and detects **no anomalous effect**. The conclusion concerns recoverability *from this chosen forward model and probe geometry*.
- Tests validate the code and its expected mathematical baseline; they are not observational evidence.

References for standard physics context:
- Newton's shell theorem: https://mathworld.wolfram.com/NewtonsShellTheorem.html
- Spherical mass distributions: https://galaxiesbook.org/chapters/I-01.-Gravitation_2-Spherical-systems%3A-Newton%27s-shell-theorems.html
- General relativistic stress-energy: https://warwick.ac.uk/fac/sci/physics/intranet/pendulum/generalrelativity/

## Execute in existing GED directory

```bash
cd ~/Desktop/Cosmology/GED
source .venv/bin/activate
unzip -n ~/Downloads/GED_V1_R10_bundle.zip -d .
python -m unittest -v test_ged_v1_r10.py
python ged_v1_r10.py
```

Requires only the **Python standard library** (Python 3.13 preferred); no old GED source module is imported and no additional dependency is required.

A fresh `runs/GED_V1_R10_.../` folder records exact input snapshots, `results.json`, `radial_samples.csv`, `report.txt`, `tests.txt`, `provenance.json`, `SHA256SUMS`, and SHA-256 input/output hashes. It will record Git HEAD *if present*, but **will not stage, commit, tag, push, or publish anything**.

For the output directory printed by your terminal:

```bash
cd ~/Desktop/Cosmology/GED/runs/GED_V1_R10_<timestamp>_<id>
shasum -a 256 -c SHA256SUMS
```

**Publication checkpoint:** V1-R10 is a **candidate for audit and a possible tagged release**. The existing V1-R4 public GitHub/Zenodo release remains unchanged. Publishing requires an R10 local reproduction, full evidence/provenance review, and the user's explicit approval. No finding of new physics has been earned that would justify an early scientific claim.

See `GED_V1_R10_PRIOR_EVIDENCE.json` for each previous scout's reported verification status. Its bundle hashes were checked against the developer-side files, **not** against every preserved local Mac run. The reported prior test total is **478 across ten separate suites**; this is not one combined fresh test execution.
