# GED-V1-R9 — Information Loss and Gravitational Reconstruction Limits

**Turnip Research Laboratory | Experimental series: Gravitational Energy Distribution | Local numerical scout**

## Why this scout

GED-V1-R8 classified three **already-known** source shapes. GED-V1-R9 advances to an inverse problem: infer **six unknown nonnegative energy fractions at six specified positions** from external gravitational measurements. This removes the three-choice menu but **does not** infer unknown source locations, chemical composition, microscopic binding, material strength or fundamental physical laws.

The goal is to distinguish *mathematical non-uniqueness* from *limited observation geometry* and *measurement noise* under conventional Newtonian gravity. Nothing here is a physical detection or novel gravity law.

## The four controls

1. **Sparse observation:** a single exterior vector probe supplies three components; knowing the total energy supplies one more equation. The six unknown shares still have a two-dimensional nullspace. The experiment **constructs two different nonnegative source distributions** that give indistinguishable measurements at that probe.
2. **Poor geometry:** six x-axis probes measuring only x-acceleration have an augmented measurement matrix of rank 3/6. Additional probes arranged the same way do not reveal the y/z antisymmetric details. Geometry matters, not merely sample count.
3. **Distributed vector measurements:** eight spatially distributed probes supply 24 Cartesian measurements. The augmented matrix is full rank **within the six-site model**; fit the unknown nonnegative fractions subject to their total summing to one. Compare noise-free and seeded Gaussian-noise recovery.
4. **Exact spherical degeneracy:** concentric ideal thin shells with different radii and allocations, but identical total mass, give exactly identical Newtonian exterior fields. This ambiguity persists under any exterior measurement precision. An interior measurement can distinguish them.

## Units and equations

The added source energy is fixed at `E = 100,000,000 J`. Its Newtonian mass equivalent is `E/c²` in kg. The gravitational acceleration at an external point `x` from source position `s_j` with fraction `w_j` is

`g(x) = -[G E/c²] Σ_j w_j (x - s_j)/|x-s_j|³`, with `Σ_j w_j = 1`, `w_j ≥ 0`.

`G`: m³/(kg·s²); `E`: J; `c`: m/s; positions: m; `g`: m/s².

For numerical stability, calculate `g_ref = G(E/c²)/r_ref²` using `r_ref = 200 m`, and divide gravitational accelerations by `g_ref`. This gives a **dimensionless** linear forward equation `y = A w`. Scaling avoids hiding information under `~10^-24 m/s²` floating-point amplitudes; it **does not** introduce new physics.

The known total-energy constraint is an additional row of the measurement matrix. The number of independent constraints is its numerical rank. A nonzero null vector satisfies `A v≈0`, meaning two different source allocations can produce identical selected measurements. For full-rank measurements, a simple finite active-set constrained least-squares solver fits nonnegative weights with fixed total. This solver enumerates 63 possible subsets of six sources; it is a small toy inverse solver, **not a general research-grade inverse package**.

## Source and instrument assumptions

- Six possible source positions, fixed on opposing x/y/z axes 50 m from the centre, are **known perfectly** before inversion. The six **energy shares** are unknown to the fitter.
- Nominal 100 m spherical object. The much stronger gravitational field of the body is **assumed removed perfectly**.
- Energy parcels are ideal point-equivalent weak gravitational sources; material transport, elastic equilibrium, pressure and support stress are **not** modeled.
- Noise is simulated independent Gaussian acceleration noise using a recorded random seed; the values are not claimed achievable by current instruments.
- Independent measurements are treated as if uncorrelated, with perfect sensor placement and timing.
- Rank is computed from a double-precision, dimensionless matrix with a fixed `1e-11` pivot threshold. It is a **numerical rank** for this chosen geometry and model; a more general inverse analysis requires singular-value decomposition and measurement uncertainty propagation.
- Constrained fits may become inaccurate or encounter simplex boundaries at higher noise; the mean L1 fraction error is `Σ_j |w_est,j-w_true,j|`, a dimensionless share error (not classification accuracy).
- Shell-theorem degeneracy is exact analytically under Newtonian ideal spherical symmetry; floating-point equality is only a numerical implementation check.

## Run in your GED Python 3.13 virtual environment

```bash
cd ~/Desktop/Cosmology/GED
source .venv/bin/activate
unzip -n ~/Downloads/GED_V1_R9_bundle.zip -d .
python -m unittest -v test_ged_v1_r9.py
python ged_v1_r9.py
```

No third-party dependencies. Python 3.13 preferred.

The experiment creates a unique `runs/GED_V1_R9_.../` directory with `results.json`, `noise_sweep.csv`, `probe_geometry.csv`, `report.txt`, `tests.txt`, `input_snapshot/`, `provenance.json` and `SHA256SUMS`. The file manifest includes SHA-256 hashes for every preserved run file except itself. Exact inputs are snapshotted and versioned by content hash; the Git commit is recorded **if available**. Nothing is committed, pushed, or submitted to Zenodo.

To check the run integrity, replace the example directory with the one printed on your terminal:

```bash
cd ~/Desktop/Cosmology/GED/runs/GED_V1_R9_.../
shasum -a 256 -c SHA256SUMS
```

## Expected interpretation

- Inverse source recovery depends both on the **assumed source basis** and the **observation geometry**.
- Rank deficiency represents mathematical non-identifiability even before adding noise.
- Full numerical rank is **not** sufficient to guarantee good noisy recovery: signal-to-noise and conditioning matter.
- Exact spherical external-field degeneracy persists even with arbitrarily accurate external observations.
- None of these facts shows that gravitational fields encode instructions for binding, fracture or material response. Those remain unproved aspects of the original GED hypothesis.

**Publication policy:** Preserve all local receipts and failures. Do not push or publish V1-R9 to GitHub or Zenodo before the V1-R10 milestone, absent an independently validated definite finding and explicit approval.

## Prior controls

The R7 asymmetric Newtonian multipole scout and R8 known-template classification scout are preserved as separate experiments. R9 does not import or change their frozen source files.
