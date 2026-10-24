# GED-V1-R8 — Gravitational Information Recoverability

**Turnip Research Laboratory | Local research scout | 2026-10-24**

## Objective

Test how much of an internal mass-energy arrangement can be inferred from its **external gravitational acceleration** under deliberately restricted, conventional Newtonian assumptions. This is a synthetic inversion control, **not** a measurement, a derivation of new gravity, or a confirmation of an information-field hypothesis.

This test follows GED-V1-R7 (multipoles), but is **standalone Python**: it does not change or import prior frozen code. It must not be published/pushed to GitHub or Zenodo before the planned V1-R10 checkpoint, absent an independently vetted definite finding and user approval.

## Three scientific controls

1. **Restricted recoverability:** Assume a source is one of three known possibilities: two opposing equal-energy parcels on X (`pair_x`), two on Y (`pair_y`), or eight at cube corners (`cube`). Every case contains **100,000,000 J** of added energy (mass-equivalent `E/c²`). Generate vector acceleration at 24 exterior probes on six signed Cartesian axes, at radii 200, 300, 500 and 800 m. Given *exact knowledge* of this menu, positions, total energy, and perfectly removed background, select the source with minimum total squared gravitational-field residual. This is **classification under restrictive priors**, not recovery of arbitrary interior matter.
2. **Exact non-uniqueness:** Compare perfectly uniform spherical thin shells with identical mass but radii 30 m and 70 m. All observations outside the 100 m nominal body are mathematically identical by Newton's shell theorem, **even with no noise and complete exterior coverage**. A probe at 50 m inside the outer shell but outside the inner shell separates the two models; it is intentionally outside the external-only measurement regime.
3. **Practical recoverability:** Inject reproducible independent Gaussian noise into every simulated acceleration component, at specified amplitudes in m/s². Repeat each source class 200 times at each noise level, and preserve full confusion matrices. These are **toy numerical noise levels**, not real instrument noise or achievable measurement sensitivities.

## Physical equations and units

The extra 100 MJ per body supplies the notional Newtonian mass equivalent

`M_add [kg] = E_add [J] / c² [m²/s²]`.

For each parcel at position `s` [m], measured at position `x` [m]:

`g(x) [m/s²] = -G [m³/(kg s²)] * m [kg] * (x-s) [m] / |x-s|³ [m³]`.

Vector contributions add linearly, without a novel coupling. The common spherical body's very much larger background field is **assumed to have been subtracted exactly**; its mass/shape is not inferred. A squared distance between two complete field signatures is `Σ_i (g_measured,i - g_template,i)²`, across 72 Cartesian components. Given independent equal-variance Gaussian noise, selecting the smallest residual is equivalent to known-template maximum likelihood. The algorithm **does not fit** unknown energy, positions, geometry, or calibration errors.

For a Newtonian uniform spherical shell, the exterior acceleration is that of a central point mass; acceleration everywhere inside the ideal shell is zero. Distinct shell radii therefore cannot be distinguished by exterior acceleration alone.

## Run in ~/Desktop/Cosmology/GED

```bash
cd ~/Desktop/Cosmology/GED
source .venv/bin/activate
unzip -n ~/Downloads/GED_V1_R8_bundle.zip -d .
python -m unittest -v test_ged_v1_r8.py
python ged_v1_r8.py
```

No third-party dependencies. Tested in Python 3.13.5; intended for your Python 3.13.17 venv.

The run creates a **new**, time-stamped directory under `runs/`, containing `results.json`, `exterior_signatures.csv`, `report.txt`, `tests.txt`, `provenance.json`, `input_snapshot/` and `SHA256SUMS`. The manifest hashes all recorded files except itself, as a self-hash is impossible. The provenance records Python/runtime, Git HEAD if available, configuration/source hashes, test command and publication hold. It does not alter Git or published releases.

Check a returned run's hashes from inside its directory:

```bash
cd ~/Desktop/Cosmology/GED/runs/GED_V1_R8_...   # replace ... with actual run directory
shasum -a 256 -c SHA256SUMS
```

## Scientific boundaries

- A **successful template match** only recovers an identifier from a three-element source menu known in advance. It does not reconstruct general material composition, atomic binding energy, stress or entropy.
- The **shell counterexample** proves that exterior Newtonian gravimetry alone does not uniquely determine arbitrary internal spherical distributions. Non-uniqueness is a mathematical property; precision improvements cannot resolve it without additional constraints or measurements.
- No physical model of parcel creation/transport, pressure, support stresses, elasticity, full GR, Einstein field equations or microscopic binding physics is solved.
- All noise is independent Gaussian numerical noise, not actual gravimeter performance. The common background is perfectly removed by assumption, which is not realistic in a precision measurement campaign.
- No material-instruction field, entropy coupling, extra gravitational information channel, experimental measurement or new physics is inferred.
- Reported recovery rates are **finite seeded simulations** and may change with noise settings, source menu, sampling, and seed; they are not statistical evidence of new physics.

## Sources for conventional physics

- [Newton's shell theorem, Wolfram MathWorld](https://mathworld.wolfram.com/NewtonsShellTheorem.html)
- [MIT classical mechanics: gravitational field of a spherical shell, Physics LibreTexts](https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Classical_Mechanics_%28Dourmashkin%29/09%3A_Circular_Motion_Dynamics/9.04%3A_Appendix_9A_The_Gravitational_Field_of_a_Spherical_Shell_of_Matter_.)
- [Acquisition and Analysis of Terrestrial Gravity Data, Chapter 7: The inversion of gravity data (Cambridge University Press)](https://www.cambridge.org/core/books/abs/acquisition-and-analysis-of-terrestrial-gravity-data/inversion-of-gravity-data/7AA217A5761455072D1D69D7A41F3390)

## Status

**DEVELOPER VERIFICATION ONLY, PENDING MAC REPRODUCTION AND INDEPENDENT RECEIPT AUDIT.** GED-V1-R7 and earlier release contents are not modified. No GitHub/Zenodo publish before V1-R10, per user's publication cadence.
