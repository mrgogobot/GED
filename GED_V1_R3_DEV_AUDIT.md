# GED-V1-R3 — Development audit

Date: 2026-10-24. Status: **development stage numerical scout**, subject to independent execution and local provenance audit by user.

## Inputs and architecture

- Reused the exact `ged_v1_r2.py` from the earlier GED-V1-R2 package rather than duplicating gravitational self-binding equations.
- SHA-256 of required GED-V1-R2 source: `c1c8bd268fcfad2f8dea23bba703ef04284ced2590680584ccfa263cd9010813`.
- The V1-R3 experiment adds a thermal-energy and entropy accounting layer with idealised constant heat capacities and heat bath.
- Input values and assumptions are frozen in `GED_V1_R3_config.json`.

## Preserved development failures (not science failures)

1. **Initial 46-test run: 3 FAIL**, caused by test arithmetic inadvertently using Python Decimal's default 28-digit context while the science model used 60 digits. Failing checks: exported-entropy identity, inverse-square `delta_g`, and half-export mass identity. The complete first run's output is preserved in `GED_V1_R3_PRECISION_FAILURE_A.log`.
2. **First test correction: 1 FAIL**: the inverse-square check used an unrealistic exact-equality assertion. The last few digits differed at 60-digit precision because two rounded division pathways differ. This was **a test precision issue, not a failed physical identity**.
3. **Final correction:** the inverse-square comparison uses a relative tolerance of `1e-55`, far below any practical measurement uncertainty. All remaining math equality checks use the explicit 60-digit context when doing arithmetic.
4. **Input-gate improvement:** enforce inclusion of the 50% escaped-energy case, since the generated report always includes that case.

**No physical equations were changed to make the tests pass.** The initially observed failures remain part of the development record.

## Explicit scientific limitations

- Release energy is computed from the difference in **Newtonian self-gravitational potential**, not nuclear/chemical binding or structural breakage.
- The B→A internal rearrangement is prescribed, with no work input or dynamics modeled; no claim that rearrangement is spontaneous.
- Temperature-dependent heat capacity, phase transitions, pressure, configurational entropy, free energy and irreversible transport are omitted.
- `Q_bath/T_bath` applies only to an *equilibrated thermal bath*, not directly to free-streaming photons.
- Total **thermal** entropy is computed; the **full** entropy change of the reconfiguration is not known from this model.
- A stable system need not continuously radiate energy to stay bound. Such radiation would require an energy source.
- No additional gravitational information coupling was fitted, derived, measured or tested.

## Verification

- Developer Python runtime: Python 3.13.x, standard library only.
- **46 tests pass** after corrections and before clean-room packaging.
- Independent local execution on the user's Python 3.13.17 remains pending.
- Provenance acceptance requires comparing receipt hashes with the actual local run.
