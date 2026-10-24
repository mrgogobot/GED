# GED-V1-R7: development audit, 2026-10-24

## Scope and verification

Conventional Newtonian multipole model with fixed 100-m spherical iron/ice background and equal 100,000,000-J additional energy in three constructed arrangements. Python 3.13 standard library; 70-digit Decimal. No new gravity or material-information mechanism demonstrated.

The equivalent added mass of each is 1.112650056053618432174090e-9 kg. The common background mass is 7.481702904524092597398785e9 kg. In all three configurations the energy monopole is identical and the source-centred dipole vanishes to numerical precision. Cubic symmetry forces the eight-parcel quadrupole to zero, yet leaves higher multipoles.

Independent formulas checked by separate unit tests:

* Pair along its axis: gx = -GM/2 [1/(r-a)^2 + 1/(r+a)^2].
* Pair transverse: gy = -GM*r/(r^2+a^2)^(3/2).
* Cube eight corners, b=a/sqrt(3): gx = -GM/2[(r-b)/((r-b)^2+2b^2)^(3/2)+(r+b)/((r+b)^2+2b^2)^(3/2)].
* Axial quadrupole for x pair: -3GMa^2/r^4. Transverse: +3GMa^2/(2r^4).

At +x, r=200m: exact pair_x minus cube acceleration is -6.141600042717311254966278e-25 m/s^2. Leading quadrupole estimate is -5.012658181655099213755182e-25 m/s^2, relative discrepancy 18.3818851%. The relative discrepancy falls with increasing radius: 5.036% at 400m, 1.291% at 800m, 0.3248% at 1600m, 0.08134% at 3200m. This is conventional multipole convergence.

## Failure preservation

GED_V1_R7_PRECISION_FAILURE_A.log preserves the first 61-test development run with four failures: exact-zero assertions for dipole/quadrupole quantities encountered Decimal round-off near 1e-78 kg*m and 1e-75 kg*m^2. Tests now use explicitly specified tolerances much smaller than meaningful physical scales. Physical equations unchanged.

GED_V1_R7_REPORT_FORMAT_FAILURE_A.log preserves a malformed 67-times repeated report header discovered after the first green 61-test run. The issue was Python string concatenation around a repeated separator; corrected with explicit + operations. Physical equations unchanged.

Final development suite: 62 tests passing. The subprocess test run sets PYTHONDONTWRITEBYTECODE=1, keeping regenerateable bytecode caches out of the frozen input snapshots. Actual clean-extraction run and SHA-256 audit are separate completion gates.

## Interpretation / limits

External Newtonian gravity distinguishes asymmetric mass-energy distributions through conventional higher multipoles. The construction does not solve elastic equilibrium, complete stress-energy, matter binding, source motion, supporting forces, or Einstein's equations; the parcels are idealized point-source proxies. No new physical law, observation, or information-field coupling is identified.

Local research only. GED version R10 is the next planned release milestone, unless a definite finding is separately reviewed and approved by the user. No git push or Zenodo publication occurs here.
