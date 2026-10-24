#!/usr/bin/env python3
"""GED-V1-R6 — radial elastic-energy distribution (numerical baseline).

WHAT:
    Compare two SPHERICALLY SYMMETRIC, idealized iron-core/ice-shell
    configurations with the same constituent mass, the same radius and
    the same TOTAL reversible elastic energy, but different radial allocations
    of that energy. Calculate a Newtonian mass-energy monopole proxy for
    gravitational acceleration both inside and outside the body.

WHY:
    See whether ordinary gravitational bookkeeping already accounts for
    changing the POSITION of stored elastic energy, without assuming a
    hypothetical material-instruction or entropy-specific gravity field.

LIMITATION:
    Uniform assigned elastic energy in radial zones is a constructed scalar
    bookkeeping profile, not a solved self-consistent static stress field.
    Elastic compatibility, force equilibrium, boundary traction, shape changes,
    gravitational self-energy changes, anisotropic stress and a full GR solution
    are NOT included. The energy correction is a weak-field Newtonian PROXY.
    Nothing here derives or detects a new physical coupling.

REPRODUCIBILITY:
    Standard-library Python; exact SHA-256-pinned GED-V1-R2 gravity helpers
    and GED-V1-R5 elasticity helpers; snapshots, tests and SHA-256 receipts.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import uuid

# WHAT: Import the existing, authoritative physics functions.
# WHY: Avoid maintaining second, subtly different implementations of
# Newtonian acceleration, mass-energy conversion, or Hookean strain energy.
import ged_v1_r2 as r2
import ged_v1_r5 as r5

ROOT = Path(__file__).resolve().parent
PRECISION = 70
ZERO, ONE, TWO, THREE = map(Decimal, ('0', '1', '2', '3'))
R2_SHA256 = 'c1c8bd268fcfad2f8dea23bba703ef04284ced2590680584ccfa263cd9010813'
R5_SHA256 = '608b67e91756e26667818743cdb906d487eeae65ca3d3f1de350bd87af50749d'
INPUT_FILES = (
    'ged_v1_r6.py', 'test_ged_v1_r6.py', 'GED_V1_R6_config.json',
    'GED_V1_R6_README.md', 'GED_V1_R6_DEV_AUDIT.md',
    'GED_V1_R6_FAILURE_A.log', 'ged_v1_r2.py', 'ged_v1_r5.py',
)


def D(value):
    """WHAT: Validate a finite decimal. WHY: retain tiny gravity differences."""
    return r2.D(value)


def sha256(path):
    """WHAT: SHA-256 file's exact bytes. WHY: identify an exact version."""
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(part)
    return digest.hexdigest()


def check_dependencies(root=ROOT):
    """WHAT: Verify BOTH inherited helper files. WHY: prevent version drift."""
    root = Path(root)
    checks = {'ged_v1_r2.py': R2_SHA256, 'ged_v1_r5.py': R5_SHA256}
    for filename, expected in checks.items():
        actual = sha256(root / filename)
        if actual != expected:
            raise ValueError(f'Pinned dependency mismatch {filename}: {actual}')
    return checks


def load_config(path=ROOT / 'GED_V1_R6_config.json'):
    """WHAT: Read immutable experiment inputs. WHY: no hidden run defaults."""
    cfg = json.loads(Path(path).read_text(encoding='utf-8'))
    if cfg.get('experiment_id') != 'GED-V1-R6':
        raise ValueError('Incorrect experiment identifier')
    if cfg.get('r2_dependency_sha256') != R2_SHA256:
        raise ValueError('Configuration pins incorrect R2 source')
    if cfg.get('r5_dependency_sha256') != R5_SHA256:
        raise ValueError('Configuration pins incorrect R5 source')
    return cfg


def validate(cfg):
    """WHAT: Check geometry, energy, material constants, probes and cases.

    WHY: No negative radii, nonfinite physical inputs, invisible missing
    probes, or unmatched total energy may silently enter a scientific result.
    """
    if cfg.get('experiment_id') != 'GED-V1-R6':
        raise ValueError('Wrong experiment')
    if cfg.get('r2_dependency_sha256') != R2_SHA256 or cfg.get('r5_dependency_sha256') != R5_SHA256:
        raise ValueError('Wrong pinned dependencies')
    p = {k: D(cfg[k]) for k in (
        'G_m3_kg_s2', 'c_m_s', 'outer_radius_m', 'core_radius_m',
        'total_elastic_energy_J')}
    if min(p.values()) <= ZERO:
        raise ValueError('All experiment constants must be positive')
    if p['core_radius_m'] >= p['outer_radius_m']:
        raise ValueError('Core radius must be strictly smaller than outer radius')
    materials = cfg.get('materials')
    if not isinstance(materials, dict) or set(materials) != {'core', 'shell'}:
        raise ValueError('Exactly core and shell material definitions required')
    mat = {}
    for region in ('core', 'shell'):
        entry = materials[region]
        E, rho = D(entry['young_modulus_Pa']), D(entry['density_kg_m3'])
        if min(E, rho) <= ZERO:
            raise ValueError('Young modulus and mass density must be positive')
        mat[region] = {'name': str(entry['name']), 'E_Pa': E, 'rho_kg_m3': rho}
    cases = cfg.get('configurations')
    if not isinstance(cases, dict) or set(cases) != {'core_weighted', 'shell_weighted'}:
        raise ValueError('Must have two independently named loading configurations')
    fractions = {}
    for key, value in cases.items():
        fraction = D(value['core_energy_fraction'])
        if not ZERO < fraction < ONE:
            raise ValueError('Each core fraction must be strictly between zero and one')
        fractions[key] = fraction
    if fractions['core_weighted'] <= fractions['shell_weighted']:
        raise ValueError('Core-weighted case must store more energy in core')
    samples = cfg.get('sample_radii_m')
    if not isinstance(samples, list) or len(samples) < 4:
        raise ValueError('At least four sample radii required')
    radii = [D(v) for v in samples]
    if radii[0] != ZERO or any(r < ZERO for r in radii):
        raise ValueError('Samples must include center and have nonnegative radii')
    if radii != sorted(set(radii)):
        raise ValueError('Sample radii must be strictly increasing')
    if p['core_radius_m'] not in radii or p['outer_radius_m'] not in radii:
        raise ValueError('Both material boundaries must be sampled')
    if radii[-1] <= p['outer_radius_m']:
        raise ValueError('At least one exterior point required')
    return p, mat, fractions, radii


def sphere_volume(radius_m):
    """WHAT: Volume of a spherical region. WHY: density * volume = mass.

    Here PI is defined by the R2 Decimal helper, rather than another
    package's approximation. The same common factor cancels in fractions.
    """
    r = D(radius_m)
    if r < ZERO:
        raise ValueError('Negative radius')
    with localcontext() as ctx:
        ctx.prec = PRECISION
        return +(Decimal(4) * r2.PI * r**3 / THREE)


def geometry(p, materials):
    """WHAT: Compute constituent volumes and rest-mass inventory.

    WHY: Every comparison must preserve the same iron and ice inventory.
    This is a conventional mass distribution, NOT a new information field.
    """
    with localcontext() as ctx:
        ctx.prec = PRECISION
        core_vol = sphere_volume(p['core_radius_m'])
        all_vol = sphere_volume(p['outer_radius_m'])
        shell_vol = all_vol - core_vol
        core_mass = core_vol * materials['core']['rho_kg_m3']
        shell_mass = shell_vol * materials['shell']['rho_kg_m3']
        return {
            'core_volume_m3': +core_vol, 'shell_volume_m3': +shell_vol,
            'core_rest_mass_kg': +core_mass, 'shell_rest_mass_kg': +shell_mass,
            'total_rest_mass_kg': +(core_mass + shell_mass),
        }


def make_case(p, materials, geom, fraction):
    """WHAT: Assign the SAME total energy, split it between two radial zones.

    WHY: Isolate WHERE energy sits while keeping total energy fixed.
    Compute associated ideal elastic strains using the EXACT R5 helper.

    Caveat: independently specified zone strains need not satisfy physical
    displacement compatibility or stress-equilibrium across the boundary.
    Treat these as constructed energy profiles, not stationary asteroids.
    """
    with localcontext() as ctx:
        ctx.prec = PRECISION
        total = p['total_elastic_energy_J']
        energies = {'core': +(total * fraction), 'shell': +(total * (ONE - fraction))}
        data = {}
        for region in ('core', 'shell'):
            volume = geom[f'{region}_volume_m3']
            mechanical = r5.state_for_target_energy(
                materials[region]['E_Pa'], volume, energies[region])
            data[region] = {
                'allocated_energy_J': energies[region],
                'energy_density_J_m3': +(energies[region] / volume),
                'stress_Pa': mechanical['stress_Pa'],
                'strain_dimensionless': mechanical['strain_dimensionless'],
            }
            if mechanical['strain_dimensionless'] >= Decimal('0.01'):
                raise ValueError('Proposed loading outside the small-strain regime')
        if data['core']['allocated_energy_J'] + data['shell']['allocated_energy_J'] != total:
            raise AssertionError('Energy allocation failed conservation')
        return data


def enclosed_fraction(r, outer_region_radius, inner_region_radius=ZERO):
    """WHAT: Fraction of a constant-volume-density spherical zone enclosed.

    WHY: Newton's shell theorem says ONLY enclosed mass contributes to
    radial gravitational acceleration in a spherically symmetric model.
    For a spherical shell, fraction=(r³-a³)/(b³-a³) in a<r<b.
    """
    r, a, b = D(r), D(inner_region_radius), D(outer_region_radius)
    if a < ZERO or b <= a or r < ZERO:
        raise ValueError('Invalid shell bounds or evaluation radius')
    with localcontext() as ctx:
        ctx.prec = PRECISION
        if r <= a:
            return ZERO
        if r >= b:
            return ONE
        return +((r**3 - a**3) / (b**3 - a**3))


def sample_field(p, materials, geom, case, r):
    """WHAT: Baseline plus monopole-equivalent elastic-energy gravity.

    WHY: Compare spatial arrangements at radius r, explicitly retaining
    the same baseline constituent mass inventory for both configurations.

    This is a NEWTONIAN toy energy-density correction to enclosed mass.
    It is NOT a tensor-gravity solution or exact stressed-solid prediction.
    """
    r = D(r)
    if r < ZERO:
        raise ValueError('Negative radius')
    with localcontext() as ctx:
        ctx.prec = PRECISION
        a, b = p['core_radius_m'], p['outer_radius_m']
        frac_core = enclosed_fraction(r, a)
        frac_shell = enclosed_fraction(r, b, a)
        enclosed_rest_mass = +(geom['core_rest_mass_kg'] * frac_core
                               + geom['shell_rest_mass_kg'] * frac_shell)
        enclosed_elastic_energy = +(case['core']['allocated_energy_J'] * frac_core
                                    + case['shell']['allocated_energy_J'] * frac_shell)
        # At the mathematical center spherical symmetry gives g=0. A
        # naive G*M/r² evaluation would incorrectly divide by zero.
        if r == ZERO:
            baseline_g = ZERO
            correction_g = ZERO
        else:
            baseline_g = +r2.exterior_g(enclosed_rest_mass, r, p['G_m3_kg_s2'])
            correction_g = +r2.delta_g(
                enclosed_elastic_energy, r, p['G_m3_kg_s2'], p['c_m_s'])
        return {
            'radius_m': r,
            'enclosed_rest_mass_kg': enclosed_rest_mass,
            'enclosed_elastic_energy_J': enclosed_elastic_energy,
            'mass_equivalent_enclosed_kg': +r2.mass_equivalent(
                enclosed_elastic_energy, p['c_m_s']),
            'baseline_g_m_s2': baseline_g,
            'elastic_delta_g_proxy_m_s2': correction_g,
            'total_g_proxy_m_s2': +(baseline_g + correction_g),
        }


def calculate(cfg):
    """WHAT: Form two loading states and compare their radial fields.

    WHY: Equal exterior field need not imply the same interior energy
    configuration; this is the conventional null baseline for GED.
    """
    p, materials, fractions, radii = validate(cfg)
    with localcontext() as ctx:
        ctx.prec = PRECISION
        geom = geometry(p, materials)
        cases = {name: make_case(p, materials, geom, f) for name, f in fractions.items()}
        rows = []
        for r in radii:
            a = sample_field(p, materials, geom, cases['core_weighted'], r)
            b = sample_field(p, materials, geom, cases['shell_weighted'], r)
            rows.append({
                'radius_m': r,
                'core_weighted': a, 'shell_weighted': b,
                'difference_delta_g_proxy_m_s2': +(a['elastic_delta_g_proxy_m_s2']
                                                   - b['elastic_delta_g_proxy_m_s2']),
            })
        # The two zone energy changes are opposites; the total sum is ZERO.
        core_transfer = +(cases['shell_weighted']['core']['allocated_energy_J']
                          - cases['core_weighted']['core']['allocated_energy_J'])
        shell_transfer = +(cases['shell_weighted']['shell']['allocated_energy_J']
                           - cases['core_weighted']['shell']['allocated_energy_J'])
        if core_transfer + shell_transfer != ZERO:
            raise AssertionError('Energy redistribution ledger not balanced')
        return {
            'experiment_id': 'GED-V1-R6',
            'verdict': 'CONVENTIONAL_NEWTONIAN_DISTRIBUTION_CONTROL_ONLY',
            'observational_data': False,
            'novel_gravitational_information_channel': 'NOT_DERIVED_OR_TESTED',
            'geometry': geom,
            'fixed_total_elastic_energy_J': p['total_elastic_energy_J'],
            'configs': cases,
            'radial_samples': rows,
            'repartition_ledger_J': {'core': core_transfer, 'shell': shell_transfer,
                                     'whole_body': +(core_transfer + shell_transfer)},
            'limitations': cfg['limitations'],
        }


def json_safe(value):
    """WHAT: Convert Decimal to decimal strings. WHY: prevent rounding."""
    if isinstance(value, Decimal):
        return str(value)
    if isinstance(value, dict):
        return {key: json_safe(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def git_state(root):
    """WHAT: Record local Git provenance if any. WHY: hashes != commits."""
    def run_git(*args):
        try:
            result = subprocess.run(['git', '-C', str(root), *args],
                                    capture_output=True, text=True, check=False)
            return result.stdout.strip() if result.returncode == 0 else None
        except FileNotFoundError:
            return None
    commit = run_git('rev-parse', 'HEAD')
    if not commit:
        return {'commit': None, 'exact_tag': None, 'dirty': None}
    status = run_git('status', '--porcelain', '--untracked-files=all')
    return {'commit': commit, 'exact_tag': run_git('describe', '--tags', '--exact-match'),
            'dirty': bool(status) if status is not None else None}


def run_experiment(root=ROOT, config_name='GED_V1_R6_config.json'):
    """WHAT: Execute one uniquely recorded run with frozen source and tests.

    WHY: Local R6 through R9 must remain reproducible WITHOUT uploading or
    publishing new versions to GitHub/Zenodo before the R10 milestone.
    """
    root = Path(root).resolve()
    check_dependencies(root)
    cfg = load_config(root / config_name)
    result = calculate(cfg)

    now = datetime.now(timezone.utc)
    run_id = f'GED_V1_R6_{now.strftime("%Y%m%dT%H%M%S%fZ")}_{uuid.uuid4().hex[:8]}'
    out = root / 'runs' / run_id
    snapshot = out / 'input_snapshot'
    snapshot.mkdir(parents=True, exist_ok=False)
    for filename in INPUT_FILES:
        src = root / filename
        if not src.is_file():
            raise FileNotFoundError(f'Missing required source: {filename}')
        shutil.copy2(src, snapshot / filename)
    hashes = {filename: sha256(snapshot / filename) for filename in INPUT_FILES}

    # WHAT: Test the saved snapshot itself, not a later edited working file.
    # WHY: This ensures the test receipts belong to this exact source version.
    cmd = [sys.executable, '-B', '-m', 'unittest', '-v', 'test_ged_v1_r6.py']
    check = subprocess.run(cmd, cwd=snapshot, capture_output=True, text=True,
                           env={**os.environ, 'GED_R6_SNAPSHOT_TESTS': '1'})
    (out / 'tests.txt').write_text(check.stdout + check.stderr, encoding='utf-8')
    if check.returncode != 0:
        (out / 'FAILED_TESTS.txt').write_text('Snapshot tests failed. DO NOT ACCEPT.\n', encoding='utf-8')
        raise RuntimeError(f'Frozen run tests failed; failure evidence kept: {out}')

    (out / 'results.json').write_text(
        json.dumps(json_safe(result), indent=2, sort_keys=True) + '\n',
        encoding='utf-8')
    with (out / 'radial_energy_gravity.csv').open('w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(('radius_m', 'core_enclosed_energy_J', 'shell_enclosed_energy_J',
                         'core_delta_g_m_s2', 'shell_delta_g_m_s2', 'delta_g_difference_m_s2'))
        for row in result['radial_samples']:
            writer.writerow(tuple(str(v) for v in (
                row['radius_m'], row['core_weighted']['enclosed_elastic_energy_J'],
                row['shell_weighted']['enclosed_elastic_energy_J'],
                row['core_weighted']['elastic_delta_g_proxy_m_s2'],
                row['shell_weighted']['elastic_delta_g_proxy_m_s2'],
                row['difference_delta_g_proxy_m_s2'],
            )))

    lines = [
        'GED-V1-R6 — RADIAL STRESS-ENERGY DISTRIBUTION SCOUT',
        '=' * 69,
        'CONVENTIONAL WEAK-FIELD NEWTONIAN MASS-ENERGY PROXY ONLY',
        'NOT A SOLUTION OF THE FULL STRESS-ENERGY TENSOR / EINSTEIN EQUATIONS',
        f"Core rest mass: {result['geometry']['core_rest_mass_kg']:.12E} kg",
        f"Shell rest mass: {result['geometry']['shell_rest_mass_kg']:.12E} kg",
        f"Total retained elastic energy per configuration: {result['fixed_total_elastic_energy_J']:.12E} J",
    ]
    for case_name, case in result['configs'].items():
        lines.append(f'\n{case_name}:')
        for region, data in case.items():
            lines.append(f'  {region}: U={data["allocated_energy_J"]:.12E} J; '
                         f'energy density={data["energy_density_J_m3"]:.12E} J/m^3; '
                         f'strain={data["strain_dimensionless"]:.12E}')
    lines.append('\nRADIUS / DIFFERENCE IN MASS-ENERGY GRAVITY PROXY:')
    for row in result['radial_samples']:
        lines.append(f'  r={row["radius_m"]} m: '
                     f'delta_g(core-weighted - shell-weighted)='
                     f'{row["difference_delta_g_proxy_m_s2"]:.12E} m/s^2')
    lines += [
        '\nRepartition ledger: core + shell = 0 J; no energy is created.',
        'Conventional prediction: DIFFERENT INTERNAL correction; IDENTICAL EXTERIOR correction.',
        'Model assumes spherical energy density without solving elastic equilibrium.',
        'No novel interaction field derived, observed or confirmed.',
        'Tests: PASS; scientific claim restricted to the specified numerical proxy.',
    ]
    (out / 'report.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    outputs = ['results.json', 'radial_energy_gravity.csv', 'report.txt', 'tests.txt']
    provenance = {
        'experiment_id': 'GED-V1-R6', 'run_id': run_id, 'utc_start': now.isoformat(),
        'python_executable': sys.executable, 'python_version': sys.version,
        'platform': platform.platform(), 'machine': platform.machine(),
        'dependencies': 'Python standard library, pinned GED-V1-R2 and GED-V1-R5',
        'command': [sys.executable, str(root / 'ged_v1_r6.py'), '--config', config_name],
        'test_command': cmd, 'test_exit_status': check.returncode,
        'git': git_state(root),
        'inputs_SHA256': hashes,
        'outputs_SHA256': {item: sha256(out / item) for item in outputs},
        'pinned_r2_expected_SHA256': R2_SHA256,
        'pinned_r5_expected_SHA256': R5_SHA256,
        'limitations': cfg['limitations'],
        'acceptance_status': 'LOCAL_NUMERICAL_TEST_PASS_PENDING_USER_REPRODUCTION_AND_AUDIT',
        'publication_policy': 'LOCAL_ONLY; next planned GitHub/Zenodo release at GED-V1-R10 or definite finding',
    }
    (out / 'provenance.json').write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    tracked = [f'input_snapshot/{name}' for name in INPUT_FILES] + outputs + ['provenance.json']
    (out / 'SHA256SUMS').write_text(''.join(
        f'{sha256(out / item)}  {item}\n' for item in sorted(tracked)), encoding='utf-8')
    print((out / 'report.txt').read_text(encoding='utf-8'))
    print('Results and receipts:', out)
    return out


def main():
    """WHAT: CLI runner. WHY: explicit local project config, no hidden I/O."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default='GED_V1_R6_config.json')
    args = parser.parse_args()
    if Path(args.config).name != args.config:
        parser.error('--config must be a filename within the GED folder')
    run_experiment(ROOT, args.config)


if __name__ == '__main__':
    main()
