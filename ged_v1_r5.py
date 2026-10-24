#!/usr/bin/env python3
"""GED-V1-R5 — stress-energy and material response.

WHAT: Compare how idealized iron and ice store elastic energy under
(a) the same applied stress, (b) the same imposed strain and
(c) the same stored elastic energy. Audit energy conservation when work
is stored, recovered, converted to heat or exported.

WHY: Separate material response (Young's modulus) from the conventional
mass-energy equivalent of retained deformation energy. This does NOT
implement or confirm a new gravitational information field.

MODELLING LIMITS: Linear elastic uniaxial Hooke-law material samples,
very weak fields, compact-source far-field point-monopole proxy, fixed
nominal sample volume. Full tensor stress-energy/GR and boundary supports
are outside this numerical scout. All numerical quantities are in SI units.

PYTHON: Standard library only, developed for Python 3.13.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import hashlib
import io
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import uuid

# WHAT: Reuse the *exact pinned* GED-V1-R2 mass-energy helpers.
# WHY: Avoid creating a second authoritative implementation of E/c^2.
# The source hash is checked before any scientific calculation proceeds.
import ged_v1_r2 as r2

PRECISION = 60
ZERO = Decimal('0')
ONE = Decimal('1')
TWO = Decimal('2')
R2_SHA256 = 'c1c8bd268fcfad2f8dea23bba703ef04284ced2590680584ccfa263cd9010813'
ROOT = Path(__file__).resolve().parent
INPUT_FILES = (
    'ged_v1_r5.py', 'test_ged_v1_r5.py', 'GED_V1_R5_config.json',
    'GED_V1_R5_README.md', 'GED_V1_R5_DEV_AUDIT.md',
    'GED_V1_R5_PRECISION_FAILURE_A.log', 'ged_v1_r2.py',
)


def D(value):
    """WHAT: Parse finite Decimal. WHY: preserve tiny E/c^2 corrections."""
    if isinstance(value, bool):
        raise ValueError('boolean is not a physical measurement')
    return r2.D(value)


def file_hash(path):
    """WHAT: SHA-256 exact bytes. WHY: immutable source/input receipts."""
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def check_dependency(root=ROOT):
    """WHAT: Verify the existing gravitational helper. WHY: no silent drift."""
    actual = file_hash(Path(root) / 'ged_v1_r2.py')
    if actual != R2_SHA256:
        raise ValueError('GED-V1-R2 dependency SHA-256 mismatch: ' + actual)
    return actual


def load_config(path=ROOT / 'GED_V1_R5_config.json'):
    """WHAT: Load fixed inputs. WHY: science runs must be reproducible."""
    cfg = json.loads(Path(path).read_text(encoding='utf-8'))
    if cfg.get('experiment_id') != 'GED-V1-R5':
        raise ValueError('incorrect experiment identifier')
    if cfg.get('r2_dependency_sha256') != R2_SHA256:
        raise ValueError('configuration requires wrong R2 version')
    return cfg


def validate(cfg):
    """WHAT: Parse and validate positive model parameters and SI values.

    WHY: Incorrect inputs should stop the experiment rather than silently
    create data that appears to describe physical behavior.
    """
    if cfg.get('experiment_id') != 'GED-V1-R5':
        raise ValueError('wrong experiment')
    names = ('G_m3_kg_s2', 'c_m_s', 'sample_volume_m3',
             'far_field_probe_m', 'equal_stress_Pa',
             'equal_strain_dimensionless', 'equal_elastic_energy_J')
    p = {name: D(cfg[name]) for name in names}
    if any(x <= ZERO for x in p.values()):
        raise ValueError('all model parameters must be positive')
    if p['equal_strain_dimensionless'] >= ONE / Decimal('100'):
        raise ValueError('small-strain scout restricts imposed strain < 0.01')
    # WHAT: Limit source size versus radius for a far-field illustration.
    # WHY: Near a deformed body, a point-monopole approximation is weak.
    with localcontext() as ctx:
        ctx.prec = PRECISION
        side_m = p['sample_volume_m3'] ** (ONE / Decimal(3))
        if p['far_field_probe_m'] < 10 * side_m:
            raise ValueError('probe must be at least 10 nominal cube side lengths away')
        if not isinstance(cfg.get('materials'), dict) or len(cfg['materials']) != 2:
            raise ValueError('exactly two materials required')
        materials = {}
        for key, data in cfg['materials'].items():
            modulus = D(data['young_modulus_Pa'])
            density = D(data['density_kg_m3'])
            if modulus <= ZERO or density <= ZERO:
                raise ValueError('material properties must be positive')
            materials[key] = {'E_Pa': modulus, 'rho_kg_m3': density}
        return p, materials


def elastic_state(E_Pa, volume_m3, *, stress_Pa=None, strain=None):
    """WHAT: Compute one reversible uniaxial Hookean deformation.

    WHY: sigma=E*epsilon; elastic energy density u=1/2*sigma*epsilon,
    equivalently sigma^2/(2E) or E*epsilon^2/2.

    Require EXACTLY ONE independent input: stress or dimensionless strain.
    Energy is J because Pa*m^3 = (N/m^2)*m^3 = N*m = J.
    """
    E = D(E_Pa)
    V = D(volume_m3)
    if E <= ZERO or V <= ZERO:
        raise ValueError('positive modulus and volume required')
    if (stress_Pa is None) == (strain is None):
        raise ValueError('provide exactly one: stress_Pa or strain')
    with localcontext() as ctx:
        ctx.prec = PRECISION
        if stress_Pa is not None:
            sigma = D(stress_Pa)
            if sigma < ZERO:
                raise ValueError('scout accepts tensile magnitude only')
            epsilon = sigma / E
        else:
            epsilon = D(strain)
            if epsilon < ZERO:
                raise ValueError('scout accepts tensile magnitude only')
            sigma = E * epsilon
        if epsilon >= Decimal('0.01'):
            raise ValueError('small-strain scope exceeded')
        energy_density = sigma * epsilon / TWO
        stored_J = energy_density * V
        return {
            'stress_Pa': +sigma,
            'strain_dimensionless': +epsilon,
            'elastic_energy_density_J_m3': +energy_density,
            'elastic_energy_J': +stored_J,
        }


def state_for_target_energy(E_Pa, volume_m3, energy_J):
    """WHAT: Find the strain required to store a specified energy.

    WHY: Provides a critical control: two different materials respond
    differently to mechanical loading but gain the SAME energy.
    Formula: epsilon=sqrt(2*U/(E*V)), for U=E*V*epsilon^2/2.
    """
    E, V, U = map(D, (E_Pa, volume_m3, energy_J))
    if min(E, V, U) <= ZERO:
        raise ValueError('E, V and target energy must be positive')
    with localcontext() as ctx:
        ctx.prec = PRECISION
        epsilon = (TWO * U / (E * V)).sqrt()
        return elastic_state(E, V, strain=epsilon)


def gravitational_proxy(energy_J, base_mass_kg, r_m, G, c):
    """WHAT: Compute an E/c^2 mass and distant monopole g correction.

    WHY: Ordinary mass-energy equivalence provides a BASELINE without
    adding a made-up GED information coupling. This correction should be
    compared ONLY at fixed external probe distance, as a bookkeeping proxy.

    WARNING: Uniaxial stress is anisotropic, loading changes geometry and
    needs supporting forces; this is NOT a GR tensor solution or a
    real laboratory signal predicted with complete boundary conditions.
    """
    U, M, r, g_const, speed = map(D, (energy_J, base_mass_kg, r_m, G, c))
    if U < ZERO or min(M, r, g_const, speed) <= ZERO:
        raise ValueError('invalid mass-energy or gravitational inputs')
    with localcontext() as ctx:
        ctx.prec = PRECISION
        dm = +r2.mass_equivalent(U, speed)
        dg = +r2.delta_g(U, r, g_const, speed)
        baseline_g = +r2.exterior_g(M, r, g_const)
        return {
            'base_mass_kg': +M,
            'point_monopole_baseline_g_m_s2': baseline_g,
            'delta_mass_if_energy_retained_kg': dm,
            'delta_g_point_monopole_proxy_m_s2': dg,
            'fractional_mass_increment': +(dm / M),
        }


def energy_ledgers(U):
    """WHAT: Track sample, loading reservoir and outside energy changes.

    WHY: Stored energy is transferred from elsewhere. A loaded sample is
    not an isolated creator of mass-energy. The SYSTEM ledger sums to 0.

    The thermalized state is a bookkeeping endpoint, not a rate or entropy
    model; the released energy stays within the sample, with no calculated
    entropy change because no temperature/path data are assumed.
    """
    energy = D(U)
    if energy < ZERO:
        raise ValueError('elastic energy cannot be negative')
    with localcontext() as ctx:
        ctx.prec = PRECISION
        cases = {
            'initial': (ZERO, ZERO, ZERO, ZERO),
            'loaded_energy_stored': (energy, ZERO, -energy, ZERO),
            'unloaded_reversibly': (ZERO, ZERO, ZERO, ZERO),
            'dissipated_as_retained_heat': (ZERO, energy, -energy, ZERO),
            'heat_exported_to_environment': (ZERO, ZERO, -energy, energy),
        }
        result = {}
        for name, (elastic, heat, reservoir, outside) in cases.items():
            total = elastic + heat + reservoir + outside
            if total != ZERO:
                raise AssertionError('non-conserving energy ledger')
            result[name] = {
                'sample_elastic_J': +elastic,
                'sample_heat_J': +heat,
                'loading_reservoir_J': +reservoir,
                'environment_J': +outside,
                'whole_system_balance_J': +total,
                'sample_total_change_J': +(elastic + heat),
            }
        return result


def calculate(cfg):
    """WHAT: Build all three controlled comparisons and energy ledgers.

    WHY: Comparing both loading protocols avoids interpreting material
    stiffness as a new force, and the fixed-energy comparison isolates
    mass-energy from the different material responses.
    """
    p, materials = validate(cfg)
    with localcontext() as ctx:
        ctx.prec = PRECISION
        probes = p['far_field_probe_m']
        V = p['sample_volume_m3']
        G, c = p['G_m3_kg_s2'], p['c_m_s']
        cases = {'equal_applied_stress': {}, 'equal_imposed_strain': {},
                 'equal_stored_elastic_energy': {}}
        for label, m in materials.items():
            E, rho = m['E_Pa'], m['rho_kg_m3']
            material_cases = {
                'equal_applied_stress': elastic_state(
                    E, V, stress_Pa=p['equal_stress_Pa']),
                'equal_imposed_strain': elastic_state(
                    E, V, strain=p['equal_strain_dimensionless']),
                'equal_stored_elastic_energy': state_for_target_energy(
                    E, V, p['equal_elastic_energy_J']),
            }
            for case_name, mechanical in material_cases.items():
                record = dict(mechanical)
                record.update(gravitational_proxy(
                    mechanical['elastic_energy_J'], rho * V, probes, G, c))
                record['energy_ledger'] = energy_ledgers(mechanical['elastic_energy_J'])
                cases[case_name][label] = record

        # WHAT: Compact machine-readable summary of the structural result.
        # WHY: Stiffness changes energy at fixed stress/strain; but matching
        # energy gives matching delta-g monopole proxies (despite different
        # stress/strain). This is not a new field's discovery.
        fixed = cases['equal_stored_elastic_energy']
        members = list(fixed)
        gravity_difference = (fixed[members[0]]['delta_g_point_monopole_proxy_m_s2']
                              - fixed[members[1]]['delta_g_point_monopole_proxy_m_s2'])
        return {
            'experiment_id': 'GED-V1-R5',
            'verdict': 'CONVENTIONAL_NUMERICAL_CONTROL_ONLY',
            'observational_data': False,
            'novel_information_channel': 'NOT_DERIVED_OR_TESTED',
            'energy_control_delta_g_difference_m_s2': gravity_difference,
            'cases': cases,
            'limits': list(cfg['limitations']),
        }


def json_safe(obj):
    """WHAT: Convert Decimals to exact decimal strings, recursively."""
    if isinstance(obj, Decimal):
        return str(obj)
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (tuple, list)):
        return [json_safe(v) for v in obj]
    return obj


def git_info(root):
    """WHAT: Record any local Git identity. WHY: hashes alone are not tags."""
    def git(*args):
        run = subprocess.run(['git', '-C', str(root), *args],
                             text=True, capture_output=True, check=False)
        return run.stdout.strip() if run.returncode == 0 else None
    commit = git('rev-parse', 'HEAD')
    if not commit:
        return {'commit': None, 'exact_tag': None, 'dirty': None}
    status = git('status', '--porcelain', '--untracked-files=all')
    return {
        'commit': commit,
        'exact_tag': git('describe', '--tags', '--exact-match'),
        'dirty': bool(status) if status is not None else None,
    }


def run_experiment(root=ROOT, cfg_name='GED_V1_R5_config.json'):
    """WHAT: Execute, retest saved inputs and write SHA-256 run receipts.

    WHY: An experiment is not accepted merely because a console summary
    looks plausible. Each run captures its precise source and environment.
    Old runs are never overwritten; archived experiments stay untouched.
    """
    root = Path(root).resolve()
    check_dependency(root)
    cfg = load_config(root / cfg_name)
    result = calculate(cfg)

    # WHAT: Assign a unique UTC execution identifier.
    # WHY: Avoid overwriting and distinguish separate repeat runs.
    now = datetime.now(timezone.utc)
    run_id = f'GED_V1_R5_{now.strftime("%Y%m%dT%H%M%S%fZ")}_{uuid.uuid4().hex[:8]}'
    out = root / 'runs' / run_id
    snapshots = out / 'input_snapshot'
    snapshots.mkdir(parents=True, exist_ok=False)

    # WHAT: Freeze all actually-used source and documentation bytes.
    # WHY: Future edits to working files must not change this run.
    for filename in INPUT_FILES:
        source = root / filename
        if not source.is_file():
            raise FileNotFoundError(f'missing necessary input: {source}')
        shutil.copy2(source, snapshots / filename)
    hashes = {f: file_hash(snapshots / f) for f in INPUT_FILES}

    # WHAT: Rerun the tests *from the frozen snapshot*.
    # WHY: Test the exact code that belongs to this run, not a later copy.
    cmd = [sys.executable, '-B', '-m', 'unittest', '-v', 'test_ged_v1_r5.py']
    check = subprocess.run(cmd, cwd=snapshots, capture_output=True,
                           text=True, check=False)
    (out / 'tests.txt').write_text(check.stdout + check.stderr, encoding='utf-8')
    if check.returncode != 0:
        (out / 'FAILED_TESTS.txt').write_text('Test command failed; run retained.\n',
                                               encoding='utf-8')
        raise RuntimeError(f'snapshot tests failed (run preserved): {out}')

    # WHAT: Serialize precise Decimal numbers as strings.
    # WHY: JSON floating numbers would lose tiny mass/gravity corrections.
    (out / 'results.json').write_text(
        json.dumps(json_safe(result), indent=2, sort_keys=True) + '\n',
        encoding='utf-8')

    # WHAT: Export a rectangular table useful for independent plotting.
    # WHY: Each case can be compared without importing our functions.
    fields = ['case', 'material', 'stress_Pa', 'strain_dimensionless',
              'elastic_energy_density_J_m3', 'elastic_energy_J',
              'delta_mass_if_energy_retained_kg',
              'delta_g_point_monopole_proxy_m_s2']
    with (out / 'elastic_cases.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for case_name, materials in result['cases'].items():
            for name, values in materials.items():
                writer.writerow({k: str(values.get(k, '')) if k not in ('case', 'material')
                                 else (case_name if k == 'case' else name)
                                 for k in fields})

    # WHAT: Write a short, human-readable summary.
    # WHY: Terminal output shouldn't be the only place the result survives.
    lines = [
        'GED-V1-R5 — STRESS-ENERGY AND MATERIAL RESPONSE',
        '=' * 64,
        'CONVENTIONAL ELASTICITY / MASS-ENERGY NUMERICAL CONTROL',
        'NO NEW GRAVITY DETECTION OR INFORMATION-FIELD COUPLING',
    ]
    for case, data in result['cases'].items():
        lines.append(f'\n{case}:')
        for material, value in data.items():
            lines.append(f'  {material}: stress={value["stress_Pa"]} Pa; '
                         f'strain={value["strain_dimensionless"]}; '
                         f'elastic energy={value["elastic_energy_J"]} J; '
                         f'delta g proxy={value["delta_g_point_monopole_proxy_m_s2"]} m/s^2')
    lines += [
        '\nEqual-energy control monopole delta-g difference: '
        + str(result['energy_control_delta_g_difference_m_s2']) + ' m/s^2',
        'Conservation: all sample/reservoir/environment energy ledgers sum to zero.',
        'Interpretation: different elastic responses do not establish extra gravity data.',
        'Full stress-energy tensor and Einstein equations are not solved.',
        'Tests: PASS',
        'Verdict: numerical scout; independently observed effect: NONE.',
    ]
    (out / 'report.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')

    outputs = ['results.json', 'elastic_cases.csv', 'report.txt', 'tests.txt']
    provenance = {
        'experiment_id': 'GED-V1-R5', 'run_id': run_id,
        'utc_start': now.isoformat(),
        'python_executable': sys.executable,
        'python_version': sys.version,
        'platform': platform.platform(), 'machine': platform.machine(),
        'dependencies': 'Python standard library; exact SHA-pinned GED-V1-R2 source',
        'command': [sys.executable, str(root / 'ged_v1_r5.py'), '--config', cfg_name],
        'test_command': cmd, 'test_exit_status': check.returncode,
        'git': git_info(root),
        'input_SHA256': hashes,
        'output_SHA256': {f: file_hash(out / f) for f in outputs},
        'r2_dependency_expected_SHA256': R2_SHA256,
        'r2_dependency_actual_SHA256': hashes['ged_v1_r2.py'],
        'limitations': cfg['limitations'],
        'acceptance_status': 'LOCAL_NUMERICAL_TEST_PASS_PENDING_EXTERNAL_AUDIT',
    }
    (out / 'provenance.json').write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + '\n', encoding='utf-8')

    # WHAT: Manifest covers the frozen inputs AND all computed outputs.
    # WHY: An auditor can detect any later accidental or deliberate change.
    tracked = [f'input_snapshot/{name}' for name in INPUT_FILES] + outputs + ['provenance.json']
    (out / 'SHA256SUMS').write_text(
        ''.join(f'{file_hash(out / name)}  {name}\n' for name in sorted(tracked)),
        encoding='utf-8')
    print((out / 'report.txt').read_text(encoding='utf-8'))
    print('Results and receipts:', out)
    return out


def main():
    """WHAT: CLI entry point. WHY: deterministic inputs, explicit command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default='GED_V1_R5_config.json',
                        help='configuration filename within the GED directory')
    args = parser.parse_args()
    # WHAT: Refuse arbitrary external paths for official experiment inputs.
    # WHY: Preserve one unambiguous project-local configuration per run.
    if Path(args.config).name != args.config:
        parser.error('--config must be a filename in the GED folder')
    run_experiment(ROOT, args.config)


if __name__ == '__main__':
    main()
