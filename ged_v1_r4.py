#!/usr/bin/env python3
"""GED-V1-R4 — equal energy, different entropy (numerical controls).

This is an established-physics *weak-field* exercise, not an observation or a
new gravity law. SI units throughout. Python standard library only.

WHAT this asks: Can entropy differ without a change in an object's external
Newtonian gravitational field? And can internal gravity change without a
change in the total energy of the idealised object?

WHY two controls: An equal-energy exterior null does not by itself isolate
entropy; a redistribution can alter the *interior* energy profile. We
therefore include a separate degenerate-state control in which the local
stress-energy description is assumed identical, and an equal-energy thermal
redistribution control whose radial distribution really does change.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import hashlib
import io
import json
from math import comb
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import uuid

# WHAT: Import the project's pinned, authoritative gravitational helpers.
# WHY: Do not duplicate constants or formulas from earlier GED scouts.
# This imported file is SHA-pinned and captured in every run's receipts.
import ged_v1_r2 as r2

PRECISION = 60
ZERO = Decimal(0)
ONE = Decimal(1)
TWO = Decimal(2)
K_B = Decimal('1.380649e-23')  # Exact SI Boltzmann constant, J/K.


def D(value):
    """WHAT: Convert to a finite Decimal. WHY: keep tiny E/c² effects."""
    return r2.D(value)


def sha256(path):
    """WHAT: Hash exact file bytes. WHY: identify the code/data actually run."""
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def load_config(path):
    """WHAT: Read the fixed experimental plan. WHY: reproducible parameters."""
    obj = json.loads(Path(path).read_text(encoding='utf-8'))
    if obj.get('experiment_id') != 'GED-V1-R4':
        raise ValueError('expected GED-V1-R4 config')
    return obj


def read_parameters(cfg):
    """WHAT: Validate dimensions, temperatures, and macrostate controls.

    WHY: An invalid parameter must fail before any result is interpreted.
    Temperatures are absolute kelvins; constant heat capacity is assumed.
    """
    keys = ('G_m3_kg_s2', 'c_m_s', 'radius_m', 'density_kg_m3',
            'specific_heat_J_kg_K', 'temperature_reference_K',
            'temperature_offset_K', 'probe_interior_fraction_R',
            'probe_exterior_multiple_R')
    p = {k: D(cfg[k]) for k in keys}
    if min(p[k] for k in keys if k not in
           ('temperature_offset_K', 'probe_interior_fraction_R')) <= 0:
        raise ValueError('physical constants and dimensions must be positive')
    if p['temperature_offset_K'] <= 0:
        raise ValueError('offset must be positive to test a nontrivial gradient')
    if p['temperature_reference_K'] <= p['temperature_offset_K']:
        raise ValueError('all absolute temperatures must remain positive')
    if not (ZERO < p['probe_interior_fraction_R'] < ONE):
        raise ValueError('internal probe must be inside sphere, away from origin')
    if p['probe_exterior_multiple_R'] < ONE:
        raise ValueError('external probe must be at/above surface')

    # WHAT: Two equal-volume zones. WHY: equal density and cp then imply
    # equal thermal capacities, so +/- temperature offsets conserve heat.
    if D(cfg['core_volume_fraction']) != Decimal('0.5'):
        raise ValueError('equal-capacity zones require half-volume core')

    # WHAT: Count microscopic binary labels for the entropy-only control.
    # WHY: Two macrostates can have different multiplicities with equal energy.
    # This does NOT assert a specific real atomic spin system has this property.
    n, k_a, k_b = (cfg[x] for x in
                     ('degenerate_site_count', 'macrostate_A_label_one_count',
                      'macrostate_B_label_one_count'))
    if any(isinstance(x, bool) or not isinstance(x, int) for x in (n, k_a, k_b)):
        raise ValueError('state counts must be integers, not booleans')
    if n < 2 or n > 100000 or not (0 <= k_a <= n and 0 <= k_b <= n):
        raise ValueError('invalid site counts')
    if k_a == k_b:
        raise ValueError('macrostate counts must differ')
    p['n'] = n
    p['k_a'] = k_a
    p['k_b'] = k_b
    return p


def microstate_entropy(n, k):
    """WHAT: S=k_B*ln(C(n,k)), Boltzmann entropy in J/K.

    WHY: A constrained macrostate with k occupied binary labels among n
    sites admits binomial(n,k) equally energetic micro-configurations.
    This is a *coarse-grained macrostate-counting* model, NOT the entropy
    of an individually specified pure microstate or a device's readout.

    Assumption: labels have exactly degenerate energy and contribute
    identical local stress-energy; interactions and external fields are off.
    """
    if not (isinstance(n, int) and isinstance(k, int) and
            not isinstance(n, bool) and not isinstance(k, bool) and
            0 <= k <= n):
        raise ValueError('invalid macrostate arguments')
    with localcontext() as ctx:
        ctx.prec = PRECISION
        multiplicity = comb(n, k)
        return +(K_B * Decimal(multiplicity).ln())


def spherical_model(p):
    """WHAT: Construct the idealised two-zone body at reference temperature.

    WHY: We need *identical rest-material geometry* in both thermal cases.
    Material density is stipulated fixed, without expansion or phase change.
    Walls/pressure/gravitational self-energy changes are outside this scout.
    """
    with localcontext() as ctx:
        ctx.prec = PRECISION
        R = p['radius_m']
        f = Decimal('0.5')
        rho = p['density_kg_m3']
        # WHAT: Reuse pinned spherical geometry from GED-V1-R2.
        # WHY: Identical volume and density -> identical base mass per zone.
        volume = +r2.sphere_volume(R)
        core_radius = +(R * f ** (ONE / Decimal(3)))
        mass_zone = +(rho * volume * f)
        heat_capacity = +(mass_zone * p['specific_heat_J_kg_K'])
        return {'R': R, 'a': core_radius, 'volume': volume,
                'mass_zone': mass_zone, 'C': heat_capacity,
                'reference_mass': +(TWO * mass_zone)}


def thermal_state(model, p, T_core, T_shell):
    """WHAT: Translate an internal temperature pattern into mass-energy.

    WHY: Ordinary relativity already associates a mass equivalent with
    retained thermal energy. Our *fixed reference* mass includes energy at
    T_reference; we add ONLY deviations C*(T-T_reference)/c^2.

    CAUTION: This is a Newtonian mass-equivalent bookkeeping approximation;
    no pressure/stress tensor, hydrostatic balance or full GR field solution.
    """
    T_c, T_s = D(T_core), D(T_shell)
    T_0 = p['temperature_reference_K']
    if T_c <= 0 or T_s <= 0:
        raise ValueError('temperature must remain above absolute zero')
    with localcontext() as ctx:
        ctx.prec = PRECISION
        C = model['C']
        dE_c = +(C * (T_c - T_0))
        dE_s = +(C * (T_s - T_0))
        dM_c = +r2.mass_equivalent(dE_c, p['c_m_s'])
        dM_s = +r2.mass_equivalent(dE_s, p['c_m_s'])
        M_c = +(model['mass_zone'] + dM_c)
        M_s = +(model['mass_zone'] + dM_s)
        if M_c <= 0 or M_s <= 0:
            raise ValueError('model mass not positive')
        # WHAT: Exact reversible entropy difference for constant heat capacity.
        # WHY: Integrate dS=C dT/T independently of any gravity rule.
        entropy_change = +(C * ((T_c/T_0).ln() + (T_s/T_0).ln()))
        return {
            'T_core_K': T_c, 'T_shell_K': T_s,
            'dE_core_J': dE_c, 'dE_shell_J': dE_s,
            'dE_total_J': +(dE_c + dE_s),
            'dM_core_kg': dM_c, 'dM_shell_kg': dM_s,
            'M_core_kg': M_c, 'M_shell_kg': M_s,
            'M_total_kg': +(M_c + M_s),
            'dS_thermal_J_K': entropy_change,
        }


def enclosed_mass(state, model, radius):
    """WHAT: Sum mass-energy interior to spherical radius r (kg).

    WHY: Newton's shell theorem says spherical shells outside r contribute
    no gravitational acceleration at r; they still influence potential.
    Within each zone, mass-energy density is treated as uniform.
    """
    r = D(radius)
    R, a = model['R'], model['a']
    if r < 0:
        raise ValueError('probe radius must be nonnegative')
    with localcontext() as ctx:
        ctx.prec = PRECISION
        if r <= a:
            return +(state['M_core_kg'] * (r/a) ** 3)
        if r < R:
            fraction_shell = +(r**3 - a**3) / (R**3 - a**3)
            return +(state['M_core_kg'] + state['M_shell_kg']*fraction_shell)
        return state['M_total_kg']


def radial_g(state, model, radius, G):
    """WHAT: Calculate inward Newtonian g=G*M_enc/r^2 in m/s².

    WHY: We compare a *radially internal* probe with an *exterior* probe.
    This is not a dynamical geodesic calculation or exact GR metric.
    """
    r = D(radius)
    g_const = D(G)
    if g_const <= 0:
        raise ValueError('G must be positive')
    if r < 0:
        raise ValueError('radius must be nonnegative')
    if r == ZERO:
        return ZERO  # Symmetry forces centre acceleration to vanish.
    with localcontext() as ctx:
        ctx.prec = PRECISION
        return +(g_const * enclosed_mass(state, model, r) / (r*r))


def calculate(cfg):
    """WHAT: Evaluate two independent controls, with no file side-effects.

    Control A: Different combinatorial macrostates but stipulated identical
    local energy/stress distribution. Gravity must be an exact null.

    Control B: Hot core/cold shell -> uniform temperature via heat exchange.
    Total internal thermal energy is conserved in this approximate ledger,
    entropy increases, internal gravity changes, and exterior gravity does
    not change. No energy escapes and no new coupling is assumed.

    Control C: Reversing hot and cold zones gives equal entropy and equal
    total energy, yet reverses the interior gravitational difference.
    """
    with localcontext() as ctx:
        ctx.prec = PRECISION
        p = read_parameters(cfg)
        body = spherical_model(p)
        T0, dT = p['temperature_reference_K'], p['temperature_offset_K']
        gradient = thermal_state(body, p, T0+dT, T0-dT)
        reverse = thermal_state(body, p, T0-dT, T0+dT)
        uniform = thermal_state(body, p, T0, T0)
        interior_r = +(body['R'] * p['probe_interior_fraction_R'])
        exterior_r = +(body['R'] * p['probe_exterior_multiple_R'])
        radii = [interior_r, body['a'], body['R'], exterior_r]
        # WHAT: Compare g before and after thermal relaxation at each radius.
        # WHY: Find exactly where energy distribution, not entropy alone,
        # changes the gravity result.
        samples = []
        for radius in radii:
            g_grad = radial_g(gradient, body, radius, p['G_m3_kg_s2'])
            g_uniform = radial_g(uniform, body, radius, p['G_m3_kg_s2'])
            samples.append({
                'radius_m': str(+radius),
                'radius_over_R': str(+(radius/body['R'])),
                'g_gradient_m_s2': str(g_grad),
                'g_uniform_m_s2': str(g_uniform),
                'delta_g_uniform_minus_gradient_m_s2': str(+(g_uniform-g_grad)),
            })

        # WHAT: Compare the two opposite temperature gradients too.
        # WHY: They have *equal entropy* (same pair of temperatures) but
        # different radial mass-energy profiles, a useful counter-control.
        reverse_probe_g = radial_g(reverse, body, body['a'], p['G_m3_kg_s2'])
        gradient_probe_g = radial_g(gradient, body, body['a'], p['G_m3_kg_s2'])
        S_a = microstate_entropy(p['n'], p['k_a'])
        S_b = microstate_entropy(p['n'], p['k_b'])
        delta_S_relax = +(uniform['dS_thermal_J_K'] - gradient['dS_thermal_J_K'])
        # WHAT: Draw a clean energy ledger for the internal redistribution.
        # WHY: Transfer from hot core to cold shell does not export energy.
        core_transfer = +(uniform['dE_core_J']-gradient['dE_core_J'])
        shell_transfer = +(uniform['dE_shell_J']-gradient['dE_shell_J'])
        result = {
            'experiment_id': 'GED-V1-R4',
            'classification': 'CONVENTIONAL_WEAK_FIELD_ENTROPY_CONTROL',
            'observation_performed': False,
            'new_gravity_channel_detected': False,
            'material': {
                'radius_m': str(body['R']),
                'core_radius_m': str(body['a']),
                'reference_mass_kg': str(body['reference_mass']),
                'heat_capacity_each_zone_J_K': str(body['C']),
                'reference_temperature_K': str(T0),
            },
            'control_A_equal_local_stress_energy': {
                'assumptions': 'perfectly degenerate noninteracting binary labels; identical local stress-energy by construction',
                'site_count': p['n'],
                'macrostate_A_count': p['k_a'],
                'macrostate_B_count': p['k_b'],
                'multiplicity_A': str(comb(p['n'], p['k_a'])),
                'multiplicity_B': str(comb(p['n'], p['k_b'])),
                'entropy_A_J_K': str(S_a),
                'entropy_B_J_K': str(S_b),
                'entropy_B_minus_A_J_K': str(+(S_b-S_a)),
                'energy_difference_J': '0',
                'gravitational_field_difference_m_s2': '0',
            },
            'control_B_thermal_redistribution': {
                'initial': {k: str(v) for k, v in gradient.items()},
                'relaxed': {k: str(v) for k, v in uniform.items()},
                'energy_exported_J': '0',
                'thermal_energy_difference_final_minus_initial_J': str(+(uniform['dE_total_J']-gradient['dE_total_J'])),
                'core_thermal_energy_transfer_J': str(core_transfer),
                'shell_thermal_energy_transfer_J': str(shell_transfer),
                'transfer_energy_balance_J': str(+(core_transfer+shell_transfer)),
                'entropy_produced_in_relaxation_J_K': str(delta_S_relax),
                'model_exterior_mass_difference_kg': str(+(uniform['M_total_kg']-gradient['M_total_kg'])),
                'interior_g_delta_at_core_boundary_m_s2': next(s['delta_g_uniform_minus_gradient_m_s2'] for s in samples if s['radius_m'] == str(body['a'])),
                'exterior_g_delta_at_probe_m_s2': samples[-1]['delta_g_uniform_minus_gradient_m_s2'],
            },
            'control_C_equal_entropy_different_interior_distribution': {
                'reverse_gradient': {k: str(v) for k, v in reverse.items()},
                'entropy_difference_reverse_minus_initial_J_K': str(+(reverse['dS_thermal_J_K']-gradient['dS_thermal_J_K'])),
                'energy_difference_reverse_minus_initial_J': str(+(reverse['dE_total_J']-gradient['dE_total_J'])),
                'model_exterior_mass_difference_kg': str(+(reverse['M_total_kg']-gradient['M_total_kg'])),
                'core_boundary_g_difference_reverse_minus_initial_m_s2': str(+(reverse_probe_g-gradient_probe_g)),
            },
            'sampling': samples,
            'limitations': [
                'No observational or laboratory data; this is a controlled numerical thought experiment.',
                'Control A uses stipulated exactly degenerate binary internal labels; it does not model an actual spin material or spontaneous dynamics.',
                'Boltzmann entropy uses macrostate multiplicity, not entropy assigned to a fully specified pure microstate.',
                'Control B assumes fixed geometry, equal constant heat capacities, and uniform energy density inside each radial zone.',
                'Heat flow is prescribed; rate, conductivity, transport stresses and equilibration dynamics are not simulated.',
                'Thermal expansion, pressure/stresses, gravitational potential-energy work and exact GR equilibrium are ignored.',
                'The fixed reference mass includes reference thermal energy; only deviations C*(T-Tref)/c^2 are added.',
                'Equality of exterior fields refers to the Newtonian spherical monopole and fixed total mass-energy in this approximate model.',
                'In complete GR, total gravitating mass includes the full stress-energy and boundary contributions; this is not a solved Einstein/TOV model.',
                'An entropy-only gravitational coupling is neither derived nor tested; no parameter is invented to force an effect.',
            ],
        }
        return result, samples


def write_json(path, value):
    """WHAT: Persist stable JSON. WHY: machine-auditable experiment receipt."""
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True,
                                     ensure_ascii=False, allow_nan=False) + '\n',
                          encoding='utf-8')


def git_provenance(project):
    """WHAT: Capture existing Git commit/tag/dirty status without modifying Git."""
    def git(*args):
        try:
            return subprocess.run(['git', *args], cwd=project, capture_output=True,
                                  text=True, check=True, timeout=6).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return None
    root = git('rev-parse', '--show-toplevel')
    return {'root': root,
            'commit': git('rev-parse', 'HEAD') if root else None,
            'exact_tag': git('describe', '--tags', '--exact-match') if root else None,
            'dirty': bool(git('status', '--porcelain')) if root else None}


def run(project, config_path):
    """WHAT: Freeze exact inputs, run tests, then write results and SHA receipts.

    WHY: The laboratory requires reproducibility and preserved failures.
    A failing test blocks science output but is still retained as evidence.
    This is the experiment's local receipt generator, NOT a claim that the
    external Turnip Run wrapper was executed.
    """
    cfg_path = Path(config_path).resolve()
    cfg = load_config(cfg_path)
    pinned = cfg['r2_dependency_sha256']
    dependency = project / 'ged_v1_r2.py'
    actual = sha256(dependency)
    if actual != pinned:
        raise RuntimeError(f'pinned V1-R2 source mismatch: expected {pinned}, got {actual}')
    # WHAT: Calculate before opening a run directory, catching bad input early.
    result, samples = calculate(cfg)
    start = datetime.now(timezone.utc)
    rid = 'GED_V1_R4_' + start.strftime('%Y%m%dT%H%M%S%fZ') + '_' + uuid.uuid4().hex[:8]
    folder = project / 'runs' / rid
    folder.mkdir(parents=True, exist_ok=False)
    snapshot = folder / 'input_snapshot'
    snapshot.mkdir()
    files = ['ged_v1_r2.py', 'ged_v1_r4.py', 'test_ged_v1_r4.py',
             'GED_V1_R4_README.md', 'GED_V1_R4_DEV_AUDIT.md',
             'GED_V1_R4_PRECISION_FAILURE_A.log']
    for filename in files:
        shutil.copy2(project / filename, snapshot / filename)
    shutil.copy2(cfg_path, snapshot / 'GED_V1_R4_config.json')

    # WHAT: Test only the frozen snapshot, not editable working files.
    # WHY: Exactly these bytes must earn the result acceptance.
    command = [sys.executable, '-B', '-m', 'unittest', '-v', 'test_ged_v1_r4.py']
    test = subprocess.run(command, cwd=snapshot, capture_output=True, text=True)
    (folder / 'tests.txt').write_text(test.stdout + test.stderr, encoding='utf-8')

    if test.returncode == 0:
        write_json(folder / 'results.json', result)
        csv_out = io.StringIO(newline='')
        writer = csv.DictWriter(csv_out, fieldnames=list(samples[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(samples)
        (folder / 'radial_gravity.csv').write_text(csv_out.getvalue(), encoding='utf-8')
        a = result['control_A_equal_local_stress_energy']
        b = result['control_B_thermal_redistribution']
        c = result['control_C_equal_entropy_different_interior_distribution']
        report = '\n'.join([
            'GED-V1-R4 — EQUAL ENERGY / DIFFERENT ENTROPY',
            '='*64,
            'CONVENTIONAL NUMERICAL CONTROL; NOT NEW GRAVITY OR AN OBSERVATION',
            f"Reference mass: {D(result['material']['reference_mass_kg']):.12E} kg",
            f"Control A macrostate entropy difference: {D(a['entropy_B_minus_A_J_K']):.12E} J/K",
            'Control A mass-energy difference: 0 J; gravity difference: 0 m/s^2 (assumed degenerate model)',
            f"Control B entropy produced on relaxation: {D(b['entropy_produced_in_relaxation_J_K']):.12E} J/K",
            f"Control B energy transfer out of hot core: {-D(b['core_thermal_energy_transfer_J']):.12E} J",
            f"Control B interior g difference at zone boundary: {D(b['interior_g_delta_at_core_boundary_m_s2']):.12E} m/s^2",
            f"Control B exterior g difference: {b['exterior_g_delta_at_probe_m_s2']} m/s^2",
            f"Control C equal-entropy reverse-gradient internal delta g: {D(c['core_boundary_g_difference_reverse_minus_initial_m_s2']):.12E} m/s^2",
            'Control C entropy difference: 0 J/K; energy difference: 0 J; exterior g difference: 0 (model)',
            'No thermal energy exported; global energy ledger is zero in this approximate control.',
            'Pressure, gravitational self-energy shifts and full GR effects are outside model.',
            'No independent entropy-coupled gravitational information channel derived or detected.',
            'Tests: PASS',
        ]) + '\n'
    else:
        report = 'GED-V1-R4 — TESTS FAILED; scientific outputs quarantined\n'
    (folder / 'report.txt').write_text(report, encoding='utf-8')

    # WHAT: Capture all runtime versions, config/code fingerprints and outputs.
    # WHY: Manifest of Destiny requires an unambiguous provenance chain.
    hashes_before = {p.name: sha256(p) for p in sorted(snapshot.iterdir())}
    output_names = sorted(p.name for p in folder.iterdir() if p.is_file())
    provenance = {
        'experiment_id': 'GED-V1-R4', 'run_id': rid,
        'utc_start': start.isoformat(),
        'status': 'CONVENTIONAL_NUMERICAL_CONTROLS_PASS' if test.returncode == 0 else 'FAILED_TESTS_QUARANTINED',
        'test_command': command, 'test_exit_status': test.returncode,
        'launch_arguments': sys.argv,
        'python_version': sys.version, 'python_executable': sys.executable,
        'platform': platform.platform(), 'machine': platform.machine(),
        'dependencies': 'Python stdlib plus exact pinned GED-V1-R2 source',
        'r2_dependency_expected_sha256': pinned,
        'r2_dependency_actual_sha256': actual,
        'git': git_provenance(project),
        'input_SHA256': hashes_before,
        'output_SHA256': {name: sha256(folder/name) for name in output_names},
        'limitations': result['limitations'],
    }
    write_json(folder / 'provenance.json', provenance)
    # WHAT: Hash all files in this particular run (other than self-manifest).
    # WHY: A future verifier can detect any modification of recorded evidence.
    receipt_lines = []
    for path in sorted(folder.rglob('*')):
        if path.is_file() and path.name != 'SHA256SUMS':
            receipt_lines.append(f'{sha256(path)}  {path.relative_to(folder).as_posix()}\n')
    (folder / 'SHA256SUMS').write_text(''.join(receipt_lines), encoding='utf-8')
    print(report, end='')
    print('Results and receipts:', folder)
    if test.returncode:
        raise RuntimeError('Tests failed; quarantined GED-V1-R4 run')
    return folder


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default='GED_V1_R4_config.json')
    args = parser.parse_args()
    project = Path(__file__).resolve().parent
    run(project, Path(args.config))


if __name__ == '__main__':
    main()
