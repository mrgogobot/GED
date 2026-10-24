#!/usr/bin/env python3
"""GED-V1-R3: entropic exhaust and structural persistence (SI units).

ESTABLISHED-PHYSICS MODEL; NOT a test of an extra gravitational information channel.
Requires the pinned GED-V1-R2 calculation module in the same directory.
Python standard library only; designed for Python 3.13+.
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

# WHAT: Reuse the previously audited binding-energy and gravitational formulas.
# WHY: GED has one authoritative implementation of those formulas rather than
# copying them into every experiment and risking silent numerical divergence.
import ged_v1_r2 as r2

# WHAT: Use Decimal at high precision for small mass/field differences.
# WHY: E/c^2 for this experiment is ~10^-10 kg, tiny beside a ~10^10 kg body.
PRECISION = 60
ZERO = Decimal("0")
ONE = Decimal("1")


def dec(value):
    """WHAT: Parse a finite number. WHY: block NaN, Infinity and float surprises."""
    return r2.D(value)


def sha256(path):
    """WHAT: Hash exact file bytes. WHY: reproduce the code and records later."""
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_config(path):
    """WHAT: Load a fixed run plan. WHY: inputs are separate from the source."""
    cfg = json.loads(Path(path).read_text(encoding="utf-8"))
    if cfg.get("experiment_id") != "GED-V1-R3":
        raise ValueError("wrong experiment_id")
    return cfg


def binding_transition(cfg):
    """WHAT: Compare equal-inventory spheres before and after layer swapping.

    WHY: The previous scout demonstrated that potential self-binding depends on
    the internal density arrangement even when total constituent mass is fixed.

    LIMITATION: We prescribe B (ice core) -> A (iron core); we do not simulate
    the mechanical process, barriers, pressure, chemistry or external work.
    Its potential-energy difference is merely an AVAILABLE energy scale.
    """
    G = dec(cfg["G_m3_kg_s2"])
    c = dec(cfg["c_m_s"])
    R = dec(cfg["radius_m"])
    fraction = dec(cfg["core_volume_fraction"])
    rho_fe = dec(cfg["density_iron_kg_m3"])
    rho_ice = dec(cfg["density_ice_kg_m3"])
    cp_fe = dec(cfg["heat_capacity_iron_J_kg_K"])
    cp_ice = dec(cfg["heat_capacity_ice_J_kg_K"])
    T = dec(cfg["initial_temperature_K"])
    bath_T = dec(cfg["environment_temperature_K"])
    probe_multiple = dec(cfg["report_probe_R"])

    # WHAT: Require half-volume layers; swapping them then preserves quantities.
    # WHY: Any other volume fraction generally changes each material's mass.
    if fraction != Decimal("0.5"):
        raise ValueError("equal material inventory requires core_volume_fraction=0.5")
    if min(G, c, R, rho_fe, rho_ice, cp_fe, cp_ice, T, bath_T) <= 0:
        raise ValueError("all physical inputs must be positive and finite")
    if probe_multiple < ONE:
        raise ValueError("probe must lie at or outside the sphere")

    with localcontext() as ctx:
        ctx.prec = PRECISION
        # WHAT: Build the two arrangements through pinned V1-R2 source.
        # WHY: Preserve one authoritative spherical binding-energy calculation.
        A = r2.make_body(R, fraction, rho_fe, rho_ice)
        B = r2.make_body(R, fraction, rho_ice, rho_fe)
        if A["M0"] != B["M0"]:
            raise AssertionError("unequal material inventory")
        U_A = r2.binding_energy(A, G)
        U_B = r2.binding_energy(B, G)
        release = +(U_B - U_A)
        if release <= 0:
            raise AssertionError("chosen B->A energy scale must be positive")

        # WHAT: C=Σ(m_i*c_p,i), total effective heat capacity in J/K.
        # WHY: For a retained amount of heat Q, temperature rises by Q/C.
        # Heat capacities are illustrative constants and independent of T.
        thermal_capacity = +(A["mc"] * cp_fe + A["ms"] * cp_ice)
        probe_r = +(R * probe_multiple)
        baseline_g = +(G * A["M0"] / probe_r**2)
        return {
            "U_A_J": +U_A, "U_B_J": +U_B,
            "release_J": release, "constituent_mass_kg": +A["M0"],
            "heat_capacity_J_K": thermal_capacity,
            "initial_temperature_K": T, "environment_temperature_K": bath_T,
            "probe_radius_m": probe_r, "baseline_g_m_s2": baseline_g,
            "G_m3_kg_s2": G, "c_m_s": c,
        }


def scenario(core, escaped_fraction, *, stable=False):
    """WHAT: Keep energy and entropy as TWO SEPARATE ledgers.

    WHY: Entropy production is not synonymous with energy leaving an object.
    A stable bound state need not produce a continuous energy 'exhaust'.

    Assumed process (non-stable cases):
      - Potential energy falls by Q=U_before-U_after (>0).
      - f*Q enters an ideal constant-temperature environmental heat bath.
      - (1-f)*Q thermalises inside an ideal constant-heat-capacity body.
      - Radiation/recoil/work/geometry are not dynamically simulated.
      - Additional structural/configurational entropy is OMITTED, so the
        computed entropy is a THERMAL ledger, not the full physical entropy.
    """
    f = dec(escaped_fraction)
    if f < ZERO or f > ONE:
        raise ValueError("escaped_fraction must lie between zero and one")
    with localcontext() as ctx:
        ctx.prec = PRECISION
        release = ZERO if stable else core["release_J"]
        Q_env = +(release * f)
        Q_body = +(release - Q_env)
        C = core["heat_capacity_J_K"]
        T0 = core["initial_temperature_K"]
        Tb = core["environment_temperature_K"]
        c = core["c_m_s"]
        G = core["G_m3_kg_s2"]
        r = core["probe_radius_m"]

        # WHAT: Thermal body gets warmer in proportion to retained heat.
        # WHY: Internal gravitational potential energy can become heat.
        T1 = +(T0 + Q_body / C)

        # WHAT: Integrate dS=C*dT/T for constant C and equilibrium states.
        # WHY: Entropy depends on energy distribution, not simply mass loss.
        delta_S_body = +(C * (T1 / T0).ln()) if Q_body else ZERO

        # WHAT: Energy deposited into an ideal, equilibrium external bath.
        # WHY: A constant-temperature bath gains ΔS=Q/T_bath.
        # IMPORTANT: This is NOT an entropy formula for unthermalised photons.
        delta_S_bath = +(Q_env / Tb)
        delta_S_thermal = +(delta_S_body + delta_S_bath)

        # WHAT: Gravitational potential decreases; heat and export increase.
        # WHY: Confirm ΔU + ΔE_internal + ΔE_outside = 0.
        delta_U = -release
        delta_internal_heat = Q_body
        delta_outside = Q_env
        energy_residual = +(delta_U + delta_internal_heat + delta_outside)

        # WHAT: The body alone loses exported energy; surroundings gain it.
        # WHY: Its weak-field gravitational mass changes by -Q_out/c^2.
        delta_mass_body = +(-Q_env / (c*c))
        delta_g = +(G * delta_mass_body / (r*r))
        fraction_g = +(delta_mass_body / core["constituent_mass_kg"])
        global_mass_energy_residual = +(delta_mass_body + Q_env/(c*c))

        return {
            "kind": "stable_reference" if stable else "prescribed_B_to_A_rearrangement",
            "escaped_fraction": str(f),
            "potential_energy_change_J": str(delta_U),
            "thermal_energy_retained_J": str(Q_body),
            "energy_to_bath_J": str(Q_env),
            "energy_balance_residual_J": str(energy_residual),
            "initial_temperature_K": str(T0),
            "final_body_temperature_K": str(T1),
            "body_entropy_change_J_K": str(delta_S_body),
            "bath_entropy_change_J_K": str(delta_S_bath),
            "modeled_thermal_entropy_change_J_K": str(delta_S_thermal),
            "body_mass_change_kg": str(delta_mass_body),
            "outside_mass_energy_equivalent_kg": str(+(Q_env/(c*c))),
            "global_mass_energy_residual_kg": str(global_mass_energy_residual),
            "exterior_delta_g_m_s2": str(delta_g),
            "exterior_fractional_delta_g": str(fraction_g),
        }


def calculate(cfg):
    """WHAT: Pure model with fixed scenarios and a predeclared fraction sweep.

    WHY: All calculations can be independently tested without writing files.
    """
    with localcontext() as ctx:
        ctx.prec = PRECISION
        core = binding_transition(cfg)
        sweep = [dec(v) for v in cfg["escaped_fractions"]]
        if not sweep or len(sweep) != len(set(sweep)):
            raise ValueError("escaped_fractions must be nonempty and unique")
        if sweep != sorted(sweep) or sweep[0] != ZERO or sweep[-1] != ONE:
            raise ValueError("sweep must be sorted and include 0 and 1")
        # WHY: The report explicitly contains a 50% retained/exported comparison.
        if Decimal("0.5") not in sweep:
            raise ValueError("sweep must include 0.5 for the fixed half-export report")
        rows = [scenario(core, ZERO, stable=True)]
        rows.extend(scenario(core, f) for f in sweep)
        out = {
            "experiment_id": "GED-V1-R3",
            "classification": "ESTABLISHED_PHYSICS_THERMAL_LEDGER_SCOUT",
            "observation_performed": False,
            "binding_energy_model": "Newtonian self-gravitational binding U, not chemical/nuclear material binding",
            "body_mass_kg": str(core["constituent_mass_kg"]),
            "U_iron_core_J": str(core["U_A_J"]),
            "U_ice_core_J": str(core["U_B_J"]),
            "available_potential_release_J": str(core["release_J"]),
            "total_thermal_capacity_J_K": str(core["heat_capacity_J_K"]),
            "initial_temperature_K": str(core["initial_temperature_K"]),
            "environment_temperature_K": str(core["environment_temperature_K"]),
            "gravitational_probe_radius_m": str(core["probe_radius_m"]),
            "baseline_gravitational_acceleration_m_s2": str(core["baseline_g_m_s2"]),
            "scenarios": rows,
            "limitations": [
                "A stable, bound system needs no continuous entropy/energy exhaust in this model.",
                "The B->A layer swap is imposed, not dynamically simulated or demonstrated to occur spontaneously.",
                "Potential energy is only Newtonian self-gravity; chemical and nuclear binding are not computed.",
                "No mechanical work, fracture, pressure, elastic strain, latent heat or thermal radiation dynamics.",
                "Heat capacities are constant illustrative values; no temperature-dependent material properties.",
                "The external bath is at constant temperature and fully thermalises exported energy.",
                "Reported entropy is an idealised THERMAL entropy ledger only; reconfiguration/mixing entropy is omitted.",
                "A direct escaping photon entropy must not be identified with Q/T_bath before thermalisation.",
                "Far-field mass/acceleration are weak-field energy bookkeeping, not a complete Einstein solution.",
                "No new gravity coupling, no information field signal, no observations, no measurement sensitivity claim.",
            ],
        }
        return out, rows


def json_write(path, obj):
    """WHAT: Write canonical readable data without invalid numeric tokens."""
    Path(path).write_text(
        json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def git_provenance(directory):
    """WHAT: Collect Git identity when available. WHY: version pinning."""
    def git(*args):
        try:
            return subprocess.run(
                ["git", *args], cwd=directory, text=True,
                capture_output=True, check=True, timeout=5,
            ).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return None
    root = git("rev-parse", "--show-toplevel")
    return {
        "root": root,
        "commit": git("rev-parse", "HEAD") if root else None,
        "exact_tag": git("describe", "--tags", "--exact-match") if root else None,
        "dirty": bool(git("status", "--porcelain")) if root else None,
    }


def run(project, config_path):
    """WHAT: Execute frozen-source tests, calculate, and preserve SHA-256 receipts.

    WHY: An official science run should record exactly which code and inputs
    produced every result. A failed test blocks a passing interpretation.
    """
    config_path = Path(config_path).resolve()
    cfg = load_config(config_path)
    dependency_path = project / "ged_v1_r2.py"
    observed_dependency_hash = sha256(dependency_path)
    pinned_dependency_hash = cfg["r2_dependency_sha256"]
    if observed_dependency_hash != pinned_dependency_hash:
        raise RuntimeError(
            "GED-V1-R2 dependency changed: expected " + pinned_dependency_hash
            + "; found " + observed_dependency_hash
        )
    # WHAT: Validate all science calculations before opening a new run folder.
    result, rows = calculate(cfg)
    start = datetime.now(timezone.utc)
    run_id = "GED_V1_R3_" + start.strftime("%Y%m%dT%H%M%S%fZ") + "_" + uuid.uuid4().hex[:8]
    folder = project / "runs" / run_id
    folder.mkdir(parents=True, exist_ok=False)
    snapshot = folder / "input_snapshot"
    snapshot.mkdir()
    filenames = [
        "ged_v1_r2.py", "ged_v1_r3.py", "test_ged_v1_r3.py",
        "GED_V1_R3_README.md", "GED_V1_R3_DEV_AUDIT.md",
        "GED_V1_R3_PRECISION_FAILURE_A.log",
    ]
    for filename in filenames:
        shutil.copy2(project / filename, snapshot / filename)
    shutil.copy2(config_path, snapshot / "GED_V1_R3_config.json")

    # WHAT: Re-run tests using frozen input copies and actual active Python.
    # WHY: Testing an editable working copy instead would not prove this run.
    test_command = [sys.executable, "-B", "-m", "unittest", "-v", "test_ged_v1_r3.py"]
    test_run = subprocess.run(test_command, cwd=snapshot, text=True, capture_output=True)
    (folder / "tests.txt").write_text(test_run.stdout + test_run.stderr, encoding="utf-8")

    if test_run.returncode == 0:
        json_write(folder / "results.json", result)
        csv_buffer = io.StringIO(newline="")
        writer = csv.DictWriter(csv_buffer, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        (folder / "entropy_energy_ledger.csv").write_text(csv_buffer.getvalue(), encoding="utf-8")
        stable, closed, *rest = rows
        open_case = rows[-1]
        half = next(r for r in rows if dec(r["escaped_fraction"]) == Decimal("0.5") and r["kind"] != "stable_reference")
        report = "\n".join([
            "GED-V1-R3 — ENTROPIC EXHAUST / STRUCTURAL PERSISTENCE",
            "="*69,
            "NUMERICAL THERMODYNAMIC CONTROL — NOT A NEW GRAVITY DETECTION",
            f"Potential release (B -> A): {dec(result['available_potential_release_J']):.12E} J",
            f"Stable case: thermal entropy ΔS={stable['modeled_thermal_entropy_change_J_K']} J/K; exported energy=0 J",
            f"Isolated case: exported energy=0 J; ΔS_body={dec(closed['body_entropy_change_J_K']):.12E} J/K",
            f"50% export: ΔS_body={dec(half['body_entropy_change_J_K']):.12E} J/K; ΔS_bath={dec(half['bath_entropy_change_J_K']):.12E} J/K",
            f"100% export: body ΔM={dec(open_case['body_mass_change_kg']):.12E} kg",
            f"100% export: body Δg at probe={dec(open_case['exterior_delta_g_m_s2']):.12E} m/s^2",
            "Every scenario: energy balance=0; modeled thermal entropy >=0.",
            "Interpretation: entropy can increase without energy leaving the body.",
            "Entropic exhaust is not a required cost of a stable bound state.",
            "No new information-carrying gravity channel derived or detected.",
            "Tests: PASS",
        ]) + "\n"
        (folder / "report.txt").write_text(report, encoding="utf-8")
    else:
        report = "GED-V1-R3 — TEST FAILURE: science results blocked\n"
        (folder / "report.txt").write_text(report, encoding="utf-8")

    outputs = sorted(p.name for p in folder.iterdir() if p.is_file())
    provenance = {
        "experiment_id": "GED-V1-R3", "run_id": run_id,
        "utc_start": start.isoformat(),
        "status": "NUMERICAL_BASELINE_TESTED_NEW_HYPOTHESIS_NOT_TESTED" if test_run.returncode == 0 else "FAILED_TESTS_QUARANTINED",
        "invocation": sys.argv, "test_command": test_command,
        "test_returncode": test_run.returncode,
        "dependency_r2_sha256_expected": pinned_dependency_hash,
        "dependency_r2_sha256_observed": observed_dependency_hash,
        "python_version": sys.version, "python_executable": sys.executable,
        "platform": platform.platform(), "machine": platform.machine(),
        "dependencies": "Python standard library and pinned GED-V1-R2 source",
        "git": git_provenance(project),
        "input_SHA256": {p.name: sha256(p) for p in sorted(snapshot.iterdir())},
        "output_SHA256": {name: sha256(folder / name) for name in outputs},
        "limitations": result["limitations"],
    }
    json_write(folder / "provenance.json", provenance)
    all_files = sorted(
        p for p in folder.rglob("*")
        if p.is_file() and p.name != "SHA256SUMS"
    )
    (folder / "SHA256SUMS").write_text(
        "".join(f"{sha256(p)}  {p.relative_to(folder).as_posix()}\n" for p in all_files),
        encoding="utf-8",
    )
    print(report, end="")
    print("Results and receipts:", folder)
    if test_run.returncode:
        raise RuntimeError("GED-V1-R3: tests failed; run remains quarantined")
    return folder


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="GED_V1_R3_config.json")
    args = parser.parse_args()
    project = Path(__file__).resolve().parent
    run(project, Path(args.config))


if __name__ == "__main__":
    main()
