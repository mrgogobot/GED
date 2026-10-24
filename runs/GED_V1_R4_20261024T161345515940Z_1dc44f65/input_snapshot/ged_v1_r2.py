#!/usr/bin/env python3
"""GED-V1-R2: Newtonian gravitational binding + energy/mass accounting.

Established-physics NUMERICAL SCOUT, not an observation, new law, or complete GR.
All units SI. No third-party packages; developed for Python 3.13.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext, InvalidOperation
import hashlib
import io
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import uuid

# WHAT: High-precision decimal arithmetic and pi.
# WHY: A gravitational binding mass correction may be ~10^-10 kg compared
# with an asteroid mass ~10^10 kg. Binary float subtraction would lose it.
PRECISION = 60
PI = Decimal("3.14159265358979323846264338327950288419716939937510582097494")


def D(value):
    """WHAT: Parse SI quantities. WHY: avoid surprise binary float rounding."""
    try:
        number = Decimal(str(value))
    except (InvalidOperation, ValueError):
        raise ValueError(f"invalid number: {value!r}") from None
    if not number.is_finite():
        raise ValueError("nonfinite numbers rejected")
    return number


def sphere_volume(radius):
    """WHAT: V=4*pi*R^3/3 in m^3. WHY: convert density to mass."""
    r = D(radius)
    if r < 0:
        raise ValueError("negative radius")
    return 4 * PI * r**3 / 3


def make_body(radius, fraction, core_density, shell_density):
    """WHAT: Build a two-layer ideal sphere with fixed constituent inventory.

    WHY: A and B swap material positions, NOT total material quantities.
    IMPORTANT: Densities here represent constituent mass per volume BEFORE
    adding gravitational self-energy. A real measured density may already
    include energy contributions; our bookkeeping intentionally separates them.
    """
    R, f, rc, rs = map(D, (radius, fraction, core_density, shell_density))
    if R <= 0 or rc <= 0 or rs <= 0 or not 0 < f < 1:
        raise ValueError("invalid geometry or density")
    with localcontext() as ctx:
        ctx.prec = PRECISION
        # Cube root: because the core occupies fraction f of total volume.
        a = R * f ** (D(1) / 3)
        mc = rc * sphere_volume(R) * f
        ms = rs * sphere_volume(R) * (1 - f)
        return dict(R=R, a=+a, rc=rc, rs=rs, mc=+mc, ms=+ms, M0=+(mc + ms))


def binding_energy(body, G):
    """WHAT: Newtonian gravitational self-potential U<0, in joules.

    WHY: The gravitational binding energy depends on internal density profile.
    Derived from U=-integral_0^R G*M_inside(r)/r dM(r) for two uniform zones.
    NOT chemical or nuclear binding energy, nor a full GR correction.
    """
    g = D(G)
    if g <= 0:
        raise ValueError("G must be positive")
    with localcontext() as ctx:
        ctx.prec = PRECISION
        a, R, rc, rs = (body[k] for k in ("a", "R", "rc", "rs"))
        pieces = (rc**2 * a**5 / 5
                  + rs**2 * (R**5 - a**5) / 5
                  + rs * (rc-rs) * a**3 * (R**2 - a**2) / 2)
        return +(-16 * PI**2 * g * pieces / 3)


def homogeneous_binding(M, radius, G):
    """Independent uniform-sphere law U=-3GM^2/(5R); useful audit check."""
    m, r, g = map(D, (M, radius, G))
    if min(m, r, g) <= 0:
        raise ValueError("M,R,G must be positive")
    return -3 * g * m**2 / (5 * r)


def mass_equivalent(delta_E, c):
    """WHAT: delta_M=delta_E/c^2. WHY: energy gravitates as mass-energy.

    Signed energy refers to energy inside the chosen BODY boundary. For
    emission, negative delta_E. Quasistatic symmetric transfer assumed;
    no net momentum or recoil. Do NOT discard energy outside the boundary.
    """
    e, speed = D(delta_E), D(c)
    if speed <= 0:
        raise ValueError("speed of light must be positive")
    return e / speed**2


def exterior_g(mass, r, G):
    """WHAT: Newtonian exterior acceleration G*M/r^2 (m/s^2)."""
    m, radius, g = map(D, (mass, r, G))
    if min(m, radius, g) <= 0:
        raise ValueError("invalid exterior gravity inputs")
    return g * m / radius**2


def delta_g(delta_E, r, G, c):
    """WHAT: Small external delta-g for fixed r: G*delta_E/(c^2*r^2).

    WHY: Directly calculate difference; subtraction of large float field
    values destroys tiny corrections. Weak-field illustration, NOT exact GR.
    """
    radius, g = D(r), D(G)
    if min(radius, g) <= 0:
        raise ValueError("invalid radius or G")
    return g * mass_equivalent(delta_E, c) / radius**2


def calculate(config):
    """WHAT: Pure calculation. WHY: allow isolated independent tests.

    A = iron core / ice shell; B = ice core / iron shell.
    Scenarios contrast energy leaving a body with heat remaining inside it.
    No GED-specific information coupling is inserted into these controls.
    """
    with localcontext() as ctx:
        ctx.prec = PRECISION
        G, c, R, f = [D(config[k]) for k in
                       ("G_m3_kg_s2", "c_m_s", "radius_m", "core_volume_fraction")]
        iron, ice = D(config["density_iron_kg_m3"]), D(config["density_ice_kg_m3"])
        heat_kg = D(config["thermal_input_J_kg"])
        probes = [D(v) for v in config["exterior_probes_R"]]
        probe_r = D(config["report_probe_R"])
        if min(G, c, R, iron, ice, heat_kg) <= 0 or not 0 < f < 1:
            raise ValueError("invalid positive input or volume fraction")
        if not probes or min(probes) < 1 or probe_r < 1:
            raise ValueError("all exterior probes must be at least 1R")
        A, B = make_body(R, f, iron, ice), make_body(R, f, ice, iron)
        if A["M0"] != B["M0"]:
            raise AssertionError("not equal constituent inventories")
        ua, ub = binding_energy(A, G), binding_energy(B, G)
        du = ua - ub     # Ice-core B -> iron-core A changes potential energy.
        if du >= 0:
            raise AssertionError("iron core should have more negative self-energy")
        Q = heat_kg * A["M0"]
        # WHAT: Four explicitly bounded energy transactions.
        # WHY: binding energy alone is NOT a command to alter gravity.
        # The total energy in the defined body matters.
        cases = {
            "absorbed_heat_kept": +Q,
            "same_heat_radiated": -Q,
            "B_to_A_binding_release_escapes": du,
            "B_to_A_binding_release_retained_as_heat": D(0),
        }
        rows = []
        for name, de in cases.items():
            for ratio in probes:
                r = ratio*R
                g0 = exterior_g(A["M0"], r, G)
                dg = delta_g(de, r, G, c)
                rows.append(dict(scenario=name, distance_over_R=str(ratio),
                                 distance_m=str(r), delta_E_J=str(de),
                                 delta_M_kg=str(mass_equivalent(de, c)),
                                 baseline_g_m_s2=str(g0), delta_g_m_s2=str(dg),
                                 fractional_delta_g=str(dg/g0)))
        at_probe = {}
        for name, de in cases.items():
            r = probe_r * R
            g0 = exterior_g(A["M0"], r, G)
            at_probe[name] = dict(delta_E_J=str(de),
                                  delta_M_kg=str(mass_equivalent(de, c)),
                                  delta_g_m_s2=str(delta_g(de, r, G, c)),
                                  fractional_delta_g=str(delta_g(de, r, G, c)/g0))
        result = {
            "experiment_id": "GED-V1-R2", "observation_performed": False,
            "classification": "ESTABLISHED_PHYSICS_WEAK_FIELD_SCOUT",
            "source_v1_r1_run": config["source_v1_r1_run"],
            "constituent_inventory_mass_kg": str(A["M0"]),
            "radius_m": str(R), "core_radius_m": str(A["a"]),
            "U_iron_core_J": str(ua), "U_ice_core_J": str(ub),
            "delta_U_B_to_A_J": str(du),
            "binding_energy_release_B_to_A_J": str(-du),
            "mass_change_if_release_escapes_kg": str(mass_equivalent(du, c)),
            "mass_change_if_release_retained_kg": "0",
            "at_distance_R": str(probe_r),
            "baseline_exterior_g_m_s2": str(exterior_g(A["M0"], probe_r*R, G)),
            "scenarios": at_probe,
            "model_notes": [
                "Only Newtonian gravitational self-binding U is computed; not nuclear, chemical or fracture energy.",
                "Binding energy is already part of ordinary mass-energy accounting; no GED extra field term.",
                "In closed rearrangement, reduced potential energy is offset by increased heat.",
                "In open rearrangement, symmetric outward radiation exports binding energy.",
                "Delta-g is fixed-radius point-mass weak-field bookkeeping, not a full GR solution.",
                "Constituent mass inventory is distinguished from already-bound measured mass.",
                "High-precision arithmetic is not evidence of measurable gravitational sensitivity.",
                "No dynamic rearrangement, external mechanical work, material failure, recoil or field propagation modeled.",
                "No observation; cannot confirm or disprove an additional information field.",
            ],
        }
        return result, rows

# WHAT: The following procedures record a reproducible official experiment.
# WHY: Source, configuration, runtime, tests and every output need receipts.

def digest(path):
    """Return SHA-256 digest of an exact on-disk file's bytes."""
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for piece in iter(lambda: stream.read(1024*1024), b""):
            h.update(piece)
    return h.hexdigest()


def write_json(path, obj):
    """Write readable JSON with stable key order and no NaN/Infinity."""
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False)+"\n",
                    encoding="utf-8")


def git_version(directory):
    """Read available Git identity; never initialize/change a repository."""
    def call(*args):
        try:
            return subprocess.run(["git", *args], cwd=directory, text=True,
                                  check=True, capture_output=True, timeout=5).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return None
    root = call("rev-parse", "--show-toplevel")
    return dict(root=root, commit=call("rev-parse", "HEAD") if root else None,
                exact_tag=call("describe", "--tags", "--exact-match") if root else None,
                dirty=bool(call("status", "--porcelain")) if root else None)


def run(project, config_path):
    """Snapshot exact inputs; test; then emit results with SHA receipts.

    WHY: A failed test must be visible and must block a passing science run.
    This standalone runner is not claiming the separate Turnip Run wrapper ran.
    """
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if config.get("experiment_id") != "GED-V1-R2":
        raise ValueError("wrong experiment ID")
    # Validate inputs before allocating a run directory.
    calculate(config)
    start = datetime.now(timezone.utc)
    run_id = "GED_V1_R2_"+start.strftime("%Y%m%dT%H%M%S%fZ")+"_"+uuid.uuid4().hex[:8]
    folder = project / "runs" / run_id
    folder.mkdir(parents=True, exist_ok=False)
    snapshot = folder / "input_snapshot"
    snapshot.mkdir()
    input_names = ["ged_v1_r2.py", "test_ged_v1_r2.py",
                   "GED_V1_R2_README.md", "GED_V1_R2_DEV_AUDIT.md"]
    for name in input_names:
        shutil.copy2(project / name, snapshot / name)
    shutil.copy2(config_path, snapshot / "GED_V1_R2_config.json")

    # Tests run against FROZEN source, not against a separately edited copy.
    command = [sys.executable, "-B", "-m", "unittest", "-v", "test_ged_v1_r2.py"]
    tested = subprocess.run(command, cwd=snapshot, text=True, capture_output=True)
    (folder / "tests.txt").write_text(tested.stdout+tested.stderr, encoding="utf-8")
    if tested.returncode != 0:
        print("UNACCEPTED — failed tests recorded at", folder / "tests.txt")
        raise RuntimeError("tests failed, no results emitted")

    results, rows = calculate(config)
    write_json(folder / "results.json", results)
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    (folder / "energy_field_samples.csv").write_text(stream.getvalue(), encoding="utf-8")
    report = "\n".join([
        "GED-V1-R2 — GRAVITATIONAL BINDING / MASS-ENERGY ACCOUNTING",
        "="*64,
        "ESTABLISHED-PHYSICS NUMERICAL MODEL — NOT A NEW GRAVITY DETECTION",
        f"Each body constituent mass: {D(results['constituent_inventory_mass_kg']):.12E} kg",
        f"Iron core U: {D(results['U_iron_core_J']):.12E} J",
        f"Ice core U:  {D(results['U_ice_core_J']):.12E} J",
        f"B -> A binding energy release: {D(results['binding_energy_release_B_to_A_J']):.12E} J",
        f"If energy ESCAPES, delta mass: {D(results['mass_change_if_release_escapes_kg']):.12E} kg",
        "If energy STAYS as heat, delta mass: 0 kg",
        f"Probe radius: {results['at_distance_R']} R",
        f"Reference exterior g: {D(results['baseline_exterior_g_m_s2']):.12E} m/s^2",
        "Verdict: numerical baseline only; no additional information channel tested.",
        "Tests: PASS",
    ])+"\n"
    (folder / "report.txt").write_text(report, encoding="utf-8")
    output_names = ["tests.txt", "results.json", "energy_field_samples.csv", "report.txt"]
    provenance = {
        "run_id": run_id, "utc_start": start.isoformat(),
        "outcome": "NUMERICAL_BASELINE_VERIFIED_NEW_HYPOTHESIS_UNTESTED",
        "command_line": sys.argv, "test_command": command,
        "test_return_code": tested.returncode,
        "python_version": sys.version, "python_executable": sys.executable,
        "platform": platform.platform(), "architecture": platform.machine(),
        "dependencies": "Python standard library only", "git": git_version(project),
        "input_snapshot_SHA256": {p.name: digest(p) for p in sorted(snapshot.iterdir())},
        "output_SHA256": {name: digest(folder/name) for name in output_names},
        "limitations": results["model_notes"],
    }
    write_json(folder / "provenance.json", provenance)
    all_files = sorted(p for p in folder.rglob("*") if p.is_file() and p.name != "SHA256SUMS")
    (folder / "SHA256SUMS").write_text(
        "".join(f"{digest(p)}  {p.relative_to(folder).as_posix()}\n" for p in all_files),
        encoding="utf-8")
    print(report, end="")
    print("Results and receipts:", folder)
    return folder


def main():
    """Run one configured model and preserve its complete receipts."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="GED_V1_R2_config.json")
    args = parser.parse_args()
    project = Path(__file__).resolve().parent
    run(project, Path(args.config).resolve())


if __name__ == "__main__":
    main()
