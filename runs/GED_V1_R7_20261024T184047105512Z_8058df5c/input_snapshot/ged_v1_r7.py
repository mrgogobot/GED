#!/usr/bin/env python3
"""GED-V1-R7: asymmetric equal-energy sources, exact Newtonian field and quadrupole.

READ FIRST:
    This script is a conventional, explicitly idealized gravitational calculation.
    It does not establish a gravity-mediated material-information mechanism.

PYTHON: 3.13 stdlib; no third-party dependencies.

WHAT IS MODELED:
    The *same* 100-m spherical iron/ice asteroid background is present in each
    experiment. On top of that common background, equal total energy is placed
    into either two opposing parcels or eight cubically placed parcels.

WHY:
    Earlier GED scouts found equal exterior fields for spherical distributions.
    Breaking spherical symmetry tests whether ordinary multipoles already
    reveal structural differences outside the object.

COMMENTING POLICY:
    Each computational section explains WHAT it calculates and WHY it matters.
    Assumptions are repeated beside any potentially overinterpreted quantity.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import platform
import subprocess
import sys
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from pathlib import Path


# WHAT: Standard decimal representation of SI quantities and a high-precision
# arithmetic context. WHY: m_energy = E/c^2 is ~1.1e-9 kg, much smaller than
# the background asteroid mass. Separately evaluating corrections avoids
# losing significant digits when subtracting full external accelerations.
D = Decimal
PRECISION = 70
FILE_SET = (
    "ged_v1_r7.py",
    "test_ged_v1_r7.py",
    "GED_V1_R7_config.json",
    "GED_V1_R7_README.md",
    "GED_V1_R7_DEV_AUDIT.md",
    "GED_V1_R7_PRECISION_FAILURE_A.log",
    "GED_V1_R7_REPORT_FORMAT_FAILURE_A.log",
)
AXES = {"x": (D(1), D(0), D(0)),
        "y": (D(0), D(1), D(0)),
        "z": (D(0), D(0), D(1))}


def dec(value: object, label: str) -> Decimal:
    """WHAT: Convert text values to finite Decimal. WHY: reject invalid inputs."""
    if isinstance(value, bool):
        raise ValueError(f"{label}: booleans are not physical numbers")
    try:
        number = D(str(value))
    except Exception as exc:
        raise ValueError(f"{label}: invalid number") from exc
    if not number.is_finite():
        raise ValueError(f"{label}: must be finite")
    return number


def vector_add(left, right):
    """WHAT: Componentwise vector addition. WHY: sum gravitational contributions."""
    return tuple(left[i] + right[i] for i in range(3))


def vector_sub(left, right):
    """WHAT: Componentwise subtraction. WHY: measure displacement to a source."""
    return tuple(left[i] - right[i] for i in range(3))


def dot(left, right):
    """WHAT: Euclidean scalar product. WHY: distances and multipoles."""
    return sum((left[i] * right[i] for i in range(3)), D(0))


def norm(vector):
    """WHAT: Length in metres or acceleration units. WHY: inverse-square law."""
    with localcontext() as ctx:
        ctx.prec = PRECISION
        return dot(vector, vector).sqrt()


def multiply(vector, scalar):
    """WHAT: Scale a vector. WHY: Newtonian acceleration along displacement."""
    return tuple(scalar * item for item in vector)


@dataclass(frozen=True)
class Parcel:
    """WHAT: One localized energy parcel. WHY: retain position and energy separately.

    IMPORTANT: A point parcel is a mathematical approximation. It does not mean
    physical energy exists at infinitely small size or that stresses vanish.
    """
    position_m: tuple[Decimal, Decimal, Decimal]
    energy_J: Decimal


@dataclass(frozen=True)
class Model:
    """WHAT: Fully validated frozen experiment inputs. WHY: prevent silent drift."""
    G: Decimal
    c: Decimal
    R: Decimal
    core_R: Decimal
    rho_core: Decimal
    rho_shell: Decimal
    E: Decimal
    a: Decimal
    radii: tuple[Decimal, ...]
    axes: tuple[str, ...]
    cases: tuple[str, ...]
    limitations: tuple[str, ...]


def model_from_config(data: dict) -> Model:
    """WHAT: Parse and validate a config. WHY: reject unphysical scout requests."""
    if data.get("experiment_id") != "GED-V1-R7":
        raise ValueError("incorrect experiment identifier")
    kwargs = {name: dec(data[field], field) for name, field in (
        ("G", "G_m3_kg_s2"), ("c", "c_m_s"),
        ("R", "outer_radius_m"), ("core_R", "core_radius_m"),
        ("rho_core", "iron_core_density_kg_m3"),
        ("rho_shell", "ice_shell_density_kg_m3"),
        ("E", "total_added_energy_J"), ("a", "parcel_offset_m"),
    )}
    if any(kwargs[k] <= 0 for k in ("G", "c", "R", "core_R", "rho_core", "rho_shell", "E", "a")):
        raise ValueError("constants, densities, dimensions, and energy must be positive")
    if not (kwargs["core_R"] < kwargs["a"] < kwargs["R"]):
        raise ValueError("parcels must sit inside shell and outside iron core")
    radii = tuple(dec(item, "sample radius") for item in data["sample_radii_m"])
    if not radii or any(r <= kwargs["R"] for r in radii):
        raise ValueError("all probes must be outside the outer radius")
    if len(set(radii)) != len(radii) or tuple(sorted(radii)) != radii:
        raise ValueError("radii must be strictly increasing")
    axes = tuple(data["sample_axes"])
    if not axes or len(set(axes)) != len(axes) or any(axis not in AXES for axis in axes):
        raise ValueError("sample axes must be distinct x/y/z entries")
    cases = tuple(data["comparison"])
    if set(cases) != {"pair_x", "pair_y", "cube_eight"} or len(cases) != 3:
        raise ValueError("comparison must include exactly the three defined cases")
    limits = tuple(data["limitations"])
    if not limits or any(not isinstance(item, str) or not item.strip() for item in limits):
        raise ValueError("scientific limitations must be present")
    return Model(**kwargs, radii=radii, axes=axes, cases=cases, limitations=limits)


def spherical_background_mass(model: Model) -> Decimal:
    """WHAT: Compute fixed iron-core + ice-shell mass [kg].

    WHY: Keep GED R6's material inventory while studying *energy placement*.
    The background is spherical and identical in every comparison, so its
    exterior Newtonian field is the same and cancels in pairwise differences.
    """
    with localcontext() as ctx:
        ctx.prec = PRECISION
        pi = D("3.1415926535897932384626433832795028841971693993751058209749445923078")
        volume_coeff = D(4) * pi / D(3)
        return volume_coeff * (model.rho_core * model.core_R**3 +
                               model.rho_shell * (model.R**3 - model.core_R**3))


def create_parcels(model: Model, case: str) -> tuple[Parcel, ...]:
    """WHAT: Construct three *equal-total-energy* configurations.

    WHY:
      pair_x: two parcels at ±a x. Nonzero quadrupole.
      pair_y: same geometry rotated 90 degrees. Rotates the quadrupole.
      cube_eight: parcels at the eight corners (±a/√3,±a/√3,±a/√3).
      The quadrupole cancels but higher multipoles remain; a cube of
      point parcels is not a continuous spherical source.

    In all cases total E is fixed and the energy centroid is at origin.
    """
    if case == "pair_x":
        positions = ((model.a, D(0), D(0)), (-model.a, D(0), D(0)))
    elif case == "pair_y":
        positions = ((D(0), model.a, D(0)), (D(0), -model.a, D(0)))
    elif case == "cube_eight":
        # WHAT: Eight equally weighted energy parcels at cube corners.
        # WHY: Exact 1/8 energy shares avoid fractional-rounding asymmetry,
        # while cubic symmetry cancels dipole and quadrupole moments.
        # Their radial distance is approximately 'a'; Decimal precision
        # limits the last digits of the square-root coordinate.
        with localcontext() as ctx:
            ctx.prec = PRECISION
            b = model.a / D(3).sqrt()
            positions = tuple((D(sx)*b, D(sy)*b, D(sz)*b)
                              for sx in (-1, 1) for sy in (-1, 1)
                              for sz in (-1, 1))
    else:
        raise ValueError(f"undefined energy-placement case: {case}")
    with localcontext() as ctx:
        ctx.prec = PRECISION
        # Allocate the same E exactly, aside from Decimal's finite precision.
        return tuple(Parcel(position, model.E / D(len(positions)))
                     for position in positions)


def energy_mass(model: Model, parcels: tuple[Parcel, ...]) -> Decimal:
    """WHAT: Additional monopole mass = ΣE/c² [kg]. WHY: standard E=mc²."""
    with localcontext() as ctx:
        ctx.prec = PRECISION
        return sum((p.energy_J for p in parcels), D(0)) / model.c**2


def dipole_kg_m(model: Model, parcels: tuple[Parcel, ...]):
    """WHAT: First mass moment Σm_i*x_i [kg m]. WHY: control centre of energy."""
    with localcontext() as ctx:
        ctx.prec = PRECISION
        result = (D(0), D(0), D(0))
        for p in parcels:
            result = vector_add(result, multiply(p.position_m, p.energy_J / model.c**2))
        return result


def quadrupole_kg_m2(model: Model, parcels: tuple[Parcel, ...]):
    """WHAT: Symmetric trace-free Q_ij = Σm(3x_i x_j - r² δ_ij).

    WHY: Q is the first shape-sensitive external moment for our zero-dipole
    sources. Units: kg m². It does *not* encode chemistry or binding rules.
    """
    with localcontext() as ctx:
        ctx.prec = PRECISION
        matrix = [[D(0) for _ in range(3)] for _ in range(3)]
        for p in parcels:
            mass = p.energy_J / model.c**2
            squared = dot(p.position_m, p.position_m)
            for i in range(3):
                for j in range(3):
                    matrix[i][j] += mass * (
                        D(3) * p.position_m[i] * p.position_m[j] -
                        (squared if i == j else D(0)))
        return tuple(tuple(row) for row in matrix)


def gravity_exact(model: Model, parcels: tuple[Parcel, ...], probe):
    """WHAT: Vector sum of exact Newtonian *energy-parcel* accelerations.

    Formula: g(r) = -G Σ (E_i/c²)(r-r_i)/|r-r_i|³  [m/s²].
    WHY: Acts as an independent target for a truncated multipole expansion.
    This is the added-energy correction, not the entire asteroid's gravity.
    """
    with localcontext() as ctx:
        ctx.prec = PRECISION
        if norm(probe) <= model.R:
            raise ValueError("probe must be outside the reference asteroid")
        g = (D(0), D(0), D(0))
        for parcel in parcels:
            displacement = vector_sub(probe, parcel.position_m)
            d = norm(displacement)
            factor = -model.G * parcel.energy_J / (model.c**2 * d**3)
            g = vector_add(g, multiply(displacement, factor))
        return g


def gravity_monopole(model: Model, parcels: tuple[Parcel, ...], probe):
    """WHAT: Point-mass field of the same *total* added energy [m/s²].

    WHY: Equal-energy sources must have identical monopole fields, whatever
    their orientations or arrangements. This omits all higher multipoles.
    """
    with localcontext() as ctx:
        ctx.prec = PRECISION
        r = norm(probe)
        if r <= model.R:
            raise ValueError("probe must be external")
        return multiply(probe, -model.G * energy_mass(model, parcels) / r**3)


def gravity_quadrupole(model: Model, parcels: tuple[Parcel, ...], probe):
    """WHAT: First non-monopole exterior-field correction from Q_ij [m/s²].

    The potential is Φ_Q = -G(r^T Q r)/(2r^5).
    Taking minus its gradient yields:
        g_Q = G[(Q r)/r^5 - (5/2)(r^T Q r)r/r^7].

    WHY: Compare this analytic approximation to the exact parcel sum.
    Valid only as a far-field expansion; a=60m, so r≫a is preferable.
    """
    with localcontext() as ctx:
        ctx.prec = PRECISION
        r = norm(probe)
        if r <= model.R:
            raise ValueError("probe must be external")
        q = quadrupole_kg_m2(model, parcels)
        qr = tuple(sum((q[i][j] * probe[j] for j in range(3)), D(0))
                   for i in range(3))
        scalar = dot(probe, qr)
        return tuple(model.G * (qr[i] / r**5 - D(5) * scalar * probe[i] /
                                (D(2) * r**7)) for i in range(3))


def gravity_through_quadrupole(model: Model, parcels, probe):
    """WHAT: Monopole plus quadrupole. WHY: test leading shape signal."""
    return vector_add(gravity_monopole(model, parcels, probe),
                      gravity_quadrupole(model, parcels, probe))


def radial_component(vector, axis):
    """WHAT: Signed projection onto a unit radial direction [m/s²]."""
    return dot(vector, AXES[axis])


def relative_error(actual: Decimal, prediction: Decimal) -> Decimal:
    """WHAT: |prediction-actual|/|actual|. WHY: avoid hidden absolute scales."""
    if actual == 0:
        raise ValueError("relative error requires nonzero exact value")
    return abs((prediction - actual) / actual)


def fmt(value: Decimal) -> str:
    """WHAT: Serializable precision. WHY: preserve tiny nonzero signals."""
    return format(value, ".24E")


def json_write(path: Path, obj):
    """WHAT: Stable JSON text output. WHY: reproducible human/machine records."""
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    """WHAT: Exact byte digest. WHY: identify code, data and results."""
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for piece in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(piece)
    return digest.hexdigest()


def git_details(folder: Path) -> dict:
    """WHAT: Record optional local Git state. WHY: receipts, never auto-commit."""
    data = {"repository_commit": None, "working_tree_dirty": None,
            "release_tag": None, "note": "Git not required or modified."}
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=folder,
                              text=True, capture_output=True, timeout=10)
        if head.returncode:
            return data
        data["repository_commit"] = head.stdout.strip()
        changes = subprocess.run(["git", "status", "--porcelain"], cwd=folder,
                                 text=True, capture_output=True, timeout=10)
        if changes.returncode == 0:
            data["working_tree_dirty"] = bool(changes.stdout.strip())
        tags = subprocess.run(["git", "tag", "--points-at", "HEAD"], cwd=folder,
                              text=True, capture_output=True, timeout=10)
        if tags.returncode == 0:
            data["release_tag"] = tags.stdout.strip().splitlines() or None
        return data
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return data


def compute(model: Model):
    """WHAT: Evaluate exact and approximate fields at declared external probes.

    WHY: isolate conventional asymmetric gravity, holding mass-energy fixed.
    We *never* subtract the huge total-body field to find tiny corrections.
    """
    with localcontext() as ctx:
        ctx.prec = PRECISION
        groups = {case: create_parcels(model, case) for case in model.cases}
        background = spherical_background_mass(model)
        energy_equiv = model.E / model.c**2
        summary = {}
        for name, parts in groups.items():
            tensor = quadrupole_kg_m2(model, parts)
            summary[name] = {
                "number_of_energy_parcels": len(parts),
                "total_energy_J": fmt(sum((p.energy_J for p in parts), D(0))),
                "energy_mass_kg": fmt(energy_mass(model, parts)),
                "dipole_kg_m": [fmt(x) for x in dipole_kg_m(model, parts)],
                "quadrupole_kg_m2": [[fmt(x) for x in row] for row in tensor],
                "quadrupole_trace_kg_m2": fmt(sum((tensor[i][i] for i in range(3)), D(0))),
            }
        rows = []
        for case in model.cases:
            parts = groups[case]
            for axis in model.axes:
                for r in model.radii:
                    probe = multiply(AXES[axis], r)
                    exact = gravity_exact(model, parts, probe)
                    point = gravity_monopole(model, parts, probe)
                    quad = gravity_quadrupole(model, parts, probe)
                    pred = vector_add(point, quad)
                    exact_radial = radial_component(exact, axis)
                    pred_radial = radial_component(pred, axis)
                    rows.append({
                        "case": case, "direction": axis, "radius_m": str(r),
                        "exact_added_radial_g_m_s2": fmt(exact_radial),
                        "monopole_added_radial_g_m_s2": fmt(radial_component(point, axis)),
                        "quadrupole_radial_g_m_s2": fmt(radial_component(quad, axis)),
                        "truncated_added_radial_g_m_s2": fmt(pred_radial),
                        "abs_radial_residual_m_s2": fmt(abs(pred_radial-exact_radial)),
                        "relative_truncation_error": fmt(relative_error(exact_radial, pred_radial)),
                        "common_background_inward_g_m_s2": fmt(model.G * background / r**2),
                    })
        at = model.radii[0]
        px = multiply(AXES["x"], at)
        g_pair = radial_component(gravity_exact(model, groups["pair_x"], px), "x")
        g_cube = radial_component(gravity_exact(model, groups["cube_eight"], px), "x")
        g_pair_Q = radial_component(gravity_quadrupole(model, groups["pair_x"], px), "x")
        g_cube_Q = radial_component(gravity_quadrupole(model, groups["cube_eight"], px), "x")
        delta_exact = g_pair - g_cube
        delta_pred = g_pair_Q - g_cube_Q
        # The ratio is conditional on the signal being nonzero; this is a
        # property of our stipulated source arrangement, not new physics.
        comparison_error = relative_error(delta_exact, delta_pred)
        return {
            "experiment_id": "GED-V1-R7",
            "classification": "CONVENTIONAL NUMERICAL MULTIPOLE SCOUT — NO NEW GRAVITY",
            "assumptions": {
                "background": "unchanged spherical iron core and ice shell",
                "parcels": "point-source representation of fixed total additional energy",
                "physics": "Newtonian inverse-square sum; E/c^2 as a weak-field source proxy",
                "higher_order": "through traceless Cartesian quadrupole only",
            },
            "outer_radius_m": str(model.R),
            "offset_m": str(model.a),
            "common_background_rest_mass_kg": fmt(background),
            "added_energy_J": str(model.E),
            "added_mass_equivalent_kg": fmt(energy_equiv),
            "cases": summary,
            "comparison_at_first_radius_on_x_axis": {
                "radius_m": str(at),
                "exact_pair_x_minus_cube_radial_g_m_s2": fmt(delta_exact),
                "quadrupole_pair_x_minus_cube_radial_g_m_s2": fmt(delta_pred),
                "quadrupole_relative_error": fmt(comparison_error),
            },
            "sample_count": len(rows),
            "scientific_verdict": "Asymmetric conventional sources can have different external Newtonian fields at matched total energy and dipole. No information-field coupling is derived or detected.",
            "limits": list(model.limitations),
        }, rows


def write_csv(path: Path, rows: list[dict]):
    """WHAT: Write one row per configuration and probe. WHY: external audits."""
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_manifest(run_folder: Path):
    """WHAT: SHA-256 every material run file except the manifest itself.

    WHY: SHA256SUMS cannot hash itself without a circular definition.
    Use relative paths so the archive can be validated anywhere.
    """
    entries = sorted(p for p in run_folder.rglob("*")
                     if p.is_file() and p.name != "SHA256SUMS")
    lines = [f"{sha256_file(p)}  {p.relative_to(run_folder).as_posix()}"
             for p in entries]
    (run_folder / "SHA256SUMS").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(project_dir: Path) -> Path:
    """WHAT: Execute frozen config and preserve source/tests/results receipts.

    WHY: Manifest of Destiny: experiments are not accepted just because a
    display says PASS. Immutable run snapshots and checksums enable audit.
    """
    with localcontext() as ctx:
        ctx.prec = PRECISION
        for name in FILE_SET:
            if not (project_dir / name).is_file():
                raise FileNotFoundError(f"required GED-V1-R7 input missing: {name}")
        data = json.loads((project_dir / "GED_V1_R7_config.json").read_text(encoding="utf-8"))
        model = model_from_config(data)
        created = datetime.now(timezone.utc)
        run_id = f"GED_V1_R7_{created.strftime('%Y%m%dT%H%M%S%fZ')}_{uuid.uuid4().hex[:8]}"
        folder = project_dir / "runs" / run_id
        folder.mkdir(parents=True, exist_ok=False)
        snapshot = folder / "input_snapshot"
        snapshot.mkdir()
        for name in FILE_SET:
            (snapshot / name).write_bytes((project_dir / name).read_bytes())
        # WHAT: Execute the preserved tests in a separate interpreter.
        # WHY: Confirm the recorded run's *own* source, not an unrelated edit.
        command = [sys.executable, "-m", "unittest", "-v", "test_ged_v1_r7.py"]
        # Do not create regenerable .pyc caches in the preserved run.
        # Provenance should track scientific inputs and evidence, not Python caches.
        test_env = dict(os.environ)
        test_env["PYTHONDONTWRITEBYTECODE"] = "1"
        check = subprocess.run(command, cwd=snapshot, capture_output=True,
                               text=True, timeout=90, env=test_env)
        (folder / "tests.txt").write_text(
            "COMMAND: " + " ".join(command) + "\nRETURN CODE: " +
            str(check.returncode) + "\n\nSTDOUT:\n" + check.stdout +
            "\nSTDERR:\n" + check.stderr, encoding="utf-8")
        if check.returncode:
            # Preserve the failure rather than producing a false PASS.
            json_write(folder / "FAILED_RUN.json", {
                "experiment_id": "GED-V1-R7",
                "reason": "frozen source test execution failed",
                "test_return_code": check.returncode,
            })
            write_manifest(folder)
            raise RuntimeError(f"Frozen-snapshot tests failed; evidence preserved: {folder}")
        results, rows = compute(model)
        json_write(folder / "results.json", results)
        write_csv(folder / "field_samples.csv", rows)
        target = results["comparison_at_first_radius_on_x_axis"]
        report = (
            "GED-V1-R7 — ASYMMETRIC ENERGY / NEWTONIAN MULTIPOLE SCOUT\n"
            + "=" * 67 + "\n"
            + "NO NEW GRAVITY OR MATERIAL-INSTRUCTION INFORMATION CHANNEL DETECTED\n"
            f"Common background mass: {results['common_background_rest_mass_kg']} kg\n"
            f"Added energy each case: {results['added_energy_J']} J\n"
            f"Added equivalent mass each case: {results['added_mass_equivalent_kg']} kg\n"
            "Monopole matches, dipoles vanish; quadrupoles differ.\n"
            "Cubic quadrupole vanishes but higher moments remain.\n"
            f"At r={target['radius_m']}m along +x, exact pair_x - cube g_r: "
            f"{target['exact_pair_x_minus_cube_radial_g_m_s2']} m/s^2\n"
            f"Quadrupole-only prediction for the difference: "
            f"{target['quadrupole_pair_x_minus_cube_radial_g_m_s2']} m/s^2\n"
            f"Relative truncation error: {target['quadrupole_relative_error']}\n"
            "Full source and all calculated samples are in the preserved files.\n"
            "No pressure/stress equilibrium, source transport, or full GR modeled.\n"
            "Scientific verdict: expected Newtonian multipole structure only.\n"
            "Tests: PASS\n"
            f"Results and receipts: {folder}\n"
        )
        (folder / "report.txt").write_text(report, encoding="utf-8")
        provenance = {
            "experiment_id": "GED-V1-R7",
            "run_id": run_id,
            "timestamp_utc": created.isoformat(),
            "project_directory": str(project_dir.resolve()),
            "snapshot_directory": "input_snapshot/",
            "python_version": sys.version,
            "python_executable": sys.executable,
            "platform": platform.platform(),
            "environment": {
                "PYTHONHASHSEED": os.environ.get("PYTHONHASHSEED"),
                "PYTHONPATH": os.environ.get("PYTHONPATH"),
                "frozen_test_PYTHONDONTWRITEBYTECODE": "1",
                "decimal_precision_digits": PRECISION,
                "third_party_dependencies": [],
            },
            "working_directory": os.getcwd(),
            "entry_command": " ".join(sys.argv),
            "frozen_test_command": " ".join(command),
            "frozen_test_return_code": check.returncode,
            "git": git_details(project_dir),
            "publication_status": "local_only; no GitHub or Zenodo upload before R10 without separate consent",
            "input_sha256": {name: sha256_file(snapshot / name) for name in FILE_SET},
            "output_sha256_before_manifest": {
                name: sha256_file(folder / name) for name in
                ("results.json", "field_samples.csv", "report.txt", "tests.txt")},
            "scientific_status": "numerical control, not a new gravitational finding",
        }
        json_write(folder / "provenance.json", provenance)
        write_manifest(folder)
        print(report, end="")
        print(f"Run-file checksum manifest: {folder / 'SHA256SUMS'}")
        return folder


if __name__ == "__main__":
    try:
        run(Path(__file__).resolve().parent)
    except (ValueError, OSError, RuntimeError, KeyError, subprocess.TimeoutExpired) as error:
        print(f"GED-V1-R7 ERROR: {error}", file=sys.stderr)
        sys.exit(1)
