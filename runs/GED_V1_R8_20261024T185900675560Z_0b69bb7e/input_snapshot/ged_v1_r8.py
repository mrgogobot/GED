#!/usr/bin/env python3
"""GED-V1-R8: Newtonian gravitational information recoverability.

WHAT: Compare gravitational measurements from candidate internal energy patterns.
WHY:  Test whether external gravity can identify source structure, and establish
      non-uniqueness and noise limitations before proposing any new physics.

This is a conventional, synthetic, restricted inverse problem—not an observation.
Python standard library only. Units appear in variable and configuration names.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import subprocess
import sys
from datetime import datetime, timezone
from typing import Any
from uuid import uuid4

EXPERIMENT_ID = "GED-V1-R8"
CONFIG_NAME = "GED_V1_R8_config.json"
# WHAT: Fixed source menu. WHY: We can test identifiability only after specifying
#       what the hypothetical observer already knows about possible sources.
CLASSES = ("pair_x", "pair_y", "cube")
AXES = ((1., 0., 0.), (-1., 0., 0.), (0., 1., 0.),
        (0., -1., 0.), (0., 0., 1.), (0., 0., -1.))
INPUT_NAMES = ("ged_v1_r8.py", "test_ged_v1_r8.py", CONFIG_NAME,
               "GED_V1_R8_README.md", "GED_V1_R8_DEV_AUDIT.md")


def validate(config: dict[str, Any]) -> None:
    """WHAT: Reject ill-posed numerical inputs. WHY: Guard the run evidence."""
    if config.get("experiment_id") != EXPERIMENT_ID:
        raise ValueError("wrong experiment identifier")
    if tuple(config.get("candidate_classes", ())) != CLASSES:
        raise ValueError("candidate classes must be the frozen source menu")
    # WHAT: Python bool is also an int. WHY: Reject true/false masquerading as
    #       physically meaningful numeric inputs.
    numeric = ("gravitational_constant_m3_kg_s2", "speed_of_light_m_s",
               "added_energy_j", "nominal_body_radius_m",
               "parcel_offset_radius_m", "shell_radius_a_m",
               "shell_radius_b_m", "inside_shell_probe_radius_m")
    for key in numeric:
        value = config.get(key)
        if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value <= 0:
            raise ValueError(f"invalid positive value: {key}")
    outer = config["nominal_body_radius_m"]
    if not (config["parcel_offset_radius_m"] < outer):
        raise ValueError("parcels must be inside the body")
    sa, sb = config["shell_radius_a_m"], config["shell_radius_b_m"]
    si = config["inside_shell_probe_radius_m"]
    if not (0 < sa < si < sb < outer):
        raise ValueError("inside shell probe must separate two shells")
    radii = config.get("exterior_probe_radii_m")
    if not isinstance(radii, list) or not radii or any(
        isinstance(r, bool) or not isinstance(r, (int, float)) or
        not math.isfinite(r) or r <= outer for r in radii
    ) or any(a >= b for a, b in zip(radii, radii[1:])):
        raise ValueError("exterior radii must be strictly increasing beyond body")
    levels = config.get("gaussian_noise_std_m_s2")
    if not isinstance(levels, list) or not levels or levels[0] != 0 or any(
        isinstance(s, bool) or not isinstance(s, (int, float)) or
        not math.isfinite(s) or s < 0 for s in levels
    ) or any(a >= b for a, b in zip(levels, levels[1:])):
        raise ValueError("noise levels must begin at zero and increase")
    trials = config.get("trials_per_class_per_noise_level")
    if isinstance(trials, bool) or not isinstance(trials, int) or not 1 <= trials <= 10000:
        raise ValueError("invalid number of trials")
    seed = config.get("random_seed")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("random seed must be a nonnegative integer")
    if config.get("claims") != {
        "new_gravity_detected": False, "observational_dataset": False,
        "source_menu_assumed_known": True, "baseline_subtracted_perfectly": True
    }:
        raise ValueError("required scientific disclosure is missing")


def vec_length(vector: tuple[float, float, float]) -> float:
    """WHAT: Euclidean norm. WHY: Newton's gravitational force uses distance."""
    return math.sqrt(sum(v * v for v in vector))


def probe_locations(config: dict[str, Any]) -> tuple[tuple[float, float, float], ...]:
    """WHAT: Generate 24 exterior probes. WHY: Multiple viewing angles matter."""
    return tuple(tuple(radius * unit for unit in axis) for radius in config["exterior_probe_radii_m"] for axis in AXES)


def parcel_positions(kind: str, offset_m: float) -> tuple[tuple[float, float, float], ...]:
    """WHAT: Build equal-weight point-energy source arrangements.

    WHY:  We vary only the spatial placement, holding total added energy fixed.
    Pair_x and pair_y have a quadrupole; the cube's quadrupole vanishes, but
    higher-order moments remain. Points represent energies, not real iron lumps.
    """
    if kind == "pair_x":
        return ((offset_m, 0., 0.), (-offset_m, 0., 0.))
    if kind == "pair_y":
        return ((0., offset_m, 0.), (0., -offset_m, 0.))
    if kind == "cube":
        a = offset_m / math.sqrt(3.)
        return tuple((sx * a, sy * a, sz * a) for sx in (-1., 1.) for sy in (-1., 1.) for sz in (-1., 1.))
    raise ValueError(f"unrecognised source arrangement: {kind}")


def point_mass_field(point: tuple[float, float, float],
                     source: tuple[float, float, float],
                     mass_kg: float, G: float) -> tuple[float, float, float]:
    """WHAT: Exact Newtonian acceleration vector from one point mass.

    WHY:  The measurement model needs the vector, not merely its magnitude.
    Equation: g(x) = -G*m*(x-source)/|x-source|^3 [m/s^2].
    """
    sep = tuple(p - s for p, s in zip(point, source))
    r = vec_length(sep)
    if r == 0:
        raise ValueError("probe coincides with a point source")
    prefactor = -G * mass_kg / (r ** 3)
    return tuple(prefactor * value for value in sep)


def pattern_field(point: tuple[float, float, float], kind: str,
                  config: dict[str, Any]) -> tuple[float, float, float]:
    """WHAT: Sum the exact field of each source energy parcel.

    WHY:  Gravity is linear in Newtonian mass density. Here energy supplies
    equivalent mass E/c^2 and each parcel gets an equal share.
    """
    if vec_length(point) <= config["nominal_body_radius_m"]:
        raise ValueError("only exterior pattern probes are supported")
    sources = parcel_positions(kind, config["parcel_offset_radius_m"])
    parcel_mass = config["added_energy_j"] / config["speed_of_light_m_s"]**2 / len(sources)
    g = [0., 0., 0.]
    for source in sources:
        part = point_mass_field(point, source, parcel_mass, config["gravitational_constant_m3_kg_s2"])
        for i in range(3):
            g[i] += part[i]
    return tuple(g)


def spherical_shell_field(point: tuple[float, float, float], radius_m: float,
                          mass_kg: float, G: float) -> tuple[float, float, float]:
    """WHAT: Newton's exact theorem for an ideal, uniform thin shell.

    WHY:  Two different shell radii with the same mass have identical exterior
    fields at *every* exterior location, even with perfect observations.
    Inside: zero field. Outside: same as a central point mass.
    """
    r = vec_length(point)
    if r == radius_m:
        raise ValueError("ideal thin shell has a discontinuity at its surface")
    if r < radius_m:
        return (0., 0., 0.)
    return point_mass_field(point, (0., 0., 0.), mass_kg, G)


def flatten_vectors(vectors: tuple[tuple[float, float, float], ...]) -> tuple[float, ...]:
    """WHAT: Convert measured x,y,z accelerations to one ordered feature vector."""
    return tuple(component for vector in vectors for component in vector)


def signatures(config: dict[str, Any]) -> dict[str, tuple[float, ...]]:
    """WHAT: Produce the known three-class library used by the inverse scout."""
    points = probe_locations(config)
    return {kind: flatten_vectors(tuple(pattern_field(p, kind, config) for p in points)) for kind in CLASSES}


def squared_distance(first: tuple[float, ...], second: tuple[float, ...]) -> float:
    """WHAT: Sum squared measurement residuals. WHY: Gaussian-noise fit."""
    if len(first) != len(second):
        raise ValueError("measurement feature counts differ")
    return math.fsum((a - b)**2 for a, b in zip(first, second))


def classify(measured: tuple[float, ...], library: dict[str, tuple[float, ...]]) -> tuple[str, dict[str, float]]:
    """WHAT: Pick the smallest least-squares residual among 3 known templates.

    WHY:  This tests *restricted classification*, not full source reconstruction.
    Known: class menu, total energy, source offset, measurement locations, and
    an exactly subtracted common background. NO free physical parameters fit.
    """
    residuals = {kind: squared_distance(measured, library[kind]) for kind in CLASSES}
    return min(CLASSES, key=lambda kind: residuals[kind]), residuals


def add_noise(clean: tuple[float, ...], sigma: float, rng: random.Random) -> tuple[float, ...]:
    """WHAT: Add independent Gaussian noise to each Cartesian component.

    WHY:  Quantify practical recoverability against numerical perturbations.
    This is *synthetic noise*, not calibrated noise from a real gravimeter.
    """
    return tuple(value + rng.gauss(0., sigma) for value in clean)


def monte_carlo(config: dict[str, Any], library: dict[str, tuple[float, ...]]) -> list[dict[str, Any]]:
    """WHAT: Repeated, seeded classification at increasing synthetic noise.

    WHY:  A single lucky draw is insufficient to assess practical decoding.
    Every truth class gets exactly N trials at each noise setting.
    """
    rng = random.Random(config["random_seed"])
    trials = config["trials_per_class_per_noise_level"]
    result = []
    for sigma in config["gaussian_noise_std_m_s2"]:
        counts = {true: {pred: 0 for pred in CLASSES} for true in CLASSES}
        for true in CLASSES:
            for _ in range(trials):
                measured = add_noise(library[true], sigma, rng)
                predicted, _ = classify(measured, library)
                counts[true][predicted] += 1
        correct = sum(counts[kind][kind] for kind in CLASSES)
        result.append({"noise_std_m_s2": sigma,
                       "correct": correct, "total": trials * len(CLASSES),
                       "accuracy": correct / (trials * len(CLASSES)),
                       "confusion_matrix": counts})
    return result


def shell_nonuniqueness(config: dict[str, Any]) -> dict[str, Any]:
    """WHAT: Compare two different uniform spherical shells.

    WHY:  This is an EXACT null-space counterexample, unlike noisy failure:
    even complete perfect measurements outside the body cannot tell them apart.
    An inside measurement distinguishes them, but penetrates our observation
    boundary and changes the available information.
    """
    mass = config["added_energy_j"] / config["speed_of_light_m_s"]**2
    points = probe_locations(config)
    s1 = flatten_vectors(tuple(spherical_shell_field(p, config["shell_radius_a_m"], mass, config["gravitational_constant_m3_kg_s2"]) for p in points))
    s2 = flatten_vectors(tuple(spherical_shell_field(p, config["shell_radius_b_m"], mass, config["gravitational_constant_m3_kg_s2"]) for p in points))
    inside = (config["inside_shell_probe_radius_m"], 0., 0.)
    i1 = spherical_shell_field(inside, config["shell_radius_a_m"], mass, config["gravitational_constant_m3_kg_s2"])
    i2 = spherical_shell_field(inside, config["shell_radius_b_m"], mass, config["gravitational_constant_m3_kg_s2"])
    return {"shell_radii_m": [config["shell_radius_a_m"], config["shell_radius_b_m"]],
            "external_signatures_equal": s1 == s2,
            "external_max_absolute_difference_m_s2": max(abs(a - b) for a, b in zip(s1, s2)),
            "internal_probe_radius_m": inside[0],
            "inside_field_a_x_m_s2": i1[0], "inside_field_b_x_m_s2": i2[0]}


def compute(config: dict[str, Any]) -> dict[str, Any]:
    """WHAT: Assemble all measured and inverse-problem calculations."""
    validate(config)
    library = signatures(config)
    noiseless = {k: classify(library[k], library)[0] for k in CLASSES}
    pair_separation = {
        f"{a}_vs_{b}": math.sqrt(squared_distance(library[a], library[b]))
        for i, a in enumerate(CLASSES) for b in CLASSES[i+1:]
    }
    return {
        "experiment_id": EXPERIMENT_ID,
        "physics": "Newtonian synthetic gravity; no additional coupling",
        "observed": False,
        "recoverability_scope": "restricted known-template classification; known added energy and exact baseline subtraction",
        "added_energy_j": config["added_energy_j"],
        "equivalent_mass_kg": config["added_energy_j"] / config["speed_of_light_m_s"]**2,
        "probe_count": len(probe_locations(config)),
        "components_per_probe": 3,
        "noiseless_classifications": noiseless,
        "template_separation_l2_m_s2": pair_separation,
        "noise_sweep": monte_carlo(config, library),
        "shell_null_control": shell_nonuniqueness(config),
        "limits": [
            "Exact uniform shell theorem nonuniqueness; classification restricted to three preselected patterns",
            "Perfectly known/subtracted background and source energy; no unknown source amplitudes or positions",
            "Independent Gaussian synthetic noise only; no real gravimeter calibration, systematic errors, or laboratory data",
            "Point energy parcels are mathematical approximations; no stress equilibrium or source transport",
            "Newtonian weak-field treatment; no complete stress-energy tensor, Einstein equations, or binding microphysics",
            "No new gravitational information channel derived, measured, or confirmed"
        ]
    }


def sha256(path: Path) -> str:
    """WHAT: Hash file bytes. WHY: Pin exact evidence content across transfers."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_state(directory: Path) -> dict[str, Any]:
    """WHAT: Capture existing Git state without creating commits or publishing."""
    if not shutil_which_git():
        return {"available": False}
    def cmd(*args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(["git", "-C", str(directory), *args], text=True, capture_output=True, check=False)
    root = cmd("rev-parse", "--show-toplevel")
    if root.returncode != 0:
        return {"available": True, "repository": False}
    commit = cmd("rev-parse", "HEAD")
    status = cmd("status", "--porcelain", "--untracked-files=no")
    return {"available": True, "repository": True, "root": root.stdout.strip(),
            "head_commit": commit.stdout.strip() if commit.returncode == 0 else None,
            "tracked_files_dirty": bool(status.stdout.strip()) if status.returncode == 0 else None,
            "new_experiment_untracked_expected": True,
            "git_writes_performed": False}


def shutil_which_git() -> bool:
    """WHAT: Probe PATH only; no dependency on an additional Python package."""
    import shutil
    return shutil.which("git") is not None


def write_report(result: dict[str, Any]) -> str:
    """WHAT: Short human-readable summary without overstating evidence."""
    lines = ["GED-V1-R8 — GRAVITATIONAL INFORMATION RECOVERABILITY",
             "=" * 67,
             "SYNTHETIC NEWTONIAN INVERSE SCOUT; NOT AN OBSERVATION OR NEW PHYSICS",
             f"Energy per arrangement: {result['added_energy_j']:.7g} J",
             f"Equivalent mass: {result['equivalent_mass_kg']:.12e} kg",
             f"Exterior probes: {result['probe_count']} (three components each)",
             "Noiseless classification among pair_x, pair_y, cube:"]
    for name, guess in result["noiseless_classifications"].items():
        lines.append(f"  {name} -> {guess}")
    lines.append("\nKNOWN-TEMPLATE SYNTHETIC NOISE SWEEP:")
    for sweep in result["noise_sweep"]:
        lines.append(f"  sigma={sweep['noise_std_m_s2']:.3e} m/s^2 : {sweep['correct']}/{sweep['total']} = {100*sweep['accuracy']:.2f}%")
    shell = result["shell_null_control"]
    lines += ["\nEXACT EXTERIOR NONUNIQUENESS (UNIFORM SPHERICAL SHELLS):",
              f"  Shell radii {shell['shell_radii_m']} m, equal total added mass",
              f"  All exterior field vectors equal: {shell['external_signatures_equal']}",
              f"  Max exterior component difference: {shell['external_max_absolute_difference_m_s2']:.3e} m/s^2",
              f"  At inside r={shell['internal_probe_radius_m']} m: shell A gx={shell['inside_field_a_x_m_s2']:.12e}, shell B gx={shell['inside_field_b_x_m_s2']:.12e} m/s^2",
              "\nVERDICT: Restricted synthetic classification only. Exterior gravity cannot",
              "uniquely identify a general internal energy distribution; added noise",
              "reduces practical recoverability. No new gravity channel derived."]
    return "\n".join(lines) + "\n"


def snapshot_and_write_run(project: Path, config: dict[str, Any], result: dict[str, Any]) -> Path:
    """WHAT: Freeze sources, results, tests, environment, command and hashes.

    WHY:  Manifest of Destiny: each run preserves reproducible inputs and
    immutable receipts. The run directory is new; older GED runs stay untouched.
    """
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    run_dir = project / "runs" / f"GED_V1_R8_{stamp}_{uuid4().hex[:8]}"
    run_dir.mkdir(parents=True, exist_ok=False)
    inputs = run_dir / "input_snapshot"
    inputs.mkdir()
    import shutil
    for filename in INPUT_NAMES:
        shutil.copy2(project / filename, inputs / filename)
    tests = subprocess.run([sys.executable, "-m", "unittest", "-v", "test_ged_v1_r8.py"],
                           cwd=project, text=True, capture_output=True, check=False)
    (run_dir / "tests.txt").write_text(tests.stdout + tests.stderr +
        f"\nRETURN_CODE={tests.returncode}\n", encoding="utf-8")
    (run_dir / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    (run_dir / "report.txt").write_text(write_report(result), encoding="utf-8")
    library = signatures(config)
    locations = probe_locations(config)
    with (run_dir / "exterior_signatures.csv").open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["probe_index", "x_m", "y_m", "z_m"] +
                        [f"{kind}_{axis}_m_s2" for kind in CLASSES for axis in "xyz"])
        for idx, point in enumerate(locations):
            writer.writerow([idx, *(format(value, ".17g") for value in point),
                             *(format(library[kind][3*idx+j], ".17e") for kind in CLASSES for j in range(3))])
    provenance = {
        "experiment_id": EXPERIMENT_ID,
        "start_utc_from_run_dir": stamp,
        "python_executable": sys.executable,
        "python_version": sys.version,
        "platform": platform.platform(),
        "implementation": platform.python_implementation(),
        "dependencies": {"third_party": [], "stdlib": True},
        "config_file": CONFIG_NAME,
        "input_sha256": {name: sha256(inputs / name) for name in INPUT_NAMES},
        "run_command": "python ged_v1_r8.py",
        "verification_command": "python -m unittest -v test_ged_v1_r8.py",
        "test_return_code": tests.returncode,
        "test_count_expected": 71,
        "random_seed": config["random_seed"],
        "git": git_state(project),
        "publication_status": "LOCAL_ONLY_NO_GITHUB_NO_ZENODO_UNTIL_V1_R10",
        "scientific_scope": result["recoverability_scope"],
        "limitation": "Synthetic Newtonian inverse control; not observational evidence"
    }
    (run_dir / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    # WHAT: Hash every run file except the manifest itself. WHY: A manifest
    #       cannot include its own hash without a self-reference problem.
    paths = sorted(p for p in run_dir.rglob("*") if p.is_file() and p.name != "SHA256SUMS")
    manifest = "".join(f"{sha256(p)}  {p.relative_to(run_dir).as_posix()}\n" for p in paths)
    (run_dir / "SHA256SUMS").write_text(manifest, encoding="utf-8")
    if tests.returncode != 0:
        raise RuntimeError(f"Run retained at {run_dir}, but test suite failed: {tests.returncode}")
    return run_dir


def main() -> None:
    project = Path(__file__).resolve().parent
    config = json.loads((project / CONFIG_NAME).read_text(encoding="utf-8"))
    result = compute(config)
    print(write_report(result), end="")
    run_dir = snapshot_and_write_run(project, config, result)
    print(f"Results and provenance: {run_dir}")
    print(f"Run-file checksums: {run_dir / 'SHA256SUMS'}")


if __name__ == "__main__":
    main()
