#!/usr/bin/env python3
"""GED-V1-R1: conventional spherical null and deliberately hypothetical toy law.

PURPOSE: Compare spherical iron-core/ice-shell with ice-core/iron-shell bodies
that have equal radius and total Newtonian mass. The external fields are
identical under established spherical gravity; the internal fields and thermal
behaviours differ. Then INVENT a free material-dependent external coupling to
show precisely what a nonstandard measurable prediction could look like.

IMPORTANT: The toy term is not a derived, validated, or observed physical law.
Dependencies: Python 3.13 standard library; units: SI throughout.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import uuid


# WHAT: Store one material's bulk density and thermal characteristics.
# WHY: This prevents mixing up values and units when comparing iron and ice.
# LIMIT: Heat of fusion is not the same as total atomic/nuclear binding energy.
@dataclass(frozen=True)
class Material:
    name: str
    density: float          # kg/m^3
    melt_temp: float        # kelvin (K)
    heat_capacity: float    # J/(kg*K)
    fusion_energy: float    # J/kg


# WHAT: Store a spherical core plus a shell, each with uniform density.
# WHY: Swapping their locations isolates structure from total composition.
@dataclass(frozen=True)
class LayeredSphere:
    radius: float                 # outer radius, m
    core_fraction: float          # fraction of total VOLUME occupied by core
    core: Material
    shell: Material

    @property
    def core_radius(self) -> float:
        # Radius scales as the cube root of volume fraction.
        return self.radius * self.core_fraction ** (1 / 3)

    @property
    def mass(self) -> float:
        # Total Newtonian mass is the sum of density*volume for each layer.
        v_core = volume(self.core_radius)
        return self.core.density * v_core + self.shell.density * (volume(self.radius) - v_core)


def volume(radius: float) -> float:
    # WHAT: Volume of a sphere; WHY: Needed to convert densities to masses.
    if not math.isfinite(radius) or radius < 0:
        raise ValueError("radius must be finite and nonnegative")
    return 4 * math.pi * radius**3 / 3


def enclosed_mass(body: LayeredSphere, distance: float) -> float:
    # WHAT: Calculate the mass enclosed inside distance r from the centre.
    # WHY: Newton's shell theorem says outer spherical shells cancel in the
    # interior ACCELERATION (not necessarily in gravitational potential).
    if not math.isfinite(distance) or distance < 0:
        raise ValueError("distance must be finite and nonnegative")
    r = min(distance, body.radius)
    core_v = volume(min(r, body.core_radius))
    return body.core.density * core_v + body.shell.density * (volume(r) - core_v)


def gravity(body: LayeredSphere, distance: float, G: float) -> float:
    # WHAT: Newtonian acceleration magnitude in m/s^2.
    # WHY: This is the accepted conventional null comparator (alpha=0).
    if not math.isfinite(G) or G <= 0 or not math.isfinite(distance) or distance < 0:
        raise ValueError("invalid gravity inputs")
    if distance == 0:
        return 0.0  # Centre has zero acceleration from perfect spherical symmetry.
    return G * enclosed_mass(body, distance) / distance**2


def point_mass_gravity(mass: float, distance: float, G: float) -> float:
    # WHAT: Independent point-mass calculation of the external field.
    # WHY: Cross-check the layered-sphere implementation at r >= R.
    if not all(math.isfinite(x) and x > 0 for x in (mass, distance, G)):
        raise ValueError("point mass requires positive finite inputs")
    return G * mass / distance**2


def melt_threshold(material: Material, T0: float) -> float:
    # WHAT: Heat per kg to reach melting point AND melt entirely, in J/kg.
    # WHY: One measurable material-response proxy for the hypothesis scout.
    # NOT atomic binding energy! Ice may sublimate in vacuum instead of melt.
    if not math.isfinite(T0) or T0 < 0 or T0 >= material.melt_temp:
        raise ValueError("initial temperature must be below melting")
    return (material.melt_temp - T0) * material.heat_capacity + material.fusion_energy


def melt_fraction(material: Material, T0: float, heat_per_kg: float) -> float:
    # WHAT: Simplified fraction of the material melted by a specified heat input.
    # WHY: The same thermal stress produces a different response in different
    # materials. Assumes constant heat capacity, no loss, proper pressure,
    # no fracture, no chemistry, no sublimation, no change in gravitational mass.
    if not math.isfinite(heat_per_kg) or heat_per_kg < 0:
        raise ValueError("invalid heat input")
    heat_to_start_melt = melt_threshold(material, T0) - material.fusion_energy
    heat_for_melting = max(0.0, heat_per_kg - heat_to_start_melt)
    return min(1.0, heat_for_melting / material.fusion_energy)


def shell_contrast(shell: Material, iron: Material, ice: Material, T0: float) -> float:
    # WHAT: Normalize thermal-response thresholds so iron=+1 and ice=-1.
    # WHY: The arbitrary toy gravity formula needs a DIMENSIONLESS input.
    # This chosen material proxy is NOT a physical information measurement.
    i = melt_threshold(iron, T0)
    w = melt_threshold(ice, T0)
    if math.isclose(i, w):
        raise ValueError("materials have indistinguishable melting thresholds")
    return (2 * melt_threshold(shell, T0) - (i + w)) / (i - w)


def toy_gravity(body: LayeredSphere, r: float, G: float,
                alpha: float, contrast: float, p: int) -> float:
    # WHAT: Introduce an INVENTED, externally measurable correction:
    #        g_TOY = (G*M/r^2) * [1 + alpha * contrast * (R/r)^p]
    # WHY: To distinguish a zero-extra-information null from a specific
    # example of how a new material-dependent coupling might behave.
    # WARNING: Not derived from GR; conservation, covariance and viability
    # have NOT been established. A nonzero simulated output proves nothing.
    if r < body.radius:
        raise ValueError("toy formula defined only outside sphere")
    if not math.isfinite(alpha) or not math.isfinite(contrast):
        raise ValueError("toy parameters must be finite")
    if type(p) is not int or p < 0:
        raise ValueError("radial exponent must be a nonnegative integer")
    correction = 1 + alpha * contrast * (body.radius / r)**p
    if correction < 0:
        raise ValueError("negative magnitude outside toy model domain")
    return point_mass_gravity(body.mass, r, G) * correction


def alpha_sensitivity(uncertainty_ppm: float, r_over_R: float, p: int) -> float:
    # WHAT: Translate an assumed measurement uncertainty into a toy-alpha
    # scale, given shell contrasts of exactly +1 and -1.
    # WHY: A null experiment can constrain alpha within THIS arbitrary model.
    # It is not a real instrument error budget or a formal confidence limit.
    if not math.isfinite(uncertainty_ppm) or uncertainty_ppm <= 0 or r_over_R < 1 or type(p) is not int or p < 0:
        raise ValueError("invalid toy sensitivity inputs")
    return uncertainty_ppm * 1e-6 * r_over_R**p / 2


def load_config(path: Path) -> dict:
    # WHAT: Read declared model parameters instead of hardcoding experiments.
    # WHY: Every repeat should use an auditable input configuration.
    with path.open(encoding="utf-8") as f:
        c = json.load(f)
    if c.get("experiment_id") != "GED-V1-R1":
        raise ValueError("wrong experiment identifier")
    if c["radius_m"] <= 0 or c["G_m3_kg_s2"] <= 0:
        raise ValueError("negative radius or G")
    if not math.isclose(c["core_volume_fraction"], 0.5, rel_tol=0, abs_tol=1e-14):
        raise ValueError("this swap control requires a 50% core volume")
    if not 0 < c["interior_probe_R"] < 1 or any(r < 1 for r in c["exterior_probes_R"]):
        raise ValueError("invalid probe radii")
    if c["toy"]["probe_R"] < 1 or c["toy"]["difference_uncertainty_ppm"] <= 0:
        raise ValueError("invalid toy probe or uncertainty")
    if type(c["toy"]["radial_power"]) is not int or c["toy"]["radial_power"] < 0:
        raise ValueError("invalid radial power")
    if c["temperature_K"] < 0 or c["heat_input_J_kg"] < 0:
        raise ValueError("temperature or energy cannot be negative")
    for material in c["materials"].values():
        if not all(math.isfinite(x) and x > 0 for x in material.values()):
            raise ValueError("all material parameters must be finite and positive")
        if c["temperature_K"] >= material["melt_K"]:
            raise ValueError("initial temperature must be below melting")
    return c


def material(name: str, c: dict) -> Material:
    # WHAT: Convert user-controlled JSON material constants into typed fields.
    # WHY: Keep mass and heat formulas free of JSON field-name clutter.
    m = c["materials"][name]
    return Material(name, m["density_kg_m3"], m["melt_K"],
                    m["heat_capacity_J_kg_K"], m["latent_fusion_J_kg"])


def calculate(c: dict) -> tuple[dict, list[dict]]:
    # WHAT: Perform the numerical experiment without reading or writing files.
    # WHY: Pure calculations can be audited and tested independent of receipts.
    iron, ice = material("iron", c), material("ice", c)
    R, f, G = c["radius_m"], c["core_volume_fraction"], c["G_m3_kg_s2"]
    A = LayeredSphere(R, f, iron, ice)  # Same material quantities as B.
    B = LayeredSphere(R, f, ice, iron)  # Only radial ordering changes.
    T0, q = c["temperature_K"], c["heat_input_J_kg"]
    t = c["toy"]
    chi_A = shell_contrast(A.shell, iron, ice, T0)
    chi_B = shell_contrast(B.shell, iron, ice, T0)
    if not math.isclose(A.mass, B.mass, rel_tol=1e-14):
        raise AssertionError("matched-control masses differ")

    # WHAT: Sample both models at five external radii for independent plotting.
    # WHY: Do not rely on one handpicked point to characterize a radial field.
    rows = []
    for r_over_R in c["exterior_probes_R"]:
        r = r_over_R * R
        g_A, g_B = gravity(A, r, G), gravity(B, r, G)
        if not math.isclose(g_A, g_B, rel_tol=1e-13):
            raise AssertionError("external Newtonian field mismatch")
        toy_A = toy_gravity(A, r, G, t["alpha"], chi_A, t["radial_power"])
        toy_B = toy_gravity(B, r, G, t["alpha"], chi_B, t["radial_power"])
        rows.append({
            "radius_in_R": r_over_R,
            "distance_m": r,
            "baseline_A_m_s2": g_A,
            "baseline_B_m_s2": g_B,
            "toy_A_m_s2": toy_A,
            "toy_B_m_s2": toy_B,
            "toy_difference_m_s2": toy_B - toy_A,
            "toy_difference_ppm_baseline": 1e6 * (toy_B - toy_A) / g_A,
        })
    probe_r = t["probe_R"] * R
    g_ref = gravity(A, probe_r, G)
    diff = (toy_gravity(B, probe_r, G, t["alpha"], chi_B, t["radial_power"])
            - toy_gravity(A, probe_r, G, t["alpha"], chi_A, t["radial_power"]))
    results = {
        "experiment_id": "GED-V1-R1",
        "verdict": "CONVENTIONAL_NULL_PLUS_HYPOTHETICAL_TOY_SCOUT",
        "observation_performed": False,
        "model": {
            "baseline": "Newtonian spherical g=G*M_enclosed(r)/r^2",
            "hypothetical": "g_TOY=(G*M/r^2)*(1+alpha*chi_shell*(R/r)^p), r>=R",
            "toy_is_derived_or_validated": False,
        },
        "geometry_and_mass": {
            "radius_m": R,
            "core_volume_fraction": f,
            "core_radius_m": A.core_radius,
            "mass_A_kg": A.mass,
            "mass_B_kg": B.mass,
            "A_core_shell": "iron/ice",
            "B_core_shell": "ice/iron",
        },
        "baseline": {
            "interior_probe_radius_m": c["interior_probe_R"] * R,
            "interior_g_A_m_s2": gravity(A, c["interior_probe_R"] * R, G),
            "interior_g_B_m_s2": gravity(B, c["interior_probe_R"] * R, G),
            "surface_g_A_m_s2": gravity(A, R, G),
            "surface_g_B_m_s2": gravity(B, R, G),
            "exterior_fields_equal": all(math.isclose(row["baseline_A_m_s2"], row["baseline_B_m_s2"], rel_tol=1e-13) for row in rows),
            "shell_melt_fraction_A": melt_fraction(A.shell, T0, q),
            "shell_melt_fraction_B": melt_fraction(B.shell, T0, q),
            "shell_heat_to_full_melt_A_J_kg": melt_threshold(A.shell, T0),
            "shell_heat_to_full_melt_B_J_kg": melt_threshold(B.shell, T0),
        },
        "hypothetical_gate": {
            "alpha_dimensionless": t["alpha"],
            "p_dimensionless": t["radial_power"],
            "chi_A": chi_A,
            "chi_B": chi_B,
            "probe_r_over_R": t["probe_R"],
            "probe_g_baseline_m_s2": g_ref,
            "toy_difference_m_s2": diff,
            "toy_difference_ppm_baseline": 1e6 * diff / g_ref,
            "hypothetical_measurement_difference_uncertainty_ppm": t["difference_uncertainty_ppm"],
            "alpha_sensitivity_at_assumed_uncertainty": alpha_sensitivity(
                t["difference_uncertainty_ppm"], t["probe_R"], t["radial_power"]),
        },
        "limitations": [
            "No gravitational measurement has been performed.",
            "Toy coupling is stipulated, not derived; it may be inconsistent with conserved/covariant theories.",
            "Heat-to-melt threshold is not nuclear, atomic, or fracture binding energy.",
            "Bulk densities fixed; no E/c^2 corrections, phase-change mass redistribution, or full GR.",
            "Static spherical idealization excludes multipoles, rotation, environment and tides.",
            "Nominal measurement uncertainty is illustrative; no actual device error budget.",
        ],
    }
    return results, rows


def sha256(path: Path) -> str:
    # WHAT: SHA-256 of byte-for-byte files, read in small chunks.
    # WHY: Preserve identity of input programs, config, and observed outputs.
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, data: dict) -> None:
    # WHAT: Stable human-readable JSON; WHY: Make comparison straightforward.
    path.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def git_info(project: Path) -> dict:
    # WHAT: Read (never modify) optional Git commit/tag/dirty status.
    # WHY: One authoritative code version must be uniquely identified.
    def git(*args: str) -> str | None:
        try:
            return subprocess.run(["git", *args], cwd=project, check=True,
                                  text=True, capture_output=True, timeout=5).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return None
    root = git("rev-parse", "--show-toplevel")
    return {"repository_root": root,
            "commit": git("rev-parse", "HEAD") if root else None,
            "exact_tag": git("describe", "--tags", "--exact-match") if root else None,
            "dirty": bool(git("status", "--porcelain")) if root else None}


def execute(project: Path, config_path: Path) -> Path:
    # WHAT: Run tests, calculate the model, and write an immutable per-run copy.
    # WHY: The Manifest of Destiny requires verifiable inputs/outputs and both
    # accepted and failed test records, not just screen output.
    c = load_config(config_path)
    timestamp = datetime.now(timezone.utc)
    run_id = "GED_V1_R1_" + timestamp.strftime("%Y%m%dT%H%M%S%fZ") + "_" + uuid.uuid4().hex[:8]
    out = project / "runs" / run_id
    out.mkdir(parents=True, exist_ok=False)
    snapshot = out / "input_snapshot"
    snapshot.mkdir()

    # Freeze the exact inputs before execution, guarding against later edits.
    input_names = ["ged_v1_r1.py", "test_ged_v1_r1.py", "GED_V1_R1_README.md", "GED_V1_R1_DEV_AUDIT.md"]
    for name in input_names:
        shutil.copy2(project / name, snapshot / name)
    shutil.copy2(config_path, snapshot / "GED_V1_R1_config.json")

    # Run the test suite against the frozen input snapshot, not another copy.
    test_cmd = [sys.executable, "-B", "-m", "unittest", "-v", "test_ged_v1_r1.py"]
    test = subprocess.run(test_cmd, cwd=snapshot, text=True, capture_output=True, check=False)
    (out / "tests.txt").write_text(test.stdout + test.stderr, encoding="utf-8")
    if test.returncode != 0:
        print("TEST FAILURE — saved unaccepted tests at", out / "tests.txt")
        raise RuntimeError("GED-V1-R1 verification failed; no science results emitted")

    results, rows = calculate(c)
    write_json(out / "results.json", results)
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    (out / "exterior_samples.csv").write_text(buffer.getvalue(), encoding="utf-8")

    s, b, t = results["geometry_and_mass"], results["baseline"], results["hypothetical_gate"]
    report = "\n".join([
        "GED-V1-R1 — SPHERICAL NULL / TOY COUPLING SCOUT",
        "=" * 64,
        f"Radius: {s['radius_m']:.9g} m; mass of each: {s['mass_A_kg']:.12g} kg",
        f"Matched conventional surface g: {b['surface_g_A_m_s2']:.12g} m/s^2",
        f"Interior g at {b['interior_probe_radius_m']} m: A {b['interior_g_A_m_s2']:.12g}; B {b['interior_g_B_m_s2']:.12g} m/s^2",
        f"Ice shell melted (A): {100*b['shell_melt_fraction_A']:.6f}%",
        f"Iron shell melted (B): {100*b['shell_melt_fraction_B']:.6f}%",
        "CONVENTIONAL: external gravitational fields identical.",
        "",
        "HYPOTHETICAL TOY ONLY — NOT DETECTED OR PHYSICALLY DERIVED:",
        f"alpha={t['alpha_dimensionless']}; p={t['p_dimensionless']}",
        f"At {t['probe_r_over_R']} R: toy delta={t['toy_difference_ppm_baseline']:.9g} ppm of baseline",
        f"At {t['probe_r_over_R']} R: toy delta={t['toy_difference_m_s2']:.12g} m/s^2",
        f"Illustrative alpha sensitivity: {t['alpha_sensitivity_at_assumed_uncertainty']:.9g}",
        "",
        "VERDICT: Conventional baseline verified numerically. No evidence of new gravity.",
        "Heat threshold is NOT full binding energy. No laboratory observation.",
        "Tests: PASS",
    ]) + "\n"
    (out / "report.txt").write_text(report, encoding="utf-8")

    # Record executable, toolchain and all file hashes. Provenance hashes the
    # result files, then SHA256SUMS hashes provenance itself (avoids cycles).
    outputs = ["tests.txt", "results.json", "exterior_samples.csv", "report.txt"]
    provenance = {
        "run_id": run_id, "utc_start": timestamp.isoformat(),
        "outcome": "NUMERICAL_BASELINE_VERIFIED_TOY_UNTESTED",
        "command_line": sys.argv, "test_command": test_cmd, "test_return_code": test.returncode,
        "python_version": sys.version, "python_executable": sys.executable,
        "platform": platform.platform(), "machine": platform.machine(),
        "dependencies": "Python standard library only", "git": git_info(project),
        "source_config_snapshot_SHA256": {
            p.name: sha256(p) for p in snapshot.iterdir() if p.is_file()
        },
        "output_SHA256": {name: sha256(out / name) for name in outputs},
        "limitations": results["limitations"],
    }
    write_json(out / "provenance.json", provenance)
    all_files = sorted(p for p in out.rglob("*") if p.is_file() and p.name != "SHA256SUMS")
    (out / "SHA256SUMS").write_text(
        "".join(f"{sha256(p)}  {p.relative_to(out).as_posix()}\n" for p in all_files),
        encoding="utf-8",
    )
    print(report, end="")
    print("Results and receipts:", out)
    return out


def main() -> None:
    # WHAT: Choose one explicit input config. WHY: Keep runs reproducible.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="GED_V1_R1_config.json")
    args = parser.parse_args()
    project = Path(__file__).resolve().parent
    execute(project, Path(args.config).resolve())


if __name__ == "__main__":
    main()
