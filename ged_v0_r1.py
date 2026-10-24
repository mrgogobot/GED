#!/usr/bin/env python3
"""GED-V0-R1: reproducible gravitational and material-response baseline.

WHAT THIS EXPERIMENT TESTS
--------------------------
1. Equal-radius iron and ice spheres: different densities -> different mass/gravity.
2. Spherically symmetric layered objects: equal mass/radius -> identical EXTERNAL
   Newtonian field despite different material placements, while INTERNAL fields differ.
3. An ice sphere divided into non-overlapping pieces: total mass conserved;
   each piece has its own field; the spatial distribution matters nearby.
4. Simple uniform-heating response: materials have different energy thresholds.

WHAT IT DOES NOT TEST
---------------------
It does not prove that gravity stores chemical, material, or quantum instructions.
The word "information" here is a proposed interpretation to be investigated,
not an additional force or a result of this calculation.

Dependencies: Python 3.13 standard library only. No numpy/pip installs required.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import uuid
from typing import Any


# WHAT: A named container for properties used by the simplified heating model.
# WHY: We must distinguish ordinary material response from gravity. Keeping the
# material properties together helps prevent mixing units or materials.
@dataclass(frozen=True)
class Material:
    name: str
    density_kg_m3: float
    melting_point_K: float
    specific_heat_J_kg_K: float
    latent_heat_fusion_J_kg: float


# WHAT: A basic model of an approximately spherical, uniform-density body.
# WHY: At its surface and outside it, the Newtonian spherical-field equations
# are simple and provide an auditable reference calculation.
@dataclass(frozen=True)
class Sphere:
    radius_m: float
    material: Material

    @property
    def volume_m3(self) -> float:
        # Volume = (4/3)*pi*R^3, in cubic metres.
        return (4.0 / 3.0) * math.pi * self.radius_m**3

    @property
    def mass_kg(self) -> float:
        # Density [kg/m^3] times volume [m^3] gives mass [kg].
        return self.material.density_kg_m3 * self.volume_m3


# WHAT: Represent a possible fragmentation result, not the collision itself.
# WHY: This prevents an imposed change of geometry from being confused with a
# physics simulation that predicts the collision or breaking mechanism.
@dataclass(frozen=True)
class Fragments:
    count: int
    mass_each_kg: float
    radius_each_m: float
    centres_m: tuple[tuple[float, float, float], ...]

    @property
    def total_mass_kg(self) -> float:
        return self.count * self.mass_each_kg


# WHAT: Convert the project's explicit JSON configuration into numbers.
# WHY: Researchers must be able to audit constants, inputs, and units without
# modifying executable source code. We intentionally reject impossible inputs.
def load_config(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        cfg = json.load(handle)

    if cfg.get("experiment_id") != "GED-V0-R1":
        raise ValueError("Configuration experiment_id must be GED-V0-R1")
    if cfg["sphere_radius_m"] <= 0:
        raise ValueError("sphere_radius_m must be positive")
    if cfg["constants"]["G_m3_kg-1_s-2"] <= 0:
        raise ValueError("G must be positive")
    if cfg["constants"]["c_m_s-1"] <= 0:
        raise ValueError("c must be positive")
    if not 0 < cfg["layered_bodies"]["iron_volume_fraction"] < 1:
        raise ValueError("Layered iron volume fraction must be between 0 and 1")
    count_axis = cfg["fragmentation"]["fragments_per_axis"]
    if not isinstance(count_axis, int) or isinstance(count_axis, bool) or count_axis < 2:
        raise ValueError("fragments_per_axis must be an integer >= 2")
    if cfg["fragmentation"]["surface_gap_fraction_of_diameter"] < 0:
        raise ValueError("The fragment separation gap must be nonnegative")
    if not 0 < cfg["layered_bodies"]["internal_probe_radius_fraction"] < 1:
        raise ValueError("Internal probe must lie strictly inside the layered sphere")
    if not cfg["layered_bodies"]["external_probe_radius_fraction"] > 1:
        raise ValueError("External probe must lie outside the layered sphere")
    if cfg["heating"]["deposited_energy_J_kg-1"] < 0 or cfg["heating"]["initial_temperature_K"] < 0:
        raise ValueError("Heating inputs must be nonnegative")
    for name, spec in cfg["materials"].items():
        if any(float(spec[key]) <= 0 for key in (
            "density_kg_m-3", "melting_point_K", "specific_heat_J_kg-1_K-1",
            "latent_heat_of_fusion_J_kg-1"
        )):
            raise ValueError(f"Nonpositive material property for {name}")
    if any(cfg["heating"]["initial_temperature_K"] >= spec["melting_point_K"]
           for spec in cfg["materials"].values()):
        raise ValueError("Simplified melt model needs solid materials initially below melting")
    return cfg


def read_material(cfg: dict[str, Any], name: str) -> Material:
    # Translate explicitly unit-labelled JSON values into the typed container.
    spec = cfg["materials"][name]
    return Material(
        name=name,
        density_kg_m3=float(spec["density_kg_m-3"]),
        melting_point_K=float(spec["melting_point_K"]),
        specific_heat_J_kg_K=float(spec["specific_heat_J_kg-1_K-1"]),
        latent_heat_fusion_J_kg=float(spec["latent_heat_of_fusion_J_kg-1"]),
    )


def surface_gravity_m_s2(G: float, mass_kg: float, radius_m: float) -> float:
    # WHAT: g = G*M/R^2, in metres per second squared.
    # WHY: A direct gravitational observable for a spherical source.
    if mass_kg <= 0 or radius_m <= 0:
        raise ValueError("mass and radius must be positive")
    return G * mass_kg / (radius_m**2)


def enclosed_mass_kg(
    radial_distance_m: float,
    outer_radius_m: float,
    core_radius_m: float,
    core_density_kg_m3: float,
    shell_density_kg_m3: float,
) -> float:
    # WHAT: Calculate mass contained INSIDE a chosen radial distance.
    # WHY: By Newton's shell theorem, spherical matter OUTSIDE that radius
    # contributes no NET Newtonian gravitational force at that location.
    # This is not the same as the gravitational potential, which can include
    # contributions from outer shells.
    if radial_distance_m < 0:
        raise ValueError("radius cannot be negative")
    radius_inside_body = min(radial_distance_m, outer_radius_m)
    core_part_radius = min(radius_inside_body, core_radius_m)
    core_vol = (4.0 / 3.0) * math.pi * core_part_radius**3
    total_enclosed_vol = (4.0 / 3.0) * math.pi * radius_inside_body**3
    shell_vol = total_enclosed_vol - core_vol
    return core_density_kg_m3 * core_vol + shell_density_kg_m3 * shell_vol


def layered_acceleration_m_s2(
    G: float,
    radial_distance_m: float,
    outer_radius_m: float,
    core_radius_m: float,
    core_density_kg_m3: float,
    shell_density_kg_m3: float,
) -> float:
    # WHAT: Magnitude of radial acceleration at a distance from the centre.
    # WHY: Evaluate both an interior and an exterior probe using one model.
    # A spherical body behaves like a point mass OUTSIDE it, but its interior
    # field depends on the actual distribution of enclosed mass.
    if radial_distance_m == 0:
        return 0.0  # Exact spherical symmetry gives zero field at the centre.
    mass_inside = enclosed_mass_kg(
        radial_distance_m, outer_radius_m, core_radius_m,
        core_density_kg_m3, shell_density_kg_m3,
    )
    return G * mass_inside / radial_distance_m**2


def melt_fraction(material: Material, T0_K: float, energy_J_kg: float) -> float:
    # WHAT: Approximate liquid fraction after a fixed energy input per kilogram.
    # WHY: Material-dependent thermal thresholds provide an observable response
    # that is not independently encoded in external Newtonian gravity.
    #
    # ASSUMPTIONS: Uniform temperature; constant specific heat; no heat loss;
    # appropriate pressure for a solid-liquid phase transition; no chemistry,
    # sublimation, fracture, radiation or gravitational work. For ice in a
    # near-vacuum, sublimation may occur instead, so this is NOT a space-weather
    # or realistic asteroid impact simulation.
    if T0_K < 0 or energy_J_kg < 0:
        raise ValueError("Temperature and energy input must be nonnegative")
    # Energy first heats the solid to its melting temperature.
    sensible_energy_J_kg = max(0.0, material.melting_point_K - T0_K) * (
        material.specific_heat_J_kg_K
    )
    # Only remaining energy is available for the phase transition.
    leftover_J_kg = max(0.0, energy_J_kg - sensible_energy_J_kg)
    return min(1.0, leftover_J_kg / material.latent_heat_fusion_J_kg)


def make_fragment_grid(sphere: Sphere, count_per_axis: int, gap_fraction: float) -> Fragments:
    # WHAT: Divide the sphere's mass into n^3 EQUAL-MASS spherical pieces.
    # WHY: Test how changing the geometry of a mass distribution modifies the
    # gravitational field. No new material is created and none disappears.
    count = count_per_axis**3
    each_mass = sphere.mass_kg / count

    # Constant density -> radius scales with the cube root of relative mass.
    each_radius = sphere.radius_m / count_per_axis

    # Put fragment centres on a centred cubic grid, with a small gap so that
    # fragments DO NOT overlap. This is a prescribed geometry, not dynamics.
    spacing_m = 2.0 * each_radius * (1.0 + gap_fraction)
    midpoint = (count_per_axis - 1) / 2.0
    centres = tuple(
        ((i - midpoint) * spacing_m,
         (j - midpoint) * spacing_m,
         (k - midpoint) * spacing_m)
        for i in range(count_per_axis)
        for j in range(count_per_axis)
        for k in range(count_per_axis)
    )
    return Fragments(count, each_mass, each_radius, centres)


def acceleration_x_at_probe_m_s2(G: float, pieces: Fragments, probe_x_m: float) -> float:
    # WHAT: Sum the x-components of gravitational acceleration at a probe.
    # WHY: This exposes differences between a concentrated sphere and the same
    # mass distributed among separated fragments, without new physical laws.
    #
    # Newtonian vector law for each fragment:
    # a_vector = G*m*(r_fragment-r_probe) / |r_fragment-r_probe|^3
    # Units: [m^3/(kg*s^2)]*[kg]*[m]/[m^3] = m/s^2.
    # The probe must lie OUTSIDE every fragment. In that regime each spherical
    # fragment's external gravitational field equals a point mass at its centre.
    components = []
    for x, y, z in pieces.centres_m:
        dx = x - probe_x_m
        distance = math.sqrt(dx*dx + y*y + z*z)
        if distance <= pieces.radius_each_m:
            raise ValueError("Probe intersects a fragment; point-mass law invalid there")
        components.append(G * pieces.mass_each_kg * dx / distance**3)
    # math.fsum improves numerical precision when adding many small values.
    return math.fsum(components)


def calculate(cfg: dict[str, Any]) -> dict[str, Any]:
    # WHAT: Run all controlled comparisons from a single immutable config.
    # WHY: Deterministic input -> deterministic physics results for audits.
    G = float(cfg["constants"]["G_m3_kg-1_s-2"])
    c = float(cfg["constants"]["c_m_s-1"])
    R = float(cfg["sphere_radius_m"])
    iron = read_material(cfg, "iron")
    ice = read_material(cfg, "ice")
    iron_sphere = Sphere(R, iron)
    ice_sphere = Sphere(R, ice)

    equal_radius = {
        name: {
            "mass_kg": body.mass_kg,
            "surface_gravity_m_s2": surface_gravity_m_s2(G, body.mass_kg, R),
            "radius_m": R,
        }
        for name, body in (("iron", iron_sphere), ("ice", ice_sphere))
    }

    # Equal total mass and outer radius: the interior arrangement is inverted.
    # Both arrangements have the SAME VOLUME of iron and of ice, so the external
    # Newtonian fields must be equal under exact spherical symmetry.
    f_iron = float(cfg["layered_bodies"]["iron_volume_fraction"])
    iron_core_radius = R * f_iron**(1.0/3.0)
    ice_core_radius = R * (1.0 - f_iron)**(1.0/3.0)
    r_internal = R * float(cfg["layered_bodies"]["internal_probe_radius_fraction"])
    r_external = R * float(cfg["layered_bodies"]["external_probe_radius_fraction"])
    iron_core = {
        "core_radius_m": iron_core_radius,
        "core_material": "iron",
        "shell_material": "ice",
        "interior_g_m_s2": layered_acceleration_m_s2(
            G, r_internal, R, iron_core_radius, iron.density_kg_m3, ice.density_kg_m3
        ),
        "exterior_g_m_s2": layered_acceleration_m_s2(
            G, r_external, R, iron_core_radius, iron.density_kg_m3, ice.density_kg_m3
        ),
    }
    ice_core = {
        "core_radius_m": ice_core_radius,
        "core_material": "ice",
        "shell_material": "iron",
        "interior_g_m_s2": layered_acceleration_m_s2(
            G, r_internal, R, ice_core_radius, ice.density_kg_m3, iron.density_kg_m3
        ),
        "exterior_g_m_s2": layered_acceleration_m_s2(
            G, r_external, R, ice_core_radius, ice.density_kg_m3, iron.density_kg_m3
        ),
    }
    # Analytically, each has this exact mass irrespective of layer placement.
    same_total_mass = iron_sphere.mass_kg*f_iron + ice_sphere.mass_kg*(1.0-f_iron)

    # Apply the same specific heat input to either material. This is NOT a
    # microscopic nuclear binding-energy calculation, nor a fracture prediction.
    T0 = float(cfg["heating"]["initial_temperature_K"])
    q = float(cfg["heating"]["deposited_energy_J_kg-1"])
    thermal = {
        "initial_temperature_K": T0,
        "deposited_energy_J_kg-1": q,
        "melt_fraction_iron": melt_fraction(iron, T0, q),
        "melt_fraction_ice": melt_fraction(ice, T0, q),
        # If a system retains added energy, that energy also contributes to
        # gravitational mass-energy. This is Δm=Q/c²; in this example tiny.
        "equivalent_added_mass_per_kg_kg": q / c**2,
    }

    # Fragment the ice sphere into 10^3 parts by default. Their combined mass
    # is unchanged; only the allocation of that mass in space has changed.
    fr_cfg = cfg["fragmentation"]
    fragments = make_fragment_grid(
        ice_sphere, int(fr_cfg["fragments_per_axis"]),
        float(fr_cfg["surface_gap_fraction_of_diameter"])
    )
    probes: dict[str, Any] = {}
    for label, probe_x in (
        ("near", float(fr_cfg["near_probe_x_m"])),
        ("far", float(fr_cfg["far_probe_x_m"])),
    ):
        original_x = -G * ice_sphere.mass_kg / probe_x**2
        fragmented_x = acceleration_x_at_probe_m_s2(G, fragments, probe_x)
        probes[label] = {
            "probe_x_m": probe_x,
            "original_acceleration_x_m_s2": original_x,
            "fragmented_acceleration_x_m_s2": fragmented_x,
            "fractional_difference_magnitude": abs(fragmented_x-original_x)/abs(original_x),
        }

    return {
        "experiment_id": "GED-V0-R1",
        "scientific_status": "SCOUT_BASELINE_NOT_EVIDENCE_FOR_INFORMATION_FIELD",
        "model": "Newtonian gravity; uniform-density spheres; prescribed fragmentation; idealized thermal model",
        "equal_radius_comparison": equal_radius,
        "layered_same_mass_radius": {
            "outer_radius_m": R,
            "total_mass_each_kg": same_total_mass,
            "interior_probe_radius_m": r_internal,
            "exterior_probe_radius_m": r_external,
            "iron_core_ice_shell": iron_core,
            "ice_core_iron_shell": ice_core,
            "external_fields_equal_by_spherical_symmetry": math.isclose(
                iron_core["exterior_g_m_s2"], ice_core["exterior_g_m_s2"],
                rel_tol=1e-12, abs_tol=0.0,
            ),
        },
        "thermal_response": thermal,
        "fragmentation_ice": {
            "count": fragments.count,
            "mass_per_fragment_kg": fragments.mass_each_kg,
            "radius_per_fragment_m": fragments.radius_each_m,
            "sum_fragment_mass_kg": fragments.total_mass_kg,
            "original_ice_mass_kg": ice_sphere.mass_kg,
            "mass_conserved": math.isclose(
                fragments.total_mass_kg, ice_sphere.mass_kg, rel_tol=1e-14
            ),
            "fragment_surface_g_m_s2": surface_gravity_m_s2(
                G, fragments.mass_each_kg, fragments.radius_each_m
            ),
            "probes": probes,
        },
    }


def file_sha256(path: Path) -> str:
    # WHAT: Hash file bytes with SHA-256.
    # WHY: A receipt can later establish whether a config, source or output
    # was altered after the experiment ran.
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_commit(directory: Path) -> str | None:
    # WHAT: Read the current Git commit if this is a Git repository.
    # WHY: A commit, together with hashes, precisely identifies code revision.
    # Missing Git is recorded as UNKNOWN rather than inventing provenance.
    try:
        result = subprocess.run(
            ["git", "-C", str(directory), "rev-parse", "HEAD"],
            check=True, capture_output=True, text=True, timeout=3,
        )
        return result.stdout.strip() or None
    except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return None


def pretty_report(result: dict[str, Any]) -> str:
    # WHAT: Produce a simple human-readable report alongside machine JSON.
    # WHY: The user can assess main physical outcomes without reading source.
    equal = result["equal_radius_comparison"]
    layered = result["layered_same_mass_radius"]
    heat = result["thermal_response"]
    frag = result["fragmentation_ice"]
    near = frag["probes"]["near"]
    far = frag["probes"]["far"]
    return "\n".join([
        "GED-V0-R1 — CONVENTIONAL PHYSICS BASELINE (NOT A HYPOTHESIS CONFIRMATION)",
        "="*70,
        f"Equal radius: {equal['iron']['radius_m']:.1f} m",
        f"Iron mass: {equal['iron']['mass_kg']:.10g} kg",
        f"Ice mass:  {equal['ice']['mass_kg']:.10g} kg",
        f"Iron surface gravity: {equal['iron']['surface_gravity_m_s2']:.10g} m/s^2",
        f"Ice surface gravity:  {equal['ice']['surface_gravity_m_s2']:.10g} m/s^2",
        "",
        "LAYERED SPHERES: identical outer radius and total mass",
        f"Each total mass: {layered['total_mass_each_kg']:.10g} kg",
        f"External fields equal: {layered['external_fields_equal_by_spherical_symmetry']}",
        f"Interior g (iron core): {layered['iron_core_ice_shell']['interior_g_m_s2']:.10g} m/s^2",
        f"Interior g (ice core):  {layered['ice_core_iron_shell']['interior_g_m_s2']:.10g} m/s^2",
        "",
        "IDEALIZED HEATING (not vacuum asteroid physics)",
        f"Iron fraction melted: {100*heat['melt_fraction_iron']:.6f}%",
        f"Ice fraction melted:  {100*heat['melt_fraction_ice']:.6f}%",
        "",
        f"FRAGMENTATION: {frag['count']} separated ice fragments",
        f"Mass conserved: {frag['mass_conserved']}",
        f"Fragment radius: {frag['radius_per_fragment_m']:.6g} m",
        f"Individual fragment surface g: {frag['fragment_surface_g_m_s2']:.10g} m/s^2",
        f"Near x={near['probe_x_m']} m relative field change: {near['fractional_difference_magnitude']:.9g}",
        f"Far x={far['probe_x_m']} m relative field change:  {far['fractional_difference_magnitude']:.9g}",
        "",
        "SCIENTIFIC INTERPRETATION:",
        "Physical materials differ. Gravitational fields respond to total mass",
        "and its distribution. This experiment does not detect an additional",
        "gravitational information channel or a material-instruction mechanism.",
        "",
    ])


def main(argv: list[str] | None = None) -> int:
    # WHAT: Parse user choices, run calculations, and write isolated receipts.
    # WHY: This keeps an experiment run auditable and avoids overwriting an
    # earlier run's science files. No external dependencies or network access.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path,
        default=Path(__file__).with_name("GED_V0_R1_config.json"),
        help="Input JSON with SI units (default: bundled config)",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=None,
        help="New output directory; must not already exist",
    )
    args = parser.parse_args(argv)
    config_path = args.config.expanduser().resolve()
    cfg = load_config(config_path)
    result = calculate(cfg)

    # A unique directory preserves previous evidence instead of overwriting it.
    utc_now = datetime.now(timezone.utc)
    if args.output_dir is None:
        stamp = utc_now.strftime("%Y%m%dT%H%M%S%fZ")
        output_dir = Path(__file__).resolve().parent / "runs" / (
            f"GED_V0_R1_{stamp}_{uuid.uuid4().hex[:8]}"
        )
    else:
        output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=False)

    result_file = output_dir / "results.json"
    report_file = output_dir / "report.txt"
    result_file.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    report_file.write_text(pretty_report(result), encoding="utf-8")

    project_dir = Path(__file__).resolve().parent
    # Every present declared input/code file is individually hashed. Missing
    # optional project documents are marked absent, not silently substituted.
    sources = {}
    for name in ("ged_v0_r1.py", "test_ged_v0_r1.py", "GED_V0_R1_config.json", "GED_V0_R1_README.md"):
        path = project_dir / name
        sources[name] = file_sha256(path) if path.is_file() else None
    sources[f"config_used::{config_path}"] = file_sha256(config_path)

    commit = git_commit(project_dir)
    manifest = {
        "experiment_id": "GED-V0-R1",
        "run_utc": utc_now.isoformat(),
        "classification": "UNACCEPTED_SCOUT_BASELINE",
        "git_commit": commit,
        "git_commit_missing": commit is None,
        "command_arguments": list(sys.argv if argv is None else [sys.argv[0], *argv]),
        "working_directory": str(Path.cwd()),
        "python_version": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "pip_dependencies_required": [],
        "source_and_input_sha256": sources,
        "output_sha256": {
            "results.json": file_sha256(result_file),
            "report.txt": file_sha256(report_file),
        },
        "notice": "Provenance supplement only; use official Turnip Run wrapper and an exact Git revision for formal acceptance.",
    }
    (output_dir / "provenance.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True)+"\n", encoding="utf-8"
    )

    print(pretty_report(result))
    print(f"Results and receipts: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
