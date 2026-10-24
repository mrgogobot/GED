#!/usr/bin/env python3
"""GED-V1-R10 — CCMEI / layered gravitational information retention.

WHAT: A three-shell, spherical Newtonian source model plus an exact-rational
      inverse calculation and constituent-energy bookkeeping.
WHY:  Separate (i) a complete constituent mass-energy inventory, (ii) the
      information surviving in a gravitational field, and (iii) any additional
      material-response information that GED has NOT demonstrated.

This program intentionally does not simulate particle binding or solve Einstein's
field equations. A synthetic binding-energy ledger is not a microscopic model.
All physical quantities include explicit SI units; no third-party packages.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from uuid import uuid4

EXPERIMENT_ID = "GED-V1-R10"
CONFIG_NAME = "GED_V1_R10_config.json"
ROOT = Path(__file__).resolve().parent
INPUT_NAMES = (
    "ged_v1_r10.py", "test_ged_v1_r10.py", CONFIG_NAME,
    "GED_V1_R10_README.md", "GED_V1_R10_DEV_AUDIT.md",
    "GED_V1_R10_PRIOR_EVIDENCE.json", "GED_V1_R10_FAILURE_A.log",
)
# WHY: Fixed precision protects tiny mass-energy contributions against rounding
# when they are added to billion-kilogram source masses.
PRECISION = 60


def D(x: object) -> Decimal:
    """WHAT: Convert a stored literal to high-precision decimal (not through float)."""
    if isinstance(x, bool) or x is None:
        raise ValueError("Expected a real-valued decimal literal, not bool/null")
    try:
        n = Decimal(str(x))
    except Exception as exc:
        raise ValueError(f"Not a decimal number: {x!r}") from exc
    if not n.is_finite():
        raise ValueError("Nonfinite number is not permitted")
    return n


def load_configuration(path: Path | None = None) -> dict:
    """WHAT: Parse the frozen JSON. WHY: No hidden constants or silent defaults."""
    config = json.loads((path or ROOT / CONFIG_NAME).read_text(encoding="utf-8"))
    validate(config)
    return config


def validate(cfg: dict) -> None:
    """WHAT: Fail closed on altered physics, invalid input, or missing limitations."""
    if cfg.get("experiment_id") != EXPERIMENT_ID:
        raise ValueError("Wrong experiment ID")
    G, c = D(cfg.get("G_m3_kg_s2")), D(cfg.get("c_m_s"))
    if G <= 0 or c <= 0:
        raise ValueError("G and c must be positive")
    boundaries = [D(x) for x in cfg.get("shell_boundaries_m", [])]
    if len(boundaries) != 4 or boundaries[0] != 0 or any(
        a >= b for a, b in zip(boundaries, boundaries[1:])
    ):
        raise ValueError("Require four strictly ordered shell boundaries starting at zero")
    for name in ("inventory_A", "inventory_B", "inventory_C"):
        obj = cfg.get(name, {})
        if obj.get("label") != name or not obj.get("description"):
            raise ValueError(f"Invalid {name} identity")
        masses = [D(x) for x in obj.get("shell_total_mass_energy_kg", [])]
        binding = [D(x) for x in obj.get("shell_binding_energy_j", [])]
        if len(masses) != 3 or len(binding) != 3 or any(m <= 0 for m in masses):
            raise ValueError(f"Invalid {name} source inventory")
        if any(b > 0 for b in binding):
            raise ValueError("Synthetic binding-energy contribution must be nonpositive")
        if any(m - b / (c*c) <= 0 for m,b in zip(masses,binding)):
            raise ValueError("Invalid constituent rest-energy ledger")
    A, B, C = [[D(x) for x in cfg[name]["shell_total_mass_energy_kg"]]
               for name in ("inventory_A","inventory_B","inventory_C")]
    if sum(A) != sum(B) or A == B or A != C:
        raise ValueError("Source contrasts must hold equal-total / same-shell controls")
    bindA = cfg["inventory_A"]["shell_binding_energy_j"]
    bindC = cfg["inventory_C"]["shell_binding_energy_j"]
    if bindA == bindC:
        raise ValueError("Same-field different-binding control must be different")
    for key in ("probe_radii_m", "inverse_exterior_radii_m", "inverse_interior_radii_m"):
        radii = [D(x) for x in cfg.get(key, [])]
        if not radii or radii != sorted(set(radii)) or any(r <= 0 for r in radii):
            raise ValueError(f"Invalid ordered radii: {key}")
    if any(D(x) <= boundaries[-1] for x in cfg["inverse_exterior_radii_m"]):
        raise ValueError("Exterior probes must be outside source")
    if any(D(x) >= boundaries[-1] for x in cfg["inverse_interior_radii_m"]):
        raise ValueError("Interior probes must be inside source")
    if set(cfg.get("scientific_scope", {}).items()) != {
        ("numerical_only", True), ("new_gravity_detected", False),
        ("general_relativity_solved", False), ("atomic_binding_measured", False),
        ("real_material_microphysics_modeled", False),
        ("equilibrium_or_thermal_pressure_modeled", False),
        ("source_geometry_assumed", True)
    }:
        raise ValueError("Scientific limitations must remain explicit")


def shell_fraction(r: Decimal, low: Decimal, high: Decimal) -> Decimal:
    """WHAT: Fraction of a UNIFORM-density spherical shell enclosed by radius r.

    WHY: For spherical Newtonian gravity only enclosed mass contributes to g(r).
    A shell's volume element is 4*pi*r^2 dr; pi cancels in its volume ratio.
    Assumes spherical symmetry, fixed shell boundaries and uniform shell density.
    """
    r, low, high = D(r), D(low), D(high)
    if r < 0 or low < 0 or high <= low:
        raise ValueError("Invalid shell geometry or radius")
    if r <= low:
        return Decimal(0)
    if r >= high:
        return Decimal(1)
    return (r**3 - low**3) / (high**3 - low**3)


def shell_fraction_exact(r: int, low: int, high: int) -> Fraction:
    """WHAT: Independent EXACT rational version for discrete inverse checks."""
    if not (r >= 0 and 0 <= low < high):
        raise ValueError("Invalid exact shell geometry")
    if r <= low:
        return Fraction(0)
    if r >= high:
        return Fraction(1)
    return Fraction(r**3-low**3, high**3-low**3)


def enclosed_mass(r: Decimal, boundaries: list[Decimal], masses: list[Decimal]) -> Decimal:
    """WHAT: Add the gravitationally sourced mass-equivalent inside radius r.

    WHY: Each shell contributes its mass times the enclosed volume fraction.
    Binding effects are already inside each total shell mass. Never add twice.
    """
    if len(boundaries) != 4 or len(masses) != 3:
        raise ValueError("Expect three shells")
    return sum((m*shell_fraction(r, a, b) for m,a,b in
                zip(masses, boundaries, boundaries[1:])), Decimal(0))


def gravity_magnitude(r: Decimal, boundaries: list[Decimal], masses: list[Decimal], G: Decimal) -> Decimal:
    """WHAT: Newtonian inward radial acceleration magnitude g=G*M(<r)/r².

    WHY: Distinguishes gravitational source distribution from atomic identity.
    The mathematical centre of a regular spherical body has g=0 by symmetry.
    """
    r = D(r)
    if r < 0:
        raise ValueError("Negative radius")
    if r == 0:
        return Decimal(0)
    return G*enclosed_mass(r, boundaries, masses)/(r*r)


def constituent_energy_ledger(masses: list[Decimal], binding: list[Decimal], c: Decimal) -> list[dict]:
    """WHAT: Split each TOTAL shell mass-energy into rest-equivalent + binding.

    WHY: Avoid treating binding energy as a new gravitational charge. A negative
    binding-energy term lowers total energy; compensating constituent energy
    permits two SYNTHETIC inventories with the same shell source energy.

    NOTE: This is a mathematical bookkeeping construction, NOT a claim that
    two realizable materials have the same complete GR stress-energy tensor.
    """
    if len(masses) != len(binding):
        raise ValueError("Inventory size mismatch")
    ledger = []
    for i, (total, bind) in enumerate(zip(masses, binding)):
        rest_equiv = total - bind/(c*c)
        recovered_total = rest_equiv + bind/(c*c)
        ledger.append({
            "shell": i,
            "total_mass_equivalent_kg": str(total),
            "constituent_rest_equivalent_kg": str(rest_equiv),
            "binding_energy_j": str(bind),
            "binding_mass_equivalent_kg": str(bind/(c*c)),
            "recovered_total_mass_equivalent_kg": str(recovered_total),
        })
    return ledger


def response_matrix_exact(probes: list[int], boundaries: list[int]) -> list[list[Fraction]]:
    """WHAT: Linear map from three shell masses to enclosed mass measurements.

    WHY: The measured g(r) can be converted to M_enc = g(r)*r²/G.
    This rescaling changes units, NOT information or matrix rank.
    Each column shows the radial fraction of one known shell's mass.
    """
    if len(boundaries) != 4:
        raise ValueError("Expected three shells")
    return [[shell_fraction_exact(r, a, b) for a,b in zip(boundaries, boundaries[1:])]
            for r in probes]


def matrix_rank_exact(matrix: list[list[Fraction]]) -> int:
    """WHAT: Gaussian elimination over fractions; rank without tolerance choices."""
    if not matrix or not matrix[0] or any(len(row)!=len(matrix[0]) for row in matrix):
        raise ValueError("Must be rectangular nonempty matrix")
    A = [[Fraction(v) for v in row] for row in matrix]
    nr,nc = len(A),len(A[0]); pivot_row=0
    for col in range(nc):
        chosen = next((k for k in range(pivot_row,nr) if A[k][col]), None)
        if chosen is None:
            continue
        A[pivot_row],A[chosen] = A[chosen],A[pivot_row]
        scale=A[pivot_row][col]
        A[pivot_row]=[v/scale for v in A[pivot_row]]
        for row in range(nr):
            if row != pivot_row and A[row][col]:
                factor=A[row][col]
                A[row]=[v-factor*w for v,w in zip(A[row], A[pivot_row])]
        pivot_row+=1
        if pivot_row == nr:
            break
    return pivot_row


def solve_square_exact(matrix: list[list[Fraction]], measurements: list[Fraction]) -> list[Fraction]:
    """WHAT: Recover three unknown shell masses only when observations uniquely identify them.

    WHY: This is an analytic inverse demonstration, not an arbitrary optimiser.
    Fractions keep the noiseless integer-mass source result exact.
    """
    n=len(matrix)
    if n == 0 or len(measurements)!=n or any(len(row)!=n for row in matrix):
        raise ValueError("Exact solver requires a square system")
    A = [[Fraction(v) for v in row]+[Fraction(b)] for row,b in zip(matrix, measurements)]
    for i in range(n):
        pivot=next((r for r in range(i,n) if A[r][i]),None)
        if pivot is None:
            raise ValueError("Singular inverse system")
        A[i],A[pivot]=A[pivot],A[i]
        scale=A[i][i]
        A[i]=[v/scale for v in A[i]]
        for k in range(n):
            if k!=i and A[k][i]:
                f=A[k][i]
                A[k]=[v-f*w for v,w in zip(A[k],A[i])]
    return [row[-1] for row in A]


def run_science(cfg: dict) -> tuple[dict, list[dict]]:
    """WHAT: Compute all pre-registered controls, not just selected successes."""
    validate(cfg)
    with localcontext() as ctx:
        ctx.prec=PRECISION
        G,c = D(cfg["G_m3_kg_s2"]),D(cfg["c_m_s"])
        boundaries=[D(x) for x in cfg["shell_boundaries_m"]]
        cases={key:[D(x) for x in cfg[key]["shell_total_mass_energy_kg"]]
               for key in ("inventory_A","inventory_B","inventory_C")}
        bindings={key:[D(x) for x in cfg[key]["shell_binding_energy_j"]]
                  for key in cases}
        radius_samples=[]
        for s in cfg["probe_radii_m"]:
            r=D(s)
            vals={key:gravity_magnitude(r,boundaries,masses,G) for key,masses in cases.items()}
            radius_samples.append({"radius_m":s,
                **{f"g_{key}_m_s2":str(v) for key,v in vals.items()},
                "delta_A_minus_B_m_s2":str(vals["inventory_A"]-vals["inventory_B"]),
                "delta_A_minus_C_m_s2":str(vals["inventory_A"]-vals["inventory_C"]),
            })
        bounds_int=[int(D(x)) for x in cfg["shell_boundaries_m"]]
        exterior=[int(D(x)) for x in cfg["inverse_exterior_radii_m"]]
        interior=[int(D(x)) for x in cfg["inverse_interior_radii_m"]]
        out_matrix=response_matrix_exact(exterior,bounds_int)
        in_matrix=response_matrix_exact(interior,bounds_int)
        masses=[Fraction(D(x)) for x in cfg["inventory_A"]["shell_total_mass_energy_kg"]]
        data=[sum((a*b for a,b in zip(row,masses)),Fraction(0)) for row in in_matrix]
        recovered=solve_square_exact(in_matrix,data)
        n_source=3
        r_ext=matrix_rank_exact(out_matrix)
        r_int=matrix_rank_exact(in_matrix)
        result={
            "experiment_id":EXPERIMENT_ID,
            "scientific_scope":cfg["scientific_scope"],
            "model":"Spherical Newtonian static uniform-density shells; synthetic CCMEI bookkeeping",
            "shell_geometry_m":[str(x) for x in boundaries],
            "total_mass_equivalent_kg":str(sum(cases["inventory_A"])),
            "case_inventory_ledgers":{
                key:constituent_energy_ledger(cases[key],bindings[key],c) for key in cases},
            "exterior_rank":r_ext,"exterior_nullity":n_source-r_ext,
            "interior_rank":r_int,"interior_nullity":n_source-r_int,
            "interior_recovered_shell_mass_kg":[str(Decimal(v.numerator)/Decimal(v.denominator)) for v in recovered],
            "reconstruction_exact":recovered==masses,
            "same_total_mass_A_B":sum(cases["inventory_A"])==sum(cases["inventory_B"]),
            "same_radial_mass_A_C":cases["inventory_A"]==cases["inventory_C"],
            "different_binding_ledger_A_C":bindings["inventory_A"]!=bindings["inventory_C"],
            "external_difference_exactly_zero":all(
                D(row["delta_A_minus_B_m_s2"])==0 for row in radius_samples if D(row["radius_m"])>=boundaries[-1]),
            "different_internal_field_A_B":any(
                D(row["delta_A_minus_B_m_s2"])!=0 for row in radius_samples if D(row["radius_m"])<boundaries[-1]),
            "same_field_A_C_at_all_probes":all(D(row["delta_A_minus_C_m_s2"])==0 for row in radius_samples),
            "hypothesis_assessment":{
                "established":"Radial mass-energy contributions obey the Newtonian shell theorem; selected internal probes can recover assumed shell totals.",
                "unresolved":"Neither atomic identities, binding rules, nor a new material-response channel are gravitationally recovered by this toy model.",
                "new_gravity_evidence":False,
                "falsifiable_next_step":"Specify and derive a source-dependent extra observable, independently constrained by physical theory and measurement, without double-counting stress-energy."
            },
            "probe_results":radius_samples,
        }
        return result, radius_samples


def sha256(path: Path) -> str:
    """WHAT: Hash exact bytes. WHY: Preserve reproducibility and detect modifications."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_context() -> dict:
    """WHAT: Capture Git identity if this local directory happens to be a repository.

    WHY: The SHA-256 snapshots stand alone, but a Git commit is useful provenance.
    This function NEVER stages, commits, pushes, tags or changes Git state.
    """
    def cmd(*args: str) -> str | None:
        p=subprocess.run(["git",*args],cwd=ROOT,capture_output=True,text=True,check=False)
        return p.stdout.strip() if p.returncode==0 else None
    try:
        commit=cmd("rev-parse","HEAD")
        dirty=cmd("status","--porcelain") if commit else None
        return {"commit":commit,"working_tree_dirty":bool(dirty) if dirty is not None else None,
                "git_commit_covers_this_experiment":False}
    except FileNotFoundError:
        return {"commit":None,"working_tree_dirty":None,"git_commit_covers_this_experiment":False}


def execute() -> Path:
    """WHAT: Run science and tests and freeze every input/output with hashes.

    WHY: A convenient console result is NOT enough for the Manifest of Destiny.
    An immutable-style timestamped run folder retains what was actually executed.
    """
    cfg=load_configuration()
    results, samples = run_science(cfg)
    timestamp=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    run=ROOT/"runs"/(EXPERIMENT_ID.replace("-","_")+"_"+timestamp+"_"+uuid4().hex[:8])
    run.mkdir(parents=True,exist_ok=False)
    snap=run/"input_snapshot";snap.mkdir()
    for name in INPUT_NAMES:
        shutil.copy2(ROOT/name,snap/name)
    # WHY: Unit test output is retained even if a test unexpectedly fails;
    # the run is marked failed rather than falsely producing an accepted result.
    proc=subprocess.run([sys.executable,"-m","unittest","-v","test_ged_v1_r10.py"],
                        cwd=ROOT,capture_output=True,text=True,check=False)
    (run/"tests.txt").write_text("COMMAND: python -m unittest -v test_ged_v1_r10.py\n"
                                  +proc.stdout+proc.stderr+f"\nEXIT CODE: {proc.returncode}\n",encoding="utf-8")
    (run/"results.json").write_text(json.dumps(results,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    with (run/"radial_samples.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(samples[0]));w.writeheader();w.writerows(samples)
    report=(
        "GED-V1-R10 — CCMEI / GRAVITATIONAL LAYER INFORMATION\n"
        "========================================================\n"
        "CONVENTIONAL NEWTONIAN MODEL; NO NOVEL GRAVITY OR ATOMIC INFORMATION CHANNEL DETECTED\n"
        f"Each object's total gravitational source mass-equivalent: {results['total_mass_equivalent_kg']} kg\n"
        f"Three-shell exterior-only rank: {results['exterior_rank']}/3; nullity {results['exterior_nullity']}\n"
        f"Selected interior-probe rank: {results['interior_rank']}/3; nullity {results['interior_nullity']}\n"
        f"Interior recovered source masses: {results['interior_recovered_shell_mass_kg']} kg\n"
        f"Exact noiseless recovery under assumed shell boundaries: {results['reconstruction_exact']}\n"
        f"A vs B differ internally: {results['different_internal_field_A_B']}; exterior equal: {results['external_difference_exactly_zero']}\n"
        f"A vs C different synthetic binding ledger, same field at all probes: {results['same_field_A_C_at_all_probes']}\n"
        "NOTE: Cases A and C are bookkeeping-equivalent gravitational source models, not\n"
        "demonstrated realizable materials with identical complete GR stress-energy.\n"
        "LIMITS: static Newtonian spheres; no atomic microphysics, stress pressure, dynamics or full GR.\n"
        f"Tests: {'PASS' if proc.returncode==0 else 'FAIL'}\n"
        "SCIENTIFIC VERDICT: No new gravitational information-channel mechanism established.\n"
        f"Results and receipts: {run}\n"
    )
    (run/"report.txt").write_text(report,encoding="utf-8")
    record={"experiment_id":EXPERIMENT_ID,"run_utc":timestamp,
            "python":platform.python_version(),"python_executable":sys.executable,
            "platform":platform.platform(),"git":git_context(),
            "command":"python ged_v1_r10.py", "test_command":"python -m unittest -v test_ged_v1_r10.py",
            "test_exit_code":proc.returncode,
            "input_sha256":{name:sha256(snap/name) for name in INPUT_NAMES},
            "output_sha256":{name:sha256(run/name) for name in
                             ("tests.txt","results.json","radial_samples.csv","report.txt")},
            "publication":"LOCAL_ONLY; do not push or submit to Zenodo without explicit approval",
            "scientific_verdict":"Numerical Newtonian control; no new physics established"}
    (run/"provenance.json").write_text(json.dumps(record,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    entries=sorted(p for p in run.rglob("*") if p.is_file() and p.name!="SHA256SUMS")
    (run/"SHA256SUMS").write_text("".join(f"{sha256(p)}  {p.relative_to(run)}\n" for p in entries),encoding="utf-8")
    print(report,end="")
    if proc.returncode:
        raise RuntimeError("R10 tests FAILED; partial run and receipts have been preserved")
    return run


if __name__=="__main__":
    execute()
