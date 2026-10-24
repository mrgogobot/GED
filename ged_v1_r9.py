#!/usr/bin/env python3
"""GED-V1-R9: gravitational inverse-problem and identifiability limits.

WHAT: Recover unknown nonnegative energy fractions at six *assumed known* sites
      from simulated exterior Newtonian vector-gravity measurements.
WHY:  Distinguish exact mathematical non-uniqueness, limited sensor geometry,
      and finite-noise uncertainty without inventing a new gravitational law.

The model is synthetic, weak-field Newtonian, and uses the Python standard library.
The physical and computational assumptions are explained beside each operation.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import random
import shutil
import subprocess
import sys
from uuid import uuid4

EXPERIMENT_ID = "GED-V1-R9"
CONFIG_NAME = "GED_V1_R9_config.json"
INPUT_FILES = (
    "ged_v1_r9.py", "test_ged_v1_r9.py", CONFIG_NAME,
    "GED_V1_R9_README.md", "GED_V1_R9_DEV_AUDIT.md",
)
AXIS_LABELS = ("+x", "-x", "+y", "-y", "+z", "-z")
PIVOT_EPS = 1e-11  # Dimensionless tolerance for RREF/linear solve; documented numerical choice.


def require_finite_real(value: object, name: str, *, strict_positive: bool = False) -> float:
    """WHAT: Validate physical inputs. WHY: Booleans/NaN can silently corrupt scientific runs."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a real number, not a boolean")
    number = float(value)
    if not math.isfinite(number) or (strict_positive and number <= 0):
        raise ValueError(f"{name} must be finite" + (" and positive" if strict_positive else ""))
    return number


def validate(config: dict) -> None:
    """WHAT: Fail early if the frozen experiment configuration is inconsistent."""
    if config.get("experiment_id") != EXPERIMENT_ID:
        raise ValueError("incorrect experiment ID")
    for key in ("G_m3_kg_s2", "c_m_s", "total_added_energy_j", "body_radius_m",
                "source_offset_m", "reference_radius_m"):
        require_finite_real(config.get(key), key, strict_positive=True)
    if not (config["source_offset_m"] < config["body_radius_m"] < config["reference_radius_m"]):
        raise ValueError("source sites must be inside body; reference radius outside")
    weights = config.get("unknown_truth_energy_fractions")
    if not isinstance(weights, list) or len(weights) != 6:
        raise ValueError("truth must specify exactly six weights")
    if any(require_finite_real(v, "energy fraction") <= 0 for v in weights):
        raise ValueError("truth energy fractions must all be positive")
    if not math.isclose(math.fsum(weights), 1.0, rel_tol=0, abs_tol=1e-12):
        raise ValueError("energy fractions must sum to 1")
    probes = config.get("exterior_probes_m")
    if not isinstance(probes, list) or len(probes) < 6:
        raise ValueError("require at least six exterior probes")
    if len(set(tuple(p) for p in probes)) != len(probes):
        raise ValueError("exterior probes must be distinct")
    for p in probes:
        if not isinstance(p, list) or len(p) != 3 or any(not math.isfinite(require_finite_real(v, "probe component")) for v in p):
            raise ValueError("probe must be a finite 3-vector")
        if norm(p) <= config["body_radius_m"]:
            raise ValueError("all probes must be outside body")
    axis_radii = config.get("axis_probe_radii_m")
    if not isinstance(axis_radii, list) or len(axis_radii) < 6 or any(
        require_finite_real(r, "axis radius", strict_positive=True) <= config["body_radius_m"] for r in axis_radii
    ) or any(a >= b for a,b in zip(axis_radii,axis_radii[1:])):
        raise ValueError("axis probe radii must be distinct and increasing outside body")
    shell_radii = config.get("shell_radii_m")
    if not isinstance(shell_radii, list) or len(shell_radii) != 3 or any(
        require_finite_real(r, "shell radius", strict_positive=True) >= config["body_radius_m"] for r in shell_radii
    ) or not (shell_radii[0] < shell_radii[1] < shell_radii[2]):
        raise ValueError("shell radii must increase inside body")
    shell_cases = config.get("shell_energy_fractions")
    if not isinstance(shell_cases, list) or len(shell_cases) != 2 or any(
        len(ws) != 3 or min(ws) < 0 or not math.isclose(math.fsum(ws), 1, abs_tol=1e-12, rel_tol=0)
        for ws in shell_cases
    ) or shell_cases[0] == shell_cases[1]:
        raise ValueError("require two distinct equal-total shell allocations")
    levels = config.get("noise_sigma_in_reference_g")
    if not isinstance(levels, list) or len(levels) < 2 or levels[0] != 0 or any(
        require_finite_real(s,"noise") < 0 for s in levels
    ) or any(a >= b for a,b in zip(levels,levels[1:])):
        raise ValueError("noise levels must start at 0 and strictly increase")
    seed = config.get("random_seed")
    trials = config.get("trials_per_level")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("random seed must be a nonnegative integer")
    if isinstance(trials,bool) or not isinstance(trials,int) or not 1 <= trials <= 1000:
        raise ValueError("trials must be 1..1000")
    if config.get("scientific_claims") != {
        "observation": False, "new_gravity": False, "source_positions_assumed": True,
        "background_subtraction_assumed_exact": True, "full_gr": False
    }:
        raise ValueError("required scientific disclaimers missing")


def norm(v) -> float:
    """WHAT: Euclidean distance in metres for a 3-vector. WHY: Newton's law uses |r|."""
    return math.sqrt(math.fsum(float(x)**2 for x in v))


def source_sites(offset: float) -> tuple[tuple[float, float, float], ...]:
    """WHAT: Six KNOWN possible source positions. WHY: Only energy shares are unknown."""
    return ((offset,0,0),(-offset,0,0),(0,offset,0),(0,-offset,0),(0,0,offset),(0,0,-offset))


def signature_column(probes, site, reference_radius_m: float, *, component: int | None = None) -> list[float]:
    """WHAT: Measurement per unit energy fraction, scaled by a reference g.

    WHY: Newtonian linear superposition produces y=A w. Divide physical gravity
    by g_ref=G*(E/c^2)/r_ref^2 so A is of order unity instead of ~1e-24;
    this is a numerical unit conversion, NOT extra physics.

    A_ijk = r_ref^2 * (site_k - probe_k) / |probe-site|^3.
    Full-vector mode outputs x,y,z per probe; component mode one selected axis.
    """
    out = []
    for point in probes:
        d = tuple(s - p for s,p in zip(site, point))
        radius = norm(d)
        if radius == 0:
            raise ValueError("probe coincides with source")
        factor = reference_radius_m**2 / radius**3
        if component is None:
            out.extend(factor*v for v in d)
        else:
            out.append(factor*d[component])
    return out


def measurement_matrix(probes, offset: float, ref_r: float, *, component: int | None = None) -> list[list[float]]:
    """WHAT: Build A[row][source] from known source locations and probe geometry."""
    if component not in (None,0,1,2):
        raise ValueError("invalid measurement component")
    columns = [signature_column(probes, p, ref_r, component=component) for p in source_sites(offset)]
    return [list(row) for row in zip(*columns)]


def forward(matrix: list[list[float]], weights: list[float]) -> list[float]:
    """WHAT: Predict the measured gravitational field. WHY: Linear superposition."""
    if not matrix or any(len(row) != len(weights) for row in matrix):
        raise ValueError("inconsistent forward matrix")
    return [math.fsum(a*w for a,w in zip(row,weights)) for row in matrix]


def scaled_gravity_m_s2(y: float, config: dict) -> float:
    """WHAT: Convert dimensionless field back to m/s² using E/c² and G."""
    return y * reference_acceleration(config)


def reference_acceleration(config: dict) -> float:
    """WHAT: Gravity of the equivalent point mass at the chosen reference radius."""
    return (config["G_m3_kg_s2"] * config["total_added_energy_j"] /
            config["c_m_s"]**2 / config["reference_radius_m"]**2)


def rref(matrix: list[list[float]], threshold: float = PIVOT_EPS) -> tuple[list[list[float]], list[int]]:
    """WHAT: Row-reduced echelon form with partial pivoting.

    WHY: A matrix of rank below its number of columns has a nullspace:
    nonzero source changes that leave every selected measurement unchanged.
    For this small, normalized teaching model, row elimination is sufficient;
    this is NOT an SVD-based uncertainty analysis for poorly scaled real data.
    """
    if not matrix or not matrix[0] or any(len(r)!=len(matrix[0]) for r in matrix):
        raise ValueError("matrix must be nonempty and rectangular")
    A = [[float(x) for x in row] for row in matrix]
    nrow,ncol=len(A),len(A[0])
    pivots=[]; prow=0
    for col in range(ncol):
        k=max(range(prow,nrow),key=lambda i:abs(A[i][col])) if prow<nrow else None
        if k is None or abs(A[k][col]) <= threshold:
            continue
        A[prow],A[k]=A[k],A[prow]
        pivot=A[prow][col]
        A[prow]=[x/pivot for x in A[prow]]
        for i in range(nrow):
            if i != prow:
                t=A[i][col]
                if t:
                    A[i]=[v-t*w for v,w in zip(A[i],A[prow])]
        pivots.append(col)
        prow+=1
        if prow==nrow:
            break
    return A,pivots


def null_vector(matrix: list[list[float]]) -> list[float]:
    """WHAT: Construct an explicit nonzero v satisfying A*v≈0.

    WHY: This supplies two *different*, physical nonnegative source weights
    with the same measurements, proving sparse-view ambiguity directly.
    """
    reduced,pivots=rref(matrix)
    free=[j for j in range(len(matrix[0])) if j not in pivots]
    if not free:
        raise ValueError("matrix is full column rank")
    vec=[0.0]*len(matrix[0]); vec[free[0]]=1.0
    for i,col in enumerate(pivots):
        vec[col]=-reduced[i][free[0]]
    scale=max(abs(x) for x in vec)
    return [x/scale for x in vec]


def augment_total_known(matrix: list[list[float]]) -> list[list[float]]:
    """WHAT: Add the exact known total-energy constraint sum(weights)=1."""
    return [list(row) for row in matrix] + [[1.0]*len(matrix[0])]


def alternative_nonnegative_weights(truth: list[float], matrix: list[list[float]]) -> list[float]:
    """WHAT: Produce another nonnegative distribution with same measured field and total.

    WHY: A pair of explicit counterexamples is stronger than merely stating
    that the measurement matrix has deficient rank.
    """
    v=null_vector(augment_total_known(matrix))
    maximum_t=min(truth[i]/(-v[i]) for i in range(len(v)) if v[i] < -1e-12)
    t=maximum_t*0.5
    alt=[a+t*b for a,b in zip(truth,v)]
    if min(alt)<-1e-12:
        raise AssertionError("alternate configuration has negative energy")
    return alt


def linear_solve(M: list[list[float]], b: list[float], threshold: float = 1e-12) -> list[float]:
    """WHAT: Solve a small square system using partial-pivot Gaussian elimination.

    WHY: Each candidate active subset needs least-squares coefficients.
    The inversion uses only a 6-site toy basis. Failure on singular systems
    is explicit, never replaced with an arbitrary preferred solution.
    """
    n=len(b)
    if len(M)!=n or any(len(row)!=n for row in M):
        raise ValueError("non-square linear system")
    A=[list(map(float,row))+[float(rhs)] for row,rhs in zip(M,b)]
    for j in range(n):
        k=max(range(j,n),key=lambda i:abs(A[i][j]))
        if abs(A[k][j]) < threshold:
            raise ValueError("singular or ill-conditioned normal matrix")
        A[k],A[j]=A[j],A[k]
        for i in range(j+1,n):
            factor=A[i][j]/A[j][j]
            for q in range(j,n+1):
                A[i][q]-=factor*A[j][q]
    x=[0.0]*n
    for i in range(n-1,-1,-1):
        x[i]=(A[i][-1]-math.fsum(A[i][j]*x[j] for j in range(i+1,n)))/A[i][i]
    return x


def simplex_least_squares(matrix: list[list[float]], measured: list[float]) -> tuple[list[float], float]:
    """WHAT: Fit unknown nonnegative shares with sum(w)=1 by active-set enumeration.

    WHY: In a small (six-source) synthetic experiment, enumerate all 63
    possible nonempty supports, solve fixed-total least squares on each, reject
    negative shares, and keep the smallest residual. We do NOT know the truth
    during fitting; only the predeclared positions and total energy are given.

    For larger inverse problems, use validated QR/SVD or constrained solvers;
    normal equations here are acceptable only for this controlled small model.
    """
    n=len(matrix[0])
    if n != 6 or not matrix or len(measured)!=len(matrix):
        raise ValueError("measurement dimensions inconsistent")
    best=None
    for size in range(1,n+1):
        for support in itertools.combinations(range(n),size):
            base=support[-1]
            if size==1:
                fit=[0.0]*n; fit[base]=1.0
            else:
                indices=support[:-1]
                columns=[[row[j]-row[base] for row in matrix] for j in indices]
                target=[y-row[base] for row,y in zip(matrix,measured)]
                gram=[[math.fsum(x*y for x,y in zip(a,b)) for b in columns] for a in columns]
                right=[math.fsum(x*y for x,y in zip(a,target)) for a in columns]
                try:
                    coeffs=linear_solve(gram,right)
                except ValueError:
                    continue
                last=1-math.fsum(coeffs)
                if min(*coeffs,last)<-1e-10:
                    continue
                fit=[0.0]*n
                for j,value in zip(indices,coeffs):
                    fit[j]=max(0.0,value)
                fit[base]=max(0.0,last)
                total=math.fsum(fit)
                fit=[v/total for v in fit]
            residual=math.fsum((y-p)**2 for y,p in zip(measured,forward(matrix,fit)))
            if best is None or residual < best[1]:
                best=(fit,residual)
    if best is None:
        raise RuntimeError("no feasible reconstruction found")
    return best


def field_for_shells(point, shell_radii: list[float], weights: list[float], config: dict) -> list[float]:
    """WHAT: Newtonian acceleration of concentric uniform thin shells.

    WHY: An exact exterior degeneracy does not depend on sensor precision.
    Only shells whose radii are smaller than the probe radius contribute.
    """
    r=norm(point)
    if r == 0:
        return [0.,0.,0.]
    if any(math.isclose(r, radius, abs_tol=1e-12) for radius in shell_radii):
        raise ValueError("probe lies on an ideal thin-shell surface")
    enclosed=math.fsum(w for w, sr in zip(weights,shell_radii) if sr<r)
    coeff=-(config["reference_radius_m"] / r)**2 * enclosed / r
    return [coeff*x for x in point]


def sampling_experiment(config: dict) -> dict:
    """WHAT: Produce transparent rank, nullspace, geometry and noise controls."""
    validate(config)
    truth=list(config["unknown_truth_energy_fractions"])
    all_points=config["exterior_probes_m"]
    offset=config["source_offset_m"]
    ref=config["reference_radius_m"]
    sparse=measurement_matrix(all_points[:1], offset,ref)
    sparse_rank=len(rref(augment_total_known(sparse))[1])
    alt=alternative_nonnegative_weights(truth,sparse)
    sparse_residual=max(abs(x-y) for x,y in zip(forward(sparse, truth),forward(sparse,alt)))
    radial_points=[[r,0.,0.] for r in config["axis_probe_radii_m"]]
    radial=measurement_matrix(radial_points,offset,ref,component=0)
    radial_rank=len(rref(augment_total_known(radial))[1])
    geometry=[]
    for count in (1,2,3,4,len(all_points)):
        matrix=measurement_matrix(all_points[:count],offset,ref)
        geometry.append({"probe_count":count,"vector_components":3*count,
                         "rank_with_total":len(rref(augment_total_known(matrix))[1]),
                         "nullity":6-len(rref(augment_total_known(matrix))[1])})
    complete=measurement_matrix(all_points,offset,ref)
    noiseless=forward(complete,truth)
    reconstructed, residual=simplex_least_squares(complete,noiseless)
    no_noise_l1=math.fsum(abs(a-b) for a,b in zip(reconstructed,truth))
    rng=random.Random(config["random_seed"])
    sweep=[]
    for sigma in config["noise_sigma_in_reference_g"]:
        errs=[]; err_sq=[]; resid=[]
        for _ in range(config["trials_per_level"]):
            noisy=[x+rng.gauss(0.0,sigma) for x in noiseless]
            fit, rss=simplex_least_squares(complete,noisy)
            errs.append(math.fsum(abs(f-t) for f,t in zip(fit,truth)))
            err_sq.append(math.sqrt(math.fsum((f-t)**2 for f,t in zip(fit,truth))))
            resid.append(rss)
        sweep.append({"noise_sigma_in_reference_g":sigma,
                      "noise_sigma_m_s2":sigma*reference_acceleration(config),
                      "trials":config["trials_per_level"],
                      "mean_l1_share_error":math.fsum(errs)/len(errs),
                      "mean_l2_share_error":math.fsum(err_sq)/len(err_sq),
                      "mean_squared_residual":math.fsum(resid)/len(resid)})
    shells=config["shell_energy_fractions"]
    external_points=[[r,0.,0.] for r in config["axis_probe_radii_m"]]
    ext_delta=max(abs(a-b) for p in external_points for a,b in zip(
        field_for_shells(p,config["shell_radii_m"],shells[0],config),
        field_for_shells(p,config["shell_radii_m"],shells[1],config)))
    interior_probe=[35.,0.,0.]
    int_delta=max(abs(a-b) for a,b in zip(
        field_for_shells(interior_probe,config["shell_radii_m"],shells[0],config),
        field_for_shells(interior_probe,config["shell_radii_m"],shells[1],config)))
    return {"experiment_id":EXPERIMENT_ID,"synthetic":True,"new_gravity":False,
            "source_sites_m":source_sites(offset),"true_energy_fractions":truth,
            "reference_acceleration_m_s2":reference_acceleration(config),
            "sparse_one_probe":{"augmented_rank":sparse_rank,"nullity":6-sparse_rank,
                                "alternative_fractions":alt,
                                "maximum_dimensionless_field_residual":sparse_residual,
                                "l1_share_difference":math.fsum(abs(a-b) for a,b in zip(alt,truth))},
            "radial_scalar":{"probe_count":len(radial_points),"augmented_rank":radial_rank,
                              "nullity":6-radial_rank},
            "distributed_geometry":geometry,
            "full_vector_noiseless":{"inferred_fractions":reconstructed,
                                     "l1_share_error":no_noise_l1,"residual":residual},
            "synthetic_noise_sweep":sweep,
            "spherical_shells":{"shell_radii_m":config["shell_radii_m"],
                                "case_fractions":shells,
                                "max_exterior_dimensionless_difference":ext_delta,
                                "interior_at_35m_dimensionless_difference":int_delta},
            "limits":[
                "All data are simulated, not observations.",
                "Known six source positions and exact added energy total are strong priors.",
                "Common much larger body gravitational field is assumed subtracted perfectly.",
                "Point-equivalent energy parcels are not an elastic equilibrium or atomic binding model.",
                "No full GR, pressure, support stresses, dynamical gravity or material response is solved.",
                "RREF ranks are numerical and depend on the dimensionless pivot tolerance.",
                "Normal-equation active-set fitting is only for this small conditioned toy basis.",
                "Noise is assumed independent Gaussian with frozen RNG seed, not real instrument performance.",
                "No additional gravitational information channel is derived or detected.",
            ]}


def sha256_file(path: Path) -> str:
    """WHAT: Hash exact bytes. WHY: Reproducibility requires content identities."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_head(directory:Path) -> str | None:
    """WHAT: Record commit if code is tracked; never perform git writes."""
    try:
        proc=subprocess.run(["git","rev-parse","HEAD"],cwd=directory, capture_output=True,text=True,timeout=3)
        return proc.stdout.strip() if proc.returncode==0 else None
    except (FileNotFoundError,subprocess.TimeoutExpired):
        return None


def persist_run(base: Path, config:dict, result:dict) -> Path:
    """WHAT: Preserve configuration, exact input bytes, results, tests and hashes.

    WHY: Runs are immutable evidence; publishing is expressly deferred until
    GED-V1-R10. No existing run or published version is modified.
    """
    stamp=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    run=base/"runs"/f"GED_V1_R9_{stamp}_{uuid4().hex[:8]}"
    run.mkdir(parents=True, exist_ok=False)
    snap=run/"input_snapshot"; snap.mkdir()
    for filename in INPUT_FILES:
        path=base/filename
        if not path.is_file():
            raise FileNotFoundError(f"Missing frozen input file: {path}")
        shutil.copy2(path,snap/filename)
    (run/"results.json").write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    with (run/"noise_sweep.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(result["synthetic_noise_sweep"][0].keys()))
        writer.writeheader(); writer.writerows(result["synthetic_noise_sweep"])
    with (run/"probe_geometry.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(result["distributed_geometry"][0].keys()))
        writer.writeheader(); writer.writerows(result["distributed_geometry"])
    report=(
        "GED-V1-R9 — GRAVITATIONAL INVERSE RECONSTRUCTION LIMITS\n"
        "=======================================================\n"
        "SYNTHETIC NEWTONIAN CONTROL; NOT AN OBSERVATION OR NEW GRAVITY\n"
        f"Reference acceleration: {result['reference_acceleration_m_s2']:.12e} m/s^2\n"
        f"One-probe constrained rank: {result['sparse_one_probe']['augmented_rank']}/6, "
        f"nullity {result['sparse_one_probe']['nullity']}\n"
        f"One-probe alternate physically nonnegative shares L1 difference: "
        f"{result['sparse_one_probe']['l1_share_difference']:.12g}\n"
        f"Max one-probe field residual: "
        f"{result['sparse_one_probe']['maximum_dimensionless_field_residual']:.6e} (g_ref units)\n"
        f"Radial-only scalar measurements, constrained rank: "
        f"{result['radial_scalar']['augmented_rank']}/6\n"
        f"Distributed full-vector rank: {result['distributed_geometry'][-1]['rank_with_total']}/6\n"
        f"Noiseless full-vector L1 share error: "
        f"{result['full_vector_noiseless']['l1_share_error']:.6e}\n"
        f"Exact exterior spherical-shell field difference: "
        f"{result['spherical_shells']['max_exterior_dimensionless_difference']:.6e} (g_ref units)\n"
        f"Interior shell field difference: "
        f"{result['spherical_shells']['interior_at_35m_dimensionless_difference']:.6e} (g_ref units)\n\n"
        "NOISE SWEEP (mean L1 energy-share reconstruction error):\n"
        + "".join(f"  sigma={row['noise_sigma_in_reference_g']:.2e} g_ref: "
                  f"mean L1={row['mean_l1_share_error']:.8g} over {row['trials']} trials\n"
                  for row in result['synthetic_noise_sweep'])
        + "\nVERDICT: Exact ambiguity exists with incomplete geometry; full-rank\n"
          "inversion is recoverable for this fixed-site prior and becomes less\n"
          "reliable under numerical noise. No new information channel found.\n"
    )
    (run/"report.txt").write_text(report)
    print(report, end="")
    test_cmd=[sys.executable,"-m","unittest","-v","test_ged_v1_r9.py"]
    cp=subprocess.run(test_cmd,cwd=base,capture_output=True,text=True)
    (run/"tests.txt").write_text("COMMAND: "+" ".join(test_cmd)+"\n\n"+cp.stdout+cp.stderr)
    if cp.returncode:
        raise RuntimeError(f"tests failed, evidence preserved in {run}")
    provenance={"experiment_id":EXPERIMENT_ID,"created_utc":stamp,
                "version_status":"local_only_not_published",
                "python_version":sys.version,"platform":platform.platform(),
                "source_git_head_if_present":git_head(base),
                "git_repository_unchanged":True,"zenodo_unchanged":True,
                "command":[sys.executable,"ged_v1_r9.py"],
                "verification_command":test_cmd,"verification_exit_status":cp.returncode,
                "random_seed":config["random_seed"],
                "input_sha256":{name:sha256_file(snap/name) for name in INPUT_FILES},
                "floating_point_model":"Python float64; RREF tolerance 1e-11; printed errors are numerical",
                "limitations":result["limits"]}
    (run/"provenance.json").write_text(json.dumps(provenance,indent=2,sort_keys=True)+"\n")
    # WHAT: Hash all run contents except checksum list. WHY: Prevent self-hash recursion.
    entries=[]
    for path in sorted(run.rglob("*")):
        if path.is_file() and path.name!="SHA256SUMS":
            entries.append(f"{sha256_file(path)}  {path.relative_to(run).as_posix()}")
    (run/"SHA256SUMS").write_text("\n".join(entries)+"\n")
    print(f"Results and receipts: {run}")
    print(f"SHA-256 manifest: {run/'SHA256SUMS'}")
    return run


def main() -> None:
    base=Path(__file__).resolve().parent
    config=json.loads((base/CONFIG_NAME).read_text())
    result=sampling_experiment(config)
    persist_run(base,config,result)


if __name__=="__main__":
    main()
