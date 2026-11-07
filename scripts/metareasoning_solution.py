#!/usr/bin/env python3
"""
CS-023/024 — Certified Search metareasoning solution for LJ13.

Strategy:
  1. Baseline: reproduce the original uniform random initialization and
     L-BFGS-B optimization (maxiter=500, gtol=1e-6), then always verify.
  2. Proposed policy: initialize LJ13 near its known icosahedral N=13
     structure (one central atom + 12 shell vertices), add deterministic
     seed-specific perturbation, and refine with tight L-BFGS-B tolerances.
  3. Submit every accepted numerical candidate to the unchanged C++ verifier.
  4. If the structured candidate fails the numerical quality gate or the
     verifier does not certify it, fall back to the original random-start
     baseline for that seed and verify the fallback candidate.

No verifier output is used to train or tune a predictor.  The strategy is a
physics-informed change to the initial-coordinate distribution, followed by
an unchanged rigorous verifier.  The C++ verifier is never modified.

Default benchmark:
  3 paired series x 100 seeds
  series 1: 1000..1099
  series 2: 1100..1199
  series 3: 1200..1299

Run from the project root or invoke this file by absolute/relative path:
  python3 scripts/metareasoning_solution.py

Optional local optimizer-only check (does not require the C++ verifier):
  python3 scripts/metareasoning_solution.py --self-test
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

try:
    import scipy
    from scipy.optimize import minimize
except ImportError:
    print("ERROR: scipy is required. Install with: python3 -m pip install scipy", file=sys.stderr)
    raise SystemExit(1)


# ---------------------------------------------------------------------------
# Project paths and fixed experiment parameters
# ---------------------------------------------------------------------------
SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parent.parent
VERIFIER_BIN = PROJECT_ROOT / "core" / "verifier" / "build" / "krawczyk_lj13"
RESULTS_PATH = PROJECT_ROOT / "metareasoning_solution_results.json"
N_ATOMS = 13
TASK_HASH = "sha256:cs023_024_lj13_strict_local_minimum_v1"

# Baseline parameters intentionally match the supplied baseline.
BASELINE_MAXITER = 500
BASELINE_GTOL = 1e-6

# Physics-informed candidate policy.
ICOSAHEDRON_SHELL_RADIUS = 1.05
ICOSAHEDRON_NOISE_SIGMA = 0.01
POLICY_MAXITER = 120
POLICY_GTOL = 1e-10
POLICY_FTOL = 0.0
POLICY_MAXLS = 100
POLICY_MAXCOR = 30

# This is only a numerical routing gate.  It does NOT certify anything.
# Candidates outside the expected LJ13 global-minimum basin fall back to the
# original uniform random initialization and still go through the verifier.
KNOWN_LJ13_MINIMUM = -44.326801419
NUMERICAL_ENERGY_TOL = 0.02
NUMERICAL_GRAD_NORM_MAX = 1e-5

TARGET_COST_PER_CERTIFICATE = 0.0370
MIN_CERTIFICATES_PER_SERIES = 45
TARGET_GAIN_PERCENT = 20.0


# ---------------------------------------------------------------------------
# Lennard-Jones energy and analytic gradient (vectorized, no Hessian)
# ---------------------------------------------------------------------------
def lj_energy_and_grad_vectorized(coords_flat: np.ndarray, n_atoms: int = N_ATOMS):
    """LJ 12-6 energy and Cartesian gradient for a flat 3N vector.

    The expression counts each pair twice and applies the corresponding
    factor 2, matching the conventional 4 * sum_{i<j}(r^-12-r^-6).
    """
    coords = np.asarray(coords_flat, dtype=np.float64).reshape((n_atoms, 3))
    diffs = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    r2 = np.einsum("ijk,ijk->ij", diffs, diffs)
    np.fill_diagonal(r2, np.inf)

    # Collision guard.  The random and icosahedral generators should never
    # trigger this, but avoid NaN/Inf propagation if coordinates are invalid.
    finite_pair_r2 = r2[np.isfinite(r2)]
    if finite_pair_r2.size == 0 or np.any(finite_pair_r2 < 1e-16):
        return 1e100, np.zeros(coords_flat.size, dtype=np.float64)

    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        r2_inv = 1.0 / r2
        r6_inv = r2_inv ** 3
        r12_inv = r6_inv ** 2
        energy = 2.0 * np.sum(r12_inv - r6_inv)
        force_mag = 4.0 * (
            -12.0 * r12_inv * r2_inv + 6.0 * r6_inv * r2_inv
        )
        grad = np.sum(force_mag[:, :, np.newaxis] * diffs, axis=1)

    if not np.isfinite(energy) or not np.all(np.isfinite(grad)):
        return 1e100, np.zeros(coords_flat.size, dtype=np.float64)
    return float(energy), grad.reshape(-1)


def gradient_check(seed: int = 20261107, n_checks: int = 12) -> float:
    """Return max relative error between analytic and central-difference gradients."""
    rng = np.random.RandomState(seed)
    x = rng.uniform(-1.35, 1.35, size=(N_ATOMS, 3)).reshape(-1)
    energy, analytic = lj_energy_and_grad_vectorized(x)
    del energy
    # 3e-5 is a stable central-difference step for LJ forces at these random coordinates.
    eps = 3e-5
    indices = rng.choice(x.size, size=min(n_checks, x.size), replace=False)
    errors = []
    for idx in indices:
        xp = x.copy()
        xm = x.copy()
        xp[idx] += eps
        xm[idx] -= eps
        ep, _ = lj_energy_and_grad_vectorized(xp)
        em, _ = lj_energy_and_grad_vectorized(xm)
        numeric = (ep - em) / (2.0 * eps)
        denom = max(1.0, abs(numeric), abs(float(analytic[idx])))
        errors.append(abs(numeric - float(analytic[idx])) / denom)
    return float(max(errors, default=0.0))


# ---------------------------------------------------------------------------
# Reproducible initial configurations
# ---------------------------------------------------------------------------
def uniform_random_initial(seed: int, n_atoms: int = N_ATOMS) -> np.ndarray:
    """Reproduce np.random.seed(seed); np.random.uniform(-2,2,(N,3))."""
    rng = np.random.RandomState(int(seed))
    return rng.uniform(-2.0, 2.0, size=(n_atoms, 3)).astype(np.float64)


def icosahedral_initial(seed: int, shell_radius: float = ICOSAHEDRON_SHELL_RADIUS,
                        noise_sigma: float = ICOSAHEDRON_NOISE_SIGMA) -> np.ndarray:
    """Central atom plus 12 slightly perturbed normalized icosahedral vertices."""
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    vertices = []

    for a in (-1.0, 1.0):
        for b in (-1.0, 1.0):
            vertices.append((0.0, a, b * phi))
    for a in (-1.0, 1.0):
        for b in (-1.0, 1.0):
            vertices.append((a, b * phi, 0.0))
    for a in (-1.0, 1.0):
        for b in (-1.0, 1.0):
            vertices.append((a * phi, 0.0, b))

    shell = np.asarray(vertices, dtype=np.float64)
    shell /= np.linalg.norm(shell, axis=1).mean()
    shell *= float(shell_radius)

    # Seed-local RNG avoids global RNG state and is fully reproducible.
    rng = np.random.RandomState(int(seed) ^ 0x5A17C023)
    shell += rng.normal(0.0, float(noise_sigma), size=shell.shape)
    coords = np.vstack((np.zeros((1, 3), dtype=np.float64), shell))
    return coords


# ---------------------------------------------------------------------------
# Optimizers
# ---------------------------------------------------------------------------
def optimize_coordinates(initial: np.ndarray, policy: str) -> tuple[np.ndarray, dict[str, Any], float]:
    start = time.perf_counter()
    if policy == "baseline":
        options = {
            "maxiter": BASELINE_MAXITER,
            "gtol": BASELINE_GTOL,
        }
    elif policy == "icosahedral":
        options = {
            "maxiter": POLICY_MAXITER,
            "gtol": POLICY_GTOL,
            "ftol": POLICY_FTOL,
            "maxls": POLICY_MAXLS,
            "maxcor": POLICY_MAXCOR,
        }
    else:
        raise ValueError(f"Unknown optimizer policy: {policy}")

    res = minimize(
        lj_energy_and_grad_vectorized,
        np.asarray(initial, dtype=np.float64).reshape(-1),
        args=(N_ATOMS,),
        method="L-BFGS-B",
        jac=True,
        options=options,
    )
    elapsed = time.perf_counter() - start
    coords = np.asarray(res.x, dtype=np.float64).reshape((N_ATOMS, 3))
    grad_norm = float(np.linalg.norm(res.jac)) if res.jac is not None else float("inf")
    meta = {
        "optimizer": "L-BFGS-B",
        "optimizer_policy": policy,
        "optimizer_success": bool(res.success),
        "optimizer_status": int(res.status),
        "optimizer_message": str(res.message),
        "energy": float(res.fun),
        "grad_norm": grad_norm,
        "nit": int(res.nit),
        "nfev": int(res.nfev),
        "njev": int(res.njev) if getattr(res, "njev", None) is not None else None,
        "optimization_time": float(elapsed),
        "options": options,
    }
    return coords, meta, float(elapsed)


def baseline_optimization(seed: int) -> tuple[np.ndarray, dict[str, Any], float]:
    started = time.perf_counter()
    initial = uniform_random_initial(seed)
    coords, meta, optimizer_elapsed = optimize_coordinates(initial, "baseline")
    elapsed = time.perf_counter() - started
    meta["initialization_time"] = max(0.0, elapsed - optimizer_elapsed)
    return coords, meta, float(elapsed)


def candidate_optimization(seed: int) -> tuple[np.ndarray, dict[str, Any], float, bool]:
    started = time.perf_counter()
    initial = icosahedral_initial(seed)
    coords, meta, optimizer_elapsed = optimize_coordinates(initial, "icosahedral")
    elapsed = time.perf_counter() - started
    meta["initialization_time"] = max(0.0, elapsed - optimizer_elapsed)

    energy_ok = abs(meta["energy"] - KNOWN_LJ13_MINIMUM) <= NUMERICAL_ENERGY_TOL
    gradient_ok = meta["grad_norm"] <= NUMERICAL_GRAD_NORM_MAX
    quality_ok = (
        meta["optimizer_success"]
        and np.all(np.isfinite(coords))
        and math.isfinite(meta["energy"])
        and energy_ok
        and gradient_ok
    )
    meta.update({
        "initialization": "icosahedral_12_shell_plus_center",
        "shell_radius": ICOSAHEDRON_SHELL_RADIUS,
        "noise_sigma": ICOSAHEDRON_NOISE_SIGMA,
        "numerical_quality_gate": bool(quality_ok),
        "energy_gate": bool(energy_ok),
        "gradient_gate": bool(gradient_ok),
    })
    return coords, meta, float(elapsed), bool(quality_ok)


# ---------------------------------------------------------------------------
# XYZ and C++ verifier
# ---------------------------------------------------------------------------
def format_xyz(coords: np.ndarray, comment: str) -> str:
    lines = [str(len(coords)), comment]
    for x, y, z in np.asarray(coords, dtype=np.float64):
        # 12 decimal places preserve the center sufficiently for radius_exp=8.
        lines.append(f"X {x:.12f} {y:.12f} {z:.12f}")
    return "\n".join(lines) + "\n"


def parse_verifier_status(stdout: str, stderr: str, returncode: int) -> str:
    if returncode != 0:
        return f"VERIFIER_ERROR_RC_{returncode}"
    combined = f"{stdout}\n{stderr}".upper()
    negative_markers = (
        "NOT_CERTIFIED", "UNCERTIFIED", "INCONCLUSIVE", "REJECTED",
        "CERTIFIED: FALSE", "CERTIFIED=FALSE", "CERTIFIED 0",
    )
    if any(marker in combined for marker in negative_markers):
        return "INCONCLUSIVE_OR_REJECTED"
    if re.search(r"\bCERTIFIED(?:_STRICT_LOCAL_MINIMUM)?\b", combined):
        return "CERTIFIED_STRICT_LOCAL_MINIMUM"
    return "INCONCLUSIVE_OR_REJECTED"


def run_verifier(xyz_content: str, radius_exp: int = 8,
                 timeout_seconds: float = 30.0) -> tuple[str, float, dict[str, Any]]:
    """Run the existing verifier exactly once with a unique temporary XYZ file."""
    fd, temp_name = tempfile.mkstemp(prefix="cs023_lj13_", suffix=".xyz")
    os.close(fd)
    temp_path = Path(temp_name)
    start = time.perf_counter()
    stdout = ""
    stderr = ""
    returncode = None
    try:
        temp_path.write_text(xyz_content, encoding="utf-8")
        result = subprocess.run(
            [str(VERIFIER_BIN), str(temp_path), str(int(radius_exp))],
            cwd=str(PROJECT_ROOT),
            capture_output=True,
            text=True,
            timeout=float(timeout_seconds),
            check=False,
        )
        stdout = result.stdout or ""
        stderr = result.stderr or ""
        returncode = int(result.returncode)
        status = parse_verifier_status(stdout, stderr, returncode)
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
        stderr = exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or "")
        status = "TIMEOUT"
    except OSError as exc:
        stderr = repr(exc)
        status = "VERIFIER_EXECUTION_ERROR"
    finally:
        elapsed = time.perf_counter() - start
        try:
            temp_path.unlink(missing_ok=True)
        except OSError:
            pass

    info = {
        "returncode": returncode,
        "stdout_tail": stdout[-2000:],
        "stderr_tail": stderr[-2000:],
    }
    return status, float(elapsed), info


# ---------------------------------------------------------------------------
# Evidence Ledger helpers
# ---------------------------------------------------------------------------
def compute_hash(content: str) -> str:
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


def add_ledger_attempt(ledger: Any, xyz: str, status: str,
                       witness_data: dict[str, Any]) -> None:
    ledger.add_attempt(
        task_hash=TASK_HASH,
        config_hash=compute_hash(xyz),
        status=status,
        witness_data=witness_data,
    )


def verify_one_candidate(ledger: Any, coords: np.ndarray, comment: str,
                         radius_exp: int, timeout_seconds: float,
                         witness: dict[str, Any]) -> tuple[str, float, dict[str, Any]]:
    xyz = format_xyz(coords, comment)
    status, elapsed, verifier_info = run_verifier(xyz, radius_exp, timeout_seconds)
    record = dict(witness)
    record.update({
        "verifier_status": status,
        "verification_time": elapsed,
        "radius_exp": int(radius_exp),
        "verifier": str(VERIFIER_BIN),
        "verifier_info": verifier_info,
    })
    add_ledger_attempt(ledger, xyz, status, record)
    return status, elapsed, verifier_info


# ---------------------------------------------------------------------------
# Attempt runners
# ---------------------------------------------------------------------------
def run_baseline_attempt(seed: int, series_index: int, ledger: Any,
                         radius_exp: int, timeout_seconds: float) -> dict[str, Any]:
    coords, opt_meta, opt_time = baseline_optimization(seed)
    xyz = format_xyz(coords, f"Baseline_uniform_seed_{seed}")
    status, verify_time, verifier_info = run_verifier(xyz, radius_exp, timeout_seconds)
    witness = {
        "series_index": series_index,
        "seed": int(seed),
        "policy": "baseline_uniform_random",
        **opt_meta,
        "verification_time": float(verify_time),
        "verifier_status": status,
        "radius_exp": int(radius_exp),
        "verifier_info": verifier_info,
    }
    add_ledger_attempt(ledger, xyz, status, witness)
    return {
        "seed": int(seed),
        "certified": status == "CERTIFIED_STRICT_LOCAL_MINIMUM",
        "status": status,
        "optimization_time": float(opt_time),
        "verification_time": float(verify_time),
        "total_time": float(opt_time + verify_time),
        "verify_calls": 1,
        "fallback_used": False,
        "energy": float(opt_meta["energy"]),
        "grad_norm": float(opt_meta["grad_norm"]),
        "nit": int(opt_meta["nit"]),
    }


def run_candidate_attempt(seed: int, series_index: int, ledger: Any,
                          radius_exp: int, timeout_seconds: float) -> dict[str, Any]:
    """Icosahedral candidate with safe pre-verifier and post-verifier fallback."""
    total_opt_time = 0.0
    total_verify_time = 0.0
    verify_calls = 0
    subattempts: list[dict[str, Any]] = []
    fallback_used = False

    coords, opt_meta, opt_time, quality_ok = candidate_optimization(seed)
    total_opt_time += opt_time

    if not quality_ok:
        # No certificate is inferred from this gate: replace the candidate by
        # the original baseline initialization and send that result to verifier.
        fallback_used = True
        subattempts.append({
            "kind": "icosahedral_quality_gate_fallback",
            "energy": opt_meta["energy"],
            "grad_norm": opt_meta["grad_norm"],
            "optimizer_success": opt_meta["optimizer_success"],
        })
        coords, fallback_meta, fallback_opt_time = baseline_optimization(seed)
        total_opt_time += fallback_opt_time
        opt_meta = fallback_meta
        candidate_comment = f"Fallback_uniform_seed_{seed}"
        policy_meta = {
            "policy": "icosahedral_then_uniform_fallback",
            "initialization": "uniform_fallback_after_numerical_gate",
            "icosahedral_optimizer": opt_meta,
            "fallback_optimizer": fallback_meta,
        }
    else:
        candidate_comment = f"Icosahedral_LJ13_seed_{seed}"
        policy_meta = {
            "policy": "icosahedral_seeded_refinement",
            "initialization": "icosahedral_12_shell_plus_center",
            "optimizer": opt_meta,
        }

    xyz = format_xyz(coords, candidate_comment)
    status, verify_time, verifier_info = run_verifier(xyz, radius_exp, timeout_seconds)
    verify_calls += 1
    total_verify_time += verify_time
    witness = {
        "series_index": series_index,
        "seed": int(seed),
        **policy_meta,
        "verification_time": float(verify_time),
        "verifier_status": status,
        "radius_exp": int(radius_exp),
        "verifier_info": verifier_info,
        "fallback_used_before_verification": bool(fallback_used),
    }
    add_ledger_attempt(ledger, xyz, status, witness)

    # If the structured candidate itself is inconclusive, attempt the original
    # random-start route once.  The verifier remains unchanged and no previous
    # verifier result is used to train or tune any predictor.
    if status != "CERTIFIED_STRICT_LOCAL_MINIMUM" and not fallback_used:
        fallback_used = True
        subattempts.append({
            "kind": "icosahedral_verifier_fallback",
            "first_status": status,
            "first_verifier_time": float(verify_time),
        })
        fallback_coords, fallback_meta, fallback_opt_time = baseline_optimization(seed)
        total_opt_time += fallback_opt_time
        fallback_xyz = format_xyz(fallback_coords, f"Fallback_uniform_after_inconclusive_seed_{seed}")
        fallback_status, fallback_verify_time, fallback_info = run_verifier(
            fallback_xyz, radius_exp, timeout_seconds
        )
        verify_calls += 1
        total_verify_time += fallback_verify_time
        fallback_witness = {
            "series_index": series_index,
            "seed": int(seed),
            "policy": "uniform_fallback_after_icosahedral_inconclusive",
            **fallback_meta,
            "verification_time": float(fallback_verify_time),
            "verifier_status": fallback_status,
            "radius_exp": int(radius_exp),
            "verifier_info": fallback_info,
            "first_candidate_status": status,
        }
        add_ledger_attempt(ledger, fallback_xyz, fallback_status, fallback_witness)
        status = fallback_status
        opt_meta = fallback_meta

    return {
        "seed": int(seed),
        "certified": status == "CERTIFIED_STRICT_LOCAL_MINIMUM",
        "status": status,
        "optimization_time": float(total_opt_time),
        "verification_time": float(total_verify_time),
        "total_time": float(total_opt_time + total_verify_time),
        "verify_calls": int(verify_calls),
        "fallback_used": bool(fallback_used),
        "energy": float(opt_meta["energy"]),
        "grad_norm": float(opt_meta["grad_norm"]),
        "nit": int(opt_meta["nit"]),
        "subattempts": subattempts,
    }


# ---------------------------------------------------------------------------
# Metrics and statistical comparison
# ---------------------------------------------------------------------------
def summarize(records: list[dict[str, Any]]) -> dict[str, Any]:
    n = len(records)
    certified = sum(bool(r["certified"]) for r in records)
    optimization_time = float(sum(r["optimization_time"] for r in records))
    verification_time = float(sum(r["verification_time"] for r in records))
    total_time = optimization_time + verification_time
    cost = total_time / certified if certified else float("inf")
    return {
        "attempts": n,
        "certified_count": int(certified),
        "certificate_rate": certified / n if n else 0.0,
        "optimization_time": optimization_time,
        "verification_time": verification_time,
        "total_time": total_time,
        "cost_per_certificate": cost,
        "verify_calls": int(sum(r["verify_calls"] for r in records)),
        "fallback_count": int(sum(bool(r["fallback_used"]) for r in records)),
    }


def gain_percent(baseline_cost: float, policy_cost: float) -> float:
    if not math.isfinite(baseline_cost) or baseline_cost <= 0 or not math.isfinite(policy_cost):
        return float("nan")
    return 100.0 * (baseline_cost - policy_cost) / baseline_cost


def paired_bootstrap_gain(baseline_records: list[dict[str, Any]],
                          policy_records: list[dict[str, Any]],
                          n_boot: int = 10000,
                          seed: int = 20261107) -> tuple[float, float, float]:
    """Paired attempt-level bootstrap CI for gain in total-time/certificate."""
    if len(baseline_records) != len(policy_records) or not baseline_records:
        return float("nan"), float("nan"), float("nan")

    b_time = np.asarray([r["total_time"] for r in baseline_records], dtype=np.float64)
    p_time = np.asarray([r["total_time"] for r in policy_records], dtype=np.float64)
    b_cert = np.asarray([r["certified"] for r in baseline_records], dtype=np.int64)
    p_cert = np.asarray([r["certified"] for r in policy_records], dtype=np.int64)

    rng = np.random.RandomState(seed)
    n = len(baseline_records)
    gains = []
    for _ in range(n_boot):
        idx = rng.randint(0, n, size=n)
        bc = int(b_cert[idx].sum())
        pc = int(p_cert[idx].sum())
        if bc == 0 or pc == 0:
            continue
        b_cost = float(b_time[idx].sum()) / bc
        p_cost = float(p_time[idx].sum()) / pc
        gains.append(gain_percent(b_cost, p_cost))

    if not gains:
        return float("nan"), float("nan"), float("nan")
    point = gain_percent(float(b_time.sum()) / int(b_cert.sum()),
                         float(p_time.sum()) / int(p_cert.sum()))
    low, high = np.percentile(np.asarray(gains), [2.5, 97.5])
    return float(point), float(low), float(high)


def print_summary_block(name: str, baseline_records: list[dict[str, Any]],
                        policy_records: list[dict[str, Any]]) -> dict[str, Any]:
    b = summarize(baseline_records)
    p = summarize(policy_records)
    gain = gain_percent(b["cost_per_certificate"], p["cost_per_certificate"])
    print(f"\n{name}")
    print("-" * 104)
    print(f"{'Policy':<35} {'Certs':>10} {'Opt(s)':>12} {'Verify(s)':>12} {'Total(s)':>12} {'s/cert':>12}")
    print("-" * 104)
    for label, row in (("Baseline (uniform random)", b), ("Metareasoning (icosahedral)", p)):
        cpc = row["cost_per_certificate"]
        cpc_text = f"{cpc:.6f}" if math.isfinite(cpc) else "N/A"
        print(f"{label:<35} {row['certified_count']:>4}/{row['attempts']:<5} "
              f"{row['optimization_time']:>12.4f} {row['verification_time']:>12.4f} "
              f"{row['total_time']:>12.4f} {cpc_text:>12}")
    print("-" * 104)
    print(f"Metareasoning gain: {gain:.2f}% (target >= {TARGET_GAIN_PERCENT:.1f}%)")
    print(f"Fallbacks: {p['fallback_count']}/{p['attempts']}; verifier calls: {p['verify_calls']}")
    return {"baseline": b, "policy": p, "gain_percent": gain}


# ---------------------------------------------------------------------------
# Self-test without verifier
# ---------------------------------------------------------------------------
def self_test() -> int:
    print("CS-023/024 optimizer self-test (verifier is not invoked)")
    rel_error = gradient_check()
    print(f"Analytic gradient: {'PASS' if rel_error < 1e-5 else 'FAIL'} "
          f"(max relative error={rel_error:.3e})")
    if rel_error >= 1e-5:
        return 2

    energies = []
    grad_norms = []
    iterations = []
    successes = 0
    start = time.perf_counter()
    for seed in range(1000, 1100):
        _, meta, _, quality_ok = candidate_optimization(seed)
        energies.append(meta["energy"])
        grad_norms.append(meta["grad_norm"])
        iterations.append(meta["nit"])
        successes += int(quality_ok)
    elapsed = time.perf_counter() - start
    print(f"Structured candidates passing numerical gate: {successes}/100")
    print(f"Energy: min={min(energies):.10f}, median={np.median(energies):.10f}, "
          f"max={max(energies):.10f}")
    print(f"Gradient norm: median={np.median(grad_norms):.3e}, "
          f"max={max(grad_norms):.3e}")
    print(f"Iterations: median={np.median(iterations):.1f}, max={max(iterations)}")
    print(f"Optimizer-only elapsed: {elapsed:.4f} s")
    print("NOTE: this self-test does not establish C++ certification.")
    return 0 if successes == 100 else 1


# ---------------------------------------------------------------------------
# Main benchmark
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="CS-023/024 LJ13 metareasoning benchmark")
    parser.add_argument("--runs", type=int, default=100,
                        help="attempts per independent series (default: 100)")
    parser.add_argument("--series", type=int, default=3,
                        help="number of independent seed series (default: 3)")
    parser.add_argument("--base-seed", type=int, default=1000,
                        help="first seed; default first series is 1000..1099")
    parser.add_argument("--radius-exp", type=int, default=8,
                        help="verifier radius exponent (default: 8; do not weaken for acceptance)")
    parser.add_argument("--timeout", type=float, default=30.0,
                        help="timeout per verifier call in seconds")
    parser.add_argument("--output", type=Path, default=RESULTS_PATH,
                        help="JSON report path")
    parser.add_argument("--self-test", action="store_true",
                        help="test gradient and structured optimization without verifier")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.self_test:
        return self_test()

    if args.runs <= 0 or args.series <= 0:
        print("ERROR: --runs and --series must be positive", file=sys.stderr)
        return 2
    if not 1 <= args.radius_exp <= 12:
        print("ERROR: --radius-exp must be between 1 and 12", file=sys.stderr)
        return 2
    if args.timeout <= 0:
        print("ERROR: --timeout must be positive", file=sys.stderr)
        return 2

    # Ensure EvidenceLedger's relative path resolves to PROJECT_ROOT/ledger/evidence.db.
    os.chdir(PROJECT_ROOT)
    sys.path.insert(0, str(PROJECT_ROOT))
    try:
        from ledger.ledger_db import EvidenceLedger
    except ImportError as exc:
        print(f"ERROR: cannot import EvidenceLedger from {PROJECT_ROOT / 'ledger'}: {exc}", file=sys.stderr)
        return 2

    if not VERIFIER_BIN.exists() or not os.access(VERIFIER_BIN, os.X_OK):
        print(f"ERROR: verifier missing or not executable: {VERIFIER_BIN}", file=sys.stderr)
        return 2

    rel_error = gradient_check()
    print("CS-023/024 Metareasoning Solution — LJ13")
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Verifier: {VERIFIER_BIN}")
    print(f"SciPy: {scipy.__version__}; NumPy: {np.__version__}")
    print(f"Series: {args.series} x {args.runs}; base seed={args.base_seed}")
    print(f"First series seeds: {args.base_seed}..{args.base_seed + args.runs - 1}")
    print(f"Verifier radius exponent: {args.radius_exp}")
    print("Policy: perturbed icosahedral LJ13 initialization + tight L-BFGS-B")
    print("Fallback: original uniform random initialization if numerical gate or verifier fails")
    print(f"Analytic gradient check: {'PASS' if rel_error < 1e-5 else 'FAIL'} "
          f"(max relative error={rel_error:.3e})")
    print("No verifier-based predictor training; every accepted result is verified.")

    if rel_error >= 1e-5:
        print("ERROR: analytic gradient check failed; refusing to run benchmark.", file=sys.stderr)
        return 2

    ledger = EvidenceLedger()
    all_baseline: list[dict[str, Any]] = []
    all_policy: list[dict[str, Any]] = []
    per_series: list[dict[str, Any]] = []
    run_started = time.perf_counter()

    for series_idx in range(args.series):
        first_seed = args.base_seed + series_idx * args.runs
        seeds = list(range(first_seed, first_seed + args.runs))
        print(f"\nSERIES {series_idx + 1}/{args.series}: seeds {seeds[0]}..{seeds[-1]}")

        baseline_records: list[dict[str, Any]] = []
        policy_records: list[dict[str, Any]] = []
        for offset, seed in enumerate(seeds):
            # Alternate execution order to reduce systematic thermal/order bias.
            if (offset + series_idx) % 2 == 0:
                b = run_baseline_attempt(seed, series_idx + 1, ledger,
                                         args.radius_exp, args.timeout)
                p = run_candidate_attempt(seed, series_idx + 1, ledger,
                                          args.radius_exp, args.timeout)
            else:
                p = run_candidate_attempt(seed, series_idx + 1, ledger,
                                          args.radius_exp, args.timeout)
                b = run_baseline_attempt(seed, series_idx + 1, ledger,
                                         args.radius_exp, args.timeout)
            baseline_records.append(b)
            policy_records.append(p)

            done = offset + 1
            if done % 10 == 0 or done == args.runs:
                print(f"  completed paired seeds: {done}/{args.runs}")

        summary = print_summary_block(
            f"SERIES {series_idx + 1} RESULT (seeds {seeds[0]}..{seeds[-1]})",
            baseline_records,
            policy_records,
        )
        summary["series_index"] = series_idx + 1
        summary["seed_start"] = seeds[0]
        summary["seed_end"] = seeds[-1]
        summary["baseline_records"] = baseline_records
        summary["policy_records"] = policy_records
        per_series.append(summary)
        all_baseline.extend(baseline_records)
        all_policy.extend(policy_records)

    overall = print_summary_block("AGGREGATE RESULT", all_baseline, all_policy)
    point_gain, ci_low, ci_high = paired_bootstrap_gain(
        all_baseline, all_policy, seed=20261107 + args.base_seed
    )
    print("\nSTATISTICAL CHECK")
    print(f"Paired bootstrap gain: {point_gain:.2f}% "
          f"(95% CI: {ci_low:.2f}% .. {ci_high:.2f}%)")

    per_series_pass = []
    for row in per_series:
        b = row["baseline"]
        p = row["policy"]
        g = row["gain_percent"]
        per_series_pass.append(
            p["certified_count"] >= MIN_CERTIFICATES_PER_SERIES
            and p["cost_per_certificate"] <= TARGET_COST_PER_CERTIFICATE
            and g >= TARGET_GAIN_PERCENT
        )

    overall_policy = overall["policy"]
    overall_pass = (
        overall_policy["certified_count"] >= MIN_CERTIFICATES_PER_SERIES * args.series
        and overall_policy["cost_per_certificate"] <= TARGET_COST_PER_CERTIFICATE
        and overall["gain_percent"] >= TARGET_GAIN_PERCENT
    )
    statistical_positive = math.isfinite(ci_low) and ci_low > 0.0
    all_series_pass = all(per_series_pass)

    print("\n" + "=" * 88)
    print("CS-023/024 ACCEPTANCE")
    print("=" * 88)
    print(f"Certificates per series >= {MIN_CERTIFICATES_PER_SERIES}: "
          f"{'PASS' if all(s['policy']['certified_count'] >= MIN_CERTIFICATES_PER_SERIES for s in per_series) else 'FAIL'}")
    print(f"Cost/certificate <= {TARGET_COST_PER_CERTIFICATE:.4f} s (aggregate): "
          f"{'PASS' if overall_policy['cost_per_certificate'] <= TARGET_COST_PER_CERTIFICATE else 'FAIL'}")
    print(f"Gain >= {TARGET_GAIN_PERCENT:.1f}% (aggregate): "
          f"{'PASS' if overall['gain_percent'] >= TARGET_GAIN_PERCENT else 'FAIL'}")
    print(f"All {args.series} series independently pass all empirical limits: "
          f"{'PASS' if all_series_pass else 'FAIL'}")
    print(f"Bootstrap 95% CI excludes zero: {'PASS' if statistical_positive else 'FAIL'}")
    print(f"Bootstrap lower CI bound >= {TARGET_GAIN_PERCENT:.1f}%: "
          f"{'PASS' if math.isfinite(ci_low) and ci_low >= TARGET_GAIN_PERCENT else 'FAIL'}")

    stage_pass = overall_pass and all_series_pass and statistical_positive
    print(f"RESULT: {'PASS' if stage_pass else 'NOT YET'}")
    if stage_pass and ci_low < TARGET_GAIN_PERCENT:
        print("NOTE: empirical target passed, but the bootstrap CI does not establish a >=20% lower bound.")

    try:
        ledger.verify_integrity()
        print("Evidence Ledger integrity check completed.")
    except Exception as exc:
        print(f"ERROR: Evidence Ledger integrity check failed: {exc}", file=sys.stderr)
        return 3

    report = {
        "schema": "cs023_024_metareasoning_solution_v1",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "project_root": str(PROJECT_ROOT),
        "verifier": str(VERIFIER_BIN),
        "radius_exp": int(args.radius_exp),
        "num_series": int(args.series),
        "runs_per_series": int(args.runs),
        "base_seed": int(args.base_seed),
        "seed_ranges": [
            [args.base_seed + i * args.runs,
             args.base_seed + i * args.runs + args.runs - 1]
            for i in range(args.series)
        ],
        "baseline_parameters": {
            "initialization": "uniform(-2,2)",
            "optimizer": "L-BFGS-B",
            "maxiter": BASELINE_MAXITER,
            "gtol": BASELINE_GTOL,
        },
        "policy_parameters": {
            "initialization": "central_atom_plus_12_icosahedral_shell_vertices",
            "shell_radius": ICOSAHEDRON_SHELL_RADIUS,
            "noise_sigma": ICOSAHEDRON_NOISE_SIGMA,
            "optimizer": "L-BFGS-B",
            "maxiter": POLICY_MAXITER,
            "gtol": POLICY_GTOL,
            "ftol": POLICY_FTOL,
            "maxls": POLICY_MAXLS,
            "maxcor": POLICY_MAXCOR,
            "fallback": "uniform baseline if numerical gate or first verifier result fails",
            "numerical_gate_is_not_a_certificate": True,
        },
        "overall": overall,
        "paired_bootstrap_gain_percent": point_gain,
        "paired_bootstrap_95ci": [ci_low, ci_high],
        "per_series_empirical_pass": per_series_pass,
        "stage_acceptance_pass": bool(stage_pass),
        "wall_elapsed_for_both_policies": time.perf_counter() - run_started,
        "series": per_series,
    }
    out_path = Path(args.output)
    if not out_path.is_absolute():
        out_path = PROJECT_ROOT / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    def json_safe(value: Any):
        if isinstance(value, dict):
            return {str(k): json_safe(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [json_safe(v) for v in value]
        if isinstance(value, (bool, np.bool_)):
            return bool(value)
        if isinstance(value, (float, np.floating)):
            return float(value) if math.isfinite(float(value)) else None
        if isinstance(value, (int, np.integer)):
            return int(value)
        return value

    out_path.write_text(json.dumps(json_safe(report), indent=2, allow_nan=False), encoding="utf-8")
    print(f"Report saved: {out_path}")
    return 0 if stage_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
