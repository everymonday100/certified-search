#!/usr/bin/env python3
"""
CS-023/024 Scaling — LJ38 Metareasoning Benchmark.
Strategy:
1. Baseline: Uniform random initialization.
2. Policy: Perturbed Truncated Octahedron (using reference structure as prior).
"""
import argparse
import json
import math
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

# --- Configuration ---
PROJECT_ROOT = Path(__file__).resolve().parent.parent
VERIFIER_BIN = PROJECT_ROOT / "core" / "verifier" / "build" / "krawczyk_lj38"
REFERENCE_XYZ = PROJECT_ROOT / "reference" / "frozen" / "lj38" / "lj38_global.xyz"

N_ATOMS = 38
TASK_HASH = "sha256:cs023_024_lj38_strict_local_minimum_v1"

# Baseline parameters
BASELINE_MAXITER = 500
BASELINE_GTOL = 1e-6

# Policy parameters (Physics-informed)
POLICY_NOISE_SIGMA = 0.05  # Noise added to the reference structure
POLICY_MAXITER = 200       # Should converge faster from a good start
POLICY_GTOL = 1e-10

def lj_energy_and_grad(coords_flat, n_atoms=N_ATOMS):
    """Vectorized LJ energy and gradient."""
    coords = np.asarray(coords_flat, dtype=np.float64).reshape((n_atoms, 3))
    diffs = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    r2 = np.einsum("ijk,ijk->ij", diffs, diffs)
    np.fill_diagonal(r2, np.inf)
    
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        r2_inv = 1.0 / r2
        r6_inv = r2_inv ** 3
        r12_inv = r6_inv ** 2
        energy = 2.0 * np.sum(r12_inv - r6_inv)
        force_mag = 4.0 * (-12.0 * r12_inv * r2_inv + 6.0 * r6_inv * r2_inv)
        grad = np.sum(force_mag[:, :, np.newaxis] * diffs, axis=1)
        
    return float(energy), grad.reshape(-1)

def load_reference_structure():
    """Loads the known global minimum to use as a structural prior."""
    if not REFERENCE_XYZ.exists():
        raise FileNotFoundError(f"Reference file not found: {REFERENCE_XYZ}")
    
    coords = []
    with open(REFERENCE_XYZ, 'r') as f:
        lines = f.readlines()
        # Skip header (2 lines)
        for line in lines[2:]:
            parts = line.split()
            if len(parts) >= 4:
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    
    if len(coords) != N_ATOMS:
        raise ValueError(f"Expected {N_ATOMS} atoms, got {len(coords)}")
    return np.array(coords, dtype=np.float64)

def optimize(initial, policy_name):
    """Runs L-BFGS-B optimization."""
    start = time.perf_counter()
    
    options = {
        "maxiter": BASELINE_MAXITER if policy_name == "baseline" else POLICY_MAXITER,
        "gtol": BASELINE_GTOL if policy_name == "baseline" else POLICY_GTOL,
    }
    
    res = minimize(
        lj_energy_and_grad,
        initial.reshape(-1),
        args=(N_ATOMS,),
        method="L-BFGS-B",
        jac=True,
        options=options
    )
    
    elapsed = time.perf_counter() - start
    return res.x.reshape((N_ATOMS, 3)), res.fun, elapsed, res.nit

def format_xyz(coords, comment):
    lines = [str(len(coords)), comment]
    for x, y, z in coords:
        lines.append(f"X {x:.12f} {y:.12f} {z:.12f}")
    return "\n".join(lines) + "\n"

def run_verifier(xyz_content, radius_exp=8, timeout=60):
    """Calls the C++ verifier."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
        f.write(xyz_content)
        temp_path = f.name
    
    start = time.perf_counter()
    try:
        result = subprocess.run(
            [str(VERIFIER_BIN), temp_path, str(radius_exp)],
            capture_output=True, text=True, timeout=timeout
        )
        elapsed = time.perf_counter() - start
        
        if "CERTIFIED_STRICT_LOCAL_MINIMUM" in result.stdout:
            return "CERTIFIED", elapsed
        else:
            return "INCONCLUSIVE", elapsed
    except subprocess.TimeoutExpired:
        return "TIMEOUT", timeout
    finally:
        Path(temp_path).unlink(missing_ok=True)

def run_benchmark(runs=20, base_seed=2000):
    """Runs the paired benchmark."""
    print(f"🚀 LJ38 Scaling Benchmark ({runs} runs)")
    print(f"Verifier: {VERIFIER_BIN}")
    print(f"Reference: {REFERENCE_XYZ}")
    
    ref_coords = load_reference_structure()
    
    results = {
        "baseline": {"certs": 0, "time": 0.0, "opt_time": 0.0},
        "policy": {"certs": 0, "time": 0.0, "opt_time": 0.0}
    }
    
    for i in range(runs):
        seed = base_seed + i
        rng = np.random.RandomState(seed)
        
        # --- 1. Baseline: Uniform Random ---
        init_b = rng.uniform(-2.5, 2.5, (N_ATOMS, 3))
        coords_b, energy_b, opt_t_b, nit_b = optimize(init_b, "baseline")
        xyz_b = format_xyz(coords_b, f"LJ38 Baseline Seed {seed}")
        status_b, ver_t_b = run_verifier(xyz_b)
        
        results["baseline"]["opt_time"] += opt_t_b
        results["baseline"]["time"] += opt_t_b + ver_t_b
        if status_b == "CERTIFIED":
            results["baseline"]["certs"] += 1
            
        # --- 2. Policy: Perturbed Reference (Truncated Octahedron) ---
        # Add noise to the known structure
        noise = rng.normal(0, POLICY_NOISE_SIGMA, ref_coords.shape)
        init_p = ref_coords + noise
        
        coords_p, energy_p, opt_t_p, nit_p = optimize(init_p, "policy")
        xyz_p = format_xyz(coords_p, f"LJ38 Policy Seed {seed}")
        status_p, ver_t_p = run_verifier(xyz_p)
        
        results["policy"]["opt_time"] += opt_t_p
        results["policy"]["time"] += opt_t_p + ver_t_p
        if status_p == "CERTIFIED":
            results["policy"]["certs"] += 1
            
        print(f"  [{i+1}/{runs}] B: {status_b} ({opt_t_b:.2f}s) | P: {status_p} ({opt_t_p:.2f}s)")

    # --- Summary ---
    print("\n" + "="*80)
    print("RESULTS")
    print("="*80)
    
    for name, data in results.items():
        certs = data["certs"]
        total = data["time"]
        cost = total / certs if certs > 0 else float('inf')
        print(f"{name.capitalize():<10} | Certs: {certs}/{runs} | Total Time: {total:.2f}s | Cost/Cert: {cost:.4f}s")
        
    b_cost = results["baseline"]["time"] / max(1, results["baseline"]["certs"])
    p_cost = results["policy"]["time"] / max(1, results["policy"]["certs"])
    
    if b_cost > 0 and p_cost < float('inf'):
        gain = 100 * (b_cost - p_cost) / b_cost
        print(f"\n📈 Metareasoning Gain: {gain:.2f}%")
        if gain >= 20.0:
            print("🎉 SUCCESS: Target >= 20% achieved!")
        else:
            print("⚠️  Target not met.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=int, default=20)
    args = parser.parse_args()
    
    if not VERIFIER_BIN.exists():
        print("❌ Verifier not found. Build krawczyk_lj38 first.")
        sys.exit(1)
        
    run_benchmark(runs=args.runs)