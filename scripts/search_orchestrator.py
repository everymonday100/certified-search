#!/usr/bin/env python3
"""
Search Orchestrator with Adaptive Metareasoning Policy (CS-023/024).
Трёхуровневая стратегия: немедленная верификация / адаптивная дооптимизация / пропуск.
Основано на диагностическом пороге grad_norm = 2.5e-03.
"""
import sys
import time
import json
import hashlib
import subprocess
import numpy as np
from pathlib import Path

try:
    from scipy.optimize import minimize
except ImportError:
    print("❌ Требуется scipy. Установите: python3 -m pip install scipy")
    sys.exit(1)

sys.path.append(str(Path(__file__).parent.parent))
from ledger.ledger_db import EvidenceLedger

VERIFIER_BIN = "./core/verifier/build/krawczyk_lj13"
LEDGER = EvidenceLedger()

def compute_hash(content: str) -> str:
    return hashlib.sha256(content.encode('utf-8')).hexdigest()

def lj_energy_flat(coords_flat, n_atoms=13):
    coords = coords_flat.reshape((n_atoms, 3))
    e = 0.0
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            dx = coords[i, 0] - coords[j, 0]
            dy = coords[i, 1] - coords[j, 1]
            dz = coords[i, 2] - coords[j, 2]
            r2 = dx*dx + dy*dy + dz*dz
            if r2 < 1e-12:
                return 1e10
            r6 = r2**-3
            r12 = r6**2
            e += 4.0 * (r12 - r6)
    return e

def optimize_candidate(initial_coords, n_atoms=13, gtol=1e-6, maxiter=500):
    res = minimize(
        lj_energy_flat, 
        initial_coords.flatten(), 
        args=(n_atoms,), 
        method='L-BFGS-B',
        options={'maxiter': maxiter, 'gtol': gtol}
    )
    grad_norm = np.linalg.norm(res.jac) if res.jac is not None else 1e10
    return res.x.reshape((n_atoms, 3)), res.success, res.fun, grad_norm, res.nit

def format_xyz(coords, strategy_name):
    n = len(coords)
    lines = [f"{n}", f"{strategy_name}"]
    for c in coords:
        lines.append(f"X {c[0]:.10f} {c[1]:.10f} {c[2]:.10f}")
    return "\n".join(lines)

def run_verifier(xyz_content: str, radius_exp: str = "8") -> tuple[str, float]:
    temp_file = Path("/tmp/orchestrator_test.xyz")
    temp_file.write_text(xyz_content)
    start_time = time.perf_counter()
    try:
        result = subprocess.run(
            [VERIFIER_BIN, str(temp_file), radius_exp],
            capture_output=True, text=True, timeout=30
        )
        elapsed = time.perf_counter() - start_time
        if "CERTIFIED" in result.stdout:
            return "CERTIFIED_STRICT_LOCAL_MINIMUM", elapsed
        else:
            return "INCONCLUSIVE_OR_REJECTED", elapsed
    except subprocess.TimeoutExpired:
        return "TIMEOUT", 30.0
    finally:
        if temp_file.exists():
            temp_file.unlink()

# --- Политики ---

def baseline_strategy(n_atoms=13, seed=42):
    """Baseline: всегда оптимизируем с gtol=1e-6 и ВСЕГДА вызываем верификатор."""
    np.random.seed(seed)
    initial = np.random.uniform(-2.0, 2.0, (n_atoms, 3))
    opt_coords, success, energy, grad_norm, nit = optimize_candidate(initial, n_atoms, gtol=1e-6)
    
    meta = {
        "grad_norm": float(grad_norm), 
        "nit": int(nit), 
        "energy": float(energy), 
        "policy": "baseline_always_verify",
        "refine_passes": 1
    }
    return format_xyz(opt_coords, "Baseline"), meta

def metareasoning_strategy(n_atoms=13, seed=42):
    """
    Адаптивная политика на основе порога 2.5e-03:
    - grad_norm < 2.5e-03: немедленная верификация
    - 2.5e-03 <= grad_norm < 1.0e-02: дооптимизация с gtol=1e-8, затем верификация
    - grad_norm >= 1.0e-02: пропуск (заведомый провал)
    """
    np.random.seed(seed)
    initial = np.random.uniform(-2.0, 2.0, (n_atoms, 3))
    
    # Первый проход оптимизации
    opt_coords, success, energy, grad_norm, nit = optimize_candidate(initial, n_atoms, gtol=1e-6)
    
    meta = {
        "grad_norm_initial": float(grad_norm),
        "nit_initial": int(nit),
        "energy_initial": float(energy),
        "refine_passes": 1
    }
    
    # ==========================================
    # ТРЁХУРОВНЕВОЕ РЕШЕНИЕ (КРИТИЧЕСКИ ВАЖНО)
    # ==========================================
    if grad_norm < 2.5e-03:
        # Уровень 1: высокая вероятность успеха, сразу верифицируем
        meta["decision"] = "verify_immediate"
        meta["grad_norm_final"] = float(grad_norm)
        return format_xyz(opt_coords, "Metareasoning_Immediate"), meta
    
    elif grad_norm < 1.0e-02:
        # Уровень 2: дооптимизация с более строгим gtol
        opt_coords2, success2, energy2, grad_norm2, nit2 = optimize_candidate(
            opt_coords, n_atoms, gtol=1e-8, maxiter=2000
        )
        
        meta["decision"] = "refine_then_verify"
        meta["grad_norm_final"] = float(grad_norm2)
        meta["nit_total"] = int(nit + nit2)
        meta["energy_final"] = float(energy2)
        meta["refine_passes"] = 2
        
        return format_xyz(opt_coords2, "Metareasoning_Refined"), meta
    
    else:
        # Уровень 3: заведомый провал, пропускаем
        meta["decision"] = "skip"
        meta["grad_norm_final"] = float(grad_norm)
        return None, meta

def run_experiment(policy_name: str, policy_func, num_runs: int = 100):
    print(f"\n🔬 Запуск политики: {policy_name} ({num_runs} попыток)")
    
    certified_count = 0
    total_opt_time = 0.0
    total_verify_time = 0.0
    verify_calls = 0
    skip_count = 0
    refine_count = 0

    for i in range(num_runs):
        start_opt = time.perf_counter()
        xyz_content, meta = policy_func(seed=1000 + i)
        total_opt_time += (time.perf_counter() - start_opt)
        
        config_hash = compute_hash(xyz_content) if xyz_content else f"skipped_{i}"
        task_hash = "sha256:lj13_metareasoning_task"
        
        if xyz_content is None:
            # Metareasoning решил пропустить
            status = "SKIPPED_BY_POLICY"
            elapsed_verify = 0.0
            skip_count += 1
        else:
            verify_calls += 1
            if meta.get("refine_passes", 1) > 1:
                refine_count += 1
            
            start_verify = time.perf_counter()
            status, elapsed_verify = run_verifier(xyz_content)
            total_verify_time += (time.perf_counter() - start_verify)
            
            if status == "CERTIFIED_STRICT_LOCAL_MINIMUM":
                certified_count += 1

        LEDGER.add_attempt(
            task_hash=task_hash,
            config_hash=config_hash,
            status=status,
            witness_data={**meta, "elapsed_verify": round(elapsed_verify, 4)}
        )

    success_rate = certified_count / num_runs
    avg_total_time = (total_opt_time + total_verify_time) / num_runs
    
    print(f"  ✅ Сертификаций: {certified_count}/{num_runs} ({success_rate:.1%})")
    print(f"  🔍 Вызовов верификатора: {verify_calls}/{num_runs} ({verify_calls/num_runs:.1%})")
    print(f"  ⏭️  Пропущено: {skip_count}/{num_runs} ({skip_count/num_runs:.1%})")
    print(f"  🔧 Дооптимизировано: {refine_count}/{num_runs} ({refine_count/num_runs:.1%})")
    print(f"  ⏱️  Ср. общее время на попытку: {avg_total_time:.4f} сек")
    
    # Стоимость одной успешной сертификации
    if certified_count > 0:
        cost_per_cert = (total_opt_time + total_verify_time) / certified_count
        print(f"  💰 Стоимость одной сертификации: {cost_per_cert:.4f} сек")
    else:
        cost_per_cert = float('inf')
        print(f"  💰 Стоимость одной сертификации: N/A (0 успехов)")
    
    return {
        "policy": policy_name,
        "success_rate": success_rate,
        "verify_call_rate": verify_calls / num_runs,
        "skip_rate": skip_count / num_runs,
        "refine_rate": refine_count / num_runs,
        "cost_per_cert": cost_per_cert,
        "certified_count": certified_count,
        "total_time": total_opt_time + total_verify_time
    }

def main():
    print("🚀 Запуск экспериментов Metareasoning (Этап 2, CS-023/024)")
    print("Цель: снизить стоимость сертификации на ≥20% по сравнению с baseline")
    print("Порог отсечения: 2.5e-03 (на основе диагностики)")
    
    if not Path(VERIFIER_BIN).exists():
        print("❌ Бинарник верификатора не найден.")
        sys.exit(1)

    NUM_RUNS = 100
    
    results = []
    results.append(run_experiment("Baseline (Always Verify)", baseline_strategy, num_runs=NUM_RUNS))
    results.append(run_experiment("Metareasoning (Adaptive)", metareasoning_strategy, num_runs=NUM_RUNS))

    print("\n" + "="*80)
    print("📊 СВОДКА СРАВНЕНИЯ ПОЛИТИК (BASELINE vs METAREASONING)")
    print("="*80)
    
    baseline = results[0]
    meta = results[1]
    
    print(f"{'Политика':<30} | {'Успех':<10} | {'Вериф.':<8} | {'Пропуск':<8} | {'Доопт.':<8} | {'Стоимость/сертификат':<20}")
    print("-" * 90)
    print(f"{baseline['policy']:<30} | {baseline['certified_count']}/{NUM_RUNS:<5} | {baseline['verify_call_rate']:.1%} | {baseline['skip_rate']:.1%} | {baseline.get('refine_rate', 0):.1%} | {baseline['cost_per_cert']:.4f} сек")
    print(f"{meta['policy']:<30} | {meta['certified_count']}/{NUM_RUNS:<5} | {meta['verify_call_rate']:.1%} | {meta['skip_rate']:.1%} | {meta['refine_rate']:.1%} | {meta['cost_per_cert']:.4f} сек")
    
    # Расчёт улучшения
    if baseline['cost_per_cert'] > 0 and meta['cost_per_cert'] != float('inf'):
        improvement = (baseline['cost_per_cert'] - meta['cost_per_cert']) / baseline['cost_per_cert'] * 100
        print("-" * 90)
        print(f"📈 Улучшение стоимости (Metareasoning Gain): {improvement:.1f}%")
        
        if improvement >= 20.0:
            print("🎉 УСПЕХ: Цель Этапа 2 достигнута (≥20% экономии)!")
        else:
            print("⚠️  Цель не достигнута. Требуется доработка политики.")
    elif meta['certified_count'] == 0:
        print("-" * 90)
        print("❌ Metareasoning не дал ни одной сертификации. Политика слишком агрессивна.")

    print("\n✅ Данные записаны в Evidence Ledger.")
    LEDGER.verify_integrity()

if __name__ == "__main__":
    main()