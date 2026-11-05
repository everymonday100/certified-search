#!/usr/bin/env python3
"""
Red-Team Fuzzer for the Verifier (CS-018).
Генерирует патологические входные данные и прогоняет 1000+ заранее созданных негативных тестов.
"""
import sys
import subprocess
import hashlib
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent))
from ledger.ledger_db import EvidenceLedger

VERIFIER_BIN = "./core/verifier/build/krawczyk_lj13"
BASE_XYZ = Path("reference/frozen/lj13/lj13.xyz")
NEG_DIR = Path("reference/negative_tests")
LEDGER = EvidenceLedger()

def compute_hash(content: str) -> str:
    return hashlib.sha256(content.encode('utf-8')).hexdigest()

def run_verifier(file_path: str, radius_exp: str = "8") -> tuple[int, str, str]:
    try:
        result = subprocess.run(
            [VERIFIER_BIN, file_path, radius_exp],
            capture_output=True, text=True, timeout=10
        )
        return result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return -1, "", "TIMEOUT"
    except Exception as e:
        return -1, "", str(e)

def main():
    print("🚀 Запуск Red-Team Fuzzing (CS-018) на 1000+ тестов...")
    
    if not Path(VERIFIER_BIN).exists():
        print("❌ Бинарник верификатора не найден.")
        sys.exit(1)

    false_positives = 0
    crashes = 0
    total_tests = 0

    # 1. Прогон сгенерированных негативных файлов (1000 шт)
    print("\n[1/2] Прогон корпуса reference/negative_tests/...")
    if NEG_DIR.exists():
        for xyz_file in sorted(NEG_DIR.glob("*.xyz")):
            total_tests += 1
            content = xyz_file.read_text()
            config_hash = compute_hash(content)
            task_hash = "sha256:fuzzed_task_spec"
            
            returncode, stdout, stderr = run_verifier(str(xyz_file))
            
            if "CERTIFIED" in stdout:
                false_positives += 1
                print(f"❌ FALSE POSITIVE: {xyz_file.name}")
                status = "FALSE_POSITIVE_CERTIFIED"
            elif returncode < 0 or "TIMEOUT" in stderr:
                crashes += 1
                status = "CRASH_OR_TIMEOUT"
            else:
                status = "CORRECTLY_REJECTED"
                
            LEDGER.add_attempt(task_hash, config_hash, status, {"file": xyz_file.name, "rc": returncode})
    else:
        print("⚠️  Папка negative_tests не найдена. Запустите generate_negative_tests.py")

    print(f"✅ Прогон корпуса завершен. Обработано: {total_tests}")

    # 2. Финальная проверка
    print("\n" + "="*50)
    print("📊 СВОДКА ФАЗЗИНГА (PHASE 1 EXIT CRITERIA)")
    print("="*50)
    print(f"Всего протестировано случаев: {total_tests}")
    print(f"Ложные сертификации (CERTIFIED): {false_positives} (ЦЕЛЬ: 0)")
    print(f"Неперехваченные краши/таймауты: {crashes}")
    
    if false_positives == 0 and total_tests >= 1000:
        print("🎉 УСПЕХ: Критерий Phase 1 выполнен! (>=1000 тестов, 0 ложных ACCEPT)")
    else:
        print("🚨 ПРОВАЛ: Критерии не выполнены.")
        sys.exit(1)

    print("\nПроверка целостности Ledger...")
    LEDGER.verify_integrity()

if __name__ == "__main__":
    main()