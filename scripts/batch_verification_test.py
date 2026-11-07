#!/usr/bin/env python3
"""
Batch Verification Test (CS-029b).
Tests the OpenMP-parallelized batch verifier on LJ38.
"""
import json
import os
import subprocess
import tempfile
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BATCH_BIN = PROJECT_ROOT / "core/verifier/build/krawczyk_lj38_batch"
REFERENCE_DIR = PROJECT_ROOT / "reference/frozen/lj38"

def main():
    if not BATCH_BIN.exists():
        print(f"❌ Batch verifier not found: {BATCH_BIN}")
        return

    # Create a temporary directory with multiple XYZ files
    with tempfile.TemporaryDirectory() as tmpdir:
        # Copy the reference file multiple times with different names
        ref_file = REFERENCE_DIR / "lj38_global.xyz"
        if not ref_file.exists():
            print(f"❌ Reference file not found: {ref_file}")
            return
        
        num_files = 20
        for i in range(num_files):
            dst = Path(tmpdir) / f"lj38_copy_{i:03d}.xyz"
            dst.write_text(ref_file.read_text())
        
        print(f"🚀 Batch Verification Test")
        print(f"Files: {num_files}")
        print(f"Binary: {BATCH_BIN}")
        print(f"Threads: auto (all cores)")
        print()
        
        # Run batch verifier
        result = subprocess.run(
            [str(BATCH_BIN), tmpdir, "8"],
            capture_output=True,
            text=True,
            timeout=120
        )
        
        if result.returncode != 0:
            print(f"❌ Batch verifier failed:\n{result.stderr}")
            return
        
        # Parse JSON output
        try:
            data = json.loads(result.stdout)
        except json.JSONDecodeError as e:
            print(f"❌ Failed to parse JSON: {e}")
            print(f"stdout: {result.stdout[:500]}")
            return
        
        print("=" * 70)
        print("RESULTS")
        print("=" * 70)
        print(f"Total files: {data['num_files']}")
        print(f"Threads used: {data['num_threads']}")
        print(f"Batch elapsed: {data['batch_elapsed_sec']:.3f} sec")
        print(f"Avg time/file: {data['avg_time_per_file']:.4f} sec")
        print(f"Certified: {data['certified_count']}/{data['num_files']}")
        print()
        
        # Compare with single-file baseline
        single_time = 0.35  # from previous tests
        speedup = single_time / data['avg_time_per_file']
        
        print("COMPARISON")
        print("-" * 70)
        print(f"Single-file baseline: ~{single_time:.2f} sec/file")
        print(f"Batch (OpenMP): {data['avg_time_per_file']:.4f} sec/file")
        print(f"Speedup: {speedup:.2f}x")
        print()
        
        if speedup >= 4.0:
            print("🎉 SUCCESS: Achieved ≥4x speedup with OpenMP batching!")
        elif speedup >= 2.0:
            print("⚠️  Moderate speedup (2-4x). Consider increasing thread count.")
        else:
            print("❌ Low speedup (<2x). Possible causes:")
            print("   - Memory bandwidth bottleneck")
            print("   - FLINT/Arb internal locking")
            print("   - Insufficient CPU cores")

if __name__ == "__main__":
    main()