#!/usr/bin/env python3
"""
Converts raw CCD coordinate text into standard .xyz format.
Validates the parse by recomputing the LJ energy and comparing to the reference value.
"""
import re
from itertools import combinations
from pathlib import Path

RAW_TXT_DIR = Path("reference/raw_txt")
XYZ_DIR = Path("reference/raw")

# Reference energies, TZ v2.1 Table 2.2 (Cambridge Cluster Database)
ENERGIES = {13: -44.326801, 38: -173.928427, 55: -279.248470}
ENERGY_TOL = 1e-4  # coordinates have 10 decimals; a wrong parse deviates by O(1) or more

def lj_energy(coords):
    """E = 4 * sum_{i<j} [ r^-12 - r^-6 ], eps=sigma=1."""
    e = 0.0
    for (xi, yi, zi), (xj, yj, zj) in combinations(coords, 2):
        r2 = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
        r6 = r2 ** -3
        e += r6 * r6 - r6
    return 4.0 * e

def main():
    XYZ_DIR.mkdir(parents=True, exist_ok=True)

    for n in [13, 38, 55]:
        txt_file = RAW_TXT_DIR / f"{n}.txt"
        if not txt_file.exists():
            print(f"⚠️ {txt_file.name} not found, skipping.")
            continue

        text = txt_file.read_text(encoding="utf-8", errors="replace")
        numbers = re.findall(r'-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?', text)
        floats = [float(x) for x in numbers]

        expected = n * 3
        if len(floats) != expected:
            print(f"❌ {txt_file.name}: found {len(floats)} numbers, expected {expected}. Aborting this cluster.")
            continue

        coords = [floats[i:i + 3] for i in range(0, len(floats), 3)]

        # Parse validation: recompute energy
        e_calc = lj_energy(coords)
        if abs(e_calc - ENERGIES[n]) > ENERGY_TOL:
            print(f"❌ {txt_file.name}: energy mismatch. computed={e_calc:.6f}, reference={ENERGIES[n]:.6f}.")
            print("   Parse order is likely wrong. File NOT converted.")
            continue
        print(f"   energy check OK: {e_calc:.6f} vs reference {ENERGIES[n]:.6f}")

        xyz_file = XYZ_DIR / f"lj{n}.xyz"
        with xyz_file.open("w") as f:
            f.write(f"{n}\n")
            f.write(f"Energy: {ENERGIES[n]}\n")
            for x, y, z in coords:
                f.write(f"X {x:.10f} {y:.10f} {z:.10f}\n")
        print(f"✅ {xyz_file}")

if __name__ == "__main__":
    main()