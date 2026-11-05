#!/usr/bin/env python3
"""
Freezes raw XYZ reference files into the reproducible corpus.
Calculates raw SHA-256 and canonical structure hash (invariant to rotation/permutation).
"""
import hashlib
import json
import shutil
import sys
from pathlib import Path
import numpy as np
from scipy.spatial.distance import pdist

RAW_DIR = Path("reference/raw")
FROZEN_DIR = Path("reference/frozen")
CHECKSUM_FILE = Path("reference/checksums.sha256")

def parse_xyz(file_path: Path) -> tuple[int, list[list[float]], float | None]:
    """Simple parser for standard XYZ format."""
    lines = file_path.read_text().strip().splitlines()
    n_atoms = int(lines[0].strip())
    comment = lines[1].strip()
    
    # Try to extract energy from comment line if present (e.g., "Energy: -44.326801")
    energy = None
    for part in comment.split():
        try:
            energy = float(part)
            break
        except ValueError:
            continue

    coords = []
    for i in range(2, 2 + n_atoms):
        parts = lines[i].split()
        # parts[0] is element, parts[1:4] are x, y, z
        coords.append([float(x) for x in parts[1:4]])
        
    return n_atoms, coords, energy

def canonical_hash(coords: list[list[float]], energy: float | None) -> str:
    """
    TZ v2.1 requirement: sorted pairwise distances, fixed precision, plus energy.
    Invariant to translation, rotation, and atom permutation.
    """
    coords_np = np.array(coords)
    dists = pdist(coords_np)
    
    # Discretize to 6 decimal places as per Cambridge DB standard
    dists_rounded = np.round(dists, 6)
    dists_sorted = np.sort(dists_rounded)
    
    payload = {
        "N": len(coords),
        "energy": round(energy, 6) if energy is not None else None,
        "pairwise_distances": dists_sorted.tolist()
    }
    
    # Canonical JSON string (sort_keys ensures determinism)
    canonical_str = json.dumps(payload, sort_keys=True)
    return hashlib.sha256(canonical_str.encode('utf-8')).hexdigest()

def main():
    if not RAW_DIR.exists():
        print(f"❌ Directory {RAW_DIR} not found. Please place your .xyz files there first.")
        sys.exit(1)

    FROZEN_DIR.mkdir(parents=True, exist_ok=True)
    CHECKSUM_FILE.write_text("") # Reset checksums

    xyz_files = list(RAW_DIR.glob("*.xyz"))
    if not xyz_files:
        print("⚠️ No .xyz files found in reference/raw/")
        sys.exit(0)

    print("❄️ Freezing reference corpus...")
    for xyz_file in xyz_files:
        cluster_name = xyz_file.stem.lower() # e.g., lj13
        target_dir = FROZEN_DIR / cluster_name
        target_dir.mkdir(exist_ok=True)
        
        # 1. Raw SHA-256
        raw_bytes = xyz_file.read_bytes()
        raw_hash = hashlib.sha256(raw_bytes).hexdigest()
        
        # 2. Parse and Canonical Hash
        n_atoms, coords, energy = parse_xyz(xyz_file)
        canon_hash = canonical_hash(coords, energy)
        
        # 3. Copy to frozen
        target_file = target_dir / xyz_file.name
        shutil.copy2(xyz_file, target_file)
        
        # 4. Write metadata
        meta = {
            "source_file": xyz_file.name,
            "N": n_atoms,
            "energy": energy,
            "raw_sha256": raw_hash,
            "canonical_hash": canon_hash
        }
        (target_dir / "metadata.json").write_text(json.dumps(meta, indent=2))
        
        # 5. Append to global checksums
        with CHECKSUM_FILE.open("a") as f:
            f.write(f"{raw_hash}  {target_file.as_posix()}\n")
            
        print(f"✅ Froze {cluster_name}: N={n_atoms}, E={energy}, CanonHash={canon_hash[:12]}...")

    print(f"\n🔒 Corpus frozen. Checksums written to {CHECKSUM_FILE}")

if __name__ == "__main__":
    main()