#!/usr/bin/env python3
"""
Generates frozen task_spec.json files and calculates task_hash.
"""
import hashlib
import json
from pathlib import Path

SPECS_DIR = Path("specs")
SPECS_DIR.mkdir(exist_ok=True)

def generate_task_spec(n_atoms: int) -> dict:
    return {
        "schema_version": "1.0",
        "model": "LJ_12_6",
        "N": n_atoms,
        "dimension": 3,
        "epsilon": "1",
        "sigma": "1",
        "boundary": "free",
        "cutoff": "none",
        "units": "reduced"
    }

def main():
    for n in [13, 38, 55]:
        spec = generate_task_spec(n)
        
        # Canonical JSON for hashing
        spec_str = json.dumps(spec, sort_keys=True)
        task_hash = hashlib.sha256(spec_str.encode('utf-8')).hexdigest()
        
        spec["task_hash"] = f"sha256:{task_hash}"
        
        file_path = SPECS_DIR / f"lj{n}.task.json"
        file_path.write_text(json.dumps(spec, indent=2))
        print(f"📝 Generated {file_path} | Hash: {task_hash[:16]}...")

if __name__ == "__main__":
    main()