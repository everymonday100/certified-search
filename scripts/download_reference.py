#!/usr/bin/env python3
"""
Downloads raw coordinate files from Cambridge Cluster Database (/points/N endpoints)
and records provenance metadata: URL, retrieval date, SHA-256 of raw bytes.
"""
import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

BASE_URL = "https://www-wales.ch.cam.ac.uk/~wales/CCD/jon/structures/LJ/points/{n}"
RAW_TXT_DIR = Path("reference/raw_txt")
META_FILE = RAW_TXT_DIR / "source_metadata.json"
CLUSTERS = [13, 38, 55]

def main():
    RAW_TXT_DIR.mkdir(parents=True, exist_ok=True)
    meta = {
        "source": "Cambridge Cluster Database",
        "retrieval_date_utc": datetime.now(timezone.utc).isoformat(),
        "files": {}
    }

    for n in CLUSTERS:
        url = BASE_URL.format(n=n)
        print(f"⬇️  {url}")
        req = urllib.request.Request(
            url,
            headers={"User-Agent": "certified-search-research/0.1 (R&D feasibility snapshot)"}
        )
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = resp.read()

        raw_hash = hashlib.sha256(data).hexdigest()
        txt_file = RAW_TXT_DIR / f"{n}.txt"
        txt_file.write_bytes(data)

        meta["files"][str(n)] = {
            "url": url,
            "sha256": raw_hash,
            "bytes": len(data)
        }
        print(f"   -> {txt_file} ({len(data)} bytes, sha256={raw_hash[:16]}...)")

    META_FILE.write_text(json.dumps(meta, indent=2))
    print(f"📄 Provenance: {META_FILE}")

if __name__ == "__main__":
    main()