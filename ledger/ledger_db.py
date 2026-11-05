#!/usr/bin/env python3
"""
Append-only Evidence Ledger with hash-chain integrity.
Stores all verification attempts (CERTIFIED, REJECTED, INCONCLUSIVE).
"""
import sqlite3
import json
import hashlib
from pathlib import Path
from datetime import datetime, timezone

class EvidenceLedger:
    def __init__(self, db_path: str = "ledger/evidence.db"):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._init_db()

    def _init_db(self):
        with sqlite3.connect(self.db_path) as conn:
            conn.execute("""
                CREATE TABLE IF NOT EXISTS attempts (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    timestamp TEXT NOT NULL,
                    task_hash TEXT NOT NULL,
                    config_hash TEXT NOT NULL,
                    status TEXT NOT NULL,
                    witness_json TEXT,
                    prev_hash TEXT NOT NULL,
                    record_hash TEXT NOT NULL UNIQUE
                )
            """)
            conn.execute("CREATE INDEX IF NOT EXISTS idx_task_hash ON attempts(task_hash)")
            conn.execute("CREATE INDEX IF NOT EXISTS idx_status ON attempts(status)")
            conn.commit()

    def _get_last_hash(self) -> str:
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute("SELECT record_hash FROM attempts ORDER BY id DESC LIMIT 1")
            row = cursor.fetchone()
            return row[0] if row else "GENESIS"

    def _compute_record_hash(self, timestamp: str, task_hash: str, config_hash: str, 
                             status: str, witness_json: str, prev_hash: str) -> str:
        payload = f"{timestamp}|{task_hash}|{config_hash}|{status}|{witness_json}|{prev_hash}"
        return hashlib.sha256(payload.encode('utf-8')).hexdigest()

    def add_attempt(self, task_hash: str, config_hash: str, status: str, witness_data: dict = None):
        timestamp = datetime.now(timezone.utc).isoformat()
        witness_json = json.dumps(witness_data, sort_keys=True) if witness_data else ""
        prev_hash = self._get_last_hash()
        
        record_hash = self._compute_record_hash(
            timestamp, task_hash, config_hash, status, witness_json, prev_hash
        )

        with sqlite3.connect(self.db_path) as conn:
            try:
                conn.execute("""
                    INSERT INTO attempts (timestamp, task_hash, config_hash, status, witness_json, prev_hash, record_hash)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                """, (timestamp, task_hash, config_hash, status, witness_json, prev_hash, record_hash))
                conn.commit()
            except sqlite3.IntegrityError:
                # Duplicate hash (should not happen unless timestamp collision, highly unlikely)
                pass
        
        return record_hash

    def verify_integrity(self) -> bool:
        """Проверяет целостность хэш-цепочки. Возвращает True, если цепочка не нарушена."""
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute("SELECT id, timestamp, task_hash, config_hash, status, witness_json, prev_hash, record_hash FROM attempts ORDER BY id ASC")
            rows = cursor.fetchall()
            
            if not rows:
                return True # Пустой ledger считается валидным

            expected_prev = "GENESIS"
            for row in rows:
                db_id, ts, t_hash, c_hash, status, w_json, db_prev_hash, db_record_hash = row
                
                if db_prev_hash != expected_prev:
                    print(f"❌ Integrity FAIL at ID {db_id}: prev_hash mismatch")
                    return False
                
                computed_hash = self._compute_record_hash(ts, t_hash, c_hash, status, w_json, db_prev_hash)
                if computed_hash != db_record_hash:
                    print(f"❌ Integrity FAIL at ID {db_id}: record_hash mismatch")
                    return False
                
                expected_prev = db_record_hash
                
            print(f"✅ Ledger integrity verified. {len(rows)} records chained correctly.")
            return True

if __name__ == "__main__":
    # Тест Ledger
    ledger = EvidenceLedger()
    ledger.add_attempt("sha256:task123", "sha256:config456", "INCONCLUSIVE", {"reason": "kappa >= 1"})
    ledger.add_attempt("sha256:task123", "sha256:config789", "CERTIFIED_STRICT_LOCAL_MINIMUM", {"energy": "[-44.3, -44.3]"})
    ledger.verify_integrity()