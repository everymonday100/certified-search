from pathlib import Path

def test_frozen_reference_fixture_exists():
    root = Path(__file__).resolve().parents[3]
    p = root / "reference" / "frozen" / "lj13" / "lj13.xyz"
    assert p.exists(), f"missing frozen fixture: {p}"
    lines = p.read_text().strip().splitlines()
    assert int(lines[0]) == 13, "LJ13 fixture must declare 13 atoms"
