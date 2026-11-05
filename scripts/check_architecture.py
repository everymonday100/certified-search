#!/usr/bin/env python3
"""
Architectural Boundary Enforcer.
Fails the build if `core.verifier` imports anything from `core.generator`.
Uses Python AST to catch absolute and relative imports.
"""
import ast
import sys
from pathlib import Path

VERIFIER_DIR = Path("core/verifier")
FORBIDDEN_PREFIXES = ["core.generator", "generator", "core.generator."]
FORBIDDEN_RELATIVE_TARGETS = ["generator", "../generator", "../../generator"]

def check_file(filepath: Path) -> list[str]:
    errors = []
    try:
        source = filepath.read_text(encoding="utf-8")
        tree = ast.parse(source)
    except Exception as e:
        return [f"{filepath}: failed to parse AST ({e})"]

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if any(alias.name.startswith(p) or alias.name == p.rstrip('.') for p in FORBIDDEN_PREFIXES):
                    errors.append(f"{filepath}:{node.lineno} forbidden import: {alias.name}")
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            if any(module.startswith(p) or module == p.rstrip('.') for p in FORBIDDEN_PREFIXES):
                errors.append(f"{filepath}:{node.lineno} forbidden from-import: {module}")
            # Catch relative imports escaping the verifier zone
            if node.level > 0 and module in FORBIDDEN_RELATIVE_TARGETS:
                errors.append(f"{filepath}:{node.lineno} forbidden relative import escaping verifier zone")
    return errors

def main() -> int:
    if not VERIFIER_DIR.exists():
        print("Verifier directory not found, skipping arch check.")
        return 0
        
    all_errors = []
    for py_file in VERIFIER_DIR.rglob("*.py"):
        all_errors.extend(check_file(py_file))
        
    # Also check C++ includes (simple grep-based heuristic for v0)
    for cpp_file in VERIFIER_DIR.rglob("*.[ch]pp"):
        content = cpp_file.read_text(encoding="utf-8")
        for line_no, line in enumerate(content.splitlines(), 1):
            if "core/generator" in line or "generator/" in line:
                all_errors.append(f"{cpp_file}:{line_no} forbidden C++ include/mention of generator path")

    if all_errors:
        print("❌ ARCHITECTURE VIOLATION DETECTED ❌")
        print("The verifier must be strictly isolated from the generator.")
        for err in all_errors:
            print(f"  - {err}")
        return 1
        
    print("✅ Architectural boundaries verified. Verifier is isolated.")
    return 0

if __name__ == "__main__":
    sys.exit(main())