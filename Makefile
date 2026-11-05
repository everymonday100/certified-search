.PHONY: help arch-lint gen-lint gen-test ver-build ver-test docker-build

help:
	@echo "Available targets:"
	@echo "  make arch-lint    - Check architectural boundaries (CRITICAL)"
	@echo "  make gen-lint     - Lint Python generator code"
	@echo "  make gen-test     - Run generator unit tests"
	@echo "  make ver-build    - Build C++ verifier locally"
	@echo "  make ver-test     - Run C++ verifier unit tests"
	@echo "  make docker-build - Build verifier Docker image"

# Architecture
arch-lint:
	@python scripts/check_architecture.py

# Generator (Python)
gen-lint:
	@cd core/generator && ruff check .
gen-test:
	@cd core/generator && pytest tests/ -v

# Verifier (C++)
ver-build:
	@cd core/verifier && cmake -B build -S . && cmake --build build
ver-test: ver-build
	@cd core/verifier && ctest --test-dir build --output-on-failure

# Container
docker-build:
	@docker build -t certified-search/verifier:local -f core/verifier/Dockerfile core/verifier/