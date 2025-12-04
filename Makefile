# Makefile for InstaNexus
# Run 'make help' to see available commands

.PHONY: help install format lint test clean check

# Display help message
help:
	@echo "Available commands:"
	@echo "  make install   - Install the package in editable mode"
	@echo "  make format    - Auto-format code using Isort, Black, and Ruff"
	@echo "  make lint      - Check code quality without modifying files (CI mode)"
	@echo "  make test      - Run unit tests with Pytest"
	@echo "  make clean     - Remove build artifacts, cache, and temporary files"

# Install package
install:
	pip install -e .

# Format code (Modifies files)
format:
	@echo ">>> Sorting imports with Isort..."
	isort .
	@echo ">>> Formatting code with Black..."
	black .
	@echo ">>> Fixing linting errors with Ruff..."
	ruff check --fix .

# Lint code (Read-only check)
lint:
	@echo ">>> Checking imports with Isort..."
	isort --check-only .
	@echo ">>> Checking formatting with Black..."
	black --check .
	@echo ">>> Checking code quality with Ruff..."
	ruff check .

# Shortcut for linting
check: lint

# Run tests
test:
	pytest

# Clean up artifacts
clean:
	@echo ">>> Cleaning up..."
	rm -rf dist/ build/ *.egg-info/ .pytest_cache/ .coverage
	find . -name "__pycache__" -type d -exec rm -rf {} +
	find . -name "*.pyc" -delete
	find . -name ".DS_Store" -delete