# CLAUDE.md — Project Guidelines

## Build & Environment
- Package Manager: `uv`
- Environment Sync: `uv sync --all-extras`
- Pre-commit: `pre-commit run --all-files`

## Code Style
- **Formatting**: f-strings only, no blank lines after `if/for/while`, blank line before `return` (unless it's the only line).
- **Typing**: Mandatory for arguments and return values.
- **Tensors**: Annotate shapes on the line before the operation (e.g., `# (?, seq_len)`).
- **Docstrings**: Google Style (no types in docstrings).
- **Paths**: Use `pathlib.Path` (local) or `cloudpathlib.CloudPath` (remote).

## Common Commands
- Test: `uv run pytest`
- Lint: `uv run ruff check .`
- Type Check: `uv run mypy --config-file pyproject.toml`
