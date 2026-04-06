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

## Architecture & Design Principles

### Optional flags gate optional features
`--metadata-json-path` is intentionally optional, following the same pattern as `--reference`.
Neither flag is required to run the core assembly pipeline.

| Flag | Without it | With it |
|---|---|---|
| `--metadata-json-path` | Minimal mode: all peptides treated as a single pool, no protease splitting, no chain filtering, no contaminant removal | Full mode: protease assignment, chain filtering, contaminant removal |
| `--reference` | No reference mapping or coverage statistics | Reference-based validation and assembly quality metrics |

**Do NOT make `--metadata-json-path` required.** This design exists to support users
starting from any de novo sequencing tool output, not just InstaNovo.

### CSV input assumptions in minimal mode
When `--metadata-json-path` is absent, **do not assume any InstaNovo-specific columns**
(including `experiment_name`). The only guaranteed columns are the peptide sequence
and a confidence score. All protease-based and chain-based logic must be gated
behind a metadata presence check.

### Where to find the gating pattern
The `--reference` flag handling in `src/instanexus/main.py` is the canonical example
of how optional features should be gated. Use the same `if metadata is not None:`
pattern when implementing metadata-optional behavior. Do not refactor unrelated code.
```

---