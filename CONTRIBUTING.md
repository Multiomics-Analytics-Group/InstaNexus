# Contributing Guidelines

## Commit Messages
We follow Conventional Commits. Use the template provided in `.gitmessage`.

## Python Style Guide
- **Docstrings**: Follow [Google style guide](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings). Do not specify types in docstrings.
- **Typing**: Use type hints for all function arguments and return values.
- **Tensors**: Always annotate shapes in the line before the operation using `?` for batch size.
  ```python
  # (?, seq_len, d_model)
  output = model(inputs)
  ```
- **Strings**: Use `f-string` to format variables in strings instead of `%` or `.format()`.
- **Paths**: Use `pathlib.Path` to deal with local paths and `cloudpathlib.CloudPath` for remote paths.

## Releasing
Releases are published to PyPI automatically by `.github/workflows/python-publish.yml`.

1. **Bump the version**: on a branch from `main`, set `version = "X.Y.Z"` in `pyproject.toml`, run `uv lock`, and open a PR.
2. **Merge** the PR into `main`.
3. **Publish a GitHub Release** with tag `vX.Y.Z` targeting `main`, e.g.
   ```bash
   gh release create vX.Y.Z --target main --title "InstaNexus vX.Y.Z" --notes-file notes.md
   ```
   Publishing the release triggers the workflow, which fails if the tag (without the leading `v`) differs from `version` in `pyproject.toml`, then builds and uploads to PyPI.

To retry a failed publish, use "Re-run jobs" on the release's workflow run; PyPI never accepts the same version twice.

