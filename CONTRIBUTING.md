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
- **Strings**: Use `f-string` to format variables in strings instead of `%` or `.format()`.
- **Paths**: Use `pathlib.Path` to deal with local paths and `cloudpathlib.CloudPath` for remote paths.

