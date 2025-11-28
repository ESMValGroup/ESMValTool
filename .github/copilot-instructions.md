## Quick context for AI code contributors

This repository contains ESMValTool — a community climate-model diagnostics and
preprocessing framework. The goal of these notes is to give an AI agent the
minimum, high-value knowledge to be productive immediately: how the project is
structured, common developer workflows, and concrete examples to follow.

### Quick setup (developer environment)

- Prefer Conda/Mamba using the provided environment files at the repository root:

```bash
# create the dev environment (mamba recommended)
mamba env create --name esmvaltool --file environment.yml
conda activate esmvaltool
# or use environment_python.yml for a lighter python-only env
```

- Docker: the repository includes `docker/Dockerfile`. Build from the repo root:

```bash
docker build -t esmvaltool:latest . -f docker/Dockerfile
```

- Installation for local development: pip install editable mode and tests extras

```bash
pip install -e .[test]
pytest   # run tests (pytest options are configured in pyproject.toml)
```

### Where to run diagnostics / recipes

- Recipes and diagnostics are described in the README and docs. Diagnostics
  live under `esmvaltool/diag_scripts/`. Example diagnostic:
  `esmvaltool/diag_scripts/weathertyping/weathertyping.py`.

### Architecture & major components (short)

- `esmvaltool/` — main package. Diagnostics and helper code live here.
  - `diag_scripts/` — each diagnostic in its own sub-folder. Diagnostics are
    Python scripts but may call NCL/R tools.
  - `cmorizers/` — scripts to convert observational datasets to CMOR format.
- `ESMValTool_stuff/`, `preval/` — user recipes, examples, and auxiliary
  analysis code. Keep changes to these under separate branches unless they are
  core fixes.
- External integration: heavy reliance on `ESMValCore` (preprocessing), Iris
  (`scitools-iris`) for cube handling, and many scientific packages (see
  `pyproject.toml` and `environment.yml`). Changes to dataflow must respect
  ESMValCore's preprocessing contracts.

### Diagnostic implementation conventions (concrete patterns)

- Entry pattern: diagnostics obtain config via the shared context manager:

```python
from esmvaltool.diag_scripts.shared import run_diagnostic

if __name__ == '__main__':
    with run_diagnostic() as config:
        run_my_diagnostic(config)
```

Follow the `weathertyping` example: functions take the `cfg` dict returned by
`run_diagnostic()` and extract preprocessor variables through helper functions
in the same diagnostic package.

- File I/O: diagnostics generally read/write netCDF files (Iris cubes). Use the
  existing helper utilities in `diag_scripts/<name>/` (e.g., `wt_utils`,
  `calc_utils`) instead of re-implementing low-level file-handling.

### Testing, linting and formatting specifics

- Tests: configured through `pyproject.toml`. After `pip install -e .[test]` run
  `pytest`. The project uses 79-character line length.
- Linting: `ruff` + `prospector` (see `pyproject.toml`). Prefer automatic fixes
  where safe (the repo config sets `fix = true` for ruff).

### Packaging and runtime notes to keep in mind

- The package expects Python >= 3.11 (see `pyproject.toml`).
- `esmvaltool/__init__.py` raises a helpful error if the package is not installed
  editable — prefer `pip install -e .` for local imports and tests.
- Many dependencies come from conda-forge and some (e.g., ESMPy) are not on
  PyPI. For CI-like runs prefer the provided `environment.yml` or the Docker
  image.

### Helpful repository-specific search patterns for agents

- Diagnostics: `esmvaltool/diag_scripts/**/` — look for `run_diagnostic` use.
- Recipes/user examples: `ESMValTool_stuff/` and `ESMValTool_output/`.
- Environment and build: `environment.yml`, `environment_python.yml`,
  `docker/Dockerfile`, `pyproject.toml`.

### What to change and what to avoid

- Safe: small diagnostic improvements, refactors within a diagnostic folder,
  tests for new behavior, documentation updates to README/docs snippets.
- Avoid: large changes to preprocessing contracts, switching major dependency
  versions, or edits to recipes in `ESMValTool_stuff/` that alter user data
  provenance without discussing with maintainers.

### Quick examples to cite in PRs

- Point reviewers to `esmvaltool/diag_scripts/weathertyping/weathertyping.py`
  as a model for diagnostic layout and `docker/Dockerfile` and
  `environment.yml` for environment reproducibility.

### When to ask a human

- If a change affects: the public API of preprocessing (ESMValCore contracts),
  Docker/CI pipelines, or requires new conda-build recipes. Also ask if you
  need access to large datasets (CMIP/obs) that are not in the repo.

If anything here is unclear or you'd like a shorter/longer variant, tell me
which parts to expand and I'll iterate.


# Copilot Agents for VSCode

## Auto-Fix Codacy Minor Issues Agent

This agent is designed to assist developers by automatically fixing **Codacy minor issues** in code, using GitHub Copilot suggestions within VSCode.

### Overview

- **Purpose:** Automatically detect and fix minor code quality issues flagged by Codacy.
- **Scope:** Minor issues only (e.g., formatting, unused imports, variable naming, unnecessary semicolons).
- **Toolchain:** VSCode + GitHub Copilot + Codacy CLI (optional for local linting).
- **Trigger:** On file save or on-demand.

### Setup

1. **Install VSCode Extensions**
   - [GitHub Copilot](https://marketplace.visualstudio.com/items?itemName=GitHub.copilot)
   - [Codacy VSCode Extension](https://marketplace.visualstudio.com/items?itemName=Codacy.codacy)
   - Optional: Linter/formatter for your language (e.g., ESLint, Prettier, Black).

2. **Enable Copilot Suggestions**
   - Go to VSCode `Settings` → `Copilot` → Enable inline suggestions.
   - Enable `Accept Copilot suggestions on save` for automated workflow (optional, may require keybinding).

3. **Codacy CLI (optional)**
   - Install Codacy CLI:  
     ```bash
     npm install -g @codacy/cli
     ```
   - Run to check issues locally:  
     ```bash
     codacy analyze
     ```

### Workflow

1. Open your project in VSCode.
2. Open the file flagged with minor Codacy issues.
3. Accept Copilot inline suggestions to automatically fix formatting, variable names, or other minor issues.
4. Save the file to trigger any auto-formatters configured.
5. Run Codacy analysis to ensure issues are resolved.

### Example Use Cases

- **Unused import:**
  ```python
  import os  # Minor issue: unused import
