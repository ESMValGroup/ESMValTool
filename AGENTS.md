# AGENTS.md

Guidance for AI coding agents working in the ESMValTool repository.
The full contributor guide is in
[doc/sphinx/source/community/](doc/sphinx/source/community/index.rst);
this file summarizes what matters most when making changes.

## Project overview

ESMValTool is a community tool for evaluating Earth system models.
This repository contains *recipes* (YAML files that define the input data, the
preprocessing, and the diagnostic scripts to run), the diagnostic scripts that
analyze and plot the preprocessed data, and scripts that reformat
("CMORize") observational datasets.
Running recipes, finding and preprocessing data, and fixing datasets is done
by [ESMValCore](https://github.com/ESMValGroup/ESMValCore), which provides the
`esmvaltool` command. Changes to that functionality belong in ESMValCore.

Key locations:

- `esmvaltool/recipes/` – recipes; `recipes/examples/` contains examples.
- `esmvaltool/diag_scripts/` – diagnostic scripts in Python, NCL, and R.
  Code shared by many diagnostics (plotting, provenance, input/output) lives
  in `diag_scripts/shared/`.
- `esmvaltool/cmorizers/data/` – the `esmvaltool data` command, with scripts
  that download (`downloaders/`) and reformat (`formatters/`) observational
  datasets.
- `esmvaltool/config-references.yml` – authors, projects, and other tags used
  in recipes and provenance records.
- `esmvaltool/references/` – BibTeX files for the references used in recipes,
  diagnostics, and datasets.
- `esmvaltool/utils/` – tools for developers and maintainers.
- `doc/sphinx/source/` – Sphinx documentation (reStructuredText), with a page
  for each recipe in `recipes/`.
- `tests/` – tests.
- `archive/` – code that is no longer maintained; do not change it.

## Development environment

The development environment is managed with [pixi](https://pixi.sh):

```bash
pixi run test                       # run the test suite
pixi run doc --fail-on-warning      # build the documentation
pixi run pytest tests/unit/test_recipes.py   # run a single test file
pixi run pre-commit run --all-files # run the formatters and basic checks
pixi run ruff check <file>          # run all ruff checks on a file
pixi run pyrefly check <file>       # check the type hints in a file
pixi run yamllint <file>            # check a recipe or other YAML file
```

The default environment only supports Python diagnostics. Use
`pixi run -e esmvaltool-dev ...` for NCL and R diagnostics (Linux only) and
`pixi run -e esmvalcore-dev ...` to test with the development version of
ESMValCore.

Useful pytest options:

- `-n auto` runs tests in parallel (pytest-xdist).
- `pixi run pytest tests/integration/test_recipes_loading.py -k <recipe>`
  checks that a recipe is valid without running it.

Running a recipe (`esmvaltool run <recipe>`) requires the input data and an
ESMValCore configuration, which may not be available. If you cannot run a
recipe that you changed, say so and ask the user to run it.

While developing a diagnostic script, it is useful to run the preprocessing
only once and then rerun just the diagnostic scripts:

```bash
# Run the preprocessing only and keep its output
esmvaltool run <recipe> --remove_preproc_dir=False --run_diagnostic=False
# Run the diagnostic scripts on the preprocessed data from the previous run
esmvaltool run <recipe> --resume_from /path/to/esmvaltool_output/recipe_<name>_<date>_<time>
```

If the data is only available on an HPC system, the user can run the
preprocessing there and copy the output directory to their own computer, so
the diagnostic scripts can be run locally.
ESMValCore only resumes from a run of an identical recipe, so after changing
the recipe, the preprocessing must be run again (if the user is sure that the
preprocessed data has not changed, the check can be commented out in the
function `parse_resume` in `esmvalcore/_main.py` in ESMValCore).

## Design

- Keep changes small and focused. Changes to CMORizer scripts always go in a
  separate pull request, also when a new recipe uses the dataset.
- Write new diagnostics and CMORizer scripts in Python. Only use R for
  new diagnostic scripts if a specific R package is required or the user is
  more comfortable with R. Do not write new NCL code: the NCL
  interpreter is no longer maintained.
- Reuse existing code, e.g. from `esmvaltool/diag_scripts/shared/`, instead of
  writing new code. Prefer ESMValCore preprocessor functions in the recipe over
  implementing preprocessing in a diagnostic script.
- Keep functions short (no more than about 50 lines), and code and
  documentation concise.
- Do not hardcode paths: get input files and output directories from the
  configuration that ESMValCore passes to the diagnostic script.
- Do not copy code from other projects unless its license is compatible with
  Apache 2.0.

## Code style

- Python code is formatted and checked with
  [ruff](https://docs.astral.sh/ruff/) with a line length of 79; see
  `[tool.ruff]` in `pyproject.toml`. pre-commit only runs `ruff format` and a
  subset of the ruff rules, so also run `ruff check` on the Python files you
  changed and fix the issues it reports. Human-written code is allowed to
  contain some violations of the ruff rules that are not checked by pre-commit,
  but this exception does not apply to code generated by AI tools.
- Write [numpy style docstrings](https://numpydoc.readthedocs.io/en/latest/format.html)
  for functions that appear in the documentation (see
  `doc/sphinx/source/api/`). Other functions only need a short docstring.
- Type hints are not required in diagnostic scripts, but if you add them they
  must be correct; check the files you changed with
  [pyrefly](https://pyrefly.org) and fix the errors in the code you changed.
  The existing code is not free of type errors, so do not check the whole
  repository or fix errors in code unrelated to the goal of the change. Ask the
  user if they want type hints for the code they are working on.
- Do not add `# noqa` comments to silence a check unless there is no
  reasonable way to fix the underlying issue.
- Check YAML files, including recipes, with yamllint. R code should follow the
  [tidyverse style guide](https://style.tidyverse.org/). NCL code is checked
  with `nclcodestyle` by the tests.
- Match the style, naming, and comment density of the surrounding code.

## Tests

- Where reasonable, add tests for new and changed code in `tests/unit/` or
  `tests/integration/`, mirroring the package structure. This applies in
  particular to shared code, utilities, and CMORizer helper functions.
  Detailed tests are not required for diagnostic scripts or CMORizer scripts,
  as these are usually scripts that compute a single figure or reformat a
  dataset and their correctness is verified by inspecting the resulting figure
  or dataset. A single integration test for a CMORizer script may be sufficient.
  If a diagnostic script provides options that can be configured in the recipe,
  some tests that the script behaves correctly for different options are useful.
  If a diagnostic script implements a complex algorithm, tests that verify that
  the algorithm is implemented correctly are helpful. Keep the tests as
  light-touch as possible, but make sure the scientific results are reliable.
- Tests use pytest; prefer plain functions, fixtures, `pytest.mark.parametrize`,
  and `pytest-mock`'s `mocker` fixture over `unittest` style.
- Do not change existing tests to make new code pass: having to change them
  may mean that existing features break. Flag such changes to the user.
- The tests check that all recipes are valid, have a `maintainer`, and that
  a BibTeX file exists for each reference.
- `xfail_strict = true` is set, so an `xfail` test that unexpectedly passes
  counts as a failure.

## Documentation

- Update the documentation in `doc/sphinx/source/` when behavior visible to
  users changes.
- In `.rst` files, keep lines at most 80 characters long and start each
  sentence on a new line.
- Link to the ESMValCore documentation with `:ref:` and the `esmvalcore:`
  prefix, e.g. `` :ref:`esmvalcore:recipe` ``.
- For functions that compute scientific results, add comments with references
  to papers and formula numbers.
- The documentation must build without warnings; the CI builds it with
  `--fail-on-warning`.

## Common tasks

### Adding or changing a recipe and diagnostic

See
[doc/sphinx/source/community/diagnostic.rst](doc/sphinx/source/community/diagnostic.rst).
In short:

1. Start from `esmvaltool/recipes/examples/recipe_python.yml` and
   `esmvaltool/diag_scripts/examples/diagnostic.py`. Put the recipe in
   `esmvaltool/recipes/` and the diagnostic script in a directory under
   `esmvaltool/diag_scripts/`. Name the recipe `recipe_<name>.yml`: the tests
   only check files that match `recipe*.yml`.
2. Fill in the `documentation` section of the recipe (title, description,
   authors, maintainer, references, projects) and add a `description`,
   `themes`, and `realms` to each diagnostic. Use `unmaintained` if the recipe
   has no maintainer. Authors, projects, and other tags must be defined in
   `esmvaltool/config-references.yml`; each reference needs a BibTeX file
   `esmvaltool/references/<tag>.bibtex`.
3. In Python diagnostics, run the main function with
   `esmvaltool.diag_scripts.shared.run_diagnostic` and save figures and data
   with `save_figure` and `save_data` from the same module, passing a
   provenance record for each output file. Use `ProvenanceLogger` for other
   output files. Also save the data shown in each plot to a NetCDF file.
   Running the recipe should not give any warnings about provenance.
4. Add a documentation page by copying
   `doc/sphinx/source/recipes/recipe_template.rst.template` to
   `doc/sphinx/source/recipes/recipe_<name>.rst` and adding it to
   `doc/sphinx/source/recipes/index.rst`. Add an example figure (PNG, about
   150 dpi) for each type of plot to `doc/sphinx/source/recipes/figures/`.

### Adding or changing a CMORizer script

See [doc/sphinx/source/develop/dataset.rst](doc/sphinx/source/develop/dataset.rst)
and [doc/sphinx/source/community/dataset.rst](doc/sphinx/source/community/dataset.rst).
In short:

1. Add the script to `esmvaltool/cmorizers/data/formatters/datasets/`, with
   instructions on how to obtain the data and the references in its header,
   and its configuration to `esmvaltool/cmorizers/data/cmor_config/`. Use the
   helper functions in `esmvaltool/cmorizers/data/utilities.py`. A script that
   downloads the data can be added to
   `esmvaltool/cmorizers/data/downloaders/datasets/`.
2. Add the dataset to `esmvaltool/cmorizers/data/datasets.yml` and to the
   table of supported datasets in `doc/sphinx/source/input.rst`, and add a
   BibTeX file for its reference to `esmvaltool/references/`.
3. Check that the CMORizer script works by running
   `esmvaltool data format --original-data-dir <dir> <dataset>`, with the raw
   data in `<dir>/Tier<N>/<dataset>/`. Use `esmvaltool data prepare` instead
   of `format` to download the data first.
4. Check that the data has been correctly CMORized by loading it with Iris
   and passing it to the function `esmvalcore.cmor.check.cmor_check` with the
   strict check level. Any Iris warnings during load or issues reported by the
   CMOR checker must be fixed.
   The [IOOS compliance checker](https://github.com/ioos/compliance-checker)
   can also be used to check that the data is CF-compliant, and its
   [WCRP plugin](https://github.com/ESGF/cc-plugin-wcrp) to check that the
   data is compliant with the WCRP data standards, though some checks may not
   be applicable. An example of how to use these compliance checkers is
   available in
   [woodpecker-esa-cci-plugin](https://github.com/ESMValGroup/woodpecker-esa-cci-plugin/tree/main/tests/integration).
5. Add the dataset to `esmvaltool/recipes/examples/recipe_check_obs.yml`, in
   a diagnostic named after the dataset, and run that recipe for the new
   dataset only with
   `esmvaltool run examples/recipe_check_obs.yml --diagnostics=<dataset>`.
   Use `--diagnostics=<dataset>/<variable>` to check a single variable. The
   CMORized data must be in the directory configured for the `OBS6` or `OBS`
   project.

## Backward compatibility

Many people use ESMValTool recipes and diagnostics, also outside this
repository. A change that causes a recipe to no longer run (*breaking change*)
or that changes scientific results (*science change*) is subject to the
[backward compatibility policy](doc/sphinx/source/community/backward_compatibility.rst).
Avoid such changes where possible and flag them clearly to the user rather
than making them silently. Code in `esmvaltool/diag_scripts/shared/` is used by
many diagnostics, so check all its users when changing it.

## Dependencies

When adding or removing a dependency, update `pyproject.toml` and run
`pixi lock` to update `pixi.lock`:

- Dependencies of diagnostics and CMORizer scripts go in the `all` list in
  `[project.optional-dependencies]` (PyPI names) and in
  `[tool.pixi.feature.all.dependencies]` (conda-forge names). Keep the two
  lists in sync.
- `[project] dependencies` and `[tool.pixi.dependencies]` only contain what is
  needed to run `recipes/examples/recipe_python.yml`; do not add to them
  unless that is really needed.
- R and NCL dependencies go in `[tool.pixi.feature.r.dependencies]` and
  `[tool.pixi.feature.ncl.dependencies]`, and tools that are only needed for
  development, building the documentation, or running the tests go in
  `[tool.pixi.feature.dev.dependencies]`,
  `[tool.pixi.feature.doc.dependencies]`, or
  `[tool.pixi.feature.test.dependencies]`.

New dependencies must be actively maintained, be available on conda-forge,
and have a license that is compatible with Apache 2.0, as must their own
dependencies. Supported Python and dependency versions follow
[SPEC-0](https://scientific-python.org/specs/spec-0000/).

## Pull requests

- Follow the [AI policy](doc/sphinx/source/community/ai_policy.rst). In
  particular: do not open pull requests autonomously (a human must review the
  code and submit it), do not write pull request descriptions, issue
  descriptions, or comments on behalf of the contributor, and remind the
  contributor to disclose which AI tools were used and for what.
- The changelog is generated from pull request titles and labels, so use a
  short, descriptive title that is understandable to users. Do not edit
  `doc/sphinx/source/changelog.rst` by hand for regular changes.
- Do not modify code copied from other projects, such as
  `esmvaltool/diag_scripts/cvdp/cvdp/`,
  `esmvaltool/diag_scripts/extreme_events/climdex.pcic.ncdf/`, or
  `esmvaltool/utils/nclcodestyle/`, unless that is the goal of the change.
