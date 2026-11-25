"""Parse recipes run output.

Parse typical batch job output files like .out and .err to identify
recipes that have succeeded or failed; display results in a convenient
Markdown format, to be added to a GitHub issue or any other such
documentation.
"""

import datetime
import os
import re
from collections.abc import Iterator
from pathlib import Path

import fire


def parse_slurm_output(dirname: str, pattern: str) -> Iterator[Path]:
    """Parse the out dir from SLURM.

    Perform a glob on dirname/pattern where dirname is the directory
    where SLURM output is stored, and pattern is the out file pattern,
    like .out. Returns all the files in dirname that have pattern
    extension.
    """
    return Path(dirname).expanduser().glob(pattern)


def parse_output_file(slurm_out_dir: str) -> dict[str, list[str]]:
    """Parse .out and .err files in a given dir.

    Returns a tuple of lists of sorted .out files for each of these
    criteria: recipes that ran successfulltm recipes that failed with
    diagnostic errors, recipes that failed due to missing data.
    """
    categories = [
        "success",
        "diagnostic error",
        "missing data",
        "out of memory",
        "out of time",
        "unknown",
    ]
    results: dict[str, list[str]] = {k: [] for k in categories}

    files = parse_slurm_output(slurm_out_dir, "*.out")
    for file in files:
        recipe = str(Path(file.stem).with_suffix(".yml"))
        with open(file, encoding="utf-8") as outfile:
            lines = outfile.readlines()
            for line in lines:
                if "Run was successful\n" in line:
                    results["success"].append(recipe)
                    break
                if "esmvalcore._task.DiagnosticError" in line:
                    results["diagnostic error"].append(recipe)
                    break
                if "ERROR   Missing data for preprocessor" in line:
                    results["missing data"].append(recipe)
                    break
            else:
                if not file.with_suffix(".err").exists():
                    results["unknown"].append(recipe)
                else:
                    err = file.with_suffix(".err").read_text(encoding="utf-8")
                    if (
                        "killed by the cgroup out-of-memory" in err
                        or "step tasks have been OOM Killed" in err
                    ):
                        results["out of memory"].append(recipe)
                    elif re.match(".* CANCELLED AT .* DUE TO TIME LIMIT", err):
                        results["out of time"].append(recipe)
                    else:
                        results["unknown"].append(recipe)

    results = {k: sorted(v) for k, v in results.items()}

    return results


def display_in_md(
    slurm_out_dir: str = ".",
    all_recipes_file: str = "all_recipes.txt",
) -> None:
    """Print out recipes in Markdown list.

    Parameters
    ----------
    slurm_out_dir:
        Directory where SLURM output files (.out and .err) are written to.

    all_recipes_file:
        Text file containing a list of all recipes.
    """
    todaynow = datetime.datetime.now()
    print(f"## Recipe running session {todaynow}\n")
    with open(all_recipes_file, encoding="utf-8") as file:
        all_recipes = [os.path.basename(line.strip()) for line in file]
    n_recipes = len(all_recipes)

    results = parse_output_file(slurm_out_dir)
    results["no run"] = sorted(
        set(all_recipes)
        - set(recipe for v in results.values() for recipe in v),
    )
    prefix = "Recipes that"
    err_prefix = f"{prefix} failed because"
    messages = {
        "success": f"{prefix} ran successfully",
        "diagnostic error": f"{err_prefix} the diagnostic script failed",
        "missing data": f"{err_prefix} of missing data",
        "out of time": f"{err_prefix} the run took too long",
        "out of memory": f"{err_prefix} they used too much memory",
        "unknown": f"{prefix} failed of other reasons or are still running",
        "no run": f"{prefix} never ran",
    }
    for type_, msg in messages.items():
        result = results[type_]
        if result:
            print(f"### {msg} ({len(result)} out of {n_recipes})")
            for recipe in result:
                print(f"- {recipe}")
            print()


if __name__ == "__main__":
    fire.Fire(display_in_md)
