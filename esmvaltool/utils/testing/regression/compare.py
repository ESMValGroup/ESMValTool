"""Compare recipe runs to previous runs."""
from __future__ import annotations

import argparse
import difflib
import filecmp
import fnmatch
import os
import re
import sys
from pathlib import Path
from textwrap import indent
from typing import Iterator, Optional

import numpy as np
import xarray as xr
from PIL import Image

try:
    import imagehash
except ImportError:
    print("Please `pip install imagehash`")

IGNORE_FILES: tuple[str, ...] = (
    '*_citation.bibtex',
    '*_data_citation_info.txt',
    '*_info.ncl',
    '*_provenance.xml',
    'metadata.yml',
)
"""Files to ignore when comparing results."""

IGNORE_GLOBAL_ATTRIBUTES: tuple[str, ...] = (
    # see https://github.com/ESMValGroup/ESMValCore/issues/1657
    'auxiliary_data_dir',
    'creation_date',
    'history',
    'provenance',
    'software',
    # see https://github.com/ESMValGroup/ESMValCore/issues/1657
    'version',
)
"""Global NetCDF attributes to ignore when comparing."""

IGNORE_VARIABLE_ATTRIBUTES: tuple[str, ...] = IGNORE_GLOBAL_ATTRIBUTES
"""Variable NetCDF attributes to ignore when comparing."""

IGNORE_VARIABLES: tuple[str, ...] = (
    # see https://github.com/ESMValGroup/ESMValTool/issues/2714
    'temp_list',  # used by perfmetrics diagnostics to save absolute paths
)
"""Variables in NetCDF files to ignore when comparing."""

COMPARE_SUBDIRS: tuple[str, ...] = (
    'plots',
    'preproc',
    'work',
)
"""Directories of subdirectories to compare."""

RECIPE_DIR_DATETIME_PATTERN: str = r'[0-9]{8}_[0-9]{6}'
"""Regex pattern for datetime in recipe output directory."""

RECIPE_DIR_PATTERN: str = (
    r'recipe_(?P<recipe_name>[^\s]*?)_' + RECIPE_DIR_DATETIME_PATTERN
)
"""Regex pattern for recipe output directories."""


def as_txt(msg: list[str]) -> str:
    """Convert lines of text to indented text."""
    return indent('\n'.join(msg), "  ")


def diff_attrs(ref: dict, cur: dict) -> str:
    """Compare two dicts and describe the difference."""
    msg = []
    for key in ref:
        if key not in cur:
            msg.append(f"missing attribute '{key}'")
        elif not np.array_equal(ref[key], cur[key]):
            msg.append(f"value of attribute '{key}' is different: "
                       f"expected '{cur[key]}' but found '{ref[key]}'")
    for key in cur:
        if key not in ref:
            msg.append(f"extra attribute '{key}' with value '{cur[key]}'")
    msg.sort()
    txt = as_txt(msg)
    if txt:
        txt = "attributes are different:\n" + txt
    return txt


def diff_array(ref: np.ndarray, cur: np.ndarray) -> str:
    """Compare two arrays and describe the difference."""
    msg = []
    if cur.shape != ref.shape:
        msg.append("data has different shape")
    elif not np.array_equal(ref, cur):
        if np.issubdtype(ref.dtype, np.inexact) and np.issubdtype(
                cur.dtype, np.inexact) and np.allclose(ref, cur):
            msg.append("data is almost but not quite the same")
        else:
            msg.append("data is different")
    return as_txt(msg)


def diff_dataarray(ref: xr.DataArray, cur: xr.DataArray, type_: str) -> str:
    """Compare two xarray DataArrays and describe the difference."""
    msg = []
    if not cur.identical(ref):
        msg.append(f"{type_} '{cur.name}' is not identical to reference")
    if diff := diff_array(ref.values, cur.values):
        msg.append(diff)
    if diff := diff_attrs(ref.attrs, cur.attrs):
        msg.append(indent(diff, "  "))
    return as_txt(msg)


def diff_dataset(ref: xr.Dataset, cur: xr.Dataset) -> str:
    """Compare two xarray Datasets and describe the difference."""
    msg = []
    if diff := diff_attrs(ref.attrs, cur.attrs):
        msg.append(diff)

    for var in ref:
        if var not in cur:
            msg.append(f"missing variable '{var}'")
        else:
            if diff := diff_dataarray(ref[var], cur[var], "variable"):
                msg.append(diff)

            for coord in ref[var].coords:
                if coord not in cur.coords:
                    msg.append(f"missing coordinate '{coord}'")
                elif diff := diff_dataarray(ref[var].coords[coord],
                                            cur[var].coords[coord],
                                            'coordinate'):
                    msg.append(diff)
            for coord in cur[var].coords:
                if coord not in ref.coords:
                    msg.append(f"extra coordinate '{coord}'")

    for var in cur:
        if var not in ref:
            msg.append(f"extra variable '{var}'")

    return '\n'.join(msg)


def adapt_attributes(attributes: dict, ignore_attributes: tuple[str, ...],
                     recipe_name: str) -> dict:
    """Remove ignored attributes and make absolute paths relative."""
    new_attrs = {}

    for (attr, attr_val) in attributes.items():

        # Ignore attributes
        if attr in ignore_attributes:
            continue

        # Convert absolute paths to relative paths using the recipe name
        if isinstance(attr_val, str):
            recipe_dir = get_recipe_dir_from_str(attr_val, recipe_name)

            # If recipe_dir is present in attribute value, assume this
            # attribute value is a path and convert it to a relative path
            if recipe_dir is not None:
                attr_val = Path(attr_val).relative_to(recipe_dir)

        new_attrs[attr] = attr_val

    return new_attrs


def load_nc(filename: Path) -> xr.Dataset:
    """Load a NetCDF file."""
    dataset = xr.open_dataset(filename, chunks={}, decode_times=False)
    recipe_name = get_recipe_name_from_file(filename)

    # Remove ignored variables
    dataset = dataset.drop_vars(IGNORE_VARIABLES, errors='ignore')

    # Remove ignored attributes and modify attributes that contain paths
    dataset.attrs = adapt_attributes(dataset.attrs, IGNORE_GLOBAL_ATTRIBUTES,
                                     recipe_name)
    for var in dataset:
        dataset[var].attrs = adapt_attributes(dataset[var].attrs,
                                              IGNORE_VARIABLE_ATTRIBUTES,
                                              recipe_name)

    return dataset


def debug_nc(reference_file: Path, current_file: Path) -> str:
    """Find out the differences between two NetCDF files."""
    ref = load_nc(reference_file)
    cur = load_nc(current_file)

    if diff := diff_dataset(ref, cur):
        msg = diff
    else:
        msg = "Unknown difference"

    return msg


def debug_txt(reference_file: Path, current_file: Path) -> str:
    """Find out the differences between two text files."""
    with reference_file.open('rt') as file:
        ref = file.readlines()
    with current_file.open('rt') as file:
        cur = file.readlines()

    msg = difflib.unified_diff(
        cur,
        ref,
        fromfile=str(current_file),
        tofile=str(reference_file),
    )
    return "".join(msg)


debug_csv = debug_txt
debug_html = debug_txt
debug_grid = debug_txt


def compare_nc(reference_file: Path, current_file: Path) -> bool:
    """Compare two NetCDF files."""
    ref = load_nc(reference_file)
    cur = load_nc(current_file)

    return cur.identical(ref)


def compare_png(reference_file: Path, current_file: Path) -> bool:
    """Compare two PNG files."""
    # Based on:
    # https://scitools-iris.readthedocs.io/en/latest/developers_guide/contributing_graphics_tests.html

    # Perceptual hash size.
    hash_size = 16
    # Maximum perceptual hash hamming distance.
    max_distance = 2

    with Image.open(reference_file) as img:
        ref = imagehash.phash(img, hash_size=hash_size)
    with Image.open(current_file) as img:
        cur = imagehash.phash(img, hash_size=hash_size)

    distance = ref - cur
    return distance < max_distance


def debug(reference_file: Path, current_file: Path) -> str:
    """Try to find out why two files are different."""
    suffix = reference_file.suffix
    fn_name = f"debug_{suffix[1:]}"
    if fn_name in globals():
        debug_fn = globals()[fn_name]
        msg = debug_fn(reference_file, current_file)
    else:
        msg = ""
    return indent(msg, "  ")


def files_equal(reference_file: Path, current_file: Path) -> bool:
    """Compare two files."""
    suffix = reference_file.suffix[1:].lower()
    fn_name = f"compare_{suffix}"
    if fn_name in globals():
        compare_fn = globals()[fn_name]
        same = compare_fn(reference_file, current_file)
    else:
        same = filecmp.cmp(reference_file, current_file, shallow=False)
    return same


def compare_files(reference_dir: Path, current_dir: Path, files: list[Path],
                  verbose: bool) -> list[str]:
    """Compare files from the reference dir to the current dir."""
    different = []
    for file in files:
        ref_file = reference_dir / file
        cur_file = current_dir / file
        if not files_equal(ref_file, cur_file):
            msg = str(file)
            if verbose:
                if info := debug(ref_file, cur_file):
                    msg += ":\n" + info
            different.append(msg)
    return different


def compare(reference_dir: Optional[Path], current_dir: Path,
            verbose: bool) -> bool:
    """Compare a recipe run to a reference run.

    Returns True if the runs were identical, False otherwise.
    """
    recipe_name = get_recipe_name_from_dir(current_dir)
    print("")
    print(f"recipe_{recipe_name}.yml: ", end='')
    if reference_dir is None:
        print("no reference run found, unable to check")
        return False

    result = []

    reference_files = find_files(reference_dir)
    current_files = find_files(current_dir)

    if missing_files := sorted(reference_files - current_files):
        result.append("Missing files:")
        result.extend(f"  - {f}" for f in missing_files)

    if extra_files := sorted(current_files - reference_files):
        result.append("Extra files:")
        result.extend(f"  - {f}" for f in extra_files)

    if differing_files := compare_files(
            reference_dir,
            current_dir,
            sorted(reference_files & current_files),
            verbose,
    ):
        result.append("Differing files:")
        result.extend(indent(f"- {f}", "  ") for f in differing_files)

    if not result:
        print("OK")
        return True

    result.insert(0, "results differ from reference run")
    result.insert(1, f"Reference run: {reference_dir}")
    result.insert(2, f"Current run: {current_dir}")
    print("\n".join(result))
    return False


def get_recipe_name_from_dir(recipe_dir: Path) -> str:
    """Extract recipe name from output dir."""
    recipe_match = re.search(RECIPE_DIR_PATTERN, recipe_dir.stem)
    return recipe_match['recipe_name']


def get_recipe_name_from_file(filename: Path) -> str:
    """Extract recipe name from arbitrary recipe output file."""
    # Iterate starting from the root dir to avoid false matches
    for parent in list(filename.parents)[::1]:
        recipe_match = re.search(RECIPE_DIR_PATTERN, str(parent))
        if recipe_match is not None:
            return recipe_match['recipe_name']
    raise ValueError(f"Failed to extract recipe name from file {filename}")


def get_recipe_dir_from_str(str_in: str, recipe_name: str) -> Optional[Path]:
    """Try to extract recipe directory from arbitrary string."""
    recipe_dir_pattern = (
        rf'recipe_{recipe_name}_' + RECIPE_DIR_DATETIME_PATTERN
    )
    recipe_dir_match = re.search(recipe_dir_pattern, str_in)

    # If recipe directory is not found in string, return None
    if recipe_dir_match is None:
        return None

    # If recipe directory is found, return entire parent directory
    # E.g., for str_in = /root/path/recipe_test_20220202_222222/work/file.nc
    # return /root/path/recipe_test_20220202_222222
    # For this, iterate from the right (::-1) through the parents
    for parent in Path(str_in).parents[::-1]:
        if recipe_dir_match[0] in str(parent):
            return parent

    # If no valid parent is found, return str_in
    # E.g., for str_in = /root/path/recipe_test_20220202_222222 no valid
    # parents are found; thus, return str_in itself
    return Path(str_in)


def find_files(recipe_dir: Path) -> set[Path]:
    """Find all NetCDF files in a recipe run directory."""
    result: set[Path] = set()
    for subdir in COMPARE_SUBDIRS:
        for root, _, found in os.walk(recipe_dir / subdir):
            files = set(found)
            for ignore_pattern in IGNORE_FILES:
                files -= set(fnmatch.filter(files, ignore_pattern))
            parent = Path(root).relative_to(recipe_dir)
            for file in files:
                result.add(parent / file)

    return result


def find_successful_runs(dirname: Path, recipe_name: str = '*') -> list[Path]:
    """Find recipe runs in `dirname`.

    `dirname` can either be a recipe run or a directory containing
    recipe runs.
    """
    runs = []
    for recipe_file in sorted(
            list(dirname.glob(f"run/recipe_{recipe_name}.yml")) +
            list(dirname.glob(f"*/run/recipe_{recipe_name}.yml"))):
        recipe_dir = recipe_file.parent.parent
        log = recipe_dir / 'run' / 'main_log.txt'
        success = log.read_text().endswith('Run was successful\n')
        if success:
            runs.append(recipe_dir)
    return sorted(set(runs))


def find_recipes(reference: Path,
                 current: list[Path]) -> Iterator[tuple[Path, Optional[Path]]]:
    """Yield tuples of current and reference directories."""
    for current_dir in current:
        for recipe_dir in find_successful_runs(current_dir):
            recipe_name = get_recipe_name_from_dir(recipe_dir)
            reference_dirs = find_successful_runs(reference, recipe_name)
            if reference_dirs:
                reference_dir = reference_dirs[-1]
            else:
                reference_dir = None
            yield (recipe_dir, reference_dir)


def main() -> int:
    """Compare recipe runs.

    Returns 0 if all recipe runs were successfully compared and
    identical to the reference, 1 otherwise.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'reference',
        default='.',
        type=Path,
        help='Directory containing results obtained with reference version.',
    )
    parser.add_argument(
        'current',
        default='.',
        nargs='+',
        type=Path,
        help=("List of recipe results or directories containing such "
              "results, obtained with current version."),
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action="store_true",
        help="Display more information on differences.",
    )

    args = parser.parse_args()

    print("Comparing recipe run(s) in:\n{}".format('\n'.join(
        str(f) for f in args.current)))
    print(f"to reference in {args.reference}")
    fail = []
    success = []
    for current_dir, reference_dir in find_recipes(args.reference,
                                                   args.current):
        same = compare(reference_dir, current_dir, verbose=args.verbose)
        recipe = f"recipe_{get_recipe_name_from_dir(current_dir)}.yml"
        if same:
            success.append(f"{recipe}:\t{current_dir}")
        else:
            fail.append(f"{recipe}:\t{current_dir}")

    # Print summary of results to screen
    summary = ["", "Summary", "======="]
    if success:
        summary.extend([
            "", "The following recipe runs were identical to previous runs:",
            *success
        ])
    if fail:
        summary.extend([
            "", "The following recipe runs need to be inspected by a human:",
            *fail
        ])
    print("\n".join(summary))
    print("")

    if fail:
        print(f"Action required: {len(fail)} out of {len(success) + len(fail)}"
              " recipe runs need to be inspected by a human.")
    else:
        print(f"All {len(success)} recipe runs were identical.")

    return int(bool(fail))


if __name__ == '__main__':
    sys.exit(main())
