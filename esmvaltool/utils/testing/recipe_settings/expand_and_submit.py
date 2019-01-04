"""Tool for testing ESMValTool."""
import argparse
import copy
import os
import subprocess
from itertools import product
from pathlib import Path
from stat import S_IXUSR

import yaml


def linear_expand(filename):
    """Create recipes from filename using the recipe options provided.

    Uses one option at a time.
    """
    filename = Path(filename)
    yield filename

    options_file = Path(__file__).parent / 'options.yml'
    options_key = filename.name
    options = yaml.safe_load(options_file.read_text())[options_key]

    recipe = yaml.safe_load(filename.read_text())

    for key in options:
        for diag_name, diagnostic in recipe['diagnostics'].items():
            for script_name, script in diagnostic['scripts'].items():
                for value in options[key]:
                    if key in script and script[key] != value:
                        outrecipe = copy.deepcopy(recipe)
                        outrecipe['diagnostics'][diag_name]['scripts'][
                            script_name][key] = value
                        outfile = Path('{}_{}_{}.yml'.format(
                            filename.stem, key,
                            str(value).replace(os.sep, '-')))
                        print("Creating", outfile)
                        outfile.write_text(yaml.safe_dump(outrecipe))
                        yield outfile


def matrix_expand(filename, max_recipes=100):
    """Create recipes from filename using the recipe options provided.

    Tries all possible combinations of options, but stops at max_recipes.
    """
    filename = Path(filename)
    options_file = Path(__file__).parent / 'options.yml'
    options = yaml.safe_load(options_file.read_text())[filename.name]

    recipe = yaml.safe_load(filename.read_text())

    n_recipes = 0
    for diag_name, diagnostic in recipe['diagnostics'].items():
        for script_name, script in diagnostic['scripts'].items():
            outrecipe = copy.deepcopy(recipe)
            keys = list(options)
            for values in product(*[options[k] for k in keys]):
                # Avoid creating a huge number of recipes
                n_recipes += 1
                if n_recipes > max_recipes:
                    print("Warning: stopping early at", max_recipes, "recipes")
                    return

                outfile = filename.stem
                for i, key in enumerate(keys):
                    value = values[i]
                    if key in script:
                        outrecipe['diagnostics'][diag_name]['scripts'][
                            script_name][key] = value
                        outfile += '_' + str(key) + '_' + str(value).replace(
                            os.sep, '-')
                outfile = Path(outfile + '.yml')
                print("Creating", outfile)
                outfile.write_text(yaml.safe_dump(outrecipe))
                yield outfile


def submit(recipe):
    """Submit a job for recipe."""
    recipe = Path(recipe)
    job_template = Path(__file__).parent / 'job.bsub.template'
    job = job_template.read_text().format(recipe)

    jobfile = Path(recipe.stem + '_' + job_template.stem)
    jobfile.write_text(job)

    print("Submitting", jobfile)
    with jobfile.open() as stdin:
        subprocess.run('bsub', stdin=stdin, check=True)


def install(args):
    """Install ESMValTool from GitHub."""
    branch = args.branch
    script_template = Path(__file__).parent / 'install.sh.template'
    script = script_template.read_text().format(branch=branch)

    script_file = Path(branch + '_' + script_template.stem)
    script_file.write_text(script)
    script_file.chmod(script_file.stat().st_mode | S_IXUSR)
    subprocess.run('./{}'.format(script_file), check=True)


def schedule(args):
    """Create recipes with the options provided and schedule."""
    expand = matrix_expand if args.matrix else linear_expand
    for input_recipe in args.recipes:
        for recipe in expand(input_recipe):
            submit(recipe)


def main():
    """Run the program."""
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers()
    parser.set_defaults(function=lambda _: parser.print_help())

    install_parser = subparsers.add_parser('install')
    install_parser.add_argument(
        'branch', help='Name of the GitHub branch to install.')
    install_parser.set_defaults(function=install)

    schedule_parser = subparsers.add_parser('schedule')
    schedule_parser.add_argument(
        'recipes', nargs='+', help='Path to the recipe files to run.')
    schedule_parser.add_argument(
        '-m',
        '--matrix',
        action='store_true',
        help=('Use all possible combinations of options instead of a single '
              'option at a time.'))
    schedule_parser.set_defaults(function=schedule)

    args = parser.parse_args()
    args.function(args)


if __name__ == '__main__':
    main()
