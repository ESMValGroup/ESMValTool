"""Tool for testing ESMValTool."""
import argparse
import copy
import os
import subprocess
from itertools import product
from pathlib import Path

import yaml


def absolute(path):
    """Make path into an absolute Path object."""
    return Path(os.path.abspath(path))


def linear_expand(filename, cwd):
    """Create recipes from filename using the recipe options provided.

    Uses one option at a time.
    """
    filename = Path(filename)
    yield filename

    options_file = Path(__file__).parent / 'options.yml'
    options = yaml.safe_load(options_file.read_text()).get(filename.name)

    recipe = yaml.safe_load(filename.read_text())

    for key in options or {}:
        for value in options[key]:
            outrecipe = copy.deepcopy(recipe)
            write_recipe = False
            for diag_name, diagnostic in recipe['diagnostics'].items():
                for script_name, script in diagnostic['scripts'].items():
                    if key in script and script[key] != value:
                        write_recipe = True
                        outrecipe['diagnostics'][diag_name]['scripts'][
                            script_name][key] = value
            if write_recipe:
                outfile = cwd / Path('{}_{}_{}.yml'.format(
                    filename.stem, key,
                    str(value).replace(os.sep, '-')))
                print("Creating", outfile)
                outfile.write_text(yaml.safe_dump(outrecipe))
                yield outfile


def matrix_expand(filename, cwd, max_recipes=100):
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
                outfile = cwd / Path(outfile + '.yml')
                print("Creating", outfile)
                outfile.write_text(yaml.safe_dump(outrecipe))
                yield outfile


def create_script(recipe, config_file, cwd):
    """Submit a job for recipe."""
    job_template = Path(__file__).parent / 'job.sh.template'
    job = job_template.read_text().format(
        recipe=recipe,
        config=config_file,
    )

    jobfile = cwd / Path(recipe.stem + '_' + job_template.stem)
    jobfile.write_text(job)
    return jobfile


def run(script, cwd, method=''):
    """Run script in cwd using method."""
    if method == 'bsub':
        print("Submitting", script, 'in', cwd)
        with open(script) as stdin:
            subprocess.run('bsub', stdin=stdin, cwd=cwd, check=True)
    elif method == 'dry-run':
        print("Would run", script, 'in', cwd)
    else:
        print("Running", script, 'in', cwd)
        subprocess.run(['bash', str(script)], cwd=cwd, check=True)


def install(args):
    """Install ESMValTool from GitHub."""
    cwd = absolute(args.directory)
    cwd.mkdir(parents=True, exist_ok=True)
    script_template = Path(__file__).parent / 'install.sh.template'
    script = script_template.read_text().format(branch=args.branch)

    script_file = cwd / Path(args.branch + '_' + script_template.stem)
    script_file.write_text(script)
    run(script_file, cwd=cwd, method=args.run_method)


def schedule(args):
    """Create recipes with the options provided and schedule."""
    cwd = absolute(args.directory)
    expand = matrix_expand if args.matrix else linear_expand
    for input_recipe in args.recipes:
        input_recipe = absolute(input_recipe)
        for recipe in expand(input_recipe, cwd=cwd):
            script_file = create_script(
                recipe,
                config_file=absolute(args.esmvaltool_config_file),
                cwd=cwd,
            )
            run(script_file, cwd=cwd, method=args.run_method)


def main():
    """Run the program."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-d',
                        '--directory',
                        default='.',
                        help='Use as a working directory.')
    parser.add_argument('-r',
                        '--run-method',
                        default='immediate',
                        choices=['immediate', 'bsub', 'dry-run'],
                        help='Choose an execution method.')
    subparsers = parser.add_subparsers()
    parser.set_defaults(function=lambda _: parser.print_help())

    install_parser = subparsers.add_parser('install')
    install_parser.add_argument('branch',
                                help='Name of the GitHub branch to install.')
    install_parser.set_defaults(function=install)

    schedule_parser = subparsers.add_parser('schedule')
    schedule_parser.add_argument('recipes',
                                 nargs='+',
                                 help='Path to the recipe files to run.')
    schedule_parser.add_argument(
        '-c',
        '--esmvaltool-config-file',
        help='Path to the ESMValTool configuration file.')
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
