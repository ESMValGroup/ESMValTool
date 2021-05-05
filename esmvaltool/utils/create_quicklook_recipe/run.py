"""Command line tool to create quicklook recipes."""

import argparse
import logging
import time
from pathlib import Path

import yaml

logging.basicConfig(format="%(asctime)s %(levelname)-8s\t%(message)s")
logging.Formatter.converter = time.gmtime

logger = logging.getLogger('create-quicklook-recipes')


def create_recipe(args):
    """Get desired single plot diagnostics."""
    # Get Header
    templates_dir = Path(args.recipe_templates_dir)
    with (templates_dir / 'general.yml').open() as stream:
        recipe = yaml.load(stream, Loader=yaml.FullLoader)

    # Get plot script sections
    with (templates_dir / 'plot_scripts.yml').open() as stream:
        plot_scripts = yaml.load(stream, Loader=yaml.FullLoader)
    for script_body in plot_scripts.values():
        script_body['results_dir'] = args.cache_dir

    run_ids = args.simulations
    recipe['diagnostics'] = {}
    if len(run_ids) > 1:
        for script_body in plot_scripts.values():
            script_body['multi_dataset_plot'] = True
            script_body['read_all_available_datasets'] = True
            patterns = [f'{run_id}_*' for run_id in run_ids]
            script_body['patterns'] = patterns
        recipe['diagnostics']['multi_run_plots'] = {
            'description': 'Plot multiple runs in one plot',
            'scripts': plot_scripts,
        }
    else:
        for recipe_type in args.templates:
            recipe_name = f'diagnostics_{recipe_type}.yml'
            with (templates_dir / recipe_name).open() as stream:
                diagnostics = yaml.load(stream, Loader=yaml.FullLoader)
            for (diag_name, diag_content) in diagnostics.items():
                diag_name = diag_name.format(run_id=run_ids[0])
                diag_content['scripts'] = plot_scripts
                diag_content['additional_datasets'] = [
                    {
                        'dataset': run_ids[0],
                        'project': 'EMAC',
                        'start_year': args.start_year,
                        'end_year': args.end_year,
                    },
                ]
                recipe['diagnostics'][diag_name] = diag_content
    return recipe


def parse_args():
    """Define the `create_quicklook_recipe` command line options."""
    # parse command line args
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'simulations',
        nargs='+',
        type=str,
        help='Set identifier for individual simulation(s). '
        'If multiple simulations are given, produce multi-run '
        'plots for all given runs.',
    )
    parser.add_argument(
        '--cache-dir',
        required=True,
        type=str,
        help='Path to store re-usable diagnostic results.',
    )
    parser.add_argument(
        '--templates',
        nargs='*',
        choices=(
            'atmos_chem',
            'atmos_dynamics',
            'forcings',
            'ocean',
        ),
        help='The names of the templates to use to fill the recipe.',
    )
    parser.add_argument(
        '--recipe-templates-dir',
        default=Path(__file__).absolute().parent,
        type=str,
        help='Path to recipe template files',
    )
    parser.add_argument(
        '--start-year',
        type=int,
        help='Set start year. Only necessary if exactly one "simulation" is '
        'given.',
    )
    parser.add_argument(
        '--end-year',
        type=int,
        help='Set end year. Only necessary if exactly one "simulation" is '
        'given',
    )
    parser.add_argument(
        '--recipe-dir',
        default='.',
        type=str,
        help='Path to write the recipe to.',
    )
    parser.add_argument(
        '-l',
        '--log-level',
        default='info',
        help=("Set the log-level"),
        choices=['debug', 'info', 'warning', 'error'],
    )
    args = parser.parse_args()

    if args.log_level:
        logger.setLevel(args.log_level.upper())

    if len(args.simulations) < 2:
        for arg in ('start_year', 'end_year'):
            if getattr(args, arg) is None:
                raise ValueError(
                    f"Argument '--{arg}' is necessary if a single "
                    f"'simulation' is specified.")
    return args


def main():
    """Create recipe for quicklook."""
    args = parse_args()
    logger.info("Creating quicklook recipe")

    # Get diagnostics
    recipe = create_recipe(args)

    # Write recipe
    recipe_dir = Path(args.recipe_dir)
    if not recipe_dir.is_dir():
        logger.info("Creating non-existent recipe directory %s", recipe_dir)
        recipe_dir.mkdir(parents=True)

    recipe_file = recipe_dir / 'recipe_{}_quicklook.yml'.format('_'.join(
        args.simulations))
    logger.info("Writing quicklook recipe %s", recipe_file)
    with recipe_file.open('w') as stream:
        stream.write(yaml.dump(recipe))

    logger.info("Done")


if __name__ == '__main__':
    main()
