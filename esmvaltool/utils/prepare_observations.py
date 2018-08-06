"""CMORize observational datasets using CIS."""
import argparse
import logging
import os

import cis

logger = logging.getLogger('esmvaltool.utils.prepare_observations')


def convert(filenames, output_dir, suffix=''):
    """Read data from filenames and save it to output/variable_suffix.nc"""
    logger.info("Reading variables from files %s", filenames)
    for var_name in cis.get_variables(filenames):
        logger.info("Reading variable %s", var_name)
        cube = cis.read_data(filenames=filenames, variable=var_name)

        output_file = var_name
        if suffix:
            output_file += '_' + suffix
        output_file += '.nc'
        output_file = os.path.join(output_dir, output_file)
        logger.info("Saving data to %s", output_file)
        cube.save_data(output_file)

    logger.info("Done")


def _main():
    """Run as a program."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('filenames', nargs='+', help='Input data files')
    parser.add_argument(
        '-o', '--output-dir', help='Output directory', default=os.getcwd())
    parser.add_argument(
        '-s', '--file-suffix', help='Output file suffix', default='')
    parser.add_argument(
        '-l',
        '--log-level',
        default='info',
        choices=['debug', 'info', 'warning', 'error'])
    args = parser.parse_args()

    logging.basicConfig(
        format="%(asctime)s %(levelname)-8s %(name)s,%(lineno)s\t%(message)s")
    logging.getLogger().setLevel(args.log_level.upper())

    # Prepare input arguments
    filenames = [os.path.abspath(p) for p in args.filenames]
    output_dir = os.path.abspath(args.output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Run
    convert(filenames, output_dir, args.file_suffix)


if __name__ == '__main__':
    _main()
