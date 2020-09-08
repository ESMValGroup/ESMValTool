"""Rename preprocessor output files so they are named according to OBS6."""

import logging
import shutil
from pathlib import Path

import iris
from esmvalcore.cmor.table import CMOR_TABLES

from esmvaltool.diag_scripts.shared import (get_diagnostic_filename,
                                            run_diagnostic)

logger = logging.getLogger(Path(__file__).name)


def main(cfg):
    """Rename preprocessed native6 file."""
    fixed_files = cfg['input_data']

    for file, info in fixed_files.items():
        stem = Path(file).stem
        basename = stem.replace('native', 'OBS')

        if info['diagnostic'] == 'daily':
            for mip in ['day', 'Eday', 'CFday']:
                if CMOR_TABLES['CMIP6'].get_variable(mip, info['short_name']):
                    basename = basename.replace('E1hr', mip)
            basename = basename.replace('E1hr', 'day')

        cube = iris.load_cube(file)
        try:
            time = cube.coord('time')
        except iris.exceptions.CoordinateNotFoundError:
            pass
        else:
            if info['diagnostic'] == "monthly":
                start = time.cell(0).point.strftime("%Y%m")
                end = time.cell(-1).point.strftime("%Y%m")
            else:
                start = time.cell(0).point.strftime("%Y%m%d")
                end = time.cell(-1).point.strftime("%Y%m%d")
            basename = f"{basename.rstrip('0123456789-')}{start}-{end}"

        outfile = get_diagnostic_filename(basename, cfg)
        logger.info('Moving %s to %s', file, outfile)
        shutil.move(file, outfile)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
