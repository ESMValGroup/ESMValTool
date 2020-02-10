"""native6 diagnostic."""

import logging
import shutil
from pathlib import Path

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
            basename = basename.replace('E1hr', 'Eday')

        if 'fx' not in basename:
            end_year = basename[-4:]
            basename = basename.replace(end_year, f'{int(end_year) - 1}')

        outfile = get_diagnostic_filename(basename, cfg)
        logger.info('Moving %s to %s', file, outfile)
        shutil.move(file, outfile)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
