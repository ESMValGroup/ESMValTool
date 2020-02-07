"""native6 diagnostic."""

import logging
import os
import shutil

from esmvaltool.diag_scripts.shared import (get_diagnostic_filename,
                                            run_diagnostic)

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    """Rename preprocessed native6 file."""
    fixed_files = cfg['input_data']

    for file, info in fixed_files.items():
        basename, _ext = os.path.splitext(os.path.basename(file))
        basename = basename.replace('native', 'OBS')


        if info['diagnostic'] == 'daily':
            basename = basename.replace('E1hr', 'Eday')

        end_year = basename[-4:]
        basename = basename.replace(end_year, f'{int(end_year) - 1}')

        outfile = get_diagnostic_filename(basename, cfg)
        logger.info('Moving %s to %s', file, outfile)
        shutil.move(file, outfile)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
