"""native6 diagnostic."""

import os
import shutil

from esmvaltool.diag_scripts.shared import (run_diagnostic,
                                            get_diagnostic_filename)


def main(cfg):
    """Rename preprocessed native6 file."""
    fixed_files = cfg['input_data']
    
    for file, dictionary in fixed_files.items():
        basename, _ext = os.path.splitext(os.path.basename(file))
        basename = basename.replace('native', 'OBS')

        if dictionary['diagnostic'] == 'daily':
            basename = basename.replace('E1hr', 'Eday')
        outfile = get_diagnostic_filename(basename, cfg)
        shutil.move(file, outfile)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)