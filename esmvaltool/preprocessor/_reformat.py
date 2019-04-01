"""Simple interface to reformat and CMORize functions."""
import os

import iris

from ..cmor.check import cmor_check, cmor_check_data, cmor_check_metadata
from ..cmor.fix import fix_data, fix_file, fix_metadata
from ._io import save


def cmor_fix_fx(fx_files_dict, variable):
    """Fix all fx files and save them to disk."""
    fixed_fx_files = {}
    for fx_var, fx_file in fx_files_dict.items():
        fx_cube = iris.load_cube(fx_file)
        fix_metadata([fx_cube], fx_var, variable['project'],
                     variable['dataset'], variable['cmor_table'],
                     'fx', 'fx')
        cmor_check_metadata(fx_cube, variable['cmor_table'],
                            'fx', fx_var, 'fx')
        filename = os.path.basename(fx_file)
        fx_save_dir = os.path.dirname(variable['filename'])
        save([fx_cube], os.path.join(fx_save_dir, filename))
        fixed_fx_files[fx_var] = os.path.join(fx_save_dir, filename)
    return fixed_fx_files
