"""
Apply FX

Module that applies an fx-type mask:
_apply_fx(cube, fx_file, fx_option) takes cube, fx_file
and fx_option = 'land' (mask land out) or 'ocean' (mask
ocean out). Returns the masked cube
"""
import iris
import numpy as np


def _get_fx_mask(fx_data, fx_option):
    """Build a 50 percent land or ocean mask"""
    inmask = np.zeros_like(fx_data, bool)
    if fx_option == 'land':
        # Mask land out
        inmask[fx_data <= 50.] = False
        inmask[fx_data > 50.] = True
    elif fx_option == 'ocean':
        # Mask ocean out
        inmask[fx_data <= 50.] = True
        inmask[fx_data > 50.] = False
    return inmask


def _apply_fx_mask(fx_mask, var_data):
    """Apply the fx mask"""
    var_mask = np.zeros(var_data.shape, bool)

    # TIME-LAT-LON
    if len(var_data.shape) == 3.:
        for i in range(var_data.shape[0]):
            var_mask[i, :] = fx_mask

    # TIME-PLEV-LAT-LON
    elif len(var_data.shape) == 4.:
        for i in range(var_data.shape[0]):
            for j in range(var_data.shape[1]):
                var_mask[i, j, :] = fx_mask

    # Aplly mask accross
    if np.ma.is_masked(var_data):
        var_mask |= var_data.mask

    # Build the new masked data
    var_data = np.ma.array(var_data,
                           mask=var_mask,
                           fill_value=1e+20)

    return var_data


def _apply_fx(cube, fx_file, fx_option):
    """Apply a land/ocean mask"""
    # fx_option is either 'land' or 'ocean'
    fx_cube = iris.load_cube(fx_file)
    landocean_mask = _get_fx_mask(fx_cube.data, fx_option)
    cube.data = _apply_fx_mask(landocean_mask, cube.data)

    return cube
