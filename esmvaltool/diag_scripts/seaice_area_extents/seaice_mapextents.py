"""diagnostic script to plot extent and differences 

based on code from Anton Steketee's COSIMA cookbook notebook
https://cosima-recipes.readthedocs.io/en/latest/DocumentedExamples
    /SeaIce_Obs_Model_Compare.html
"""

import logging
import os
import calendar
from cartopy import crs
import xarray as xr
import xesmf
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd

from esmvaltool.diag_scripts.shared import run_diagnostic, save_figure
from esmvaltool.diag_scripts.shared._base import get_plot_filename


# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def map_diff(mod_si_ls, obs_si, months):
    "create figure mapping extents for models and months"
    # get lat max for regridding
    latmax = obs_si.lat.max().values.item()

    # fig set up, width for 2 models, check len mod_si_ls
    figure = plt.figure(figsize = (9, len(months) * 4))
    j = 0  # to iterate through positions on figure

    for mon in months:
        cdr = obs_si.siconc.sel(time=obs_si.siconc.time.dt.month.isin(mon)
                                ).mean('time')
        i = 1
        for mod_label, mod_si in mod_si_ls.items():

            lat_min = mod_si.lat.min().values.item()
            mod_si = mod_si.where(mod_si.lat < latmax, drop=True)
            regrd_out = obs_si.where(obs_si.lat > lat_min, drop=True)

            # regrid model data to observations
            regridder_access_sh = xesmf.Regridder(
                mod_si.isel(time=0).drop(['i','j']),
                regrd_out.isel(time=0),
                'bilinear',
                periodic=True,
                unmapped_to_nan=True
            )

            model_mean = mod_si.siconc.sel(time=mod_si.siconc.time.dt.month.isin(mon)
                                   ).mean('time')
            mod_regrid = regridder_access_sh(model_mean)
            diff_ds = mod_regrid - cdr

            axes = plt.subplot(len(months), 3, i + j * 3, projection=
                             crs.SouthPolarStereo(true_scale_latitude=-70))

            diffmap = axes.contourf(
                diff_ds.x, diff_ds.y, diff_ds,
                levels=np.arange(-90,91,20), cmap='RdBu'
            )
            cs_cdr = cdr.plot.contour(levels=[15], ax=axes)
            cs_mod = mod_regrid.plot.contour(levels=[15], ax=axes, 
                                             colors=['black'])

            plt.title(calendar.month_abbr[mon] + ' ' + mod_label)

            i += 1
        j += 1

    color_cdr = cs_cdr.collections[0].get_edgecolor()
    line_cdr = mlines.Line2D([], [], color=color_cdr, label="Observed Extent")

    color_mod = cs_mod.collections[0].get_edgecolor()
    line_mod = mlines.Line2D([], [], color=color_mod, label="Modelled Extent")

    plt.legend(handles=[line_cdr,line_mod], loc='center left', 
               bbox_to_anchor=(1.2,0.5))
    cax = plt.axes([0.7,0.55,0.04,0.3])
    _ = plt.colorbar(diffmap, cax=cax, 
                     label='Difference in \nSea Ice Concentration')
    return figure


def main(cfg):
    """Compute."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()
    data = []

    # Find input datasets to use
    for dataset in input_data:
        
        input_file = [dataset['filename'], dataset['dataset']]
        # drop areacello dataset for map
        if dataset['short_name'] == 'siconc':
            data.append(input_file)

    inputs_df = pd.DataFrame(data, columns=['filename', 'dataset'])

    logger.info(inputs_df)
    mod_si_dict = {}
    for filepath, data_name in inputs_df.itertuples(index=False):
        if data_name == 'NSIDC-G02202-sh':
            # Load the data
            obs_si = xr.open_dataset(filepath)
        else:
            mod_si_dict[data_name] = xr.open_dataset(filepath)

    logger.info("creating map differences")
    mapfig = map_diff(mod_si_dict, obs_si, cfg['months'])
    # Save output
    output_path = get_plot_filename('map_difference', cfg)
    # mapfig.savefig(output_path) # use esmvaltool convenience function

    provenance_record = get_provenance_record(inputs_df['filename'].to_list())
    save_figure(output_path, provenance_record, cfg, figure=mapfig)


def get_provenance_record(ancestor_files):
    "build provenance dictionary"
    record = {
        'ancestors': ancestor_files,
        'authors': [
            'chun_felicity',
        ],
        'caption': '',
        'domains': ['shpolar'],
        'plot_types': ['polar'],
        'references': [],
        'statistics': ['diff'],
            }
    return record


if __name__ == '__main__':


    with run_diagnostic() as config:
        main(config)
