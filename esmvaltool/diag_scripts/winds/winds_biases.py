import cartopy.crs as ccrs 
import esmvalcore.preprocessor as eprep
import iris
import iris.plot as iplt
import logging
import matplotlib.pyplot as plt
import os
import numpy as np
import glob

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic, select_metadata, group_metadata, get_diagnostic_filename, save_data, ProvenanceLogger
import esmvaltool.diag_scripts.shared.plot as eplot
from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import ProvenanceLogger

# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
# logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def main(cfg):

    if int(cfg['hours'][:-2])<25:
        mean_levs = np.arange(-15, 1)
        std_levs =  np.arange(0, 26)
    else:
        mean_levs = np.arange(-25, 1)
        std_levs =  np.arange(0, 36)

    input_data = cfg['input_data']

    datasets = group_metadata(input_data.values(), 'dataset', sort=True)

    for dataset in datasets.keys():
        vars_data = datasets[dataset]
        try:
            wind_from_h = select_metadata(vars_data, variable_group='derived')[0]['filename']
            wind_from_h_cb = iris.load_cube(wind_from_h)
        except:
            uas = select_metadata(vars_data, variable_group='uas')[0]['filename']
            uas_cb = iris.load_cube(uas)
            vas = select_metadata(vars_data, variable_group='vas')[0]['filename']
            vas_cb = iris.load_cube(vas)
            wind_cblst = iris.cube.CubeList([uas_cb, vas_cb])
            wind_from_h_cb = eprep._derive.sfcwind.DerivedVariable().calculate(wind_cblst)
            # wind_from_h_cb = (uas_cb**2 + vas_cb**2)**0.5
        wind_h = select_metadata(vars_data, variable_group='real')[0]['filename']
        wind_h_cb = iris.load_cube(wind_h)

        y = str(input_data[wind_h]['start_year'])

        try:
            diff_cb = (wind_from_h_cb - wind_h_cb)*100/wind_h_cb
        except:
            diff_cb = (wind_from_h_cb - wind_h_cb.data)*100/wind_h_cb.data
        
        mean_diff = diff_cb.collapsed('time', iris.analysis.MEAN)
        std_diff = diff_cb.collapsed('time', iris.analysis.STD_DEV)

        fig_mean, ax_mean = plt.subplots(ncols=1, nrows=1,
                                        subplot_kw={'projection': ccrs.Robinson()})
        fig_mean.set_size_inches(12, 6); fig_mean.set_dpi(300)
        levs = iplt.contourf(mean_diff,axes=ax_mean, levels=mean_levs, cmap='Blues_r', extend='both')
        ax_mean.coastlines(linewidth=0.5)
        plt.tight_layout()
        fig_mean.colorbar(levs)
        fig_mean.suptitle(cfg['hours']+' '+dataset+' mean wind speed difference in '+y, fontsize='xx-large')
        fig_mean.savefig(os.path.join(cfg['plot_dir'], 'mean_diff_'+dataset+'_'+y + diagtools.get_image_format(cfg)))
        
        fig_std, ax_std = plt.subplots(ncols=1, nrows=1,
                                        subplot_kw={'projection': ccrs.Robinson()})
        fig_std.set_size_inches(12, 6); fig_std.set_dpi(300)
        
        levels = iplt.contourf(std_diff, axes=ax_std, levels=std_levs, cmap='Reds', extend='both')
        ax_std.coastlines(linewidth=0.5)
        plt.tight_layout()
        fig_std.colorbar(levels, label='Difference, %')
        fig_std.suptitle(cfg['hours']+' STd of '+dataset+' wind speed difference in '+y, fontsize='xx-large')
        fig_std.savefig(os.path.join(cfg['plot_dir'], 'std_diff_'+dataset+'_'+y + diagtools.get_image_format(cfg)))

    logger.info('Success')


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)
