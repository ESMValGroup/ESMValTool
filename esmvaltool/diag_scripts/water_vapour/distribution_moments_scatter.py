import logging
import iris
import numpy as np
import matplotlib.pyplot as plt
from iris.util import equalise_attributes
import esmvalcore.preprocessor as eprep
import os

from fastkde import fastKDE

from esmvaltool.diag_scripts.water_vapour import distribution_moments_slopes as dist_slps
from esmvaltool.diag_scripts.water_vapour.distribution_moments_maps import reform_inp_data

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic, select_metadata, group_metadata, get_diagnostic_filename, save_data, ProvenanceLogger
import esmvaltool.diag_scripts.shared.plot as eplot
from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import ProvenanceLogger

# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def contour_levels_2d(mypdf,v1,v2,prob=0.95,nlevels=4):
    """
    This function determines the level values that are used 
    to plot two-dimensional contour estimates. The function
    determines the levels such that low probability noise is 
    filtered through the value in 'prob' which determines how 
    much data is encompassed by the contours (0.95 seems to work)
    relatively generically but can be set through the keys. 

    NOTE: The function should generally work for any equally spaced
          probability density estimate. However, it has been only been 
          tested extensively with the fastKDE python package.

    :mypdf:     2-d field of probability density function
    :v1:        vector of values in x direction
    :v2:        vector of values in y direction
    :prob:      fraction of data that is encompassed by the outermost
                contour level
    :nlevels:   number of contour levels
    """

    # Define the space of each grid-point for cumulative integration
    area=(v1[1]-v1[0])*(v2[1]-v2[0])
    # Retrieve the unique likelihood values and there respective counts
    unique,counts=np.unique(mypdf,return_counts=True) 
    # Flip the arrays in order to integrate from largest likelihood to lowest
    unique = np.flip(unique)
    counts = np.flip(counts)
    #Cumulative sum
    cum_sum = 0 
    # Integer to loop through unique likelihood values
    k=0
    # Loop until 'prob' percent is encompassed by the contours
    while cum_sum<prob:
        # Integrate probability density function from largest to lowest probablity
        cum_sum+=unique[k]*counts[k]*area
        k+=1

    # Determine the levels for the contours  
    clevels = np.linspace(unique[k],unique[0],nlevels+1)[0:-1]
    
    return clevels 

def plot_scatter(dtst_dic, dataset, cfg):

    st_file = eplot.get_path_to_mpl_style(cfg.get('mpl_style'))
    plt.style.use(st_file)

    nseas = len(dtst_dic.keys())
    ncols = 1

    fig_scatter, ax_scatter = plt.subplots(ncols=ncols, nrows=nseas, squeeze=False,
                                            sharex=True, sharey=True)
    fig_scatter.set_size_inches(ncols*3 +1, nseas*3)
    fig_scatter.set_dpi(300)

    for n, seas in enumerate(dtst_dic.keys()):
        mean = dtst_dic[seas]['mean'].data.ravel()
        skew = dtst_dic[seas]['skew'].data.ravel()
        ax_scatter[n, 0].scatter(mean, skew, c='silver',s=0.1)
        mypdf,axis = fastKDE.pdf(mean,skew); v1,v2=axis
        ax_scatter[n, 0].contour(v1, v2, mypdf, colors='k',linewidths=0.6,
                      levels=contour_levels_2d(mypdf,v1,v2,prob=0.9,nlevels=3))
        coords = dist_slps.pca_axes(mean.reshape(-1,1), skew.reshape(-1,1))
        ax_scatter[n, 0].plot(coords['pca1']['v1'],coords['pca1']['v2'],c='red',linewidth=1)
        ax_scatter[n, 0].set_xlim(0,70)
        ax_scatter[n, 0].set_ylim(-2.2,2.2)
        ax_scatter[n, 0].set_ylabel('skewness (kg/m$^2$)')
        ax_scatter[n, 0].text(0.075, 0.075, seas, transform=ax_scatter[n, 0].transAxes)
        if n == nseas - 1:
            ax_scatter[n, 0].set_xlabel('mean (kg/m$^2$)')


    fig_scatter.suptitle('TCWV ' +dataset)
    plt.tight_layout()
    fig_scatter.savefig(os.path.join(cfg['plot_dir'], 'scatter_'+dataset+ diagtools.get_image_format(cfg)))
    plt.close(fig_scatter)   

    return


def main(cfg):

    input_data = cfg['input_data']; reform_inp_data(input_data)

    datasets = group_metadata(input_data.values(), 'dataset', sort=True)

    for dataset in datasets.keys():
        dtst_dic = {} 
        groups = group_metadata(datasets[dataset], 'variable_group')
        for group in groups.keys():
            filepaths = list(group_metadata(groups[group], 'filename').keys())
            n_real = len(filepaths)
            if n_real == 1: 
                mod_cb = iris.load_cube(filepaths[0])
                mod_weight = np.ones(mod_cb.shape)
                dtst_dic[group] = dist_slps.calculate_single_cube_stats(mod_cb)
            else:
                mod_cblst = iris.load(filepaths)
                equalise_attributes(mod_cblst)
                [cb.add_aux_coord(iris.coords.AuxCoord(n, long_name='n_order',
                        var_name='n_order')) for n, cb in enumerate(mod_cblst)]
                mod_cb = mod_cblst.merge_cube()
                mod_weight = np.ones(mod_cb.shape)/n_real
                dtst_dic[group] = dist_slps.calculate_cube_stats(mod_cb, mod_weight)
        plot_scatter(dtst_dic, dataset, cfg)

    logger.info('Success')

if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)