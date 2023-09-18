import logging
import iris
import numpy as np
import matplotlib.pyplot as plt
from iris.util import equalise_attributes
import esmvalcore.preprocessor as eprep
import os
import cf_units

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic, select_metadata, group_metadata, get_diagnostic_filename, save_data, ProvenanceLogger
import esmvaltool.diag_scripts.shared.plot as eplot
from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import ProvenanceLogger
from esmvaltool.diag_scripts.water_vapour.distribution_moments_maps import reform_inp_data

# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))

def calculate_single_cube_stats(data_cube):

    data_dic = {}

    mean_cb = eprep.climate_statistics(data_cube, operator='mean')
    diffs = data_cube - mean_cb
    sqrt_f = iris.analysis.maths.IFunc(np.sqrt, lambda cube: cf_units.Unit('1'))
    sqrt_var_cb = sqrt_f(eprep.climate_statistics(diffs**2, operator='mean'))
    skew_cb = eprep.climate_statistics(diffs**3, operator='mean')/sqrt_var_cb**3
    skew_cb.standard_name=mean_cb.standard_name
    skew_cb.long_name=mean_cb.long_name
    skew_cb.var_name=mean_cb.var_name
    skew_cb.units=mean_cb.units
    skew_cb.attributes=mean_cb.attributes

    data_dic = { 'mean': mean_cb,
                 'skew': skew_cb}

    return data_dic


def calculate_cube_stats(data_cube, weights):

    data_dic = {}

    sqrt_f = iris.analysis.maths.IFunc(np.sqrt, lambda cube: cf_units.Unit('1'))

    mean_cb = data_cube.collapsed(['time','n_order'], iris.analysis.MEAN, weights=weights)
    diffs = data_cube - mean_cb
    sq_diffs = diffs**2
    cb_diffs = diffs**3
    std_cb = sqrt_f(sq_diffs.collapsed(['time','n_order'], iris.analysis.MEAN, weights=weights))
    skew_cb = cb_diffs.collapsed(['time','n_order'], iris.analysis.MEAN, weights=weights)/std_cb**3
    skew_cb.standard_name=mean_cb.standard_name
    std_cb.standard_name=mean_cb.standard_name
    skew_cb.long_name=mean_cb.long_name
    skew_cb.var_name=mean_cb.var_name
    skew_cb.units=mean_cb.units
    skew_cb.attributes=mean_cb.attributes
      
    data_dic = { 'mean': mean_cb,
                 'skew': skew_cb}

    return data_dic


def pca_axes(v1,v2,n_components=2,scale=2):

    """
    This routine simply returns the coordinates of the two principal components
    of two datasets that are plotted against each other
    
    :v1:            x-data to be plotted
    :v2:            y-data to be plotted
    :n_components:  the number of principle components to return
    :scale:         length of the axes to be plotted representing
                    the number of standard devitions
    """

    # Retrieve parameters to standardize the input data
    v1_sc = StandardScaler().fit(v1)
    v2_sc = StandardScaler().fit(v2)

    # Standardize the input data
    v1_std = v1_sc.transform(v1)
    v2_std = v2_sc.transform(v2)

    # Determine the eigenvalues and eigenvectors for PCA
    pca = PCA(n_components=n_components)
    pca.fit(np.hstack((v1_std,v2_std)))

    count = 1
    coords = {}

    # retrieve the eigenvectors and the variance in that component
    for length,vector in zip(pca.explained_variance_,pca.components_):
        v = vector * scale * np.sqrt(length)
        
        vec_1 = pca.mean_ + v
        vec_2 = pca.mean_ - v

        vec_x = np.array([vec_1[0], vec_2[0]])
        vec_y = np.array([vec_1[1], vec_2[1]])

        coords['pca'+str(count)] ={}
        coords['pca'+str(count)]['v1']=v1_sc.inverse_transform(vec_x.reshape(-1,1)).ravel()
        coords['pca'+str(count)]['v2']=v2_sc.inverse_transform(vec_y.reshape(-1,1)).ravel()
        count+=1

    return coords

def plot_slopes(data_dic, group, cfg):

    st_file = eplot.get_path_to_mpl_style(cfg.get('mpl_style'))
    plt.style.use(st_file)

    def_sty = eplot.get_dataset_style('default', cfg['mpl_style'] +'.yml')
    lbls = []
    fig_slopes, ax_slopes = plt.subplots(ncols=1, nrows=1)
    fig_slopes.set_size_inches(7,6)
    fig_slopes.set_dpi(300)

    for dataset in data_dic.keys():
        dtst_sty = eplot.get_dataset_style(dataset, cfg['mpl_style'] +'.yml')
        coords = pca_axes(data_dic[dataset]['mean'].data.ravel().reshape(-1,1), data_dic[dataset]['skew'].data.ravel().reshape(-1,1))
        if dtst_sty['color'] != def_sty['color']:
            lbls.append([dataset, dtst_sty['color']])
        ax_slopes.plot(coords['pca1']['v1'], coords['pca1']['v2'], c=dtst_sty['color'])
    
    for n, lbl_inf in enumerate(lbls): 
        ax_slopes.text(0.75, 0.975-n*0.06, lbl_inf[0], c=lbl_inf[1], 
                                               transform=ax_slopes.transAxes)
    ax_slopes.text(0.75, 0.975-(n+1)*0.06, 'other CMIP6', c=def_sty['color'],
                                               transform=ax_slopes.transAxes)

    ax_slopes.set_xlim(0,70) 
    ax_slopes.set_xlabel('mean TCWV, kg/m$^2$')
    ax_slopes.set_ylim(-2.2, 2.2)
    ax_slopes.set_ylabel('skewness TCWV, kg/m$^2$')   
    
    fig_slopes.suptitle('First mode of PCA in '+cfg['region']+' for '+ group+' TCWV \n ('+cfg['time_range']+')')
    plt.tight_layout()
    fig_slopes.savefig(os.path.join(cfg['plot_dir'], 'slopes_'+group+'_'+ cfg.get('region')+ diagtools.get_image_format(cfg)))
    plt.close(fig_slopes)


    return

def main(cfg):

    input_data = cfg['input_data'] ; reform_inp_data(input_data)

    groups = group_metadata(input_data.values(), 'variable_group', sort=True)
    datasets = group_metadata(input_data.values(), 'dataset', sort=True)

    for group in groups.keys():
        group_dic = {} 
        datasets = group_metadata(groups[group], 'dataset', sort=True)
        for dataset in datasets.keys():
            filepaths = list(group_metadata(datasets[dataset], 'filename').keys())
            n_real = len(filepaths)
            if n_real == 1: 
                mod_cb = iris.load_cube(filepaths[0])
                mod_weight = np.ones(mod_cb.shape)
                group_dic[dataset] = calculate_single_cube_stats(mod_cb)
            else:
                mod_cblst = iris.load(filepaths)
                equalise_attributes(mod_cblst)
                [cb.add_aux_coord(iris.coords.AuxCoord(n, long_name='n_order',
                        var_name='n_order')) for n, cb in enumerate(mod_cblst)]
                mod_cb = mod_cblst.merge_cube()
                mod_weight = np.ones(mod_cb.shape)/n_real
                group_dic[dataset] = calculate_cube_stats(mod_cb, mod_weight)
        plot_slopes(group_dic, group, cfg)

    logger.info('Success')

if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)