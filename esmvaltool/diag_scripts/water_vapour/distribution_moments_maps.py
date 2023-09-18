import logging
import cf_units
import iris
from iris.util import equalise_attributes, unify_time_units
import iris.plot as iplt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs 
import esmvalcore.preprocessor as eprep
import os

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic, select_metadata, group_metadata, get_diagnostic_filename, save_data, ProvenanceLogger
import esmvaltool.diag_scripts.shared.plot as eplot
from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import ProvenanceLogger

# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))

def calculate_single_cube_stats(data_cube):

    data_dic = {}

    mean_cb = eprep.climate_statistics(data_cube, operator='mean')
    std_cb = eprep.climate_statistics(data_cube, operator='std_dev')
    p95_cb = data_cube.collapsed('time', iris.analysis.PERCENTILE, percent=95)
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
                 'std': std_cb,
                 'skew': skew_cb,
                 'p95': p95_cb}

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
    # checking the lat/lon mask
    cropped_cb = data_cube[0,0,:,:]
    lat_lon_mask = cropped_cb.data.mask
    if not np.all(lat_lon_mask==False):
        lon_bad_idx = np.where(np.all(lat_lon_mask, axis=0))[0]
        lat_bad_idx = np.where(np.all(lat_lon_mask, axis=1))[0]
        good_lons = np.delete(cropped_cb.coord('longitude').points, 
                                lon_bad_idx)
        good_lats = np.delete(cropped_cb.coord('latitude').points, lat_bad_idx)
        data_cube = data_cube.extract(iris.Constraint(latitude=good_lats, 
                                                    longitude=good_lons))
        weights = np.delete(weights, lon_bad_idx, axis=3) 
        weights = np.delete(weights, lat_bad_idx, axis=2) 
    p95_cb = data_cube.collapsed(['time','n_order'], iris.analysis.WPERCENTILE, 
                                    percent=95, weights=weights)
    skew_cb.standard_name=mean_cb.standard_name
    std_cb.standard_name=mean_cb.standard_name
    skew_cb.long_name=mean_cb.long_name; std_cb.long_name=mean_cb.long_name
    skew_cb.var_name=mean_cb.var_name; std_cb.var_name=mean_cb.var_name
    skew_cb.units=mean_cb.units; std_cb.units=mean_cb.units
    skew_cb.attributes=mean_cb.attributes; std_cb.attributes=mean_cb.attributes
      
    data_dic = { 'mean': mean_cb,
                 'std': std_cb,
                 'skew': skew_cb,
                 'p95': p95_cb}

    return data_dic

def plot_maps(data_dic, dataset, cfg, ensemble=None, bias=False, ref_data=None):

    n_seas = len(data_dic.keys())
    if n_seas == 1: 
        ncols = 2 ; nrows = 2 ; fig_title = list(data_dic.keys())[0]
    else:
        ncols = 4 ; nrows = n_seas ; fig_title=''       

    plot_dic = {'mean': {'levels': np.linspace(0,65,25), 
                         'ticks': np.arange(0, 65, 20), 'cbar': 'GnBu'}, 
                'std': {'levels': np.linspace(0,21,21), 
                        'ticks': np.arange(0, 21, 5), 'cbar': 'GnBu'},
                'skew': {'levels': np.linspace(-2.2,2.2,24), 
                         'ticks': np.arange(-2, 2.5, 1), 'cbar': 'RdBu'},
                'p95': {'levels': np.linspace(0,70,20), 
                        'ticks': np.arange(0, 70, 20), 'cbar': 'GnBu'}}
    if bias: 
        plot_dic = {'mean': {'levels': np.linspace(-20,20,20), 
                             'ticks': np.arange(-20, 21, 10), 'cbar': 'RdBu'}, 
                    'std': {'levels': np.linspace(-8,8,15), 
                            'ticks': np.arange(-8, 9, 4), 'cbar': 'RdBu'},
                    'skew': {'levels': np.linspace(-2,2,10),
                             'ticks': np.arange(-2, 2.5, 1), 'cbar': 'RdBu'},
                    'p95': {'levels': np.linspace(-20,20,20), 
                            'ticks': np.arange(-20, 21, 10), 'cbar': 'RdBu'}}        
        
    fig_maps, ax_maps = plt.subplots(ncols=ncols, nrows=nrows,
                                      subplot_kw={'projection': ccrs.Robinson()})
    fig_maps.set_size_inches(ncols*2.25, nrows*1.25+0.25)
    fig_maps.set_dpi(300)
    ax_maps = ax_maps.flatten(); n=0

    fig_name = 'wv_' ; fig_title += 'TCWV '
    if bias: 
        fig_name += 'bias_'
        fig_title += 'bias '
    fig_name += dataset
    fig_title += dataset
    if ensemble != None:
        fig_name += '_' +ensemble
        fig_title += ' ' +ensemble
    fig_title = fig_title +' ('+cfg['time_range']+')'
    fig_maps.suptitle(fig_title)
    
    for n_s, seas in enumerate(data_dic.keys()):
        if n_seas>1:
            ax_maps[n].text(-0.06, 0.4, seas, transform=ax_maps[n].transAxes,
                                                       rotation='vertical')
        for n_st, stat in enumerate(data_dic[seas].keys()):
            if bias == True:
                # to subtract one cube from the other
                tmp_cblst= iris.cube.CubeList([data_dic[seas][stat], ref_data[seas][stat]])
                equalise_attributes(tmp_cblst)
                overlap_cblst= tmp_cblst.extract_overlapping(['latitude', 'longitude'])
                levels = iplt.contourf(iris.analysis.maths.subtract(overlap_cblst[0],
                                                         overlap_cblst[1].data), 
                              axes=ax_maps[n], levels=plot_dic[stat]['levels'],
                              extend='both', cmap=plot_dic[stat]['cbar'])
            else:
                levels = iplt.contourf(data_dic[seas][stat], axes=ax_maps[n], 
                              levels=plot_dic[stat]['levels'], extend='both', 
                              cmap=plot_dic[stat]['cbar'])
            ax_maps[n].coastlines(linewidth=0.5)
            if n_s == 0:
                ax_maps[n].set_title(stat)
            if n_seas == 1:
                ax_loc = make_axes_locatable(ax_maps[n])
                cax = ax_loc.append_axes('bottom', size='6%', pad=0.07, axes_class=plt.Axes)
                cbar = fig_maps.colorbar(levels, cax=cax, orientation='horizontal')
                cbar.set_ticks(plot_dic[stat]['ticks'], labels=plot_dic[stat]['ticks'], fontsize='small')
                cbar.set_label('kg/m$^2$', fontsize='small')
            elif (n_seas > 1)&(n_s == n_seas - 1): 
                fig_maps.subplots_adjust(left=0.015, bottom=0.08, right=0.999, top=0.92, wspace=0.01, hspace=0.02)
                posn = ax_maps[n].get_position()
                cax = fig_maps.add_axes([posn.x0, 0.06, posn.width, 0.01])
                cbar = fig_maps.colorbar(levels, cax=cax, orientation='horizontal')
                cbar.set_ticks(plot_dic[stat]['ticks'], labels=plot_dic[stat]['ticks'], fontsize='small')
                cbar.set_label('kg/m$^2$', fontsize='small', labelpad=0.01)
            n = n + 1
    
    if n_seas==1:
        plt.tight_layout()

    fig_maps.savefig(os.path.join(cfg['plot_dir'], fig_name + diagtools.get_image_format(cfg)))
    plt.close(fig_maps)

    return

def reform_inp_data(input_data):

    # this is a function to reform input_data
    # to separate CanESM5-1 p1 and p2 realisations

    canesm = select_metadata(input_data.values(), dataset='CanESM5-1')
    canesm_ens = group_metadata(canesm, 'ensemble')

    for ens in canesm_ens.keys():
        if 'p1' in ens:
            ens_mtdata = select_metadata(canesm, ensemble=ens)
            filnames = group_metadata(ens_mtdata, 'filename')
            for filname in filnames.keys():
                input_data[filname]['dataset'] = 'CanESM5-1(p1)'
        elif 'p2' in ens:
            ens_mtdata = select_metadata(canesm, ensemble=ens)
            filnames = group_metadata(ens_mtdata, 'filename')
            for filname in filnames.keys():
                input_data[filname]['dataset'] = 'CanESM5-1(p2)'

    return


def main(cfg):

    input_data = cfg['input_data'] ; reform_inp_data(input_data)

    ref_dic = {} ; mmm_dic = {} 

    data_groups = list(group_metadata(input_data.values(), 'variable_group').keys())

    ref_dataset = list(group_metadata(input_data.values(), 'reference_dataset').keys())[0]

    for group in data_groups: 
        ref_f_name = select_metadata(input_data.values(), 
                                       dataset=ref_dataset, 
                                       variable_group=group)[0]['filename']
        ref_cube = iris.load_cube(ref_f_name)
        ref_dic[group] = calculate_single_cube_stats(ref_cube)
        mmm_dic[group] = {'data': [], 'weights': []}

    plot_maps(ref_dic, ref_dataset, cfg)

    datasets = group_metadata(input_data.values(), 'dataset', sort=True)
    datasets.pop(ref_dataset)

    for dataset in datasets.keys():
        dtst_dic = {} 
        groups = group_metadata(datasets[dataset], 'variable_group')
        for group in groups.keys():
            filepaths = list(group_metadata(groups[group], 'filename').keys())
            n_real = len(filepaths)
            if n_real == 1: 
                mod_cb = iris.load_cube(filepaths[0])
                mod_weight = np.ones(mod_cb.shape)
                dtst_dic[group] = calculate_single_cube_stats(mod_cb)
            else:
                mod_cblst = iris.load(filepaths)
                equalise_attributes(mod_cblst)
                [cb.add_aux_coord(iris.coords.AuxCoord(n, long_name='n_order',
                        var_name='n_order')) for n, cb in enumerate(mod_cblst)]
                mod_cb = mod_cblst.merge_cube()
                mod_weight = np.ones(mod_cb.shape)/n_real
                dtst_dic[group] = calculate_cube_stats(mod_cb, mod_weight)
            mmm_dic[group]['data'].append(mod_cb) 
            mmm_dic[group]['weights'].append(mod_weight)
        plot_maps(dtst_dic, dataset, cfg)
        plot_maps(dtst_dic, dataset, cfg, bias=True, ref_data=ref_dic)

    logger.info('Success')

if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)