import iris
import iris.plot as iplt
import cartopy.crs as ccrs 
import cartopy.feature as cf
import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os


# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic, group_metadata
import esmvaltool.diag_scripts.shared.plot as eplot


# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
# logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

def create_map_plot(data_dic, mixns, ano_map_cb, cfg):

    b_min = np.floor(np.min([np.percentile(ano_map_cb.data, 1, method='closest_observation'), mixns['min']])) 
    b_max = np.ceil(np.max([np.percentile(ano_map_cb.data, 99, method='closest_observation'), mixns['max']]))
    prov_borders = cf.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none', edgecolor='k')

    prov_stat_dic = {'AB':{'lat': 54, 'lon': -114.5}, 'BC': {'lat': 54, 'lon': -123},
                     'MB': {'lat': 54, 'lon': -99}, 'NB': {'lat': 46.6, 'lon': -67.5},
                     'NL': {'lat': 53, 'lon': -62}, 'NS': {'lat': 44.25, 'lon': -65.7},
                     'NT': {'lat': 63, 'lon': -121}, 'NU': {'lat': 65, 'lon': -99},
                     'ON': {'lat': 50, 'lon': -86}, 'PE': {'lat': 46.4, 'lon': -63.2},
                     'QC': {'lat': 53, 'lon': -72}, 'SK':{'lat': 54, 'lon': -107}, 
                     'YT': {'lat': 63, 'lon': -136}, 'ID':{'lat': 44, 'lon': -116},
                     'ME': {'lat': 45, 'lon': -70.2}, 'MI': {'lat': 44.8, 'lon': -85.7},
                     'MN': {'lat': 46, 'lon': -95}, 'MT': {'lat': 46.5, 'lon': -110},
                     'ND': {'lat': 47, 'lon': -101}, 'NY': {'lat': 43, 'lon': -76},
                     'OR': {'lat': 43.5, 'lon': -121}, 'SD': {'lat': 44.5, 'lon': -101},  
                     'WA': {'lat': 47, 'lon': -121}, 'WI': {'lat': 44, 'lon': -90},
                     'WY': {'lat': 43, 'lon': -108}}   

    if cfg['cbar_positive']: 
        levels = np.arange(b_min, b_max+0.1, 0.1)
        cmap = 'Reds'
    else:
        abs_max = np.abs([b_min, b_max]).max()
        levels = np.arange(-1*abs_max, abs_max+0.1, 0.1)
        cmap = 'RdBu_r'

    ratio_colors = [(0, (215/256, 179/256, 119/256)), 
                    (0.4, (245/256, 240/256, 223/256)), 
                    # (0.286, (245/256, 240/256, 223/256)), 
                    (1.0, (66/256, 106/256, 90/256))]
    ratio_cmap = LinearSegmentedColormap.from_list('ratio_cmap',ratio_colors)

    for dataset in data_dic.keys(): 
        fig_map, ax_map = plt.subplots(nrows=2, ncols=1, sharex=False, sharey=False,
                                       subplot_kw={'projection': ccrs.LambertConformal(central_longitude=-91.87,
                                                                                       standard_parallels=(49,77))})
        ax_map = ax_map.flatten()

        fig_map.set_size_inches(7., 9.)
        fig_map.set_dpi(200)

        model_data = list()
        weights = list()

        for i in range(0,len(data_dic[dataset])):
            model_data.append(data_dic[dataset][i].data)
            weights.append(data_dic[dataset][i].attributes['ensemble_weight']*
                            data_dic[dataset][i].attributes['reverse_dtsts_n'])
        weights = np.array(weights)
        model_data = np.array(model_data)
        # if dataset == 'Multi-Model-Mean':
        #     mean_arr = np.average(model_data, axis=0, weights=weights)
        # else:
        mean_arr = np.average(model_data, axis=0)
        
        map_c = ax_map[0].contourf(data_dic[dataset][i].coord('longitude').points,
                                   data_dic[dataset][i].coord('latitude').points,
                                   mean_arr, levels=levels, extend='both', cmap=cmap)
        ax_map[0].set_title(dataset+' '+cfg['title_var_label'])

        iplt.contourf(ano_map_cb, axes=ax_map[1], levels=levels, extend='both', cmap=cmap)
        ax_map[1].set_title('ERA5 '+ cfg['title_var_label'])

        [a.coastlines(linewidth=0.5) for a in ax_map]
        [a.add_feature(cf.BORDERS) for a in ax_map]
        [a.add_feature(prov_borders) for a in ax_map]
        [a.set_extent([-63, -123,37,75]) for a in ax_map]
        fig_map.subplots_adjust(left=0.02, bottom=0.1, right=0.99, top=0.92, wspace=0.02, hspace=0.1)
        cax = fig_map.add_axes([0.1,0.07,0.8,0.015])
        fig_map.colorbar(map_c, cax=cax, orientation='horizontal', label=cfg['cbar_label'] +', '+ cfg['cbar_units'])

        fig_map.savefig(os.path.join(cfg['plot_dir'], 'map_anomalies_'+dataset+'.'+cfg['output_file_type']))

        if cfg['ratio']:
            fig_rat, ax_rat = plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False, 
                                           subplot_kw={'projection': ccrs.LambertConformal(central_longitude=-91.87,
                                                                                       standard_parallels=(49,77))})

            fig_rat.set_size_inches(7., 5.)
            fig_rat.set_dpi(200)

            map_rat = iplt.contourf(1/(ano_map_cb/mean_arr), axes=ax_rat, levels=np.arange(0.2, 2.1, 0.1), extend='both', cmap=ratio_cmap)

            ax_rat.set_title('Ratio of '+cfg['title_var_label']+' ('+dataset+'/ERA5)')
            ax_rat.coastlines(linewidth=0.5)
            ax_rat.add_feature(cf.BORDERS, linewidth=0.5) 
            ax_rat.add_feature(prov_borders, linewidth=0.5)
            for pr_st in prov_stat_dic.keys():
                ax_rat.text(prov_stat_dic[pr_st]['lon'],prov_stat_dic[pr_st]['lat'], pr_st, transform=ccrs.PlateCarree())

            ax_rat.set_extent([-60.5, -121.5,37.5,73.5])

            fig_rat.subplots_adjust(left=0.02, bottom=0.15, right=0.99, top=0.92, wspace=0.02, hspace=0.1)
            cax = fig_rat.add_axes([0.1,0.13,0.8,0.012])
            fig_rat.colorbar(map_rat, cax=cax, orientation='horizontal', label=cfg['cbar_label'] +' ratio')

            fig_rat.savefig(os.path.join(cfg['plot_dir'], 'ratio_anomalies_'+dataset+'.'+cfg['output_file_type']))


    return


def main(cfg):

    input_data = cfg['input_data']
    
    groups = group_metadata(input_data.values(), 'variable_group', sort=True)
    obs_info = groups.pop('obs')
    ano_map_cb = iris.load_cube(obs_info[0]['filename'])

    groups_l = list(groups.keys())

    plotting_dic = {}
    mixns = {}

    for group in groups_l:
        mins = list(); maxs = list()
        plotting_dic[group] = {}; mixns[group] = {}
        group_data = groups[group]
        datasets = group_metadata(group_data, 'dataset')
        ens_cubelist = iris.cube.CubeList()
        for dataset in datasets.keys():
            filepaths = list(group_metadata(datasets[dataset], 'filename').keys())
            n_real = len(filepaths)
            mod_cubelist = iris.cube.CubeList()
            for filepath in filepaths:
                mod_cb = iris.load_cube(filepath)
                mod_cb.attributes['ensemble_weight'] = 1 / n_real
                mod_cb.attributes['reverse_dtsts_n'] = 1/ len(datasets)
                ens_cubelist.append(mod_cb)
                mod_cubelist.append(mod_cb)
                mins.append(mod_cb.data.min()); maxs.append(mod_cb.data.max())
            plotting_dic[group][dataset] =  mod_cubelist
        mixns[group]['max'] = np.asarray(maxs).max()
        mixns[group]['min'] = np.asarray(mins).min()
        plotting_dic[group]['Multi-Model-Mean'] =  ens_cubelist

    st_file = eplot.get_path_to_mpl_style(cfg.get('mpl_style'))
    plt.style.use(st_file)

    create_map_plot(plotting_dic['txx_all'], mixns['txx_all'], ano_map_cb, cfg)

    logger.info('Success')


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)