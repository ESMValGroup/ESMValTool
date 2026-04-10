import iris
import iris.plot as iplt
import cartopy.crs as ccrs 
import cartopy.feature as cf
from shapely.geometry import shape
import fiona
import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import mpl_fontkit as fk
# from matplotlib.colors import LinearSegmentedColormap
import os


# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic, group_metadata
import esmvaltool.diag_scripts.shared.plot as eplot


# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
# logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

def plot_shape_outline(ax, shapefile, proj):

    reg_locs = {'Inuvik': (-131, 66, False),
                'Manitoba': (-100.0, 55.0, False),
                'W.Ontario': (-94.0, 50.0, False),
                'E.Ontario': (-83.0, 45.5, False),
                'Kivalliq': (-101.0, 65.0, False),
                'S.Quebec': (-79.0, 48.0, False),
                'N.Quebec': (-76.0, 58.0, False),
                'Kitikmeot': (-114.0, 66.0, True),
                'S.Qikiqtaaluk': (-85, 71.0, True),
                'N.Qikiqtaaluk': (-96.0, 75.0, True),
                'Atl.Canada': (-59.5, 49.5, True)}
    
    with fiona.open(shapefile) as src:
        for feature in src:
            name = feature["properties"]["ID"]
            geom = shape(feature["geometry"])

            if geom.geom_type == "Polygon":
                geoms = [geom]
            else:
                geoms = geom.geoms

            for poly in geoms:
                x, y = poly.exterior.xy
                ax.plot(x, y, linewidth=0.5, transform=proj, color='k')
            
            if name in reg_locs.keys():
                # label at predefined location
                if reg_locs[name][2]:
                    ax.text(reg_locs[name][0], reg_locs[name][1], name,
                             bbox={'facecolor': 'white', 'edgecolor': 'None', 
                                   'alpha': 0.9, 'pad': 0.15}, transform=proj)
                else:
                    ax.text(reg_locs[name][0], reg_locs[name][1], name, transform=proj)
            else:
                # label at centroid
                cx, cy = geom.centroid.xy
                ax.text(cx[0]-2, cy[0]-1, name, transform=proj)

def create_map_plot(data_dic, mixns, ano_map_cb, obs_name, cfg):

    b_min = np.floor(np.min([np.percentile(ano_map_cb.data, 1, method='closest_observation'), mixns['min']])) 
    b_max = np.ceil(np.max([np.percentile(ano_map_cb.data, 99, method='closest_observation'), mixns['max']]))

    shapefile = cfg['shapefile']

    if cfg['cbar_positive']: 
        levels = np.arange(b_min, b_max+0.1, 0.1)
        cmap = 'Reds'
    else:
        abs_max = np.abs([b_min, b_max]).max()
        levels = np.arange(-1*abs_max, abs_max+0.1, 0.1)
        cmap = 'RdBu_r'

    # ratio_colors = [(0, "#6a4c93"), 
    #                 (0.4, "#ece9e4"), 
    #                 # (0.286, (245/256, 240/256, 223/256)), 
    #                 (1.0, "#297430")]
    # ratio_cmap = LinearSegmentedColormap.from_list('ratio_cmap',ratio_colors)

    # norm = mcolors.TwoSlopeNorm(vmin=0.5, vcenter=1.0, vmax=2.1)
    # norm = mcolors.LogNorm(vmin=float(np.log10(0.5)), vmax=float(np.log10(2.1)))
    levels = np.logspace(np.log10(0.5), np.log10(2.0), 41)
    norm = mcolors.BoundaryNorm(boundaries=levels, ncolors=256, extend='both')

    for dataset in data_dic.keys(): 
        fig_map, ax_map = plt.subplots(nrows=2, ncols=1, sharex=False, sharey=False,
                                       subplot_kw={'projection': ccrs.LambertConformal(central_longitude=-91.87,
                                                                                       standard_parallels=(49,77))})
        ax_map = ax_map.flatten()

        fig_map.set_size_inches(7., 7.5)
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
        [a.set_extent([-63, -123,37,75]) for a in ax_map]
        fig_map.subplots_adjust(left=0.02, bottom=0.1, right=0.99, top=0.92, wspace=0.02, hspace=0.1)
        cax = fig_map.add_axes([0.1,0.07,0.8,0.015])
        fig_map.colorbar(map_c, cax=cax, orientation='horizontal', label=cfg['cbar_label'] +', '+ cfg['cbar_units'])

        fig_map.savefig(os.path.join(cfg['plot_dir'], 'map_anomalies_'+dataset+'.'+cfg['output_file_type']))

        if cfg['ratio']:
            
            plt.rcParams.update({'figure.autolayout': False})

            fig_rat, ax_rat = plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False, 
                                           subplot_kw={'projection': ccrs.LambertConformal(central_longitude=-91.87,
                                                                                       standard_parallels=(49,77))})

            fig_rat.set_size_inches(7., 6.7)
            fig_rat.set_dpi(300)

            map_rat = iplt.contourf(1/(ano_map_cb/mean_arr), axes=ax_rat, extend='both', cmap='PuOr_r', norm=norm, levels=levels)
            # map_rat = iplt.contourf(1/(ano_map_cb/mean_arr), axes=ax_rat, levels=np.arange(0.2, 2.1, 0.1), extend='both', cmap=ratio_cmap)

            ax_rat.set_title(f"Ratio of {cfg['title_var_label']} ({dataset}/{obs_name})")
            ax_rat.coastlines(linewidth=0.5)
            plot_shape_outline(ax_rat, shapefile, ccrs.PlateCarree())

            ax_rat.set_extent([-60, -119, 38, 79])

            fig_rat.subplots_adjust(left=0.02, bottom=0.12, right=0.99, top=0.95, wspace=0.02, hspace=0.1)
            cax = fig_rat.add_axes([0.05,0.08,0.9,0.012])
            fig_rat.colorbar(map_rat, cax=cax, orientation='horizontal', label=cfg['cbar_label'] +' ratio')

            fig_rat.savefig(os.path.join(cfg['plot_dir'], 'ratio_anomalies_'+dataset+'.'+cfg['output_file_type']))


    return


def main(cfg):

    input_data = cfg['input_data']
    
    groups = group_metadata(input_data.values(), 'variable_group', sort=True)
    obs_info = groups.pop('obs')
    obs_name = obs_info[0]['dataset']
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

    fk.install("Montserrat"); fk.set_font("Montserrat")

    create_map_plot(plotting_dic['txx_all'], mixns['txx_all'], ano_map_cb, obs_name, cfg)

    logger.info('Success')


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)