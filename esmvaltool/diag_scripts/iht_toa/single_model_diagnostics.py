""" Single model diagnostics
1. Solve the Poisson solver
2. Produce and save plots
"""

import logging
from pathlib import Path

import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec, rcParams
from poisson_solver import spherical_poisson

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_figure,
)

# Initialise logger
logger = logging.getLogger(Path(__file__).stem)

rcParams.update({
    'font.size': 14,
    'text.latex.preamble': [r"\usepackage{amsmath}"],
    'xtick.major.pad': 10,
    'ytick.major.pad': 10,
    'xtick.major.size': 10,
    'ytick.major.size': 10,
    'xtick.minor.size': 5,
    'ytick.minor.size': 5,
    'axes.linewidth': 2,
    'lines.markersize': 8,
    'lines.linewidth': 2
})


def get_provenance_record(plot_type, ancestor_files):
    # Create a provenance record describing the diagnostic data and plot.

    record = {
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_types': [plot_type],
        'authors': [
            'andela_bouwe',
            'righi_mattia',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': ancestor_files,
    }
    return record


def var_name_constraint(var_name):
    return iris.Constraint(cube_func=lambda c: c.var_name == var_name)


def call_poisson(flux_cube, latitude='latitude', longitude='longitude'):

    if flux_cube.coord(latitude).bounds is None:
        flux_cube.coord(latitude).guess_bounds()
    if flux_cube.coord(longitude).bounds is None:
        flux_cube.coord(longitude).guess_bounds()

    # Remove average of flux field to account for storage term
    data = flux_cube.data.copy()
    grid_areas = iris.analysis.cartography.area_weights(flux_cube)
    data_mean = flux_cube.collapsed(["longitude", "latitude"],
                                    iris.analysis.MEAN,
                                    weights=grid_areas).data
    data -= data_mean

    logger.info("Calling spherical_poisson")
    poisson, mht = spherical_poisson(logger,
                                     forcing=data * (6371e3**2.0),
                                     tolerance=2.0e-4)
    logger.info("Ending spherical_poisson")

    # Energy flux potential (P)
    p_cube = flux_cube.copy()
    p_cube.var_name = "{}_efp".format(flux_cube.var_name)
    p_cube.long_name = "energy_flux_potential_of_{}".format(flux_cube.var_name)
    p_cube.standard_name = None
    p_cube.units = 'J s-1'
    p_cube.data = poisson[1:-1, 1:-1]

    # MHT data cube
    mht_cube = flux_cube.copy()
    mht_cube = mht_cube.collapsed('longitude', iris.analysis.MEAN)
    mht_cube.var_name = "{}_mht".format(flux_cube.var_name)
    mht_cube.long_name = "meridional_heat_transport_of_{}".format(
        flux_cube.var_name)
    mht_cube.standard_name = None
    mht_cube.units = 'W'
    mht_cube.data = mht

    return p_cube, mht_cube


def compute_efp_and_mht(model_dict):
    # Sort experiments by variable
    exp_dataset = group_metadata(model_dict, 'exp', sort='variable_group')

    # Define dict to hold the data
    iht_dict = {}

    # Loop over experiments
    for exp_nme in exp_dataset:
        logger.info("Experiment: %s", exp_nme)
        p_cube = iris.cube.CubeList()
        mht_cube = iris.cube.CubeList()
        flux_cube = iris.cube.CubeList()
        files_dict = {}
        # Loop over variables
        for attributes in exp_dataset[exp_nme]:
            input_file = attributes['filename']
            variable_name = attributes['short_name']
            logger.info("Loading %s", input_file)
            cube = iris.load_cube(input_file)
            files_dict[variable_name] = input_file
            # Call the Poisson solver
            efp, mht = call_poisson(cube)
            p_cube.append(efp)
            mht_cube.append(mht)
            flux_cube.append(cube)
            del efp, mht, cube

        # Add iris cubelist to dictionary
        iht_dict[exp_nme] = implied_heat_transport(mht_cube, p_cube, flux_cube,
                                                   files_dict)
        del p_cube, mht_cube, flux_cube

    return iht_dict


class implied_heat_transport:
    def __init__(self, mht_clim, efp_clim, flx_clim, flx_files):
        self.mht_clim = mht_clim
        self.efp_clim = efp_clim
        self.flx_clim = flx_clim
        self.flx_files = flx_files

        # Grid
        vname_constraint = var_name_constraint("rtnt")
        flux = self.flx_clim.extract_cube(vname_constraint)
        self.grid = iris.analysis.cartography.area_weights(flux,
                                                           normalize=True)
        self.lat = flux.coord('latitude').points
        self.lon = flux.coord('longitude').points
        del flux

        # Compute derived fields
        self.derived_fields()
        return

    def derived_fields(self):
        for cube in self.mht_clim:
            if cube.var_name == "rlut_mht":
                dcube = cube.copy()
                dcube.data = -dcube.data
                dcube.var_name = "rlnt_mht"
                dcube.long_name = "rlnt_mht"
                self.mht_clim.append(dcube)

    def print(self):
        print("=== implied_heat_transport object ===")
        print(self.mht_clim)
        for x in self.mht_clim:
            print(x.long_name, x.var_name)
        print(self.efp_clim)
        for x in self.efp_clim:
            print(x.long_name, x.var_name)
        print(self.flx_clim)
        for x in self.flx_clim:
            print(x.long_name, x.var_name)
        print(self.flx_files)

    def hemispheric_symmetry(self, data):
        # Calculates hemispheric symmetry value
        # S = 0 is perfectly symmetrical

        nh = data[90:]
        sh = data[:90]
        sh = sh[::-1]
        grid = np.sum(self.grid, axis=1)[90:]

        diff = np.abs((nh + sh) * grid)
        hem = np.sum(diff)
        trop = np.sum(diff[:30])
        mid = np.sum(diff[30:60])
        high = np.sum(diff[60:])
        extratrop = np.sum(diff[30:90])
        return hem, trop, mid, high, extratrop

    def mht_plot(self, var_names, legend_label, ylim=(-10, 10)):
        """MHT plot Produces a single plot comparing the estimated MHT due to
        the input variables.

        MHT is presented in PW, plotted against latitude. Up to three
        variables are on each plot.
        """
        plt.figure()
        for i in range(len(var_names)):
            mht = self.mht_clim.extract_cube(var_name_constraint(
                var_names[i])).data / 1e15
            plt.plot(self.lat, mht, label=legend_label[i])
        plt.hlines(0, -90, 90, color='k', linestyles=':')
        plt.vlines(0, -10, 10, color='k', linestyles=':')
        plt.xlim(-90, 90)
        plt.ylim(ylim[0], ylim[1])
        plt.xticks(np.arange(-90, 120, 30))
        plt.xlabel('Latitude')
        plt.ylabel('MHT (PW)')
        plt.legend()
        plt.tight_layout()
        return

    def cre_mht_plot(self,
                     var_names_l,
                     legend_l,
                     var_names_r,
                     legend_r,
                     ylim=(-1.5, 1.5)):
        plt.figure(figsize=(11, 5))
        ax1 = plt.subplot(121)
        for i in range(len(var_names_l)):
            mht = self.mht_clim.extract_cube(
                var_name_constraint(var_names_l[i])).data / 1e15
            ax1.plot(self.lat, mht, label=legend_l[i])
        ax1.axhline(0, color='k', ls=':')
        ax1.axvline(0, color='k', ls=':')
        ax1.set_xlim(-90, 90)
        ax1.set_xticks(np.arange(-90, 120, 30))
        ax1.set_xlabel('Latitude')
        ax1.set_ylim(ylim[0], ylim[1])
        ax1.set_ylabel('MHT (PW)')
        ax1.annotate('(a)',
                     xy=(0.01, 0.95),
                     xycoords='axes fraction',
                     color='k')
        plt.legend()

        ax2 = plt.subplot(122)
        col = ['C3', 'C7']
        for i in range(len(var_names_r)):
            mht = self.mht_clim.extract_cube(
                var_name_constraint(var_names_r[i])).data / 1e15
            ax2.plot(self.lat, mht, label=legend_r[i], color=col[i])
        ax2.axhline(0, color='k', ls=':')
        ax2.axvline(0, color='k', ls=':')
        ax2.set_xlim(-90, 90)
        ax2.set_xticks(np.arange(-90, 120, 30))
        ax2.set_xlabel('Latitude')
        ax2.set_ylim(ylim[0], ylim[1])
        ax2.set_ylabel('MHT (PW)')
        ax2.annotate('(b)',
                     xy=(0.01, 0.95),
                     xycoords='axes fraction',
                     color='k')
        plt.legend(loc='lower right')
        plt.tight_layout()
        return

    def quiver_steps(self, n, step):
        n2 = (n - 2) // step
        start = (n - 2 - n2 * step) // 2
        return start

    def quiver_subplot(self,
                       var_name,
                       wmin,
                       wmax,
                       nwlevs,
                       wlevstep,
                       vmin,
                       vmax,
                       label=[['(a)', '(b)'], ['(c)', '(d)'], ['(e)', '(f)']],
                       xy_label=(0, 1.05),
                       title=[['', ''], ['', ''], ['', '']]):
        x, y = np.meshgrid(self.lon, self.lat)

        levels1 = np.linspace(vmin, vmax, 11)
        levels2 = np.linspace(wmin, wmax, nwlevs)
        nrows = len(var_name)
        # Calculate sampling for vector plot
        nlon = len(self.lon)
        nlat = len(self.lat)
        stepx = 20
        stepy = 20
        startx = self.quiver_steps(nlon, stepx)
        starty = self.quiver_steps(nlat, stepy)

        if nrows == 3:
            plt.figure(figsize=(10, 10))
            gs = gridspec.GridSpec(22, 2)
            gs.update(wspace=0.25, hspace=1.5)
        elif nrows == 2:
            plt.figure(figsize=(10, 6.5))
            gs = gridspec.GridSpec(15, 2)
            gs.update(wspace=0.25, hspace=1.5)
        elif nrows == 1:
            plt.figure(figsize=(12, 4))
            gs = gridspec.GridSpec(8, 2)
            gs.update(wspace=0.25, hspace=1.5)

        for i in range(nrows):
            efp = self.efp_clim.extract_cube(
                var_name_constraint(var_name[i][0])).data
            flx = self.flx_clim.extract_cube(
                var_name_constraint(var_name[i][1])).data
            efp -= np.average(efp)  # Arbitrary choice of origin
            flx -= np.average(flx)
            v, u = np.gradient(efp, 1e14, 1e14)
            u = u[1:-1, 1:-1]
            v = v[1:-1, 1:-1]

            ax1 = plt.subplot(gs[i * 7:(i * 7) + 7, 0],
                              projection=ccrs.PlateCarree())
            cb1 = ax1.contourf(self.lon,
                               self.lat,
                               efp / 1e15,
                               levels=levels1,
                               transform=ccrs.PlateCarree(central_longitude=0))
            plt.gca().coastlines()
            xq = x[starty::stepy, startx::stepx]
            yq = y[starty::stepy, startx::stepx]
            uq = u[starty::stepy, startx::stepx]
            vq = v[starty::stepy, startx::stepx]
            if i == 0:
                Q = ax1.quiver(xq,
                               yq,
                               uq,
                               vq,
                               pivot='mid',
                               color='w',
                               width=0.005)
                Q._init()
            else:
                ax1.quiver(xq,
                           yq,
                           uq,
                           vq,
                           pivot='mid',
                           scale=Q.scale,
                           color='w')
            ax1.set_xticks(np.arange(-180, 190, 60))
            ax1.set_xticklabels(
                ['180', '120W', '60W', '0', '60E', '120E', '180'])
            ax1.set_yticks(np.arange(-90, 100, 30))
            ax1.set_yticklabels(
                ['90S', '60S', '30S', 'Eq', '30N', '60N', '90N'])
            ax1.annotate(label[i][0],
                         xy=xy_label,
                         xycoords='axes fraction',
                         color='k')
            ax1.set_title(title[i][0])
            del u, v, uq, vq

            ax2 = plt.subplot(gs[i * 7:(i * 7) + 7, 1],
                              projection=ccrs.PlateCarree())
            cb2 = ax2.contourf(self.lon,
                               self.lat,
                               flx,
                               levels=levels2,
                               transform=ccrs.PlateCarree(central_longitude=0),
                               cmap='RdBu_r')
            plt.gca().coastlines()
            ax2.set_xticks(np.arange(-180, 190, 60))
            ax2.set_xticklabels(
                ['180', '120W', '60W', '0', '60E', '120E', '180'])
            ax2.set_yticks(np.arange(-90, 100, 30))
            ax2.set_yticklabels(
                ['90S', '60S', '30S', 'Eq', '30N', '60N', '90N'])
            ax2.annotate(label[i][1],
                         xy=xy_label,
                         xycoords='axes fraction',
                         color='k')
            ax2.set_title(title[i][1])

        ax1 = plt.subplot(gs[-1, 0])
        plt.colorbar(cb1,
                     cax=ax1,
                     orientation='horizontal',
                     label='Energy flux potential (PW)')

        ax2 = plt.subplot(gs[-1, 1])
        plt.colorbar(cb2,
                     cax=ax2,
                     orientation='horizontal',
                     label=r'Flux (Wm$^{-2}$)',
                     ticks=levels2[1::wlevstep])

        if nrows == 3:
            plt.subplots_adjust(left=0.1, right=0.94, top=1.0, bottom=0.11)
        elif nrows == 2:
            plt.subplots_adjust(left=0.11, right=0.9, top=1.0, bottom=0.13)
        elif nrows == 1:
            plt.subplots_adjust(left=0.11, right=0.9, top=1.0, bottom=0.20)

        return


def efp_maps(iht, model, experiment, cfg):
    # Figure 2
    iht.quiver_subplot(
        [['rtnt_efp', 'rtnt'], ['rsnt_efp', 'rsnt'], ['rlut_efp', 'rlut']],
        wmin=-180,
        wmax=180,
        nwlevs=19,
        wlevstep=4,
        vmin=-1.2,
        vmax=1.2,
        title=[['$P_{TOA}^{TOT}$', "$Delta F_{TOA}^{TOT}$"],
               ['$P_{TOA}^{SW}$', "$Delta F_{TOA}^{SW}$"],
               ['$P_{TOA}^{LW}$', "$Delta F_{TOA}^{LW}$"]])
    provenance_record = get_provenance_record(plot_type='map',
                                              ancestor_files=iht.flx_files)
    figname = "efp_and_flux_toa_net_{}_{}".format(model, experiment)
    save_figure(figname, provenance_record, cfg)
    # Figure 4
    iht.quiver_subplot(
        [['netcre_efp', 'netcre'], ['swcre_efp', 'swcre'],
         ['lwcre_efp', 'lwcre']],
        wmin=-60,
        wmax=60,
        nwlevs=13,
        wlevstep=2,
        vmin=-0.3,
        vmax=0.3,
        title=[['$P_{TOA}^{TOTCRE}$', '$Delta CRE_{TOA}^{TOT}$'],
               ['$P_{TOA}^{SWCRE}$', '$Delta CRE_{TOA}^{SW}$'],
               ['$P_{TOA}^{LWCRE}$', '$Delta CRE_{TOA}^{LW}$']])
    provenance_record = get_provenance_record(plot_type='map',
                                              ancestor_files=iht.flx_files)
    figname = "efp_and_flux_toa_cre_{}_{}".format(model, experiment)
    save_figure(figname, provenance_record, cfg)
    # Figure 5
    iht.quiver_subplot(
        [['rsutcs_efp', 'rsutcs'], ['rsut_efp', 'rsut']],
        wmin=-100,
        wmax=100,
        nwlevs=21,
        wlevstep=3,
        vmin=-0.3,
        vmax=0.3,
        title=[['$P_{TOA}^{SWup, clr}$', '$Delta CRE_{TOA}^{SWup, clr}$'],
               ['$P_{TOA}^{SWup, all}$', '$Delta CRE_{TOA}^{SWup, all}$']])
    provenance_record = get_provenance_record(plot_type='map',
                                              ancestor_files=iht.flx_files)
    figname = "efp_and_flux_toa_rsut_{}_{}".format(model, experiment)
    save_figure(figname, provenance_record, cfg)


def mht_plots(iht, model, experiment, cfg):
    # Figure 1
    iht.mht_plot(["rtnt_mht", "rsnt_mht", "rlnt_mht"], ['Net', 'SW', 'LW'])
    provenance_record = get_provenance_record(plot_type='zonal',
                                              ancestor_files=iht.flx_files)
    figname = "mht_toa_net_{}_{}.png".format(model, experiment)
    save_figure(figname, provenance_record, cfg)
    # Figure 3
    iht.cre_mht_plot(['netcre_mht', 'swcre_mht', 'lwcre_mht'],
                     ['Net CRE', 'SW CRE', 'LW CRE'],
                     ['rsut_mht', 'rsutcs_mht'],
                     ['-1 x OSR (all-sky)', '-1 x OSR (clear-sky)'])
    figname = "mht_toa_cre_and_osr_{}_{}".format(model, experiment)
    save_figure(figname, provenance_record, cfg)
    return


def plot_single_model_diagnostics(iht_dict, cfg):
    # iht_dict is a two-level dictionary: iht_dict[model][experiment]
    for model, iht_model in iht_dict.items():
        for experiment, iht_experiment in iht_model.items():
            mht_plots(iht_experiment, model, experiment, cfg)
            efp_maps(iht_experiment, model, experiment, cfg)


def main(cfg):
    """Solve the Poisson equation and estimate the meridional heat
    transport."""
    input_data = cfg['input_data'].values()

    # Dictionary of iht objects. Each entry of the dictionary
    # contains the iht object of a single model/dataset.
    iht = {}

    # Group data by dataset
    logger.info("Group input data by model, sort by experiment")
    model_dataset = group_metadata(input_data, 'dataset', sort='exp')

    # Solve Poisson equation for each dataset
    for model_name in model_dataset:
        if model_name == 'CERES-EBAF':
            # Ignore any observational data
            continue
        logger.info("Processing model data: %s", model_name)
        iht[model_name] = compute_efp_and_mht(model_dataset[model_name])

    # Produce plots
    plot_single_model_diagnostics(iht, cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
