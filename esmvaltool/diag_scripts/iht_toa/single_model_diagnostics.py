""" Single model diagnostics
1. Solve the Poisson solver
2. Produce and save plots
"""

import logging
import sys
from copy import deepcopy
from pathlib import Path

import cartopy.crs as ccrs
import iris
import iris.plot as iplt
import matplotlib.dates as mdates
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
    """Create a provenance record describing the diagnostic data and plot."""

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


def area_average(cube, latitude='latitude', longitude='longitude', mdtol=1):
    """Returns the gridbox weighted area average of a cube (and optionally a
    sub-region).

    Compulsory arguments:
        cube            Cube to be area averaged.
    Optional arguments:
        latitude, longitude  Names of coordinates to be used for latitude
                             and longitude
        mdtol           Tolerance for missing data.
    Returns:
        cube_avg        cube which has had area averaging applied across
                        latitude and longitude coordinates
    """

    if cube.coord(latitude).bounds is None:
        cube.coord(latitude).guess_bounds()
    if cube.coord(longitude).bounds is None:
        cube.coord(longitude).guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(cube)
    cube_avg = cube.collapsed([longitude, latitude],
                              iris.analysis.MEAN,
                              weights=grid_areas,
                              mdtol=mdtol)
    return cube_avg


def var_name_constraint(var_name):
    """Shortcut to create constraint for variable name."""
    return iris.Constraint(cube_func=lambda c: c.var_name == var_name)


def call_poisson(flux_cube, latitude='latitude', longitude='longitude'):
    """Top-level function that calls the Poisson solver for source cube."""

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


def symmetry_metric(data, grid):
    """Calculates hemispheric symmetry value.

    :param data: zonal mean of the input variable
    :param grid: grid weights
    """
    # S = 0 is perfectly symmetrical.
    # As coded, the calculation of the symmetry metrics needs the number of
    # latitude points to be multiple of 6, i.e. it needs 30 deg bands.
    Nlat = data.shape[0]
    if (Nlat % 6) != 0:
        logger.error("Grid not compatible with symmetry metric calculation.")
        sys.exit(1)

    Nlat_hem = Nlat // 2
    Nlat_trop = Nlat_hem // 3
    nh = data[Nlat_hem:]
    sh = data[:Nlat_hem]
    sh = sh[::-1]

    diff = np.abs((nh + sh) * grid)
    hem = np.sum(diff)
    trop = np.sum(diff[:Nlat_trop])
    extratrop = np.sum(diff[Nlat_trop:Nlat_hem])
    return hem, trop, extratrop


class implied_heat_transport:
    def __init__(self, flx_files):
        self.flx_files = flx_files

        # Create cube lists for the different datasets
        self.flx_clim = iris.cube.CubeList()
        self.mht_clim = iris.cube.CubeList()
        self.efp_clim = iris.cube.CubeList()
        self.mht_rolling_mean = iris.cube.CubeList()
        self.symmetry_metric = iris.cube.CubeList()

        # Calculate 12-month rolling means for time series.
        self.flx_rolling_mean = iris.cube.CubeList()
        for flx_file in flx_files:
            flx = iris.load_cube(flx_file)
            if len(flx.shape) == 3:
                self.flx_rolling_mean.append(
                    flx.rolling_window('time', iris.analysis.MEAN, 12))
            else:
                self.flx_clim.append(flx)

        # Save grid information
        vname_constraint = var_name_constraint("rtnt")
        flx = self.flx_clim.extract_cube(vname_constraint)
        self.grid = iris.analysis.cartography.area_weights(flx, normalize=True)
        self.lat = flx.coord('latitude').points
        self.lon = flx.coord('longitude').points

        # Compute derived fluxes
        self.derived_fluxes()

        # Calculate EFP and MHT for each flux
        self.compute_efp_and_mht()

        # Times series of MHT symmetry metric
        self.mht_symmetry_metrics()

        self.print()
        return

    def compute_efp_and_mht(self):
        """Calculation of Energy Flux Potential and Meridional heat transport
        for all the data fields."""
        # Loop over climatologies
        for flx in self.flx_clim:
            efp, mht = call_poisson(flx)
            self.efp_clim.append(efp)
            self.mht_clim.append(mht)
        # Loop over rolling means
        for flx_rm in self.flx_rolling_mean:
            mht_series = iris.cube.CubeList()
            for flx in flx_rm.slices_over('time'):
                efp, mht = call_poisson(flx)
                mht_series.append(mht)
            # Append MHT rolling mean after merging time series.
            self.mht_rolling_mean.append(mht_series.merge_cube())

    def derived_fluxes(self):
        """Calculation of the derived fluxes:

        rtntcs_clim: climatology of clear-sky net TOA
        rtntcs_rolling_mean: 12-month rolling mean of rtntcs
        """
        # Net LW (change sign to the upwelling flux)
        for cube in self.flx_clim:
            if cube.var_name == "rlut":
                dcube = cube.copy()
                dcube.data = -dcube.data
                dcube.var_name = "rlnt"
                dcube.long_name = "radiative_flux_of_rlnt"
                self.flx_clim.append(dcube)
            if cube.var_name == "rsdt":
                rsdt_clim = cube.copy()
            if cube.var_name == "rsutcs":
                rsutcs_clim = cube.copy()
            if cube.var_name == "rlutcs":
                rlutcs_clim = cube.copy()
        rtntcs_clim = rsdt_clim - rsutcs_clim - rlutcs_clim
        rtntcs_clim.var_name = "rtntcs"
        rtntcs_clim.long_name = "radiative_flux_of_rtntcs"
        self.flx_clim.append(rtntcs_clim)
        # Annual rolling means clear-sky net total TOA
        for cube in self.flx_rolling_mean:
            if cube.var_name == "rsdt":
                rsdt_rolling_mean = cube.copy()
            if cube.var_name == "rsutcs":
                rsutcs_rolling_mean = cube.copy()
            if cube.var_name == "rlutcs":
                rlutcs_rolling_mean = cube.copy()
        rtntcs_rolling_mean = rsdt_rolling_mean - rsutcs_rolling_mean - \
            rlutcs_rolling_mean
        rtntcs_rolling_mean.var_name = "rtntcs"
        rtntcs_rolling_mean.long_name = "radiative_flux_of_rtntcs"
        self.flx_rolling_mean.append(rtntcs_rolling_mean)

    def print(self):
        """Print variable names of all cubes in the IHT object."""
        logger.info("=== implied_heat_transport object ===")
        print(self.mht_clim)
        for x in self.mht_clim:
            print(x.long_name, x.var_name)

        print(self.efp_clim)
        for x in self.efp_clim:
            print(x.long_name, x.var_name)

        print(self.flx_clim)
        for x in self.flx_clim:
            print(x.long_name, x.var_name)

        print(self.mht_rolling_mean)
        for x in self.mht_rolling_mean:
            print(x.long_name, x.var_name)

        print(self.symmetry_metric)
        for x in self.symmetry_metric:
            print(x.long_name, x.var_name)

        print(self.flx_files)

    def mht_symmetry_metrics(self):
        """Calculate the symmetry metrics for all time series of 12-month
        rolling means of MHT."""
        # As coded, the calculation of the symmetry metrics needs the number of
        # latitude points to be multiple of 6, i.e. it needs 30 deg bands.
        if (self.grid.shape[0] % 6) != 0:
            logger.error(
                "Grid not compatible with symmetry metric calculation.")
            sys.exit(1)
        Nlat_2 = self.grid.shape[0] // 2
        grid = np.sum(self.grid, axis=1)[Nlat_2:]

        for mht_series in self.mht_rolling_mean:
            time_coord = mht_series.coord('time')
            Ntime = time_coord.shape[0]
            hem = np.zeros(Ntime)
            trop = np.zeros(Ntime)
            extratrop = np.zeros(Ntime)
            for i in np.arange(Ntime):
                hem[i], trop[i], extratrop[i] = symmetry_metric(
                    mht_series.data[i], grid)
            # Create the cubes for each metric
            long_name = "symmetry_hemisphere_of_{}".format(
                mht_series.long_name)
            var_name = "s_hem_{}".format(mht_series.var_name)
            cube_h = iris.cube.Cube(hem / 1.0e15,
                                    long_name=long_name,
                                    var_name=var_name,
                                    units="PW",
                                    dim_coords_and_dims=[(time_coord, 0)])
            long_name = "symmetry_tropics_of_{}".format(mht_series.long_name)
            var_name = "s_tro_{}".format(mht_series.var_name)
            cube_t = iris.cube.Cube(trop / 1.0e15,
                                    long_name=long_name,
                                    var_name=var_name,
                                    units="PW",
                                    dim_coords_and_dims=[(time_coord, 0)])
            long_name = "symmetry_extratropics_of_{}".format(
                mht_series.long_name)
            var_name = "s_ext_{}".format(mht_series.var_name)
            cube_e = iris.cube.Cube(extratrop / 1.0e15,
                                    long_name=long_name,
                                    var_name=var_name,
                                    units="PW",
                                    dim_coords_and_dims=[(time_coord, 0)])
            self.symmetry_metric.append(cube_h)
            self.symmetry_metric.append(cube_t)
            self.symmetry_metric.append(cube_e)

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
        """Plots of CRE MHTs (Figure 3)."""
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
            ax2.plot(self.lat, -mht, label=legend_r[i], color=col[i])
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

    def quiver_start(self, n, step):
        """Calculate start point for quiver plot."""
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
                       nlevs,
                       label=[['(a)', '(b)'], ['(c)', '(d)'], ['(e)', '(f)']],
                       xy_label=(0, 1.05),
                       title=[['', ''], ['', ''], ['', '']],
                       change_sign=[[False, False], [False, False],
                                    [False, False]]):
        """Plot figures with energy flux potential and gradient in the left
        hand columns and the corresponding source term in the right-hand
        column.

        Figures 2, 4, and 5.
        """
        x, y = np.meshgrid(self.lon, self.lat)
        levels1 = np.linspace(vmin, vmax, nlevs)
        levels2 = np.linspace(wmin, wmax, nwlevs)
        nrows = len(var_name)
        # Calculate sampling for vector plot
        nlon = len(self.lon)
        nlat = len(self.lat)
        stepx = nlon // 20
        stepy = nlat // 10
        startx = self.quiver_start(nlon, stepx)
        starty = self.quiver_start(nlat, stepy)

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
                var_name_constraint(var_name[i][0]))
            flx = self.flx_clim.extract_cube(
                var_name_constraint(var_name[i][1]))
            efp.data -= np.average(efp.data)  # Arbitrary choice of origin
            flx.data -= area_average(flx).data
            if change_sign[i][0]:
                efp.data = -efp.data
            if change_sign[i][1]:
                flx.data = -flx.data
            v, u = np.gradient(efp.data, 1e14, 1e14)
            u = u[1:-1, 1:-1]
            v = v[1:-1, 1:-1]
            efp.convert_units("PW")
            ax1 = plt.subplot(gs[i * 7:(i * 7) + 7, 0],
                              projection=ccrs.PlateCarree(central_longitude=0))
            cb1 = iplt.contourf(efp, levels=levels1)
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
                              projection=ccrs.PlateCarree(central_longitude=0))
            cb2 = iplt.contourf(flx, levels=levels2, cmap='RdBu_r')
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

    def plot_symmetry_time_series(self):
        """Plot the all-sky and clear-sky time series of the symmetry metrics
        (Figure 6)"""
        var_list = [["s_hem_rtnt_mht", "s_hem_rtntcs_mht"],
                    ["s_tro_rtnt_mht", "s_tro_rtntcs_mht"],
                    ["s_ext_rtnt_mht", "s_ext_rtntcs_mht"]]
        col = ['C0', 'C1']
        label = [
            r'Global: 0$^\mathrm{o}$ - 90$^\mathrm{o}$',
            r'Tropics: 0$^\mathrm{o}$ - 30$^\mathrm{o}$',
            r'Extratropics: 30$^\mathrm{o}$ - 90$^\mathrm{o}$'
        ]
        legend_label = ["TOA net all-sky", "TOA net clear-sky"]

        plt.figure(figsize=(6, 12))
        for i, var_name in enumerate(var_list):
            y0 = self.symmetry_metric.extract_cube(
                var_name_constraint(var_name[0]))
            y1 = self.symmetry_metric.extract_cube(
                var_name_constraint(var_name[1]))
            ax = plt.subplot(3, 1, i + 1)
            iplt.plot(y0, lw=4, linestyle='-', label=legend_label[0])
            iplt.plot(y1, lw=4, linestyle='-', label=legend_label[1])
            ax.annotate(r'$\sigma$: {:5.3f}'.format(np.std(y0.data)),
                        (0.05, 0.55),
                        xycoords='axes fraction',
                        color=col[0])
            ax.annotate(r'$\sigma$: {:5.3f}'.format(np.std(y1.data)),
                        (0.05, 0.45),
                        xycoords='axes fraction',
                        color=col[1])
            ax.set_ylim(0, 0.8)
            ax.set_ylabel(r'$S$ (PW)')
            ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=1))
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
            ax.set_title(label[i])
            if i == 0:
                plt.legend(loc=5)
        plt.tight_layout()
        return


def efp_maps(iht, model, experiment, cfg):
    """Wrapper function that produces Figures 2, 4, and 5."""
    # Figure 2
    iht.quiver_subplot(
        [['rtnt_efp', 'rtnt'], ['rsnt_efp', 'rsnt'], ['rlnt_efp', 'rlnt']],
        wmin=-180,
        wmax=180,
        nwlevs=19,
        wlevstep=4,
        vmin=-1.2,
        vmax=1.2,
        nlevs=11,
        title=[['$P_{TOA}^{TOT}$', r'$\Delta F_{TOA}^{TOT}$'],
               ['$P_{TOA}^{SW}$', r'$\Delta F_{TOA}^{SW}$'],
               ['$P_{TOA}^{LW}$', r'$\Delta F_{TOA}^{LW}$']])
    provenance_record = get_provenance_record(plot_type='map',
                                              ancestor_files=iht.flx_files)
    figname = "figure2_{}_{}".format(model, experiment)
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
        nlevs=11,
        title=[['$P_{TOA}^{TOTCRE}$', r'$\Delta CRE_{TOA}^{TOT}$'],
               ['$P_{TOA}^{SWCRE}$', r'$\Delta CRE_{TOA}^{SW}$'],
               ['$P_{TOA}^{LWCRE}$', r'$\Delta CRE_{TOA}^{LW}$']])
    provenance_record = get_provenance_record(plot_type='map',
                                              ancestor_files=iht.flx_files)
    figname = "figure4_{}_{}".format(model, experiment)
    save_figure(figname, provenance_record, cfg)
    # Figure 5
    iht.quiver_subplot(
        [['rsutcs_efp', 'rsutcs'], ['rsut_efp', 'rsut']],
        wmin=-100,
        wmax=100,
        nwlevs=21,
        wlevstep=3,
        vmin=-0.35,
        vmax=0.35,
        nlevs=11,
        title=[['$P_{TOA}^{SWup, clr}$', r'$\Delta F_{TOA}^{SWup, clr}$'],
               ['$P_{TOA}^{SWup, all}$', r'$\Delta F_{TOA}^{SWup, all}$']],
        change_sign=[[True, True], [True, True]])
    provenance_record = get_provenance_record(plot_type='map',
                                              ancestor_files=iht.flx_files)
    figname = "figure5_{}_{}".format(model, experiment)
    save_figure(figname, provenance_record, cfg)


def mht_plots(iht, model, experiment, cfg):
    """Wrapper function that produces Figures 1 and 3."""
    # Figure 1
    iht.mht_plot(["rtnt_mht", "rsnt_mht", "rlnt_mht"], ['Net', 'SW', 'LW'])
    provenance_record = get_provenance_record(plot_type='zonal',
                                              ancestor_files=iht.flx_files)
    figname = "figure1_{}_{}".format(model, experiment)
    save_figure(figname, provenance_record, cfg)
    # Figure 3
    iht.cre_mht_plot(['netcre_mht', 'swcre_mht', 'lwcre_mht'],
                     ['Net CRE', 'SW CRE', 'LW CRE'],
                     ['rsut_mht', 'rsutcs_mht'],
                     ['-1 x OSR (all-sky)', '-1 x OSR (clear-sky)'])
    figname = "figure3_{}_{}".format(model, experiment)
    save_figure(figname, provenance_record, cfg)
    return


def symmetry_plots(iht, model, experiment, cfg):
    """Wrapper function that produces Figure 6."""
    # Figure 6
    iht.plot_symmetry_time_series()
    provenance_record = get_provenance_record(plot_type='times',
                                              ancestor_files=iht.flx_files)
    figname = "figure6_{}_{}".format(model, experiment)
    save_figure(figname, provenance_record, cfg)
    return


def plot_single_model_diagnostics(iht_dict, cfg):
    """Wrapper function that produces all plots."""
    # iht_dict is a two-level dictionary: iht_dict[model][experiment]
    for model, iht_model in iht_dict.items():
        logger.info("Plotting model: %s", model)
        for experiment, iht_experiment in iht_model.items():
            logger.info("Plotting experiment: %s", experiment)
            mht_plots(iht_experiment, model, experiment, cfg)
            efp_maps(iht_experiment, model, experiment, cfg)
            symmetry_plots(iht_experiment, model, experiment, cfg)
    return


def main(cfg):
    """Solve the Poisson equation and estimate the meridional heat
    transport."""
    input_data = deepcopy(list(cfg['input_data'].values()))
    input_data = group_metadata(input_data, 'dataset', sort='variable_group')

    # Arrange input flux files in a 2-level dictionary [model_name][dataset]
    flx_files = {}
    for model_name, datasets in input_data.items():
        flx_files[model_name] = {}
        for dataset in datasets:
            if dataset['dataset'] in flx_files[model_name]:
                flx_files[model_name][dataset['dataset']].append(
                    dataset['filename'])
            else:
                flx_files[model_name][dataset['dataset']] = [
                    dataset['filename']
                ]

    # Create dictionary of iht objects. It's a 2-level dictionary
    # like flx_files. This is where all the calculations are done.
    iht = {}
    for model_name, datasets in flx_files.items():
        logger.info("Model %s", model_name)
        iht[model_name] = {}
        for dataset_name, files in datasets.items():
            logger.info("Dataset %s", dataset_name)
            iht[model_name][dataset_name] = implied_heat_transport(files)

    print(iht)
    # Produce plots
    plot_single_model_diagnostics(iht, cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
