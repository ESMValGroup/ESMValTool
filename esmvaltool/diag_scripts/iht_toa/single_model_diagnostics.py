"""(C) Crown Copyright 2023, the Met Office.

Single model diagnostics
1. Solve the Poisson solver
2. Produce and save plots
"""

import datetime
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
from poisson_solver import SphericalPoisson

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
    """Return area-weighted average of a cube."""
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


def weight_zm(cube, latitude=None):
    """Weight zonal-mean by normalised gridbox areas."""
    if cube.coord('latitude').bounds is None:
        cube.coord('latitude').guess_bounds()
    cube_areas = cube.copy()
    cube_areas.data = iris.analysis.cartography.area_weights(cube,
                                                             normalize=True)
    if latitude is not None:
        cube = cube.intersection(latitude=latitude)
        cube_areas = cube_areas.intersection(latitude=latitude)
    return cube.data * cube_areas.data


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
    sphpo = SphericalPoisson(logger,
                             source=data * (6371e3**2.0),
                             tolerance=2.0e-4)
    sphpo.solve()
    sphpo.calc_mht()
    logger.info("Ending spherical_poisson")

    # Energy flux potential (P)
    p_cube = flux_cube.copy()
    p_cube.var_name = f"{flux_cube.var_name}_efp"
    p_cube.long_name = f"energy_flux_potential_of_{flux_cube.var_name}"
    p_cube.standard_name = None
    p_cube.units = 'J s-1'
    p_cube.data = sphpo.efp[1:-1, 1:-1]

    # MHT data cube
    mht_cube = flux_cube.copy()
    mht_cube = mht_cube.collapsed('longitude', iris.analysis.MEAN)
    mht_cube.var_name = f"{flux_cube.var_name}_mht"
    mht_cube.long_name = f"meridional_heat_transport_of_{flux_cube.var_name}"
    mht_cube.standard_name = None
    mht_cube.units = 'W'
    mht_cube.data = sphpo.mht

    return p_cube, mht_cube


def symmetry_metric(cube):
    """Calculate symmetry metrics.

    A perfectly symmetrical latitude band gives S=0. As coded, the
    calculation of the symmetry metrics needs the number of latitude
    points to be multiple of 6, i.e. it needs 30 deg bands. It returns
    the metric for 3 regions: globe, tropics and extratropics.
    """
    # Hemisphere
    hem = np.abs(weight_zm(cube, latitude=(0, 90, False, False))[::-1] +
                 weight_zm(cube, latitude=(-90, 0, False, False))).sum()
    # Tropics
    tro = np.abs(weight_zm(cube, latitude=(0, 30, False, False))[::-1] +
                 weight_zm(cube, latitude=(-30, 0, False, False))).sum()
    # Extra-tropics
    etr = np.abs(weight_zm(cube, latitude=(30, 90, False, False))[::-1] +
                 weight_zm(cube, latitude=(-90, -30, False, False))).sum()
    return hem, tro, etr


def format_plot(axx, label, title):
    """Format plots in quiver panel."""
    axx.set_xticks(np.arange(-180, 190, 60))
    axx.set_xticklabels(['180', '120W', '60W', '0', '60E', '120E', '180'])
    axx.set_yticks(np.arange(-90, 100, 30))
    axx.set_yticklabels(['90S', '60S', '30S', 'Eq', '30N', '60N', '90N'])
    axx.annotate(label, xy=(0, 1.05), xycoords='axes fraction', color='k')
    axx.set_title(title)


class ImpliedHeatTransport:
    """Class that solves IHT for a given dataset."""

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

        # Compute derived fluxes
        self.derived_fluxes()

        # Calculate EFP and MHT for each flux
        self.compute_efp_and_mht()

        # Times series of MHT symmetry metric
        self.mht_symmetry_metrics()

        self.print()

    def compute_efp_and_mht(self):
        """Calculate Energy Flux Potential and meridional heat transport.

        Calculate EFP and MHT for the climatologies of radiative fluxes
        and the 12-month rolling means of radiative fluxes.
        """
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
        """Calculate derived radiative fluxes.

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
        for zzz in self.mht_clim:
            print(zzz.long_name, zzz.var_name)

        print(self.efp_clim)
        for zzz in self.efp_clim:
            print(zzz.long_name, zzz.var_name)

        print(self.flx_clim)
        for zzz in self.flx_clim:
            print(zzz.long_name, zzz.var_name)

        print(self.mht_rolling_mean)
        for zzz in self.mht_rolling_mean:
            print(zzz.long_name, zzz.var_name)

        print(self.symmetry_metric)
        for zzz in self.symmetry_metric:
            print(zzz.long_name, zzz.var_name)

        print(self.flx_files)

    def mht_symmetry_metrics(self):
        """Calculate symmetry metrics.

        Produce 12-month rolling means for all monthly time time series
        of MHT.
        """
        for mht_series in self.mht_rolling_mean:
            time_coord = mht_series.coord('time')
            ntime = time_coord.shape[0]
            hem = np.zeros(ntime)
            trop = np.zeros(ntime)
            extratrop = np.zeros(ntime)
            for i in np.arange(ntime):
                hem[i], trop[i], extratrop[i] = symmetry_metric(mht_series[i])
            # Create the cubes for each metric
            long_name = f"symmetry_hemisphere_of_{mht_series.long_name}"
            var_name = f"s_hem_{mht_series.var_name}"
            cube_h = iris.cube.Cube(hem / 1.0e15,
                                    long_name=long_name,
                                    var_name=var_name,
                                    units="PW",
                                    dim_coords_and_dims=[(time_coord, 0)])
            long_name = f"symmetry_tropics_of_{mht_series.long_name}"
            var_name = f"s_tro_{mht_series.var_name}"
            cube_t = iris.cube.Cube(trop / 1.0e15,
                                    long_name=long_name,
                                    var_name=var_name,
                                    units="PW",
                                    dim_coords_and_dims=[(time_coord, 0)])
            long_name = f"symmetry_extratropics_of_{mht_series.long_name}"
            var_name = f"s_ext_{mht_series.var_name}"
            cube_e = iris.cube.Cube(extratrop / 1.0e15,
                                    long_name=long_name,
                                    var_name=var_name,
                                    units="PW",
                                    dim_coords_and_dims=[(time_coord, 0)])
            self.symmetry_metric.append(cube_h)
            self.symmetry_metric.append(cube_t)
            self.symmetry_metric.append(cube_e)

    def mht_plot(self, var_names, legend_label, ylim=(-10, 10)):
        """Produce a single multi-line plot with MHT.

        MHT is presented in PW. Up to three variables are on each plot.
        """
        plt.figure()
        for i, vname in enumerate(var_names):
            mht = self.mht_clim.extract_cube(var_name_constraint(vname))
            plt.plot(mht.coord('latitude').points,
                     mht.data / 1e15,
                     label=legend_label[i])
        plt.hlines(0, -90, 90, color='k', linestyles=':')
        plt.vlines(0, -10, 10, color='k', linestyles=':')
        plt.xlim(-90, 90)
        plt.ylim(ylim[0], ylim[1])
        plt.xticks(np.arange(-90, 120, 30))
        plt.xlabel('Latitude')
        plt.ylabel('MHT (PW)')
        plt.legend()
        plt.tight_layout()

    def cre_mht_plot(self, left, right, ylim=(-1.5, 1.5)):
        """Produce two multiline plots of MHT."""
        plt.figure(figsize=(11, 5))
        ax1 = plt.subplot(121)
        for i, vname in enumerate(left['vname']):
            mht = self.mht_clim.extract_cube(var_name_constraint(vname))
            ax1.plot(mht.coord('latitude').points,
                     mht.data / 1e15,
                     label=left['legend'][i])
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
        for i, vname in enumerate(right['vname']):
            mht = self.mht_clim.extract_cube(var_name_constraint(vname))
            ax2.plot(mht.coord('latitude').points,
                     -mht.data / 1e15,
                     label=right['legend'][i],
                     color=col[i])
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

    def quiver_start(self, ntot, step):
        """Calculate start point for quiver plot."""
        nnn = (ntot - 2) // step
        start = (ntot - 2 - nnn * step) // 2
        return start

    def quiver_maps_data(self, vnames, change_sign):
        """Obtain data for one row of plots."""
        efp = self.efp_clim.extract_cube(var_name_constraint(vnames[0]))
        flx = self.flx_clim.extract_cube(var_name_constraint(vnames[1]))
        efp.data -= np.average(efp.data)  # Arbitrary choice of origin
        flx.data -= area_average(flx).data
        if change_sign[0]:
            efp.data = -efp.data
        if change_sign[1]:
            flx.data = -flx.data
        vvv, uuu = np.gradient(efp.data, 1e14, 1e14)
        uuu = uuu[1:-1, 1:-1]
        vvv = vvv[1:-1, 1:-1]
        efp.convert_units("PW")
        return {'efp': efp, 'flx': flx, 'uuu': uuu, 'vvv': vvv}

    def quiver_subplot(self, dargs, title, label, change_sign):
        """Produce panel with EFPs (left column) and fluxes (right column).

        Plot figures with energy flux potential and gradient in the left
        hand columns and the corresponding source term in the right-hand
        column.
        """
        mshgrd = np.meshgrid(self.flx_clim[0].coord('longitude').points,
                             self.flx_clim[0].coord('latitude').points)
        nrows = len(dargs['var_name'])
        # Calculate sampling for vector plot
        dxy = [mshgrd[0].shape[1] // 20, mshgrd[0].shape[0] // 10]
        startx = self.quiver_start(mshgrd[0].shape[1], dxy[0])
        starty = self.quiver_start(mshgrd[0].shape[0], dxy[1])

        if nrows == 3:
            plt.figure(figsize=(10, 10))
            grds = gridspec.GridSpec(22, 2)
            grds.update(wspace=0.25, hspace=1.5)
        elif nrows == 2:
            plt.figure(figsize=(10, 6.5))
            grds = gridspec.GridSpec(15, 2)
            grds.update(wspace=0.25, hspace=1.5)
        elif nrows == 1:
            plt.figure(figsize=(12, 4))
            grds = gridspec.GridSpec(8, 2)
            grds.update(wspace=0.25, hspace=1.5)

        cbs = []
        for i in range(nrows):
            data = self.quiver_maps_data(dargs['var_name'][i], change_sign[i])
            plt.subplot(grds[i * 7:(i * 7) + 7, 0],
                        projection=ccrs.PlateCarree(central_longitude=0))
            cbs.append(
                iplt.contourf(data['efp'],
                              levels=np.linspace(dargs['vmin'], dargs['vmax'],
                                                 dargs['nlevs'])))
            plt.gca().coastlines()
            if i == 0:
                qqq = plt.quiver(mshgrd[0][starty::dxy[1], startx::dxy[0]],
                                 mshgrd[1][starty::dxy[1], startx::dxy[0]],
                                 data['uuu'][starty::dxy[1], startx::dxy[0]],
                                 data['vvv'][starty::dxy[1], startx::dxy[0]],
                                 pivot='mid',
                                 color='w',
                                 width=0.005)
            else:
                plt.quiver(mshgrd[0][starty::dxy[1], startx::dxy[0]],
                           mshgrd[1][starty::dxy[1], startx::dxy[0]],
                           data['uuu'][starty::dxy[1], startx::dxy[0]],
                           data['vvv'][starty::dxy[1], startx::dxy[0]],
                           pivot='mid',
                           scale=qqq.scale,
                           color='w')
            format_plot(plt.gca(), label[i][0], title[i][0])

            plt.subplot(grds[i * 7:(i * 7) + 7, 1],
                        projection=ccrs.PlateCarree(central_longitude=0))
            cbs.append(
                iplt.contourf(data['flx'],
                              levels=np.linspace(dargs['wmin'], dargs['wmax'],
                                                 dargs['nwlevs']),
                              cmap='RdBu_r'))
            plt.gca().coastlines()
            format_plot(plt.gca(), label[i][1], title[i][1])

        plt.subplot(grds[-1, 0])
        plt.colorbar(cbs[0],
                     cax=plt.gca(),
                     orientation='horizontal',
                     label='Energy flux potential (PW)')
        plt.subplot(grds[-1, 1])
        plt.colorbar(cbs[1],
                     cax=plt.gca(),
                     orientation='horizontal',
                     label=r'Flux (Wm$^{-2}$)',
                     ticks=np.linspace(dargs['wmin'], dargs['wmax'],
                                       dargs['nwlevs'])[1::dargs['wlevstep']])

        if nrows == 3:
            plt.subplots_adjust(left=0.1, right=0.94, top=1.0, bottom=0.11)
        elif nrows == 2:
            plt.subplots_adjust(left=0.11, right=0.9, top=1.0, bottom=0.13)
        elif nrows == 1:
            plt.subplots_adjust(left=0.11, right=0.9, top=1.0, bottom=0.20)

    def plot_symmetry_time_series(self):
        """Produce Figure 6.

        All-sky and clear-sky time series of the symmetry metrics for
        three regions: globe, tropics and extratropics.
        """
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
            yy0 = self.symmetry_metric.extract_cube(
                var_name_constraint(var_name[0]))
            yy1 = self.symmetry_metric.extract_cube(
                var_name_constraint(var_name[1]))
            axx = plt.subplot(3, 1, i + 1)
            dtx = [
                datetime.datetime.strptime(str(cell[0]), '%Y-%m-%d %H:%M:%S')
                for cell in yy0.coord('time').cells()
            ]
            plt.plot(dtx, yy0.data, lw=4, linestyle='-', label=legend_label[0])
            plt.plot(dtx, yy1.data, lw=4, linestyle='-', label=legend_label[1])
            axx.annotate(rf'$\sigma$: {np.std(yy0.data):5.3f}', (0.05, 0.55),
                         xycoords='axes fraction',
                         color=col[0])
            axx.annotate(rf'$\sigma$: {np.std(yy1.data):5.3f}', (0.05, 0.45),
                         xycoords='axes fraction',
                         color=col[1])
            axx.set_ylim(0, 0.8)
            axx.set_ylabel(r'$S$ (PW)')
            axx.xaxis.set_major_locator(mdates.YearLocator(3, month=1, day=1))
            axx.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
            axx.xaxis.set_minor_locator(mdates.YearLocator())
            axx.set_title(label[i])
            if i == 0:
                plt.legend(loc=5)
        plt.tight_layout()


def efp_maps(iht, model, experiment, cfg):
    """Produce Figures 2, 4, and 5."""
    # Figure 2
    iht.quiver_subplot(
        {
            'var_name': [['rtnt_efp', 'rtnt'], ['rsnt_efp', 'rsnt'],
                         ['rlnt_efp', 'rlnt']],
            'wmin':
            -180,
            'wmax':
            180,
            'nwlevs':
            19,
            'wlevstep':
            4,
            'vmin':
            -1.2,
            'vmax':
            1.2,
            'nlevs':
            11
        },
        title=[['$P_{TOA}^{TOT}$', r'$\Delta F_{TOA}^{TOT}$'],
               ['$P_{TOA}^{SW}$', r'$\Delta F_{TOA}^{SW}$'],
               ['$P_{TOA}^{LW}$', r'$\Delta F_{TOA}^{LW}$']],
        label=[['(a)', '(b)'], ['(c)', '(d)'], ['(e)', '(f)']],
        change_sign=[[False, False], [False, False], [False, False]])
    provenance_record = get_provenance_record(plot_type='map',
                                              ancestor_files=iht.flx_files)
    figname = f"figure2_{model}_{experiment}"
    save_figure(figname, provenance_record, cfg)
    # Figure 4
    iht.quiver_subplot(
        {
            'var_name': [['netcre_efp', 'netcre'], ['swcre_efp', 'swcre'],
                         ['lwcre_efp', 'lwcre']],
            'wmin':
            -60,
            'wmax':
            60,
            'nwlevs':
            13,
            'wlevstep':
            2,
            'vmin':
            -0.3,
            'vmax':
            0.3,
            'nlevs':
            11
        },
        title=[['$P_{TOA}^{TOTCRE}$', r'$\Delta CRE_{TOA}^{TOT}$'],
               ['$P_{TOA}^{SWCRE}$', r'$\Delta CRE_{TOA}^{SW}$'],
               ['$P_{TOA}^{LWCRE}$', r'$\Delta CRE_{TOA}^{LW}$']],
        label=[['(a)', '(b)'], ['(c)', '(d)'], ['(e)', '(f)']],
        change_sign=[[False, False], [False, False], [False, False]])
    provenance_record = get_provenance_record(plot_type='map',
                                              ancestor_files=iht.flx_files)
    figname = f"figure4_{model}_{experiment}"
    save_figure(figname, provenance_record, cfg)
    # Figure 5
    iht.quiver_subplot(
        {
            'var_name': [['rsutcs_efp', 'rsutcs'], ['rsut_efp', 'rsut']],
            'wmin': -100,
            'wmax': 100,
            'nwlevs': 21,
            'wlevstep': 3,
            'vmin': -0.35,
            'vmax': 0.35,
            'nlevs': 11
        },
        title=[['$P_{TOA}^{SWup, clr}$', r'$\Delta F_{TOA}^{SWup, clr}$'],
               ['$P_{TOA}^{SWup, all}$', r'$\Delta F_{TOA}^{SWup, all}$']],
        label=[['(a)', '(b)'], ['(c)', '(d)']],
        change_sign=[[True, True], [True, True]])
    provenance_record = get_provenance_record(plot_type='map',
                                              ancestor_files=iht.flx_files)
    figname = f"figure5_{model}_{experiment}"
    save_figure(figname, provenance_record, cfg)


def mht_plots(iht, model, experiment, cfg):
    """Produce Figures 1 and 3."""
    # Figure 1
    iht.mht_plot(["rtnt_mht", "rsnt_mht", "rlnt_mht"], ['Net', 'SW', 'LW'])
    provenance_record = get_provenance_record(plot_type='zonal',
                                              ancestor_files=iht.flx_files)
    figname = f"figure1_{model}_{experiment}"
    save_figure(figname, provenance_record, cfg)
    # Figure 3
    iht.cre_mht_plot(
        {
            'vname': ['netcre_mht', 'swcre_mht', 'lwcre_mht'],
            'legend': ['Net CRE', 'SW CRE', 'LW CRE']
        }, {
            'vname': ['rsut_mht', 'rsutcs_mht'],
            'legend': ['-1 x OSR (all-sky)', '-1 x OSR (clear-sky)']
        })
    figname = f"figure3_{model}_{experiment}"
    save_figure(figname, provenance_record, cfg)


def symmetry_plots(iht, model, experiment, cfg):
    """Produce Figure 6."""
    iht.plot_symmetry_time_series()
    provenance_record = get_provenance_record(plot_type='times',
                                              ancestor_files=iht.flx_files)
    figname = f"figure6_{model}_{experiment}"
    save_figure(figname, provenance_record, cfg)


def plot_single_model_diagnostics(iht_dict, cfg):
    """Produce plots for a single model and experiment.

    iht_dict is a two-level dictionary: iht_dict[model][experiment]
    """
    for model, iht_model in iht_dict.items():
        logger.info("Plotting model: %s", model)
        for experiment, iht_experiment in iht_model.items():
            logger.info("Plotting experiment: %s", experiment)
            mht_plots(iht_experiment, model, experiment, cfg)
            efp_maps(iht_experiment, model, experiment, cfg)
            symmetry_plots(iht_experiment, model, experiment, cfg)


def main(cfg):
    """Produce implied heat transport plots.

    Produce Figures 1 to 6 of Pearce and Bodas-Salcedo (2023) for each
    model and dataset combination.
    """
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
            iht[model_name][dataset_name] = ImpliedHeatTransport(files)

    print(iht)
    # Produce plots
    plot_single_model_diagnostics(iht, cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
