# (C) Crown Copyright 2023, the Met Office.
""" Single model diagnostics
1. Solve the Poisson solver over the ocean
2. Produce and save plots
"""

import datetime
import logging
from copy import deepcopy

import cartopy.crs as ccrs
import iris
import iris.plot as iplt
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from iris import NameConstraint
from matplotlib import gridspec, rcParams
from poisson_solver_oht import SphericalPoisson
from scipy.ndimage import binary_fill_holes

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_figure,
)

# Initialise logger
logger = logging.getLogger(__name__)

# rcParams.update({
#     'font.size': 14,
#     'xtick.major.pad': 10,
#     'ytick.major.pad': 10,
#     'xtick.major.size': 10,
#     'ytick.major.size': 10,
#     'xtick.minor.size': 5,
#     'ytick.minor.size': 5,
#     'axes.linewidth': 2,
#     'lines.markersize': 8,
#     'lines.linewidth': 2
# })

# Figure captions
caption = {
    'F1': 'Figure 1. The implied heat transport in the global ocean'
          'and in individual basins.',
    'F2': 'Figure 2.'
}


def get_provenance_record(filenames, figure_caption):
    """Return a provenance record describing the plot.

    Parameters
    ----------
    filenames : list of strings
        The filenames containing the data used to create the plot.
    figure_caption : string
        Detailed description of the figure.

    Returns
    -------
    dictionary
        The provenance record describing the plot.
    """
    record = {
        'ancestors': filenames,
        'caption': figure_caption,
        'references': ['pearce23jclim']
    }
    return record


def matching_strings(list_of_strings, substrings):
    """Return subset of ``list_of_strings`` with matches in ``substrings``.

    Parameters
    ----------
    list_of_strings : list of strings
        List of strings to be searched.
    substrings : list of strings
        The list of search strings.

    Returns
    -------
    list
        The elements in ``list_of_strings`` that contain
        any of the substrings.
    """
    matches = []
    for element in list_of_strings:
        for var in substrings:
            if var in element:
                matches.append(element)
    return matches


def area_average(cube, latitude='latitude', longitude='longitude', mdtol=1):
    """Return area-weighted average of a cube.

    Parameters
    ----------
    cube : :class:`iris.cube.Cube`
        Input cube.
    latitude : string
        Name of latitude coordinate in ``cube``.
    longitude : string
        Name of longitude coordinate in ``cube``.
    mdtol : float
        Tolerance to missing data, between 0 and 1.


    Returns
    -------
    :class:`iris.cube.Cube`
        Collapsed cube with the weighted average.
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

def plot_array(x,figname):
    plt.figure()
    plt.pcolormesh(x)
    plt.colorbar()
    plt.savefig(figname)
    print(figname)
    print(x.min(), x.max(), x.mean())

def call_poisson(flux_cube, mask, latitude='latitude', longitude='longitude'):
    """Call the Poisson solver with the data in ``flux_cube`` as source term.

       Return the energy flux potential and implied meridional heat transport
       as cubes.

    Parameters
    ----------
    flux_cube : :class:`iris.cube.Cube`
        Input cube.
    latitude : string
        Name of latitude coordinate in ``cube``.
    longitude : string
        Name of longitude coordinate in ``cube``.

    Returns
    -------
    efp_cube: :class:`iris.cube.Cube`
        Energy flux potential cube.
    mht_cube: :class:`iris.cube.Cube`
        Implied meridional heat transport associated
        with the source flux field.
    """
    earth_radius = 6371e3  # Earth's radius in m
    if flux_cube.coord(latitude).bounds is None:
        flux_cube.coord(latitude).guess_bounds()
    if flux_cube.coord(longitude).bounds is None:
        flux_cube.coord(longitude).guess_bounds()

    # Remove average of flux field to account for storage term
    grid_areas = iris.analysis.cartography.area_weights(flux_cube, normalize=True)
    data_mean = np.average(flux_cube.data,
                           weights=np.ma.array(grid_areas, mask=(1-mask)))
    data = flux_cube.data - data_mean
    data *= mask
    plot_array(data, "data.png")

    logger.info("Calling spherical_poisson")
    sphpo = SphericalPoisson(logger,
                             source=data * (earth_radius**2.0),
                             mask=mask,
                             tolerance=2.0e-4)
    sphpo.solve()
    sphpo.calc_meridional_heat_transport()
    logger.info("Ending spherical_poisson")

    # Energy flux potential
    efp_cube = iris.cube.Cube(sphpo.energy_flux_potential[1:-1, 1:-1],
                              long_name=f"energy_flux_potential"
                                        f"_of_{flux_cube.var_name}",
                              var_name=f"{flux_cube.var_name}_efp",
                              units='J s-1',
                              dim_coords_and_dims=[(flux_cube.coords()[0], 0),
                                                   (flux_cube.coords()[1], 1)])

    # MHT data cube
    collapsed_longitude = iris.coords.AuxCoord(180.0,
                                               bounds=(0.0, 360.0),
                                               long_name='longitude',
                                               standard_name='longitude',
                                               units='degrees')
    dim_coords_and_dims = [(flux_cube.coord('latitude'), 0)]
    aux_coords_and_dims = [(flux_cube.coord('time'), None),
                           (collapsed_longitude, None)]
    mht_cube = iris.cube.Cube(sphpo.meridional_heat_transport,
                              long_name=f"meridional_heat_transport_of"
                                        f"_{flux_cube.var_name}",
                              var_name=f"{flux_cube.var_name}_mht",
                              units='W',
                              dim_coords_and_dims=dim_coords_and_dims,
                              aux_coords_and_dims=aux_coords_and_dims)
    return efp_cube, mht_cube


def format_plot(axes, label, title):
    """Format plots in quiver panel.

    Parameters
    ----------
    axes : :class:`matplotlib.axes.Axes`
        Input axes.
    label : string
        Top-left plot label.
    title : string
        Plot title.
    """
    axes.set_xticks(np.arange(-180, 190, 60))
    axes.set_xticklabels(['180', '120W', '60W', '0', '60E', '120E', '180'])
    axes.set_yticks(np.arange(-90, 100, 30))
    axes.set_yticklabels(['90S', '60S', '30S', 'Eq', '30N', '60N', '90N'])
    axes.annotate(label, xy=(0, 1.05), xycoords='axes fraction', color='k')
    axes.set_title(title)


class OceanHeatTransport:
    """Class that solves implied heat transport for an input dataset.

    These are the physical meanings of the main acronyms
    used in the variable names:
       FLX: radiative flux
       EFP: energy flux potential
       MHT: meridional heat transport
    """

    def __init__(self, flx_files, config):
        """Calculate all the diagnostics for all fluxes in ``flx_files``.

        Parameters
        ----------
        flx_files : list
            List of files with input data.
        """
        self.flx_files = flx_files

        # Create cube lists for the different datasets
        self.flx_clim = iris.cube.CubeList()
        self.mht_clim = iris.cube.CubeList()
        self.efp_clim = iris.cube.CubeList()

        # Read in climatologies
        for flx_file in flx_files:
            flx = iris.load_cube(flx_file)
            self.flx_clim.append(flx)

        # Read in mask
        self.mask = masks(config)

        # Compute derived fluxes
        self.derived_fluxes()

        # Calculate Energy Flux Potential and Meridional Heat Transport
        # for each flux component
        self.compute_efp_and_mht()

        self.print()

    def derived_fluxes(self):
        """Calculate derived radiative fluxes.

        nsf_clim: climatology of net heat surface flux
        """
        # Derived net surface flux
        rsds_clim = self.flx_clim.extract_cube(
            NameConstraint(var_name="rsds"))
        rsus_clim = self.flx_clim.extract_cube(
            NameConstraint(var_name="rsus"))
        rlds_clim = self.flx_clim.extract_cube(
            NameConstraint(var_name="rlds"))
        rlus_clim = self.flx_clim.extract_cube(
            NameConstraint(var_name="rlus"))
        hfss_clim = self.flx_clim.extract_cube(
            NameConstraint(var_name="hfss"))
        hfls_clim = self.flx_clim.extract_cube(
            NameConstraint(var_name="hfls"))
        nsf_clim = (rsds_clim - rsus_clim) + \
                   (rlds_clim - rlus_clim) - \
                   (hfss_clim + hfls_clim)
        nsf_clim = rlus_clim.copy() # Hack
        nsf_clim.var_name = "nsf"
        nsf_clim.long_name = "downward_heat_flux_in_air"
        self.flx_clim.append(nsf_clim)

    def compute_efp_and_mht(self):
        """Calculate Energy Flux Potential and meridional heat transport.

        Loop over input data and calculate EFP and MHT of the
        climatologies of radiative fluxes and the 12-month
        rolling means of radiative fluxes.
        """
        hfls_clim = self.flx_clim.extract_cube(NameConstraint(var_name="nsf"))
        efp, mht = call_poisson(hfls_clim, self.mask)
        self.efp_clim.append(efp)
        self.mht_clim.append(mht)
        # # Loop over climatologies
        # for flx in self.flx_clim:
        #     efp, mht = call_poisson(flx, self.mask)
        #     self.efp_clim.append(efp)
        #     self.mht_clim.append(mht)

    def print(self):
        """Print variable names of all cubes in an IHT object."""
        logger.info("=== implied_heat_transport object ===")
        logger.info(self.mht_clim)
        info_message = "Long name: %s; Variable: %s."
        for climatology in self.mht_clim:
            logger.info(info_message,
                        climatology.long_name,
                        climatology.var_name)

        logger.info(self.efp_clim)
        for climatology in self.efp_clim:
            logger.info(info_message,
                        climatology.long_name,
                        climatology.var_name)

        logger.info(self.flx_clim)
        for climatology in self.flx_clim:
            logger.info(info_message,
                        climatology.long_name,
                        climatology.var_name)

        logger.info(self.flx_files)

    def mht_plot(self, var_names, legend_label, ylim=(-10, 10)):
        """Produce a single multi-line plot of MHT components.

        MHT is presented in PW. Up to three variables are on each plot.

        Parameters
        ----------
        var_names : list
            Variable names to plot, e.g. ["rtnt_mht", "rsnt_mht"].
        legend_label : list
            List of labels for each line.
        ylim : tuple
            y axis limits.
        """
        plt.figure()
        for i, vname in enumerate(var_names):
            mht = self.mht_clim.extract_cube(NameConstraint(var_name=vname))
            mht.convert_units('PW')
            plt.plot(mht.coord('latitude').points,
                     mht.data,
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

    def quiver_start(self, ntot, step):
        """Calculate start point for quiver plot.

        Parameters
        ----------
        ntot : int
            Total number of points of the full vector.
        step : int
            Sampling step.
        """
        start = (ntot - 2 - ((ntot - 2) // step) * step) // 2
        return start

    def quiver_maps_data(self, vnames, change_sign):
        """Obtain data for one row of plots.

        Parameters
        ----------
        vnames : list
            Two-element list with the names of the EFP and
            flux variables.
        change_sign : list
            Two-element list of booleans to indicate if
            the variable has to be plotted with the sign changed.
        """
        efp = self.efp_clim.extract_cube(NameConstraint(var_name=vnames[0]))
        flx = self.flx_clim.extract_cube(NameConstraint(var_name=vnames[1]))
        # The choice of origin for efp is arbitrary,
        # we choose the unweighted mean.
        efp = efp - efp.collapsed(efp.coords(), iris.analysis.MEAN)
        flx = flx - area_average(flx)
        if change_sign[0]:
            efp = -efp
        if change_sign[1]:
            flx = -flx
        efp.convert_units("PW")
        v_component, u_component = np.gradient(efp.data)
        u_component = u_component[1:-1, 1:-1]
        v_component = v_component[1:-1, 1:-1]
        return {'efp': efp, 'flx': flx, 'uuu': u_component, 'vvv': v_component}

    def quiver_subplot(self, dargs):
        """Produce panel with maps of EFPs and fluxes.

        Plot figures with energy flux potential and gradient in the left-hand
        column and the corresponding source term in the right-hand column.

        Parameters
        ----------
        dargs : dictionary
            Dictionary with variable names and plot configuration
            data.
        """
        mshgrd = np.meshgrid(self.flx_clim[0].coord('longitude').points,
                             self.flx_clim[0].coord('latitude').points)
        nrows = len(dargs['var_name'])
        # Calculate sampling for vector plot
        dxy = [mshgrd[0].shape[1] // 20, mshgrd[0].shape[0] // 10]
        startx = self.quiver_start(mshgrd[0].shape[1], dxy[0])
        starty = self.quiver_start(mshgrd[0].shape[0], dxy[1])

        # Set grid layout depending on number of rows.
        # Place figures every grid_step rows in the grid.
        grid_step = 7
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
            data = self.quiver_maps_data(dargs['var_name'][i],
                                         dargs['change_sign'][i])
            plt.subplot(grds[i * grid_step:(i * grid_step) + grid_step, 0],
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
            format_plot(plt.gca(), dargs['label'][i][0], dargs['title'][i][0])

            plt.subplot(grds[i * grid_step:(i * grid_step) + grid_step, 1],
                        projection=ccrs.PlateCarree(central_longitude=0))
            cbs.append(
                iplt.contourf(data['flx'],
                              levels=np.linspace(dargs['wmin'], dargs['wmax'],
                                                 dargs['nwlevs']),
                              cmap='RdBu_r'))
            plt.gca().coastlines()
            format_plot(plt.gca(), dargs['label'][i][1], dargs['title'][i][1])

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


def masks(config):
    mask = iris.load_cube(config['mask'])

    # Close any seas not connected to the oceans
    # e.g. Mediterranean Sea at some resolutions
    (y, x) = np.shape(mask)
    wrap_mask = np.zeros([y, x + 2])
    wrap_mask[:, 1:-1] = mask.data
    wrap_mask[:, 0] = mask.data[:, -1]
    wrap_mask[:, -1] = mask.data[:, 0]
    wrap_mask = binary_fill_holes(wrap_mask)
    #wrap_mask[:,:] = 0

    plt.figure()
    plt.pcolormesh(mask.data)
    plt.colorbar()
    plt.savefig('mask.png')

    plt.figure()
    plt.pcolormesh(wrap_mask)
    plt.colorbar()
    plt.savefig('wrap_mask.png')

    wrap_mask = 1 - wrap_mask
    return wrap_mask[:, 1:-1]

def efp_maps(iht, model, experiment, config):
    """Produce Figure 2.

    Parameters
    ----------
    iht : :class: OceanHeatTransport
        Object with the recipe datasets.
    model : string
        Model name.
    experiment : string
        Experiment name.
    config : dict
        The ESMValTool configuration.
    """
    # Figure 2
    iht.quiver_subplot(
        {
            'var_name': [['nsf_efp', 'nsf']],
            'title': [['$P_{SFC}^{TOT}$', r'$\Delta F_{SFC}^{TOT}$']],
            'label': [['(a)', '(b)']],
            'change_sign': [[False, False]],
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
        })
    flx_files = matching_strings(iht.flx_files, ['nsf/'])
    provenance_record = get_provenance_record(flx_files, caption['F2'])
    figname = f"figure2_{model}_{experiment}"
    save_figure(figname, provenance_record, config)


def mht_plots(iht, model, experiment, config):
    """Produce Figures 1 and 3.

    Parameters
    ----------
    iht : :class: OceanHeatTransport
        Object with the recipe datasets.
    model : string
        Model name.
    experiment : string
        Experiment name.
    config : dict
        The ESMValTool configuration.
    """
    # Figure 1
    iht.mht_plot(["nsf_mht"], ['Net'])
    flx_files = matching_strings(iht.flx_files, ['nsf/'])
    provenance_record = get_provenance_record(flx_files, caption['F1'])
    figname = f"figure1_{model}_{experiment}"
    save_figure(figname, provenance_record, config)


def plot_single_model_diagnostics(iht_dict, config):
    """Produce plots for a single model and experiment.

    Parameters
    ----------
    iht_dict : dict
        iht_dict is a two-level dictionary: iht_dict[model][experiment]
    config : dict
        The ESMValTool configuration.
    """
    for model, iht_model in iht_dict.items():
        logger.info("Plotting model: %s", model)
        for experiment, iht_experiment in iht_model.items():
            logger.info("Plotting experiment: %s", experiment)
            mht_plots(iht_experiment, model, experiment, config)
            efp_maps(iht_experiment, model, experiment, config)


def main(config):
    """Produce all the recipe's plots.

    Parameters
    ----------
    config : dict
        The ESMValTool configuration.
    """
    input_data = deepcopy(list(config['input_data'].values()))
    input_data = group_metadata(input_data, 'dataset', sort='variable_group')

    # Arrange input flux files in a 2-level dictionary [model_name][dataset]
    flux_files = {}
    for model_name, datasets in input_data.items():
        flux_files[model_name] = {}
        for dataset in datasets:
            if dataset['dataset'] in flux_files[model_name]:
                flux_files[model_name][dataset['dataset']].append(
                    dataset['filename'])
            else:
                flux_files[model_name][dataset['dataset']] = [
                    dataset['filename']
                ]

    # Create dictionary of implied_heat_transport objects.
    # It's a 2-level dictionary like flux_files.
    # This is where all the calculations are done.
    oht = {}
    for model_name, datasets in flux_files.items():
        logger.info("Model %s", model_name)
        oht[model_name] = {}
        for dataset_name, files in datasets.items():
            logger.info("Dataset %s", dataset_name)
            oht[model_name][dataset_name] = OceanHeatTransport(files, config)

    # Produce plots
    plot_single_model_diagnostics(oht, config)


if __name__ == '__main__':

    with run_diagnostic() as configuration:
        main(configuration)
