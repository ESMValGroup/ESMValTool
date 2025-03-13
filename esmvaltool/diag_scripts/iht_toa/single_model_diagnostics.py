# (C) Crown Copyright 2023, the Met Office.
"""Single model diagnostics.

Apply Poisson solver to input fluxes and produce plots.
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

from esmvaltool.diag_scripts.iht_toa.poisson_solver import SphericalPoisson
from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_figure,
)

# Initialise logger
logger = logging.getLogger(__name__)

rcParams.update(
    {
        "font.size": 14,
        "xtick.major.pad": 10,
        "ytick.major.pad": 10,
        "xtick.major.size": 10,
        "ytick.major.size": 10,
        "xtick.minor.size": 5,
        "ytick.minor.size": 5,
        "axes.linewidth": 2,
        "lines.markersize": 8,
        "lines.linewidth": 2,
    }
)

# Figure captions
caption = {
    "F1": "Figure 1. The implied heat transport due to TOA net flux (blue), "
    "split into the contributions from SW (orange) and LW (green).",
    "F2": "Figure 2. The TOA energy flux potentials for (a) TOT, (c) "
    "SW, and (e) LW net fluxes, alongside maps of the spatial anomalies "
    "of the input fluxes [(b),(d),(f)]. The implied heat transport is "
    "the gradient of the energy flus potential, shown by the white "
    "vector arrows (with the same magnitude scale across all subplots). "
    "Heat is directed from the blue minima of the potential field to "
    "yellow maxima, with the magnitude implied by the density of "
    "contours. All maps of the same type share the same color bar at "
    "the bottom of the column so that it is possible to directly "
    "compare the results from different fluxes.",
    "F3": "Figure 3. Direct radiative effects of clouds on the meridional "
    "heat transport. (a) Contributions from TOT CRE (blue), SW CRE "
    "(orange), and LW CRE (green). (b) Contributions from all-sky and "
    "clear-sky OSR. Both curves have been multiplied by -1 such that "
    "positive heat transport is northward.",
    "F4": "Figure 4. As in Figure 2, but for cloud radiative effects.",
    "F5": "Figure 5. As in Figure 2, but for energy flux potentials and "
    "spatial radiative anomalies associated with all-sky and clear-sky "
    "outgoing shortwave radiation. ",
    "F6": "Figure 6. A measure of the symmetry between heat transport in the "
    "Northern and Southern Hemispheres, calculated for the 12-month "
    "running mean of MHT in (a) the full hemisphere, (b) from the "
    "equator to 30 deg latitude, and (c) between 30 and 90 deg "
    "latitude. Symmetry values obtained when including (blue) and "
    "excluding (orange) the effect of clouds are shown. The "
    "climatological symmetry values for the two cases are shown as "
    "black lines in each subplot. The standard deviations of the "
    "time series are shown in each subplot.",
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
        "ancestors": filenames,
        "caption": figure_caption,
        "references": ["pearce23jclim"],
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


def area_average(cube, latitude="latitude", longitude="longitude", mdtol=1):
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
    cube_avg = cube.collapsed(
        [longitude, latitude],
        iris.analysis.MEAN,
        weights=grid_areas,
        mdtol=mdtol,
    )
    return cube_avg


def weight_zm(cube, latitude=None):
    """Weight zonal-mean by normalised gridbox areas.

    Parameters
    ----------
    cube : :class:`iris.cube.Cube`
        Input cube.
    latitude : tuple
        Four-element tuple defining the latitude range.
        The last two elements must be False, e.g.
        latitude=(-90, 0, False, False).

    Returns
    -------
    :class:`numpy.array`
        Zonal-mean in the selected latitude range, weighted
        by the normalised areas.
    """
    if cube.coord("latitude").bounds is None:
        cube.coord("latitude").guess_bounds()
    areas_data = iris.analysis.cartography.area_weights(cube, normalize=True)
    cube_areas = iris.cube.Cube(
        areas_data,
        long_name="normalised_area",
        var_name="area",
        units="1",
        dim_coords_and_dims=[(cube.coords()[0], 0)],
    )
    if latitude is not None:
        cube = cube.intersection(latitude=latitude)
        cube_areas = cube_areas.intersection(latitude=latitude)
    return cube.data * cube_areas.data


def call_poisson(flux_cube, latitude="latitude", longitude="longitude"):
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
    grid_areas = iris.analysis.cartography.area_weights(flux_cube)
    data_mean = flux_cube.collapsed(
        ["longitude", "latitude"], iris.analysis.MEAN, weights=grid_areas
    ).data
    data = flux_cube.data - data_mean

    logger.info("Calling spherical_poisson")
    sphpo = SphericalPoisson(
        logger, source=data * (earth_radius**2.0), tolerance=2.0e-4
    )
    sphpo.solve()
    sphpo.calc_meridional_heat_transport()
    logger.info("Ending spherical_poisson")

    # Energy flux potential
    efp_cube = iris.cube.Cube(
        sphpo.energy_flux_potential[1:-1, 1:-1],
        long_name=f"energy_flux_potential_of_{flux_cube.var_name}",
        var_name=f"{flux_cube.var_name}_efp",
        units="J s-1",
        dim_coords_and_dims=[
            (flux_cube.coords()[0], 0),
            (flux_cube.coords()[1], 1),
        ],
    )

    # MHT data cube
    collapsed_longitude = iris.coords.AuxCoord(
        180.0,
        bounds=(0.0, 360.0),
        long_name="longitude",
        standard_name="longitude",
        units="degrees",
    )
    dim_coords_and_dims = [(flux_cube.coord("latitude"), 0)]
    aux_coords_and_dims = [
        (flux_cube.coord("time"), None),
        (collapsed_longitude, None),
    ]
    mht_cube = iris.cube.Cube(
        sphpo.meridional_heat_transport,
        long_name=f"meridional_heat_transport_of_{flux_cube.var_name}",
        var_name=f"{flux_cube.var_name}_mht",
        units="W",
        dim_coords_and_dims=dim_coords_and_dims,
        aux_coords_and_dims=aux_coords_and_dims,
    )
    return efp_cube, mht_cube


def symmetry_metric(cube):
    """Calculate symmetry metrics for a zonal-mean cube.

    It returns the symmetry metric S, as defined in Pearce and
    Bodas-Salcedo, JClim, 2023, for 3 regions: entire hemisphere,
    tropics (0 to 30 deg latitude) and extratropics
    (30 to 90 degrees latitude). Perfectly symmetrical latitude
    bands give S=0.

    Parameters
    ----------
    cube : :class:`iris.cube.Cube`
        Input cube.

    Returns
    -------
    hemisphere: float
        Metric for the whole hemisphere.
    tropics: float
        Metric for the tropics.
    extra_tropics: float
        Metric for the extra-tropics.
    """
    hemisphere = np.abs(
        weight_zm(cube, latitude=(0, 90, False, False))[::-1]
        + weight_zm(cube, latitude=(-90, 0, False, False))
    ).sum()
    tropics = np.abs(
        weight_zm(cube, latitude=(0, 30, False, False))[::-1]
        + weight_zm(cube, latitude=(-30, 0, False, False))
    ).sum()
    extra_tropics = np.abs(
        weight_zm(cube, latitude=(30, 90, False, False))[::-1]
        + weight_zm(cube, latitude=(-90, -30, False, False))
    ).sum()
    return hemisphere, tropics, extra_tropics


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
    axes.set_xticklabels(["180", "120W", "60W", "0", "60E", "120E", "180"])
    axes.set_yticks(np.arange(-90, 100, 30))
    axes.set_yticklabels(["90S", "60S", "30S", "Eq", "30N", "60N", "90N"])
    axes.annotate(label, xy=(0, 1.05), xycoords="axes fraction", color="k")
    axes.set_title(title)


class ImpliedHeatTransport:
    """Class that solves implied heat transport for an input dataset.

    These are the physical meanings of the main acronyms
    used in the variable names:
       FLX: radiative flux
       EFP: energy flux potential
       MHT: meridional heat transport
    """

    def __init__(self, flx_files):
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
        self.mht_rolling_mean = iris.cube.CubeList()
        self.symmetry_metric = iris.cube.CubeList()

        # Calculate 12-month rolling means for time series.
        self.flx_rolling_mean = iris.cube.CubeList()
        for flx_file in flx_files:
            flx = iris.load_cube(flx_file)
            if len(flx.shape) == 3:
                self.flx_rolling_mean.append(
                    flx.rolling_window("time", iris.analysis.MEAN, 12)
                )
            else:
                self.flx_clim.append(flx)

        # Compute derived fluxes
        self.derived_fluxes()

        # Calculate Energy Flux Potential and Meridional Heat Transport
        # for each flux component
        self.compute_efp_and_mht()

        # Times series of MHT symmetry metric
        self.mht_symmetry_metrics()

        self.print()

    def compute_efp_and_mht(self):
        """Calculate Energy Flux Potential and meridional heat transport.

        Loop over input data and calculate EFP and MHT of the
        climatologies of radiative fluxes and the 12-month
        rolling means of radiative fluxes.
        """
        # Loop over climatologies
        for flx in self.flx_clim:
            efp, mht = call_poisson(flx)
            self.efp_clim.append(efp)
            self.mht_clim.append(mht)
        # Loop over rolling means
        for flx_rm in self.flx_rolling_mean:
            mht_series = iris.cube.CubeList()
            for flx in flx_rm.slices_over("time"):
                efp, mht = call_poisson(flx)
                mht_series.append(mht)
            # Append MHT rolling mean after merging time series.
            self.mht_rolling_mean.append(mht_series.merge_cube())

    def derived_fluxes(self):
        """Calculate derived radiative fluxes.

        rlnt_clim: climatology of net LW TOA
        rtntcs_clim: climatology of clear-sky net TOA
        rtntcs_rolling_mean: 12-month rolling mean of rtntcs
        """
        # Derived TOA climatologies: rlnt_clim, rtntcs_clim
        rlnt_clim = -1.0 * self.flx_clim.extract_cube(
            NameConstraint(var_name="rlut")
        )
        rlnt_clim.var_name = "rlnt"
        rlnt_clim.long_name = "radiative_flux_of_rlnt"
        self.flx_clim.append(rlnt_clim)
        rtntcs_clim = (
            self.flx_clim.extract_cube(NameConstraint(var_name="rsdt"))
            - self.flx_clim.extract_cube(NameConstraint(var_name="rsutcs"))
            - self.flx_clim.extract_cube(NameConstraint(var_name="rlutcs"))
        )
        rtntcs_clim.var_name = "rtntcs"
        rtntcs_clim.long_name = "radiative_flux_of_rtntcs"
        self.flx_clim.append(rtntcs_clim)
        # Annual rolling means clear-sky net total TOA
        rtntcs_rolling_mean = (
            self.flx_rolling_mean.extract_cube(NameConstraint(var_name="rsdt"))
            - self.flx_rolling_mean.extract_cube(
                NameConstraint(var_name="rsutcs")
            )
            - self.flx_rolling_mean.extract_cube(
                NameConstraint(var_name="rlutcs")
            )
        )
        rtntcs_rolling_mean.var_name = "rtntcs"
        rtntcs_rolling_mean.long_name = "radiative_flux_of_rtntcs"
        self.flx_rolling_mean.append(rtntcs_rolling_mean)

    def print(self):
        """Print variable names of all cubes in an IHT object."""
        logger.info("=== implied_heat_transport object ===")
        logger.info(self.mht_clim)
        info_message = "Long name: %s; Variable: %s."
        for climatology in self.mht_clim:
            logger.info(
                info_message, climatology.long_name, climatology.var_name
            )

        logger.info(self.efp_clim)
        for climatology in self.efp_clim:
            logger.info(
                info_message, climatology.long_name, climatology.var_name
            )

        logger.info(self.flx_clim)
        for climatology in self.flx_clim:
            logger.info(
                info_message, climatology.long_name, climatology.var_name
            )

        logger.info(self.mht_rolling_mean)
        for rolling_mean in self.mht_rolling_mean:
            logger.info(
                info_message, rolling_mean.long_name, rolling_mean.var_name
            )

        logger.info(self.symmetry_metric)
        for metric in self.symmetry_metric:
            logger.info(info_message, metric.long_name, metric.var_name)

        logger.info(self.flx_files)

    def mht_symmetry_metrics(self):
        """Calculate symmetry metrics.

        Produce 12-month rolling means for all monthly time series
        of MHT.
        """
        petaunit = 1.0e15
        for mht_series in self.mht_rolling_mean:
            time_coord = mht_series.coord("time")
            ntime = time_coord.shape[0]
            hemisphere = np.zeros(ntime)
            tropics = np.zeros(ntime)
            extra_tropics = np.zeros(ntime)
            for i in np.arange(ntime):
                hemisphere[i], tropics[i], extra_tropics[i] = symmetry_metric(
                    mht_series[i]
                )
            # Create the cubes for each metric
            long_name = f"symmetry_hemisphere_of_{mht_series.long_name}"
            var_name = f"s_hem_{mht_series.var_name}"
            cube_h = iris.cube.Cube(
                hemisphere / petaunit,
                long_name=long_name,
                var_name=var_name,
                units="PW",
                dim_coords_and_dims=[(time_coord, 0)],
            )
            long_name = f"symmetry_tropics_of_{mht_series.long_name}"
            var_name = f"s_tro_{mht_series.var_name}"
            cube_t = iris.cube.Cube(
                tropics / petaunit,
                long_name=long_name,
                var_name=var_name,
                units="PW",
                dim_coords_and_dims=[(time_coord, 0)],
            )
            long_name = f"symmetry_extratropics_of_{mht_series.long_name}"
            var_name = f"s_ext_{mht_series.var_name}"
            cube_e = iris.cube.Cube(
                extra_tropics / petaunit,
                long_name=long_name,
                var_name=var_name,
                units="PW",
                dim_coords_and_dims=[(time_coord, 0)],
            )
            self.symmetry_metric.append(cube_h)
            self.symmetry_metric.append(cube_t)
            self.symmetry_metric.append(cube_e)

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
            mht.convert_units("PW")
            plt.plot(
                mht.coord("latitude").points, mht.data, label=legend_label[i]
            )
        plt.hlines(0, -90, 90, color="k", linestyles=":")
        plt.vlines(0, -10, 10, color="k", linestyles=":")
        plt.xlim(-90, 90)
        plt.ylim(ylim[0], ylim[1])
        plt.xticks(np.arange(-90, 120, 30))
        plt.xlabel("Latitude")
        plt.ylabel("MHT (PW)")
        plt.legend()
        plt.tight_layout()

    def cre_mht_plot(self, left, right, ylim=(-1.5, 1.5)):
        """Produce two multiline plots of MHT components.

        Parameters
        ----------
        left : dictionary
            Dictionary with variable names and labels for
            the LHS plot, e.g.
            {'vname': ['netcre_mht', 'swcre_mht', 'lwcre_mht'],
            'legend': ['Net CRE', 'SW CRE', 'LW CRE']}
        right : dictionary
            As ``left`` but for the RHS plot
        ylim : tuple
            y axis limits.
        """
        plt.figure(figsize=(11, 5))
        ax1 = plt.subplot(121)
        for i, vname in enumerate(left["vname"]):
            mht = self.mht_clim.extract_cube(NameConstraint(var_name=vname))
            mht.convert_units("PW")
            ax1.plot(
                mht.coord("latitude").points, mht.data, label=left["legend"][i]
            )
        ax1.axhline(0, color="k", ls=":")
        ax1.axvline(0, color="k", ls=":")
        ax1.set_xlim(-90, 90)
        ax1.set_xticks(np.arange(-90, 120, 30))
        ax1.set_xlabel("Latitude")
        ax1.set_ylim(ylim[0], ylim[1])
        ax1.set_ylabel("MHT (PW)")
        ax1.annotate(
            "(a)", xy=(0.01, 0.95), xycoords="axes fraction", color="k"
        )
        plt.legend()

        ax2 = plt.subplot(122)
        col = ["C3", "C7"]
        for i, vname in enumerate(right["vname"]):
            mht = self.mht_clim.extract_cube(NameConstraint(var_name=vname))
            mht.convert_units("PW")
            ax2.plot(
                mht.coord("latitude").points,
                -mht.data,
                label=right["legend"][i],
                color=col[i],
            )
        ax2.axhline(0, color="k", ls=":")
        ax2.axvline(0, color="k", ls=":")
        ax2.set_xlim(-90, 90)
        ax2.set_xticks(np.arange(-90, 120, 30))
        ax2.set_xlabel("Latitude")
        ax2.set_ylim(ylim[0], ylim[1])
        ax2.set_ylabel("MHT (PW)")
        ax2.annotate(
            "(b)", xy=(0.01, 0.95), xycoords="axes fraction", color="k"
        )
        plt.legend(loc="lower right")
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
        return {"efp": efp, "flx": flx, "uuu": u_component, "vvv": v_component}

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
        mshgrd = np.meshgrid(
            self.flx_clim[0].coord("longitude").points,
            self.flx_clim[0].coord("latitude").points,
        )
        nrows = len(dargs["var_name"])
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
            data = self.quiver_maps_data(
                dargs["var_name"][i], dargs["change_sign"][i]
            )
            plt.subplot(
                grds[i * grid_step : (i * grid_step) + grid_step, 0],
                projection=ccrs.PlateCarree(central_longitude=0),
            )
            cbs.append(
                iplt.contourf(
                    data["efp"],
                    levels=np.linspace(
                        dargs["vmin"], dargs["vmax"], dargs["nlevs"]
                    ),
                )
            )
            plt.gca().coastlines()
            if i == 0:
                qqq = plt.quiver(
                    mshgrd[0][starty :: dxy[1], startx :: dxy[0]],
                    mshgrd[1][starty :: dxy[1], startx :: dxy[0]],
                    data["uuu"][starty :: dxy[1], startx :: dxy[0]],
                    data["vvv"][starty :: dxy[1], startx :: dxy[0]],
                    pivot="mid",
                    color="w",
                    width=0.005,
                )
            else:
                plt.quiver(
                    mshgrd[0][starty :: dxy[1], startx :: dxy[0]],
                    mshgrd[1][starty :: dxy[1], startx :: dxy[0]],
                    data["uuu"][starty :: dxy[1], startx :: dxy[0]],
                    data["vvv"][starty :: dxy[1], startx :: dxy[0]],
                    pivot="mid",
                    scale=qqq.scale,
                    color="w",
                )
            format_plot(plt.gca(), dargs["label"][i][0], dargs["title"][i][0])

            plt.subplot(
                grds[i * grid_step : (i * grid_step) + grid_step, 1],
                projection=ccrs.PlateCarree(central_longitude=0),
            )
            cbs.append(
                iplt.contourf(
                    data["flx"],
                    levels=np.linspace(
                        dargs["wmin"], dargs["wmax"], dargs["nwlevs"]
                    ),
                    cmap="RdBu_r",
                )
            )
            plt.gca().coastlines()
            format_plot(plt.gca(), dargs["label"][i][1], dargs["title"][i][1])

        plt.subplot(grds[-1, 0])
        plt.colorbar(
            cbs[0],
            cax=plt.gca(),
            orientation="horizontal",
            label="Energy flux potential (PW)",
        )
        plt.subplot(grds[-1, 1])
        plt.colorbar(
            cbs[1],
            cax=plt.gca(),
            orientation="horizontal",
            label=r"Flux (Wm$^{-2}$)",
            ticks=np.linspace(dargs["wmin"], dargs["wmax"], dargs["nwlevs"])[
                1 :: dargs["wlevstep"]
            ],
        )

        if nrows == 3:
            plt.subplots_adjust(left=0.1, right=0.94, top=1.0, bottom=0.11)
        elif nrows == 2:
            plt.subplots_adjust(left=0.11, right=0.9, top=1.0, bottom=0.13)
        elif nrows == 1:
            plt.subplots_adjust(left=0.11, right=0.9, top=1.0, bottom=0.20)

    def plot_symmetry_time_series(self):
        """Produce Figure 6.

        All-sky and clear-sky time series of the symmetry metrics for
        three regions: globe, tropics and extra-tropics.
        """
        var_list = [
            ["s_hem_rtnt_mht", "s_hem_rtntcs_mht"],
            ["s_tro_rtnt_mht", "s_tro_rtntcs_mht"],
            ["s_ext_rtnt_mht", "s_ext_rtntcs_mht"],
        ]
        col = ["C0", "C1"]
        label = [
            r"Global: 0$^\mathrm{o}$ - 90$^\mathrm{o}$",
            r"Tropics: 0$^\mathrm{o}$ - 30$^\mathrm{o}$",
            r"Extratropics: 30$^\mathrm{o}$ - 90$^\mathrm{o}$",
        ]
        legend_label = ["TOA net all-sky", "TOA net clear-sky"]

        plt.figure(figsize=(6, 12))
        for count, (var_name_1, var_name_2) in enumerate(var_list):
            yy0 = self.symmetry_metric.extract_cube(
                NameConstraint(var_name=var_name_1)
            )
            yy1 = self.symmetry_metric.extract_cube(
                NameConstraint(var_name=var_name_2)
            )
            axx = plt.subplot(3, 1, count + 1)
            dtx = [
                datetime.datetime.strptime(str(cell[0]), "%Y-%m-%d %H:%M:%S")
                for cell in yy0.coord("time").cells()
            ]
            plt.plot(dtx, yy0.data, lw=4, linestyle="-", label=legend_label[0])
            plt.plot(dtx, yy1.data, lw=4, linestyle="-", label=legend_label[1])
            axx.annotate(
                rf"$\sigma$: {np.std(yy0.data):5.3f}",
                (0.05, 0.55),
                xycoords="axes fraction",
                color=col[0],
            )
            axx.annotate(
                rf"$\sigma$: {np.std(yy1.data):5.3f}",
                (0.05, 0.45),
                xycoords="axes fraction",
                color=col[1],
            )
            axx.set_ylim(0, 0.8)
            axx.set_ylabel(r"$S$ (PW)")
            axx.xaxis.set_major_locator(mdates.YearLocator(3, month=1, day=1))
            axx.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
            axx.xaxis.set_minor_locator(mdates.YearLocator())
            axx.set_title(label[count])
            if count == 0:
                plt.legend(loc=5)
        plt.tight_layout()


def efp_maps(iht, model, experiment, config):
    """Produce Figures 2, 4, and 5.

    Parameters
    ----------
    iht : :class: ImpliedHeatTransport
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
            "var_name": [
                ["rtnt_efp", "rtnt"],
                ["rsnt_efp", "rsnt"],
                ["rlnt_efp", "rlnt"],
            ],
            "title": [
                ["$P_{TOA}^{TOT}$", r"$\Delta F_{TOA}^{TOT}$"],
                ["$P_{TOA}^{SW}$", r"$\Delta F_{TOA}^{SW}$"],
                ["$P_{TOA}^{LW}$", r"$\Delta F_{TOA}^{LW}$"],
            ],
            "label": [["(a)", "(b)"], ["(c)", "(d)"], ["(e)", "(f)"]],
            "change_sign": [[False, False], [False, False], [False, False]],
            "wmin": -180,
            "wmax": 180,
            "nwlevs": 19,
            "wlevstep": 4,
            "vmin": -1.2,
            "vmax": 1.2,
            "nlevs": 11,
        }
    )
    flx_files = matching_strings(iht.flx_files, ["rtnt/", "rsut/", "rlut/"])
    provenance_record = get_provenance_record(flx_files, caption["F2"])
    figname = f"figure2_{model}_{experiment}"
    save_figure(figname, provenance_record, config)
    # Figure 4
    iht.quiver_subplot(
        {
            "var_name": [
                ["netcre_efp", "netcre"],
                ["swcre_efp", "swcre"],
                ["lwcre_efp", "lwcre"],
            ],
            "title": [
                ["$P_{TOA}^{TOTCRE}$", r"$\Delta CRE_{TOA}^{TOT}$"],
                ["$P_{TOA}^{SWCRE}$", r"$\Delta CRE_{TOA}^{SW}$"],
                ["$P_{TOA}^{LWCRE}$", r"$\Delta CRE_{TOA}^{LW}$"],
            ],
            "label": [["(a)", "(b)"], ["(c)", "(d)"], ["(e)", "(f)"]],
            "change_sign": [[False, False], [False, False], [False, False]],
            "wmin": -60,
            "wmax": 60,
            "nwlevs": 13,
            "wlevstep": 2,
            "vmin": -0.3,
            "vmax": 0.3,
            "nlevs": 11,
        }
    )
    flx_files = matching_strings(
        iht.flx_files, ["netcre/", "swcre/", "lwcre/"]
    )
    provenance_record = get_provenance_record(flx_files, caption["F4"])
    figname = f"figure4_{model}_{experiment}"
    save_figure(figname, provenance_record, config)
    # Figure 5
    iht.quiver_subplot(
        {
            "var_name": [["rsutcs_efp", "rsutcs"], ["rsut_efp", "rsut"]],
            "title": [
                ["$P_{TOA}^{SWup, clr}$", r"$\Delta F_{TOA}^{SWup, clr}$"],
                ["$P_{TOA}^{SWup, all}$", r"$\Delta F_{TOA}^{SWup, all}$"],
            ],
            "label": [["(a)", "(b)"], ["(c)", "(d)"]],
            "change_sign": [[True, True], [True, True]],
            "wmin": -100,
            "wmax": 100,
            "nwlevs": 21,
            "wlevstep": 3,
            "vmin": -0.35,
            "vmax": 0.35,
            "nlevs": 11,
        }
    )
    flx_files = matching_strings(iht.flx_files, ["rsut/", "rsutcs/"])
    provenance_record = get_provenance_record(flx_files, caption["F5"])
    figname = f"figure5_{model}_{experiment}"
    save_figure(figname, provenance_record, config)


def mht_plots(iht, model, experiment, config):
    """Produce Figures 1 and 3.

    Parameters
    ----------
    iht : :class: ImpliedHeatTransport
        Object with the recipe datasets.
    model : string
        Model name.
    experiment : string
        Experiment name.
    config : dict
        The ESMValTool configuration.
    """
    # Figure 1
    iht.mht_plot(["rtnt_mht", "rsnt_mht", "rlnt_mht"], ["Net", "SW", "LW"])
    flx_files = matching_strings(iht.flx_files, ["rtnt/", "rsut/", "rlut/"])
    provenance_record = get_provenance_record(flx_files, caption["F1"])
    figname = f"figure1_{model}_{experiment}"
    save_figure(figname, provenance_record, config)
    # Figure 3
    iht.cre_mht_plot(
        {
            "vname": ["netcre_mht", "swcre_mht", "lwcre_mht"],
            "legend": ["Net CRE", "SW CRE", "LW CRE"],
        },
        {
            "vname": ["rsut_mht", "rsutcs_mht"],
            "legend": ["-1 x OSR (all-sky)", "-1 x OSR (clear-sky)"],
        },
    )
    flx_files = matching_strings(
        iht.flx_files, ["netcre/", "swcre/", "lwcre/", "rsut/", "rsutcs/"]
    )
    provenance_record = get_provenance_record(flx_files, caption["F3"])
    figname = f"figure3_{model}_{experiment}"
    save_figure(figname, provenance_record, config)


def symmetry_plots(iht, model, experiment, config):
    """Produce Figure 6.

    Parameters
    ----------
    iht : :class: ImpliedHeatTransport
        Object with the recipe datasets.
    model : string
        Model name.
    experiment : string
        Experiment name.
    config : dict
        The ESMValTool configuration.
    """
    iht.plot_symmetry_time_series()
    flx_files = matching_strings(
        iht.flx_files,
        ["rtnt_monthly", "rsutcs_monthly", "rlutcs_monthly", "rsdt_monthly"],
    )
    provenance_record = get_provenance_record(flx_files, caption["F6"])
    figname = f"figure6_{model}_{experiment}"
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
            symmetry_plots(iht_experiment, model, experiment, config)


def main(config):
    """Produce all the recipe's plots.

    Produce Figures 1 to 6 of Pearce and Bodas-Salcedo (2023) for each
    model and dataset combination.

    Parameters
    ----------
    config : dict
        The ESMValTool configuration.
    """
    input_data = deepcopy(list(config["input_data"].values()))
    input_data = group_metadata(input_data, "dataset", sort="variable_group")

    # Arrange input flux files in a 2-level dictionary [model_name][dataset]
    flux_files = {}
    for model_name, datasets in input_data.items():
        flux_files[model_name] = {}
        for dataset in datasets:
            if dataset["dataset"] in flux_files[model_name]:
                flux_files[model_name][dataset["dataset"]].append(
                    dataset["filename"]
                )
            else:
                flux_files[model_name][dataset["dataset"]] = [
                    dataset["filename"]
                ]

    # Create dictionary of implied_heat_transport objects.
    # It's a 2-level dictionary like flux_files.
    # This is where all the calculations are done.
    iht = {}
    for model_name, datasets in flux_files.items():
        logger.info("Model %s", model_name)
        iht[model_name] = {}
        for dataset_name, files in datasets.items():
            logger.info("Dataset %s", dataset_name)
            iht[model_name][dataset_name] = ImpliedHeatTransport(files)

    # Produce plots
    plot_single_model_diagnostics(iht, config)


if __name__ == "__main__":
    with run_diagnostic() as configuration:
        main(configuration)
