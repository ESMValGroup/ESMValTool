"""
ESMValTool diagnostic for calculating 4GST.

This is doen using LAI from CMIP6 and satellite obseravtions.
"""

import logging

import iris
import iris.coord_categorisation as icc
import iris.plot as iplt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import dask.array as da
from distributed import Client
from distributed import LocalCluster
from iris.fileformats.netcdf.loader import CHUNK_CONTROL
from iris import COMBINE_POLICY
import cartopy.crs as ccrs

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
)


import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap, BoundaryNorm



logger = logging.getLogger(__name__)


def _get_input_cubes(metadata):
    """Load the data files into cubes.

    Based on the hydrology diagnostic.

    Inputs:
    metadata = List of dictionaries made from the preprocessor config

    Outputs:
    inputs = Dictionary of cubes
    ancestors = Dictionary of filename information
    """
    inputs = {}
    ancestors = {}
    for attributes in metadata:
        short_name = attributes["short_name"]
        filename = attributes["filename"]
        logger.info("Loading variable %s", short_name)
        cube = iris.load_cube(filename)
        cube.attributes.clear()
        inputs[short_name] = cube
        ancestors[short_name] = [filename]

    return inputs, ancestors



def _get_provenance_record(attributes, ancestor_files):
    """Create the provenance record dictionary.

    Inputs:
    attributes = dictionary of ensembles/models used, the region bounds
                 and years of data used.
    ancestor_files = list of data files used by the diagnostic.

    Outputs:
    record = dictionary of provenance records.
    """
    caption = (
        "Timeseries of ESA CCI LST difference to mean of "
        "model ensembles calculated over region bounded by latitude "
        "{lat_south} to {lat_north}, longitude {lon_west} to {lon_east} "
        "and for model/ensembles {ensembles}. "
        + "Shown for years {start_year} to {end_year}.".format(**attributes)
    )

    record = {
        "caption": caption,
        "statistics": ["mean", "stddev"],
        "domains": ["reg"],
        "plot_types": ["times"],
        "authors": ["king_robert"],
        # 'references': [],
        "ancestors": ancestor_files,
    }

    return record

# DASK stiff
def setup(n_workers=1, threads_per_worker=4, processes=False):
    """_summary_

    Args:
        n_workers (int, optional): _description_. Defaults to 1.
        threads_per_worker (int, optional): _description_. Defaults to 4.
        processes (bool, optional): _description_. Defaults to False.

    Returns:
        _type_: _description_
    """
    cluster = LocalCluster(n_workers = n_workers,
                           threads_per_worker = threads_per_worker,
                           processes = processes)
    client = cluster.get_client()

    return cluster, client

# definition of the basic "threshold" calculation = find threshold exceedance ignoring all after max and before prior min
def threshcalc(arr, alpha):
    """
    Calculate time-index-of-first-threshold-exceedance.

    For a data array (ny, nx, ..., nt)
    = a time sequence at each (y, x, ...) location
    Perform a *separate* time-sequence calculation at each location.

    Returns
        INTEGER array (ny, nx), of time-indexes

    For use in dask.array.map_blocks, we need to consider how it "knows" about the expected relation to the passed array.
    This requires the following:
        * the time dimension must be complete in each block -- i.e. data must NOT be chunked in the time dim (use rechunk if needed)
        * the calc doesn't support a trial call with zero-length data : must use "meta=" keyword
        * the first dim will be dropped : use "drop_dims=(0,)" keyword
        * the result always has dtype "i8" : use 'dtype' keyword

    """
    # orig_arr = arr[...]
    nt = arr.shape[-1]
    inds_shape = (1,) * (arr.ndim - 1) + (nt,)
    timeinds = np.arange(nt).reshape(inds_shape) * np.ones(arr.shape)  # time index expanded to full shape
    # find max + blank times after it, at each landpoint
    maxs = np.max(arr, axis=-1)
    maxinds = np.argmax(arr, axis=-1)
    # re-add a degenerate final dim
    # N.B. direct assignment here is problematic, because of reshape on boolean indexing
    #  - if costly, could use stack + index instead of 'where' ?
    wherefn = np.ma.where if np.ma.is_masked(arr) else np.where 
    arr = wherefn(timeinds > maxinds[..., None], maxs[..., None], arr)
    # find min + blank times before it, at each landpoint
    # NB must be done AFTER blanking out times>max-time !
    mins = np.min(arr, axis=-1)
    mininds = np.argmin(arr, axis=-1)
    arr = wherefn(timeinds < mininds[..., None], mins[..., None], arr)
    # calculate threshold values (at each landpoint)
    threshs = mins + alpha * (maxs - mins)
    # calculate "time-index of first exceedance of threshold"
    threshinds = np.argmax(arr > threshs[..., None], axis=-1)
    return threshinds

def plot_map(cube, model = 'CMCC', year=2000):

    # ============================================================
    # Create 12-month × 4-week colourblind-friendly colormap
    # ============================================================

    month_colours = [
        "#0072B2",  # Jan
        "#56B4E9",  # Feb
        "#009E73",  # Mar
        "#7CAE00",  # Apr
        "#D6B000",  # May (darker amber than F0E442)
        "#E69F00",  # Jun
        "#D55E00",  # Jul
        "#CC79A7",  # Aug
        "#882255",  # Sep
        "#332288",  # Oct
        "#44AA99",  # Nov
        "#999933",  # Dec
    ]

    week_fades = [0.55, 0.35, 0.15, 0.0]
    colours = [(1, 1, 1)]
    
    for month_colour in month_colours:
        rgb = np.array(mcolors.to_rgb(month_colour))

        for fade in week_fades:
            # Mix with white to produce week shades
            colours.append(rgb * (1 - fade) + fade)

    this_cmap = ListedColormap(colours)
    this_norm = BoundaryNorm(np.arange(50), this_cmap.N)

    # ============================================================
    # Convert DOY (1-365/366) -> month/week category (0-47)
    # ============================================================

    lookup = np.zeros(367, dtype=np.int16)

    for doy in range(1, 367):

        date = pd.Timestamp("2001-01-01") + pd.Timedelta(days=doy - 1)

        month = date.month                    # 1-12
        week = min((date.day - 1) // 7, 3)    # 0-3

        lookup[doy] = (month - 1) * 4 + week

    # Copy cube and replace DOY values with category values
    plot_cube = cube.copy()
    min_val = np.nanmin(cube.data)
    is_min = cube.data == min_val
    # Convert DOY -> category
    plot_cube.data = lookup[cube.data.astype(int)]
    # Set minima to category 0 (white)
    plot_cube.data[is_min] = 0



    # ============================================================
    # Plot
    # ============================================================

    fig = plt.figure(figsize=(16,9))
    
    pcm = iplt.pcolormesh(
        plot_cube,
        cmap=this_cmap,
        norm=this_norm
    )

    # ============================================================
    # Colorbar
    # ============================================================

    cbar = plt.colorbar(pcm,
                        ticks=[0.5] + list(np.arange(2.5, 49, 4)),
                        orientation="horizontal",
                        )

    cbar.ax.set_xticklabels(["N/A",
        "Jan", "Feb", "Mar", "Apr",
        "May", "Jun", "Jul", "Aug",
        "Sep", "Oct", "Nov", "Dec"
    ])

    cbar.set_label("Month (shade = week within month)")

    plt.gca().coastlines()

    plt.title(f"Vegetation Onset: {model} {year}", fontsize=24)

    plt.savefig(f'onset_{model}_{year}.png')
    plt.close()
    
    return None


def _diagnostic(config):
    """Perform the control for the ESA CCI LST diagnostic.

    Parameters
    ----------
    config: dict
        the preprocessor nested dictionary holding
        all the needed information.

    Returns
    -------
    figures made by make_plots.
    """
    # this loading function is based on the hydrology diagnostic
    input_metadata = config["input_data"].values()

    loaded_data = {}
    ancestor_list = []
    for dataset, metadata in group_metadata(input_metadata, "dataset").items():
        cubes, ancestors = _get_input_cubes(metadata)
        loaded_data[dataset] = cubes
        ancestor_list.append(ancestors["lai"][0])


    logger.info(f"{loaded_data}")
    # data is nested dictionaries MODEL LAI

    for MODEL in loaded_data.keys():
        if 'lai' in loaded_data[MODEL].keys():
            # follow the Dask, onset proceedure
            cluster, client = setup(n_workers=2, threads_per_worker=1, processes=True)

            lazarr = loaded_data[MODEL]['lai'].core_data()
            # this is needed for the OBS data where a lot of days are all NaNs
            # need to find a generic way to do this wit all OBS and MODELS....
            sam = lazarr[:,  0,0].compute()
            good_day_inds = np.where(~np.isnan(sam))
            logger.info(f'{good_day_inds=}')
            print(f'{good_day_inds=}')
            print('**********************************')
            # why this note on this line? this should work what ever the data NaN structure????
            good_days = lazarr[good_day_inds]  # NOTE: this is not correct, would only work if every month has 30 days
            data = good_days.transpose((1, 2, 0))
            
            # this needs a way to be generic
            data_r = data.rechunk({-1:-1, 1:20}) # 1186 was C3S LAI

            thresh_inds = da.map_blocks(
                threshcalc,
                data_r,
                alpha=0.2, # can this be passed in from the recipe???????
                dtype=int,
                drop_axis=[-1],
                meta=np.ma.array(0),
            )

            print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
            lat_coord = loaded_data[MODEL]['lai'].coord('latitude')
            lon_coord = loaded_data[MODEL]['lai'].coord('longitude')
            logger.info(f"{lat_coord=}")
            logger.info(f"{lon_coord=}")
            result_cube = iris.cube.Cube(thresh_inds,
                                         dim_coords_and_dims = ((lat_coord,0),
                                                                (lon_coord,1)),
                                         long_name = "Vegeation Onset Index"
                                         )
            try:
                icc.add_day_of_year( loaded_data[MODEL]['lai'], 'time')
            except:
                pass

           
            print('/#/#/#/#/##/#/#/##/#/#')

            doy_values = loaded_data[MODEL]['lai'].coord('day_of_year').points
            doy_data = doy_values[thresh_inds]

            doy_cube = iris.cube.Cube(doy_data,
                                      dim_coords_and_dims = ((lat_coord,0),
                                                             (lon_coord,1)),
                                      long_name = "Vegeation Onset "
                                         )
            
            print(f'{doy_cube=}')

            # change to esmvaltool save path for run

            # iris.save(result_cube, f'/home/users/robking/CMUG/ESMValTool/esmvaltool/cube_index_{MODEL}.nc')
            iris.save(doy_cube, f'/home/users/robking/CMUG/ESMValTool/esmvaltool/cube_doy_{MODEL}.nc')

           
        else:
            print('CONTINUE CONTINUE CONTINUE ---------------------------------')
            continue
 
        plot_map(doy_cube, model=MODEL, year = 2000)
    # record = _get_provenance_record(data_attributes, ancestor_list)
    # plot_file = get_plot_filename("timeseries", config)
    # with ProvenanceLogger(config) as provenance_logger:
    #     provenance_logger.log(plot_file, record)


if __name__ == "__main__":
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        _diagnostic(config)
