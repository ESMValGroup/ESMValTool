# Calculate anomalies, and then plot our model groups in a boxplot

from esmvaltool.diag_scripts.shared import (
    run_diagnostic,
    group_metadata,
    select_metadata,
    extract_variables,
)

import iris
import iris.quickplot as qplt
import cartopy
import cartopy.crs as ccrs

import os
import logging

import matplotlib.pyplot as plt

logger = logging.getLogger(os.path.basename(__file__))


def process_projections_dict(proj_dict, season):
    # recursive function to pull out data from dictionary
    out_data = {}
    for k, v in proj_dict.items():
        if isinstance(v, dict):
            vals = process_projections_dict(v, season)
            for k1, v1 in vals.items():
                out_data[f"{k} {k1}"] = v1
        else:
            if v is None:
                continue
            # extract required season
            season_con = iris.Constraint(season_number=season)
            data = v.extract(season_con)
            # if the result is a scalar cube, just store the value
            # else store the whole cube
            if data.ndim == 1:
                out_data[k] = data.data.item()
            else:
                out_data[k] = data
    return out_data


def get_anomalies(ds_list, base_clim_start, fut_clim_start, relative=False):
    # construct baseline
    base_metadata = select_metadata(ds_list, start_year=base_clim_start)
    if base_metadata == []:
        logging.warning(f"Base climatology (start {base_clim_start}) not found")
        return None
    base_file = base_metadata[0]["filename"]
    base_cube = iris.load_cube(base_file)

    # get future
    fut_metadata = select_metadata(ds_list, start_year=fut_clim_start)
    if fut_metadata == []:
        logging.warning(f"Future climatology (start {fut_clim_start}) not found")
        return None
    fut_file = fut_metadata[0]["filename"]
    fut_cube = iris.load_cube(fut_file)

    if relative:
        diff = fut_cube - base_cube
        anomaly = (diff / base_cube) * 100
        anomaly.units = "%"
    else:
        anomaly = fut_cube - base_cube

    return anomaly


def main(cfg):
    # The config object is a dict of all the metadata from the pre-processor

    # get variable processed
    var = list(extract_variables(cfg).keys())
    assert len(var) == 1
    var = var[0]

    if var == "pr":
        rel_change = True
    else:
        rel_change = False

    # establish the time periods of our datasets
    start_years = list(group_metadata(cfg["input_data"].values(), "start_year"))
    base_start = min(start_years)
    fut_start = max(start_years)

    # first group datasets by project..
    # this creates a dict of datasets keyed by project (CMIP5, CMIP6 etc.)
    projects = group_metadata(cfg["input_data"].values(), "project")
    # how to uniquely define a dataset varies by project, for CMIP it's simple, just dataset...
    # for CORDEX, combo of dataset and driver (and possibly also domain if we start adding those)
    # also gets more complex if we start adding in different ensembles..

    # This section of the code loads and organises the data to be ready for plotting
    logger.info("Loading data")
    # empty dict to store results
    projections = {}
    model_lists = {}
    cordex_drivers = []
    # loop over projects
    for proj in projects:
        # we now have a list of all the data entries..
        # for CMIPs we can just group metadata again by dataset then work with that..
        models = group_metadata(projects[proj], "dataset")

        # empty dict for results
        projections[proj] = {}
        # loop over the models
        for m in models:
            if proj[:6].upper() == "CORDEX":
                # then we need to go one deeper in the dictionary to deal with driving models
                drivers = group_metadata(models[m], "driver")
                projections[proj][m] = dict.fromkeys(drivers.keys())
                for d in drivers:
                    logging.info(f"Calculating anomalies for {proj} {m} {d}")
                    anoms = get_anomalies(drivers[d], base_start, fut_start, rel_change)
                    if anoms is None:
                        continue
                    projections[proj][m][d] = anoms
                    if proj not in model_lists:
                        model_lists[proj] = []
                    model_lists[proj].append(f"{m} {d}")
                    cordex_drivers.append(d)
            elif proj == "UKCP18":
                # go deeper to deal with ensembles and datasets
                # split UKCP into seperate GCM and RCM
                proj_key = f"UKCP18 {m}"
                ensembles = group_metadata(models[m], "ensemble")
                projections[proj_key] = dict.fromkeys(ensembles.keys())
                for ens in ensembles:
                    logging.info(f"Calculating anomalies for {proj_key} {ens}")
                    anoms = get_anomalies(
                        ensembles[ens], base_start, fut_start, rel_change
                    )
                    if anoms is None:
                        continue
                    projections[proj_key][ens] = anoms
                    if proj_key not in model_lists:
                        model_lists[proj_key] = []
                    model_lists[proj_key].append(f"{proj_key} {ens}")
            else:
                logging.info(f"Calculating anomalies for {proj} {m}")
                anoms = get_anomalies(models[m], base_start, fut_start, rel_change)
                if anoms is None:
                    continue
                projections[proj][m] = anoms
                if proj not in model_lists:
                    model_lists[proj] = []
                model_lists[proj].append(f"{m}")
        # remove any empty categories (i.e. UKCP18 which has been split into rcm and gcm)
        if projections[proj] == {}:
            del projections[proj]
    cordex_drivers = set(cordex_drivers)

    # this section of the code does the plotting..
    # we now have all the projections in the projections dictionary

    # now lets plot them
    # first we need to process the dictionary, and move the data into a list of vectors
    # the projections object is the key one that contains all our data..
    seasons = {0: "DJF", 1: "MAM", 2: "JJA", 3: "OND"}
    logger.info("Plotting")
    extent = (
        cfg["domain"]["start_longitude"] - 2,
        cfg["domain"]["end_longitude"] + 2,
        cfg["domain"]["start_latitude"] - 2,
        cfg["domain"]["end_latitude"] + 2,
    )
    for s in seasons.keys():
        # make directory
        try:
            os.mkdir(f"{cfg['plot_dir']}/{seasons[s]}")
        except FileExistsError:
            pass
        for p in projections:
            pdata = process_projections_dict(projections[p], s)

            for m in pdata:
                title = f"{p} {m} {seasons[s]} {var} change"
                plt.figure(figsize=(12.8, 9.6))
                ax = plt.axes(projection=ccrs.PlateCarree())
                ax.set_extent(extent)
                # set scales
                if var == "pr":
                    vmn = -50
                    vmx = 50
                    cmap = "brewer_RdYlBu_11"
                else:
                    vmn = 0.5
                    vmx = 8
                    cmap = "brewer_YlOrRd_09"
                # ensure longitude coordinates straddle the meridian for GCM origin data
                if pdata[m].coord("longitude").ndim == 1:
                    # TODO This will probably cause issues if it's ever run with data
                    # that straddles the dateline, so a check should be added.
                    plot_cube = pdata[m].intersection(longitude=(-180.0, 180.0))
                    plot_cube.coord("longitude").circular = False
                else:
                    plot_cube = pdata[m]
                qplt.pcolormesh(plot_cube, vmin=vmn, vmax=vmx, cmap=cmap)
                plt.title(title)
                ax.coastlines()
                ax.add_feature(cartopy.feature.BORDERS, linestyle=":")
                plt.savefig(
                    f"{cfg['plot_dir']}/{seasons[s]}/{p}_{m}_map_{seasons[s]}.png"
                )
                plt.close()

    # print all datasets used
    print("Input models for plots:")
    for p in model_lists.keys():
        print(f"{p}: {len(model_lists[p])} models")
        print(model_lists[p])
        print("")


if __name__ == "__main__":
    with run_diagnostic() as cfg:
        main(cfg)
