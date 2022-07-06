# Calculate anomalies, and then plot our model groups in a boxplot

from esmvaltool.diag_scripts.shared import (
    run_diagnostic,
    group_metadata,
    select_metadata,
    extract_variables,
)

import iris

import os
import logging
import re

import numpy as np
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
            # this should be a scalar cube, add the value to a dictionary
            out_data[k] = data.data.item()
    return out_data


def get_anomalies(ds_list, base_clim_start, fut_clim_start, relative=False):
    # construct baseline
    base_metadata = select_metadata(ds_list, start_year=base_clim_start)
    if base_metadata == []:
        logging.warning(
            f"Base climatology (start {base_clim_start}) not found")
        return None
    base_file = base_metadata[0]["filename"]
    base_cube = iris.load_cube(base_file)

    # get future
    fut_metadata = select_metadata(ds_list, start_year=fut_clim_start)
    if fut_metadata == []:
        logging.warning(
            f"Future climatology (start {fut_clim_start}) not found")
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


def save_anoms_txt(data, fname):
    # iterate over the supplied dictionary and write the data to a textfile
    # sort the data
    sorted_data = sorted(data.items(), key=lambda x: x[1])

    # open the file for writing
    with open(fname, mode="w") as f:
        for d in sorted_data:
            # write a line of data
            f.write(f"{d[0]}: {d[1]:.2f}\n")


def get_var(cfg):
    # get variable processed
    var = list(extract_variables(cfg).keys())
    assert len(var) == 1
    var = var[0]

    return var


def plot_boxplots(projections, legend_models):
    # we now have all the projections in the projections dictionary
    # check what driving models we have from CORDEX and decide on some fixed colours for them..
    # get default colours
    prop_cycle = plt.rcParams["axes.prop_cycle"]
    colours = prop_cycle.by_key()["color"]
    colours_it = enumerate(colours)
    legend_colours = {}
    for k in legend_models:
        legend_colours[k] = next(colours_it)[1]

    # special models that have an extra large symbol
    # special_models = ["RACMO22E", "HadREM3-GA7-05"]
    special_models = []

    # get variable for title later
    var = get_var(cfg)

    # plot legend type
    legend_type = cfg['cordex_legend_type']

    # now lets plot them
    # first we need to process the dictionary, and move the data into a list of vectors
    # the projections object is the key one that contains all our data..
    seasons = {0: "DJF", 1: "MAM", 2: "JJA", 3: "OND"}
    p_keys = reorder_keys(list(projections.keys()))

    logger.info("Plotting")
    for s in seasons.keys():
        proj_plotting = dict.fromkeys(p_keys)
        for p in p_keys:
            proj_plotting[p] = process_projections_dict(projections[p], s)
            save_anoms_txt(proj_plotting[p],
                           f"{cfg['work_dir']}/{p}_{seasons[s]}_values.txt")

        # eventually plotting code etc. will go in a seperate module I think.
        projects, plot_values = zip(*proj_plotting.items())
        # plots_values is a list of dictionaries of model names and associated values
        plt.figure(figsize=(12.8, 9.6))
        plt.boxplot([list(v.values()) for v in plot_values])
        plotted_models = set()
        for i, p in enumerate(proj_plotting.keys()):
            for m, v in proj_plotting[p].items():
                if p[:6].upper() == "CORDEX":
                    # extract the driving model and RCM from the string
                    rcm, driver = m.split(" ")
                    if legend_type == 'driver':
                        color = legend_colours[driver]
                        label = driver
                    else:
                        color = legend_colours[rcm]
                        label = rcm
                    alpha = 1
                    if any(i in m for i in special_models):
                        sz = 100
                    else:
                        sz = 25
                    # Check if we have already plotted this model before..
                    # This means we only use the label once, and just end up with one legend entry per model
                    if any(i in m for i in plotted_models):
                        label = None
                    else:
                        plotted_models.add(label)
                elif "CMIP" in p:
                    model = [x for x in legend_models if m in x]
                    # check if model matches any of the CORDEX ones
                    if model:
                        model = model[0]
                        color = legend_colours[model]
                        alpha = 1
                        sz = 25
                        # Check if we have already plotted this model before..
                        # This means we only use the label once
                        if any(i in model for i in plotted_models):
                            label = None
                        else:
                            plotted_models.add(model)
                    else:
                        color = "k"
                        alpha = 0.3
                        label = None
                        sz = 10
                else:
                    color = "k"
                    alpha = 0.3
                    label = None
                    sz = 10
                # Add some random "jitter" to the x-axis
                x = np.random.normal(i + 1, 0.05, size=1)
                plt.scatter(x, v, label=label, color=color, alpha=alpha, s=sz)
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
        plt.axhline(0, linestyle="dotted", color="k")
        plt.gca().set_xticklabels(projects)
        plt.title(f"{seasons[s]} {var} change")
        plt.tight_layout()
        plt.savefig(f"{cfg['plot_dir']}/boxplot_{seasons[s]}.png")
        plt.close()


def simple_dots_plot(projections, cordex_drivers):
    # get variable for title later
    var = get_var(cfg)

    # now lets plot them
    # first we need to process the dictionary, and move the data into a list of vectors
    # the projections object is the key one that contains all our data..
    seasons = {0: "DJF", 1: "MAM", 2: "JJA", 3: "OND"}
    logger.info("Plotting")
    p_keys = reorder_keys(list(projections.keys()))
    for s in seasons.keys():
        proj_plotting = dict.fromkeys(p_keys)
        for p in p_keys:
            proj_plotting[p] = process_projections_dict(projections[p], s)
            save_anoms_txt(proj_plotting[p],
                           f"{cfg['work_dir']}/{p}_{seasons[s]}_values.txt")

        # plots_values is a list of dictionaries of model names and associated values
        plt.figure(figsize=(12.8, 9.6))
        projects, values = zip(*proj_plotting.items())
        for i, p in enumerate(p_keys):
            for m, v in proj_plotting[p].items():
                if p == "CMIP5":
                    if any(m in d for d in cordex_drivers):
                        plt.plot(i + 1,
                                 v,
                                 marker="o",
                                 fillstyle="full",
                                 color="k",
                                 markersize=12)
                if 'CMIP' not in p:
                    plt.plot(
                        i + 1,
                        v,
                        marker="o",
                        fillstyle="none",
                        color="k",
                    )
            if 'CMIP' in p:
                plt.violinplot(list(values[i].values()),
                               positions=[i + 1],
                               showmedians=True)
        plt.gca().set_xticks(range(1, len(projects) + 1))
        plt.gca().set_xticklabels(projects)
        plt.title(f"{seasons[s]} {var} change")
        plt.tight_layout()
        plt.savefig(f"{cfg['plot_dir']}/violin_{seasons[s]}.png")
        plt.close()


def reorder_keys(keys):
    # order keys so that first 3 projects are CMIP6, CMIP5, CORDEX
    if 'CORDEX' in keys:
        # remove from present location and move to front
        i = keys.index('CORDEX')
        keys.insert(0, keys.pop(i))
    if 'CMIP5' in keys:
        # remove from present location and move to front
        i = keys.index('CMIP5')
        keys.insert(0, keys.pop(i))
    if 'CMIP6' in keys:
        # remove from present location and move to front
        i = keys.index('CMIP6')
        keys.insert(0, keys.pop(i))

    return keys


def main(cfg):
    # The config object is a dict of all the metadata from the pre-processor

    # get variable processed
    var = get_var(cfg)

    if var == "pr":
        rel_change = True
    else:
        rel_change = False

    # establish the time periods of our datasets
    start_years = list(group_metadata(cfg["input_data"].values(),
                                      "start_year"))
    base_start = min(start_years)
    fut_start = max(start_years)

    # first group datasets by project..
    # this creates a dict of datasets keyed by project (CMIP5, CMIP6 etc.)
    projects = group_metadata(cfg["input_data"].values(), "project")
    # how to uniquely define a dataset varies by project, for CMIP it's simple, just dataset...
    # for CORDEX, combo of dataset and driver (and possibly also domain if we start adding those)
    # also gets more complex if we start adding in different ensembles..

    # get "special" RCMS if being used
    if 'special_rcms' in cfg:
        spec_rcms = cfg['special_rcms']
    else:
        spec_rcms = None

    # This section of the code loads and organises the data to be ready for plotting
    logger.info("Loading data")
    # empty dict to store results
    projections = {}
    model_lists = {}
    cordex_drivers = []
    cordex_rcms = []
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
                    anoms = get_anomalies(drivers[d], base_start, fut_start,
                                          rel_change)
                    if anoms is None:
                        continue
                    projections[proj][m][d] = anoms
                    if proj not in model_lists:
                        model_lists[proj] = []
                    model_lists[proj].append(f"{m} {d}")
                    cordex_drivers.append(d)
                    cordex_rcms.append(m)
            elif proj == "UKCP18":
                # go deeper to deal with ensembles and datasets
                # split UKCP into seperate GCM and RCM
                proj_key = f"UKCP18 {m}"
                ensembles = group_metadata(models[m], "ensemble")
                projections[proj_key] = dict.fromkeys(ensembles.keys())
                for ens in ensembles:
                    logging.info(f"Calculating anomalies for {proj_key} {ens}")
                    anoms = get_anomalies(ensembles[ens], base_start,
                                          fut_start, rel_change)
                    if anoms is None:
                        continue
                    projections[proj_key][ens] = anoms
                    if proj_key not in model_lists:
                        model_lists[proj_key] = []
                    model_lists[proj_key].append(f"{proj_key} {ens}")
            else:
                logging.info(f"Calculating anomalies for {proj} {m}")
                anoms = get_anomalies(models[m], base_start, fut_start,
                                      rel_change)
                if anoms is None:
                    continue
                projections[proj][m] = anoms
                if proj not in model_lists:
                    model_lists[proj] = []
                model_lists[proj].append(f"{m}")

        # seperate CORDEX RCMs into special and normal if needed
        if spec_rcms and proj[:6].upper() == "CORDEX":
            # create new dictionary entry
            projections['CORDEX_aerosol'] = {}
            for m in models:
                if m in spec_rcms:
                    data = projections[proj].pop(m)
                    projections['CORDEX_aerosol'][m] = data

        # remove any empty categories (i.e. UKCP18 which has been split into rcm and gcm)
        if projections[proj] == {}:
            del projections[proj]
    cordex_drivers = set(cordex_drivers)
    cordex_rcms = set(cordex_rcms)

    if cfg['cordex_legend_type'] == 'driver':
        legend_models = cordex_drivers
    else:
        legend_models = cordex_rcms

    # this section of the code does all the plotting..
    plot_boxplots(projections, legend_models)
    simple_dots_plot(projections, cordex_drivers)

    # print all datasets used
    print("Input models for plots:")
    for p in model_lists.keys():
        print(f"{p}: {len(model_lists[p])} models")
        print(model_lists[p])
        print("")


if __name__ == "__main__":
    with run_diagnostic() as cfg:
        main(cfg)
