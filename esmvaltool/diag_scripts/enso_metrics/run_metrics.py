"""ENSO Metrics package wrapper"""

import os
import logging
from collections import defaultdict

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata


logger = logging.getLogger(os.path.basename(__file__))

from cdms2 import open as CDMS2open
from copy import deepcopy
from getpass import getuser as GETPASSgetuser
from glob import iglob as GLOBiglob
from inspect import stack as INSPECTstack
import json
from os import environ as OSenviron
from os.path import join as OSpath__join
import sys


# ENSO_metrics package
# set new path where to find programs
# sys.path.insert(0, "/home/" + user_name + "/Test001/ENSO_metrics/lib")
# sys.path.insert(0, "/home/" + user_name + "/Test001/ENSO_metrics/plots")
# sys.path.insert(1, "/home/" + user_name + "/Test001/ENSO_metrics/scripts")
from ensometrics.EnsoCollectionsLib import (
    CmipVariables,
    defCollection,
    ReferenceObservations,
)
from ensometrics.EnsoComputeMetricsLib import compute_collection



class ENSOMetrics(object):
    def __init__(self, config):
        self.cfg = config
        self.filenames = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.metrics_collection = defCollection(self.cfg['metrics_collection'])
        self.metrics_collection_name = self.cfg['metrics_collection']

    def compute(self):

        # list of variables
        required_vars = set()
        for name, info in self.metrics_collection['metrics_list'].items():
            required_vars.update((var for var in info['variables']))
        logger.info(required_vars)

        dict_var = CmipVariables()["variable_name_in_file"]
        logger.info(dict_var)

        data = group_metadata(self.cfg['input_data'].values(), n.ALIAS)
        dict_ens = {}
        VARIABLE_ALIAS = {'sst': 'tos'}

        models = defaultdict(dict)
        obs = defaultdict(dict)

        for alias in data:
            variables = group_metadata(data[alias], n.SHORT_NAME)
            for var in required_vars:
                var_alias = VARIABLE_ALIAS.get(var, var)
                if var_alias not in variables:
                    continue
                var_info = variables[var_alias][0]
                if var_info[n.DATASET] == var_info['reference_dataset']:
                    dictionary = obs[alias]
                else:
                    dictionary = models[alias]
                dictionary[var] = {
                    'path + filename': var_info[n.FILENAME],
                    'varname': var_info[n.SHORT_NAME],
                }
                if 'ocean' in var_info[n.MODELING_REALM]:
                    cell_area = 'areacello'
                    land_sea_mask = 'sftof'
                else:
                    cell_area = 'areacella'
                    land_sea_mask = 'sftlf'

                if cell_area in variables:
                    dictionary[var]["path + filename_area"] = variables[cell_area][0][n.FILENAME]
                    dictionary[var]["areaname"] = cell_area
                if land_sea_mask in variables:
                    dictionary[var]["path + filename_landmask"] = variables[land_sea_mask][0][n.FILENAME]
                    dictionary[var]["landmaskname"] = land_sea_mask

        for alias in models:

            datasets = {"model": models, "observations": obs}
            results = compute_collection(
                self.cfg['metrics_collection'],
                datasets,
                alias,
                netcdf=self.cfg[n.WRITE_NETCDF],
                debug=self.cfg[n.LOG_LEVEL] == 'debug',
                netcdf_name=os.path.join(
                    self.cfg[n.WORK_DIR],
                    alias,
                )
            )
            dict_ens[alias] = results[0]
        # save json
        json_name = os.path.join(
            self.cfg[n.WORK_DIR],
            self.metrics_collection_name,
        )
        self._save_json(
            dict_ens,
            json_name=json_name,
            metric_only=False,
        )
        with open(f"{json_name}_raw.json", "w") as outfile:
            json.dump(dict_ens, outfile, sort_keys=True)


    def _save_json(self, dict_in, json_name, metric_only=True):
        # reshape dictionary
        liste = sorted(dict_in.keys())
        listm = sorted(dict_in[liste[0]]["value"].keys())
        dict_out = dict()
        for met in listm:
            dict1 = dict()
            for ens in liste:
                # metadata (nyears)
                dict_meta = dict()
                for key1 in list(
                    dict_in[ens]["metadata"]["metrics"][met]["diagnostic"].keys()
                ):
                    if key1 not in [
                        "time_frequency",
                        "ref",
                        "method",
                        "method_nonlinearity",
                        "name",
                    ]:
                        if key1 == "units":
                            dict_meta[key1] = dict_in[ens]["metadata"]["metrics"][
                                met
                            ]["diagnostic"][key1]
                        else:
                            dict_meta[key1] = dict_in[ens]["metadata"]["metrics"][
                                met
                            ]["diagnostic"][key1]["nyears"]
                units = dict_in[ens]["metadata"]["metrics"][met]["metric"]["units"]
                if metric_only is True:
                    # metrics
                    dict2 = dict()
                    for key1 in list(dict_in[ens]["value"][met]["metric"].keys()):
                        tmp = dict_in[ens]["value"][met]["metric"][key1]["value"]
                        tmp_key = key1.replace("ref_", "")
                        dict2[tmp_key] = {
                            "metric": tmp,
                            "nyears_obs": dict_meta[tmp_key],
                            "units": units,
                        }
                        del tmp, tmp_key
                else:
                    # metrics
                    dict2 = {"metric": {}, "diagnostic": {}}
                    for key1 in list(dict_in[ens]["value"][met]["metric"].keys()):
                        tmp = dict_in[ens]["value"][met]["metric"][key1]["value"]
                        tmp_key = key1.replace("ref_", "")
                        dict2["metric"][tmp_key] = {
                            "value": tmp,
                            "nyears_obs": dict_meta[tmp_key],
                            "units": units,
                        }
                        del tmp, tmp_key
                    # dive down diagnostics
                    for key1 in list(
                        dict_in[ens]["value"][met]["diagnostic"].keys()
                    ):
                        tmp = dict_in[ens]["value"][met]["diagnostic"][key1][
                            "value"
                        ]
                        if key1 == "model":
                            dict2["diagnostic"][ens] = {
                                "value": tmp,
                                "nyears": dict_meta[key1],
                                "units": dict_meta["units"],
                            }
                        else:
                            dict2["diagnostic"][key1] = {
                                "value": tmp,
                                "nyears": dict_meta[key1],
                                "units": dict_meta["units"],
                            }
                        del tmp
                dict1[ens] = dict2
                del dict_meta, dict2
            dict_out[met] = dict1
            del dict1
        # save as json file
        with open(json_name + ".json", "w") as outfile:
            json.dump(dict_out, outfile, sort_keys=True)
        return

# # list of variables

# # list of observations
# list_obs = list()
# for metric in list_metric:
#     dict_var_obs = dict_mc["metrics_list"][metric]["obs_name"]
#     for var in list(dict_var_obs.keys()):
#         for obs in dict_var_obs[var]:
#             if obs not in list_obs:
#                 list_obs.append(obs)
# list_obs = sorted(list_obs)
# if mc_name == "MC1":
#     list_obs = ["Tropflux"]
# elif mc_name == "ENSO_perf":
#     list_obs = [
#         "ERA-Interim"
#     ]  # ['ERA-Interim', 'HadISST', 'Tropflux', 'GPCPv2.3']#['Tropflux','GPCPv2.3']#['HadISST']#
# elif mc_name == "ENSO_tel":
#     list_obs = ["ERA-Interim"]  # ['ERA-Interim', 'HadISST', 'GPCPv2.3']
# elif mc_name == "ENSO_proc":
#     list_obs = ["ERA-Interim"]  # ['ERA-Interim', 'HadISST', 'GPCPv2.3']
# elif mc_name == "ENSO_test":
#     list_obs = ["AVISO", "ERA-Interim", "HadISST", "Tropflux", "GPCPv2.3"]
# print("\033[95m" + str(list_obs) + "\033[0m")

def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        ENSOMetrics(config).compute()


if __name__ == "__main__":
    main()