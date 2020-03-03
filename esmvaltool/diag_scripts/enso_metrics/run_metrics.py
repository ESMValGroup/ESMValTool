"""ENSO Metrics package wrapper"""

import os
import logging
from collections import defaultdict
import json

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata

from ensometrics.EnsoCollectionsLib import CmipVariables, defCollection
from ensometrics.EnsoComputeMetricsLib import compute_collection

logger = logging.getLogger(os.path.basename(__file__))


class ENSOMetrics(object):
    VARIABLE_ALIAS = {'sst': 'tos'}

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

        models = defaultdict(dict)
        obs = defaultdict(dict)

        for alias in data:
            variables = group_metadata(data[alias], n.SHORT_NAME)
            for var in required_vars:
                var_alias = ENSOMetrics.VARIABLE_ALIAS.get(var, var)
                if var_alias not in variables:
                    continue
                var_info = variables[var_alias][0]
                if var_info[n.DATASET] == var_info['reference_dataset']:
                    dictionary = obs[alias]
                else:
                    dictionary = models[alias]
                dictionary[var] = self.create_var_dict(var_info, variables)

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

    def create_var_dict(self, var_info, variables):
        var_dict = {
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
            filename = variables[cell_area][0][n.FILENAME]
            var_dict["path + filename_area"] = filename
            var_dict["areaname"] = cell_area
        if land_sea_mask in variables:
            filename = variables[land_sea_mask][0][n.FILENAME]
            var_dict["path + filename_landmask"] = filename
            var_dict["landmaskname"] = land_sea_mask
        return var_dict

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
                dict_met = dict_in[ens]["metadata"]["metrics"][met]
                diag_dict = dict_met["diagnostic"]
                for key1 in diag_dict:
                    if key1 not in [
                        "time_frequency",
                        "ref",
                        "method",
                        "method_nonlinearity",
                        "name",
                    ]:
                        if key1 == "units":
                            dict_meta[key1] = diag_dict[key1]
                        else:
                            dict_meta[key1] = diag_dict[key1]["nyears"]
                units = dict_met["metric"]["units"]
                metric_dict = dict_in[ens]["value"][met]["metric"]
                if metric_only is True:
                    # metrics
                    dict2 = dict()
                    for key1 in list(metric_dict.keys()):
                        tmp = metric_dict[key1]["value"]
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
                    for key1 in metric_dict:
                        tmp = metric_dict[key1]["value"]
                        tmp_key = key1.replace("ref_", "")
                        dict2["metric"][tmp_key] = {
                            "value": tmp,
                            "nyears_obs": dict_meta[tmp_key],
                            "units": units,
                        }
                        del tmp, tmp_key
                    # dive down diagnostics
                    diag_dict = dict_in[ens]["value"][met]["diagnostic"]
                    for key1 in diag_dict:
                        tmp = diag_dict[key1]["value"]
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


def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        ENSOMetrics(config).compute()


if __name__ == "__main__":
    main()
