"""Data finder and Data selector functions."""
import copy
import os
import subprocess
from pprint import pprint
from itertools import combinations
from template_parser import extract_requirements
from esmvalcore._data_finder import get_start_end_year
class VariableCap():
    "class that handles finding data sets for one variable"
    def __init__(self,var_name,diag_name,facets):
        self.var_name=var_name
        self.diag_name=diag_name
        self.special_flag= self.detect_special_flag(facets)
        if self.check_facets(facets):
            self.facets=facets
        else:
            raise Exception("Facets for diag %s and var %s are not valid!" %(diag_name,var_name))

    def check_facets(self,facets):
        """Check that facets contain mandatory keys"""
        mandatory_keys=["project","short_name","mip","exp"]
        return all([key in facets for key in mandatory_keys])

    def detect_special_flag(self,facets):
        """Check whether dictionnary facets contains special cases : mixed_scenario,piControl,facets with no_time"""
        if facets.setdefault("special_flag",None)== "mixed_scenario":
            return "mixed_scenario"
        elif facets.setdefault("start_year",None)== "duration":
            return "piControl"
        elif "start_year" not in facets or "end_year" not in facets:
            return "no_time"
        else:
            return "standard"


    def find_datasets(self):
        """Return list of datasets available at the local replica given self.facets"""
        project= self.facets["project"]
        special_flag = self.special_flag
        kept_facets = ['dataset', 'grid', 'ensemble', 'institute']
        kept_facets_piControl = [
        'dataset', 'grid', 'ensemble', 'institute', 'start_year', 'end_year']
        requirement = self.facets
        final_var_datasets = []
        #Handles mixed_scenario special_flag
        if special_flag == "mixed_scenario":
            assert (
                len(requirement["exp"]) >= 2
                and len(requirement["start_year"]) >= 2
            ), "Special flag mixed_scenario is used but only one exp is used in the template !"
            assert (len(requirement["exp"])==len(requirement["start_year"]) and len(requirement["exp"])==len(requirement["end_year"])), "Number of experiments, start_years and end_years do not match when they should"
            new_req={}
            #Transform mixed_scenario VariableCap to a DiagCap with standard special_flag
            for i in range(len(requirement["exp"])):
                key_i = "var"+str(i)
                new_req[key_i] = copy.deepcopy(requirement)
                new_req[key_i]["exp"] = requirement["exp"][i]
                new_req[key_i]["start_year"] = requirement["start_year"][i]
                new_req[key_i]["end_year"] = requirement["end_year"][i]
                new_req[key_i].pop("special_flag",None)
            temp_diag_name=self.diag_name+"_"+self.var_name+"_"+"mixed_scenario"
            diag_req = DiagnosticCap(temp_diag_name,new_req,False)
            mixed_scenario_datasets = diag_req.find_datasets()
            #Add start_year and end_year in the found datasets
            final_var_datasets=[]
            for ds in mixed_scenario_datasets:
                for i in range(0,len(requirement["exp"])):
                    copy_ds=copy.deepcopy(ds)
                    copy_ds["start_year"]= requirement["start_year"][i]
                    copy_ds["end_year"]= requirement["end_year"][i]
                    copy_ds["exp"]= requirement["exp"][i]
                    final_var_datasets.append(copy_ds)
                    del copy_ds
        #Handles other cases : piControl, no_time, standard
        else:
            #get all datasets file available from freva/tphoon search, plus check years etc...
            if requirement["project"]=="CMIP6":
                found_datasets = [
                    typhoon2datasets(item)
                    for item in get_info_from_typhoon(**requirement)
                ]
            elif requirement["project"]=="CMIP5":
                found_datasets = [
                        freva2datasets(item)
                        for item in get_info_from_freva(**requirement)
                        ]
            else: raise Exception("Project %s for diag %s and var %s is not handled!"%(project,self.diag_name,self.var_name))
            #Handle case when unvalid (because their paths are not following template) datasets are found
            if found_datasets==[None]*len(found_datasets):
                return []
            #Take care of each dataset one by one, removing all similar files to avoid duplicates
            while len(found_datasets) > 0 :
                current_ds = found_datasets[0]
                #get all dict (i.e nc files) that have same facets than current_ds
                list_similar_ds = get_similar_datasets(
                    found_datasets,
                    current_ds,
                    skipped_keys=["start_year", "end_year"])
                #remove those files from the freva_datasets
                for similar_ds in list_similar_ds:
                    found_datasets.remove(similar_ds)
                #Handles no_time special_flag
                if special_flag=="no_time":
                    facets_to_remove = set(current_ds.keys()).difference(kept_facets)
                    for d in facets_to_remove:
                        del current_ds[d]
                    final_var_datasets.append(current_ds)
                #Handles piControl special_flag
                elif special_flag=="piControl":
                    dataset_found, start_year, end_year = years_search_for_piControl(
                        list_similar_ds, requirement["end_year"])
                    if dataset_found:
                        facets_to_remove = set(current_ds.keys()).difference(
                            kept_facets_piControl)
                        for d in facets_to_remove:
                            del current_ds[d]
                        current_ds["start_year"] = start_year
                        current_ds["end_year"] = end_year
                        if current_ds not in final_var_datasets:final_var_datasets.append(current_ds)
                #Handles standard special_flag
                elif years_compatible(list_similar_ds,
                                      start_year=requirement["start_year"],
                                      end_year=requirement["end_year"]):
                    facets_to_remove = set(
                        current_ds.keys()).difference(kept_facets)
                    for d in facets_to_remove:
                        del current_ds[d]
                    final_var_datasets.append(current_ds)
        return final_var_datasets

    def modify_recipe_var(self,recipe_dict):
        """Modify template for current diag,variable given special_flag so it is a valid recipe"""
        special_flag=self.special_flag
        diag= self.diag_name
        var= self.var_name
        if special_flag=="piControl":
            recipe_dict["diagnostics"][diag]["variables"][var].pop("start_year",None)
            recipe_dict["diagnostics"][diag]["variables"][var].pop("end_year",None)
        elif special_flag=="mixed_scenario":
            recipe_dict["diagnostics"][diag]["variables"][var].pop("start_year",None)
            recipe_dict["diagnostics"][diag]["variables"][var].pop("end_year",None)
            recipe_dict["diagnostics"][diag]["variables"][var].pop("exp",None)
            recipe_dict["diagnostics"][diag]["variables"][var].pop("special_flag",None)
        else:
            pass
        return 0

class DiagnosticCap():
    "class that handles finding data sets for one diagnostic"
    def __init__(self,diag_name,requirements,single_ensemble):
        self.diag_name=diag_name
        self.requirements=requirements
        self.single_ensemble= single_ensemble
        var_cap_list=[]
        for key,facets in requirements.items():
            v= VariableCap(key,diag_name,facets)
            var_cap_list.append(v)
        self.VariableCap_list= var_cap_list

    def find_datasets(self):
        """Get available datasets for current diagnostic after getting them for each variable"""
        single_ensemble= self.single_ensemble
        datasets_per_variable=[]
        for variable in self.VariableCap_list:
            datasets_per_variable.append(variable.find_datasets())
        #keeping only datasets available for all variables of the diagnostic
        final_diag_datasets = intersection_of_datasets_for_all_var(
        datasets_per_variable)
        #only keep on ensemble members if asked
        only_one_ensemble_datasets= []
        while len(final_diag_datasets)>0:
            current_ds=final_diag_datasets[0]
            if single_ensemble:
                list_ensemble_members = get_similar_datasets(
                    final_diag_datasets,
                    current_ds,
                    skipped_keys=["ensemble","grid"])
            else:
                list_ensemble_members = get_similar_datasets(
                    final_diag_datasets,
                    current_ds,
                    skipped_keys=["grid"])
            only_one_ensemble_datasets.append(current_ds)
            for discarded_ds in list_ensemble_members:
                final_diag_datasets.remove(discarded_ds)

        return only_one_ensemble_datasets

    def modify_recipe_diag(self,recipe_dict,diag_dataset):
        #loop over variable and modify some part of the current diagnostic
        for variable in self.VariableCap_list:
            variable.modify_recipe_var(recipe_dict)
        #add found datasets
        recipe_dict["diagnostics"][self.diag_name].setdefault("additional_datasets",[]).extend(diag_dataset)
        return 0

def intersection_of_datasets_for_all_var(diag_dataset: list) -> list:
    """
    diag_dataset is a list of list of dictionaries.
    = [[list of dict( i.e available datasets) for var1 in diag1],[list of dict for var2 in diag1],...[list of dict for var n ]]
    we output a list of dict for dicts that are in the intersection of each list of dict
    in our case : it is a list of datasets which are common to every var for diag1
    """
    intersected_datasets = []
    temp_datasets_list = diag_dataset
    if len(diag_dataset) == 1:
        return diag_dataset[0]
    if len(diag_dataset) == 0:
        return []
    #get the intersection for a pair of list of dict until there is only two lists
    while len(temp_datasets_list) > 2:
        #get intersection between first 2 lists(dict that are in both lists), remove those first 2 lists from temp_list and add the intersection found in temp_list
        temp_datasets_list = temp_datasets_list[2:] + [[
            ds for ds in temp_datasets_list[0] if ds in temp_datasets_list[1]
        ]]
    intersected_datasets = [
        ds for ds in temp_datasets_list[0] if ds in temp_datasets_list[1]
    ]
    return intersected_datasets


def get_similar_datasets(freva_ds: list, reference_dict: dict,
                         skipped_keys: list) -> list:
    """
    return a list of all dictionnaries inside freva_ds which have same values as reference dict for all keys except for the skipped keys specified 
    """
    output_list = list()
    for ds in freva_ds:
        try:
            if all([
                    ds[key] == reference_dict[key] for key in reference_dict
                    if key not in skipped_keys
            ]):
                output_list.append(ds)
        except:
            "Key error between reference dictionary and looked-up dict"
    return output_list


def years_search_for_piControl(list_ds, duration):
    year_list = sorted([ds["start_year"] for ds in list_ds],
                       key=lambda x: x[0])
    new_start_year = int(year_list[0])
    new_end_year = int(new_start_year) + int(duration)
    return years_compatible(list_ds, new_start_year,
                            new_end_year), new_start_year, new_end_year


def years_compatible(list_ds: list, start_year: int, end_year: int) -> bool:
    """
    return true if the current datasets have compatible years with start_year and end_year : i.e. do we have files matching with the requested start-year and end-year
    """
    years_to_cover= set(range(start_year,end_year+1))
    list_years_available= [list(range(int(list_ds[i]["start_year"]),int(list_ds[i]["end_year"])+1)) for i in range(len(list_ds))]
    flatten_years_available=[]
    for sublist in list_years_available:
        for item in sublist:
            flatten_years_available.append(item)
    years_available=set(flatten_years_available)
    return years_to_cover.issubset(years_available)

def typhoon2datasets(pathname):
    """Reformat the typhoon output to facet dictionary."""
    this_should_come_from_config = '[institute]/[dataset]/[exp]/[ensemble]/[mip]/[short_name]/[grid]/[latestversion]'.replace(
        '[', '').replace(']', '')
    facets = reversed(this_should_come_from_config.split('/'))
    facet_values = reversed(os.path.dirname(pathname).split(os.sep))
    out = dict(zip(facets, facet_values))
    times = get_start_end_year(os.path.basename(pathname))
    try :
        out['start_year'] = times[0]
        out['end_year'] = times[1]
    except IndexError:
        pass
    return None if pathname is None else out


def freva2datasets(pathname):
    """Reformat the freva output of CMIP5 only to facet dictionary."""
    this_should_come_from_config = '[root]/[institute]/[dataset]/[exp]/[frequency]/[modeling_realm]/[mip]/[ensemble]/[latestversion]/[short_name]'.replace(
        '[', '').replace(']', '')
    facets = reversed(this_should_come_from_config.split('/'))
    facet_values = reversed(os.path.dirname(pathname).split(os.sep))
    out = dict(zip(facets, facet_values))
    times = _get_start_end_year(os.path.basename(pathname))
    out['start_year'] = times[0]
    out['end_year'] = times[1]
    return None if pathname is None else out


def get_info_from_freva(**kwargs):
    """Get paths of available data via `freva`."""
    def _format_output(stdout):
        out = [item for item in stdout.split('\n') if len(item) != 0]
        return out

    translate = {
        'exp': 'experiment',
        'mip': 'cmor_table',
        'short_name': 'variable'
    }

    include = ['short_name', 'mip', 'exp','project']
    freva_cmd = "freva"
    module = "module load cmip6-dicad/1.0"
    facets = []
    if "project" in kwargs.keys() : kwargs["project"]= kwargs["project"].lower()
    for key, value in kwargs.items():
        if key in include:
            facets.append("{0}={1}".format(translate.get(key, key), value))
    #print(facets)
    cmd = " ".join([
        "{0} &> /dev/null; {1} ".format(module, freva_cmd), "--databrowser"
    ] + facets)
    #print("Sending freva query : %s" % cmd)
    return _format_output(subprocess.check_output(cmd, shell=True).decode())

def get_info_from_typhoon(**kwargs):
    """Get paths of available data via typhoon`."""
    typhoon_path="/work/bd0854/CMIP6_index/typhoon.txt"

    def _format_output(stdout):
        out = [item for item in stdout.split('\n') if len(item) != 0]
        return out

    cmd="grep "+ typhoon_path+" -e /"+list(kwargs.values())[0]+"/"
    facets = []
    include= ['short_name', 'mip', 'exp']
    for key, value in list(kwargs.items())[1:]:
        if key in include:
            cmd+="| grep -e /"+ value+"/"
    return _format_output(subprocess.check_output(cmd, shell=True).decode())



if __name__ == '__main__':
    PATHNAME = os.path.join(os.path.dirname(__file__), 'recipes',
                            'template_test.yml')

    LIST_OF_FACETS = extract_requirements(PATHNAME)
    for diag in LIST_OF_FACETS:
        pprint(get_available_datasets(LIST_OF_FACETS[diag]))
