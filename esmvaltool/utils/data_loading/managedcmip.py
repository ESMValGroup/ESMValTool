import os.path
import subprocess

import yaml

test_recipe = '/home/h02/ktomkins/projects/esmvaltool/ESMValTool/' \
              'esmvaltool/recipes/recipe_radiation_budget.yml'
test_recipe2 = '/home/h02/ktomkins/projects/esmvaltool/ESMValTool' \
               '/esmvaltool/recipes/' \
               'recipe_autoassess_landsurface_permafrost.yml'

managedcmip_dir = '/project/champ/data'


def read_in_recipe_yml_file(recipe):
    with open(recipe) as file:
        return yaml.safe_load(file)


def get_dataset_dicts_from_recipe(loaded_recipe):
    return loaded_recipe.get('datasets')


def get_path_items_from_recipe(datasets, path_item):
    items = list()
    for dataset in datasets:
        item_name = dataset.get(path_item)
        items.append(item_name)
    return list(dict.fromkeys(items))


def find_search_root_path(projects=None):
    root_paths = list()
    for project in projects:
        path = os.path.join(managedcmip_dir, project)
        root_paths.append(path)
    return root_paths


def find_search_paths(datasets=None, exps=None, ensembles=None):
    search_paths = list()
    for dataset in datasets:
        for exp in exps:
            for ensemble in ensembles:
                search_path = os.path.join(dataset, exp, ensemble)
                search_paths.append(search_path)
    return list(dict.fromkeys(search_paths))


def search_managed_cmip_for_data(root_paths=None, search_paths=None):
    for root_path in root_paths:
        for search_path in search_paths:
            subprocess.Popen(root_path)
            print([
                'find', root_path, '-type d', '-name *' + search_path, '-print'
            ])
            subprocess.run([
                'find', {root_path}, '-type d', '-name *' + search_path,
                '-print'
            ])


def main(recipe):
    loaded_recipe = read_in_recipe_yml_file(recipe)
    datasets = get_dataset_dicts_from_recipe(loaded_recipe)
    dataset_names = get_path_items_from_recipe(datasets, 'dataset')
    project_names = get_path_items_from_recipe(datasets, 'project')
    exp_names = get_path_items_from_recipe(datasets, 'exp')
    ensemble_names = get_path_items_from_recipe(datasets, 'ensemble')
    find_root_path = find_search_root_path(projects=project_names)
    find_search_path = find_search_paths(datasets=dataset_names,
                                         exps=exp_names,
                                         ensembles=ensemble_names)
    search_managed_cmip_for_data(root_paths=find_root_path,
                                 search_paths=find_search_path)


if __name__ == '__main__':
    main(test_recipe)
    main(test_recipe2)
