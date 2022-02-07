"""Main function for the automatization tool of ESMValTool"""
import filecmp
import argparse
import glob
import logging
import os
from pprint import pprint, pformat
import copy
import shutil
try:
    import oyaml as yaml
except ImportError:
    import yaml
import subprocess
from template_parser import extract_requirements
from dataset_selector import DiagnosticCap
import sys
_ROOT_PATH = os.path.dirname(os.path.realpath(__file__))
_DEFAULT_CONFIG_FILE = os.path.join(_ROOT_PATH, 'config.yml')
_LOG_FORMATTER = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')


def get_args():
    """Get input arguments."""
    argparser = argparse.ArgumentParser(
        description='Automated creation of CMIP analysis results')
    argparser.add_argument('-c',
                           '--configfile',
                           type=str,
                           default=_DEFAULT_CONFIG_FILE,
                           help='Specify configuration file')
    argparser.add_argument(
        '-r',
        '--runmode',
        type=str,
        default='generate',
        help='Specify run mode: can be "generate" (default) to generate recipes only,\
            or can be "esmvaltool" to run esmvaltool from recipes in \
            recipes/generated_recipes that are not in old_recipes or can be "all" to run\
            generate and then esmvaltool mode')
    argparser.add_argument(
        '-f',
        '--one_file',
        type=str,
        default=None,
        help='Specify one file to run the automatization instead of the files in recipe_dir of the config_file: only possible when runmode is generate')
    argparser.add_argument(
        '-e',
        '--one_ensemble',
        action="store_false",
        help="Specify whether we want only one ensemble members per models or not in the generated recipe")


    return argparser.parse_args()


def normalize_path(path):
    """Normalize path."""
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))


def read_config_file(config_path):
    """Read configuration file."""
    with open(config_path, 'r') as config_file:
        cfg = yaml.safe_load(config_file)

    # Check if necessary keys are available
    necessary_keys = {
        'log_file': os.path.join('~', 'cmip_results.log'),
        'recipe_dir': os.path.join(_ROOT_PATH, 'recipes/template'),
        'created_recipe_dir': os.path.join(_ROOT_PATH, 'recipes/generated_recipes/'),
        'old_created_recipe_dir': os.path.join(_ROOT_PATH, 'recipes/old_generated_recipes/'),
        'work_dir': os.path.join('~', 'cmip_results/'),
        'job_dir': os.path.join(_ROOT_PATH,'job_output/'),
        'sbatch_dir': os.path.join(_ROOT_PATH,'temp/'),

    }
    for (key, default) in necessary_keys.items():
        if key not in cfg:
            print(f"""WARNING: No '{key}' set in configuration file,
                    defaulting to '{default}' """)
            cfg[key] = default
        cfg[key] = normalize_path(cfg[key])

    # Check if paths exists if not create and use default dict
    if not os.path.exists(cfg['recipe_dir']):
        print(f"WARNING: Recipe directory '{cfg['recipe_dir']}' does not "
                f"exist, using default one instead")
    for directory in [
            'work_dir', 'created_recipe_dir', 'old_created_recipe_dir',
            'job_dir', 'sbatch_dir','recipe_dir'
    ]:
        if not os.path.exists(cfg[directory]):
            print(f"WARNING: {directory} directory '{cfg[directory]}' does not "
                     f"exist, using default one {necessary_keys[directory]} instead")
            if not os.path.exists(necessary_keys[directory]) :
                print(
                f"INFO: Created non-existent path to default directory {directory} '{necessary_keys[directory]}'"
                )
                os.makedirs(necessary_keys[directory])
            cfg[directory]= necessary_keys[directory]
        else : cfg[directory]= cfg[directory]
    return cfg


def setup_logger(cfg):
    """Setup logger."""
    levels = {
        'debug': logging.DEBUG,
        'info': logging.INFO,
        'warning': logging.WARNING,
        'error': logging.ERROR,
    }
    default = 'info'
    log_level = cfg.get('log_level', default)
    if log_level not in levels:
        print(f"WARNING: Got invalid log level '{log_level}', defaulting to "
              "{default}")
    log_level = levels.get(log_level, levels[default])
    logger = logging.getLogger()
    logger.setLevel(log_level)
    file_log_handler = logging.FileHandler(cfg['log_file'], mode='a')
    file_log_handler.setFormatter(_LOG_FORMATTER)
    logger.addHandler(file_log_handler)
    console_log_handler = logging.StreamHandler()
    console_log_handler.setFormatter(_LOG_FORMATTER)
    logger.addHandler(console_log_handler)
    return logger


def list_of_recipes_to_rerun(new_recipes_dir: str,
                             old_recipes_dir: str) -> list:
    """Give list of recipes in new_recipes_dir which have a modification compared to the previous version in old_recipes_dir """
    output_list = []
    for new_recipe_path in glob.glob(os.path.join(new_recipes_dir, '*.yml')):
        recipe_name = os.path.basename(os.path.normpath(new_recipe_path))
        old_recipe_path = os.path.join(old_recipes_dir, recipe_name)
        #check if new_recipe is not in old_recipes_folder
        if not os.path.isfile(old_recipe_path):
            output_list.append(new_recipe_path)
        else:
            #compare new_recipe and old_recipe
            if not filecmp.cmp(new_recipe_path, old_recipe_path,
                               shallow=False):
                output_list.append(new_recipe_path)
    return output_list


def write_and_run_sbatch(recipes_to_run, sbatch_folder, job_folder,
                         config_esm):
    """ Create a job taskis and send tehm for each recipes in recipes_to_run"""
    args = get_args()
    cfg = read_config_file(args.configfile)
    logger = setup_logger(cfg)
    if os.path.isdir(sbatch_folder):
        shutil.rmtree(sbatch_folder, ignore_errors=True)
    os.mkdir(sbatch_folder)
    if os.path.isdir(job_folder):
        shutil.rmtree(job_folder, ignore_errors=True)
    os.mkdir(job_folder)
    for recipe in recipes_to_run:
        current_recipe_name = os.path.basename(os.path.normpath(recipe))
        try:
            file_string = f"""#!/bin/bash
#SBATCH --partition=compute
#SBATCH --ntasks=16
#SBATCH --time=08:00:00
#SBATCH --account=bd1083
#SBATCH --output={job_folder}/job_{current_recipe_name}_%j.out.log\n\
#SBATCH --error={job_folder}/job_{current_recipe_name}_%j.err.log\n\
#SBATCH --job-name={current_recipe_name[7:-4]}\n\
###############################################################################

echo "Executing the ESMValTool."
plotdir=$(esmvaltool --skip-nonexistent -c {config_esm} {recipe} | grep 'PLOTDIR    ='|head -1|cut -d'=' -f2| tr -d ' ')

outdir="/scratch/b/${USER}/automatization_output/GLOBALOUT"

echo "Postprocessing the image files."
echo "Creating: "
echo "Processing directory $plotdir"
python {_ROOT_PATH}/provenance_v2_to_v1.py -i $plotdir -o ${outdir} -r {current_recipe_name[:-4]}
"""
            sbatch_file = "{}sbatch_{}.sh".format(sbatch_folder,
                                                  current_recipe_name[:-4])
            with open(sbatch_file, "w") as text_file:
                logger.info("Writing sbatch script %s" % sbatch_file)
                text_file.write(file_string)
                text_file.close()
            CMD = "sbatch {}".format(sbatch_file)
            proc = subprocess.call([CMD], shell=True)
            logger.info("Running sbatch script %s" % sbatch_file)
        except Exception as exp:
            print(exp)
    return (0)


def main():
    """Main function."""
    args = get_args()
    one_file = args.one_file
    one_ensemble= args.one_ensemble
    runmode = args.runmode
    cfg = read_config_file(args.configfile)
    logger = setup_logger(cfg)
    recipe_dir = cfg['recipe_dir']
    new_recipe_dir = cfg['created_recipe_dir']
    logger.info("AUTOMATIZATION TOOL getting started !")

    #Generate new recipes from templates
    if runmode == "all" or runmode == "generate":
        recipes_list= glob.glob(os.path.join(recipe_dir, '*.yml')) if one_file is None else [one_file]

        #For each template in the recipe_dir, generated a recipe from available datasets
        for recipe_path in recipes_list:

            logger.info("Generating recipe from template : %s" % recipe_path)
            with open(recipe_path, 'r') as recipe_file:
                recipe = yaml.safe_load(recipe_file)
            logger.debug("Processing recipe '%s'", recipe_path)

            #Get facets (requirements) from the template
            req = extract_requirements(recipe_path)

            #Get the datasets for each diagnostic in the template
            for diag in req:
                logger.info("Finding datasets for diagnostic %s" % diag)
                current_diag = DiagnosticCap(diag,req[diag],one_ensemble)
                #Use data finder and data selector given facets in req
                diag_datasets = current_diag.find_datasets()
                if len(diag_datasets) == 0:
                    logger.warning(
                        "0 dataset found for diagnostic %s"
                        % diag)
                #Modify diagnostics in the template so it can become a valid recipe, and add found datasets
                current_diag.modify_recipe_diag(recipe,diag_datasets)

            #Write the generated recipe in the folder specified in config.yml
            current_recipe_name = os.path.basename(
                os.path.normpath(recipe_path))
            #replace template by recipe in current_recipe_name
            current_recipe_name = current_recipe_name.replace(
                "template", "recipe")
            generated_recipe_name = os.path.join(new_recipe_dir,
                                                 current_recipe_name)
            with open(generated_recipe_name, 'w') as recipe_out:
                yaml.safe_dump(recipe, recipe_out)
                logger.info("Writing generated recipe in : %s \n" %
                            generated_recipe_name)

    #Run esmvaltool for the new recipes if asked
    if runmode == "all" or runmode == "esmvaltool":
        old_recipe_dir = cfg['old_created_recipe_dir']
        recipes_to_run = list_of_recipes_to_rerun(new_recipe_dir,
                                                  old_recipe_dir)
        with open("recipes_to_run.yml", 'w') as yamlfile:
            yaml.safe_dump(recipes_to_run, yamlfile)
        logger.info(
            "Writing which recipes have changed and need to be run in : %s \n"
            % "recipes_to_run.yml")
        logger.info("Recipes to rerun are:%s" % ("\n".join(recipes_to_run)))
        job_folder = cfg["job_dir"]
        sbatch_folder = cfg["sbatch_dir"]
        config_esm = cfg["esm_config_file"]
        write_and_run_sbatch(recipes_to_run, sbatch_folder, job_folder,
                             config_esm)
    return 0


if __name__ == '__main__':
    main()
