# automatization
Functionalities for the automated creation of CMIP results.
Goal: automatically generate ESMValTool recipes given available datasets from a *template recipe* specifying the project, variable, years, the mip, etc ... needed for the recipe. Those generated recipes can optionally be automatically sent as a job.

### Installation guidelines :
- install *ESMValGroup/ESMValTool*

### How to run :
- Activate your "ESMValTool specific" virtual environment
- Go to *esmvaltool/utils/automatization* folder of ESMValTool

- Adapt the config.yml to your need : you can specify the config file used by esmvaltool (by default it will use *config-esm.yml* which you can adapt to your needs), the directory of the template recipes, the directory of the generated recipes, the directoy of previously ran generated recipes. The job_dir is used to store the job output after the sbatch submission. The sbatch_dir contains all the sbatch files that were submitted. **IMPORTANT** : The content of those two last files (job_dir and sbatch_dir) will be **deleted every time the automation is run**. So do not store anything important in them !

- Adapt *config-esm.yml* to your need like the *output_dir* field , the *CMIP5, OBS and CMIP6 paths*. 
If you run into memory error later on try to reduce the *max_parallel_tasks* to fix it.

- Then all you have to do is : `python __main__.py [option] ... [-r all/generate/esmvaltool] [-f TEMPLATE_XXX.yml] -e [True or False]`

### Possible command options :
You have the option of running three different modes:
- Runmode "generate": `python __main__.py` or `python __main__.py -r generate`. The default run mode. It will run the first part of runmode "all" i.e generate all the recipes from the templates (so it will fill in the template with the available datasets and save them in the *created_recipe_dir*)
- Runmode "esmvaltool": `python __main__.py -r esmvaltool`. It will run the second part of runmode "all"  i.e create the sbatch files in *sbatch_dir* used to run esmvaltool for each generated recipe. And it will submit those sbatch jobs. You can then get the job err and out files in the *job_dir*. And you will find all the output of the esmvaltool in the *output_dir* specified in *esm_config_file*.
- Runmode "all": `python __main__.py -r all`.First, it will run the generate mode (see below for more about the generate mode) on all the recipes in the template folder ( so it will fill in the template with the available datasets and save them in the *created_recipe_dir*). Then, it will create the sbatch files in *sbatch_dir* used to run esmvaltool for each generated recipe. Finally it will submit those sbatch jobs. You can then get the job err and out files in the *job_dir*. And you will find all the output of the esmvaltool in the *output_dir* specified in *esm_config_file*.
-
There is also the following options to run the "generate" runmode part of the tool :
- Option "one_file": example : `python __main__.py -r generate -f template_to_use.py`. To run the generate runmode on one file only. The resulting recipe will be created in the *generated_recipe_dir*. The default name of the created recipe replace the "template" string by "recipe". It is usefull when the only goal is to create a recipe from only one template with all available data quickly, instead of creating a new folder and modifying the *config.yml* file.

- Option "one_ensemble": when using the generate/all mode to fill in the generated recipe with one ensemble members per model (default), or all members available with the `-e` option.

### Additionnal information:
- Have a look at recipe/template to see already existing templates

### How to create a template of a recipe:
COMING SOON
#### Special flag and special cases

### Improvements to make:
- Have a files which all of the available datasets which is refresh every day
- use config file of esmvaltool instead of a custom one
