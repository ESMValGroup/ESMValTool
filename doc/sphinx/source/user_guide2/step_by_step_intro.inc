
Introduction step by step
=========================

Currently this documentation shows a step by step learning on ESMValTool v2 without any background in v1.
It is supposed that the user has a installed and working version of the tool in a computer or server, 
basic background with Python would be useful.

File Structure
===============


├── conda_build_config.yaml
├── doc
│   └── sphinx
│       └── source
├── docker
│   ├── Dockerfile
│   └── docker_guide.md
├── environment.yml
├── esmvaltool ---------------------> core folder of the library/tool
│   ├── cmor
│   │   ├── check.py
│   │   ├── _fixes
│   │   ├── fix.py
│   │   ├── __init__.py
│   │   ├── table.py
│   │   └── tables
│   ├── config-developer.yml
│   ├── config-logging.yml
│   ├── _config.py
│   ├── config-references.yml
│   ├── config-user.yml
│   ├── _data_finder.py
│   ├── diag_scripts
│   │   ├── autoassess
│   │   ├── climate_metrics
│   │   ├── clouds
│   │   ├── examples
│   │   ├── __init__.py
│   │   ├── ipcc_ar5
│   │   ├── ocean
│   │   ├── perfmetrics
│   │   ├── shared
│   │   └── validation.py
│   ├── __init__.py
│   ├── interface_scripts
│   │   ├── auxiliary.ncl
│   │   ├── constants.ncl
│   │   ├── data_handling.ncl
│   │   ├── interface.ncl
│   │   ├── logging.ncl
│   │   ├── write_header.ncl
│   │   └── write_references.ncl
│   ├── _main.py
│   ├── preprocessor
│   │   ├── _area_pp.py
│   │   ├── _derive.py
│   │   ├── _download.py
│   │   ├── __init__.py
│   │   ├── _io.py
│   │   ├── _mask.py
│   │   ├── _multimodel.py
│   │   ├── ne_masks
│   │   ├── _reformat.py
│   │   ├── _regrid.py
│   │   ├── _time_area.py
│   │   └── _volume_pp.py
│   ├── _recipe.py
│   ├── recipes
│   │   ├── examples
│   │   ├── recipe_autoassess_landsurface_snow.yml
│   │   ├── recipe_autoassess_landsurface_surfrad.yml
│   │   ├── recipe_autoassess_radiation_rms_Amon_all.yml
│   │   ├── recipe_autoassess_radiation_rms_Amon_obs.yml
│   │   ├── recipe_autoassess_radiation_rms_cfMon_all.yml
│   │   ├── recipe_autoassess_stratosphere.yml
│   │   ├── recipe_clouds_bias.yml
│   │   ├── recipe_clouds_ipcc.yml
│   │   ├── recipe_flato13ipcc.yml
│   │   ├── recipe_lauer13jclim.yml
│   │   ├── recipe_my_personal_diagnostic.yml
│   │   ├── recipe_OceanBGC.yml
│   │   ├── recipe_OceanPhysics.yml
│   │   ├── recipe_OceanQuadMap.yml
│   │   ├── recipe_perfmetrics_CMIP5.yml
│   │   ├── recipe_SeaIceExtent.yml
│   │   └── recipe_validation.yml
│   ├── recipe_schema.yml
│   ├── _task.py
│   ├── utils
│   │   ├── batch-jobs
│   │   ├── editor-enhancements
│   │   ├── __init__.py
│   │   ├── nclcodestyle
│   │   └── xml2yml
│   └── _version.py
├── LICENSE
├── meta.yaml
├── NOTICE
├── README.md
├── readthedocs.yml
├── RELEASE_NOTES
├── setup.cfg
├── setup.py ---------------------------> setup.py for manual installation
└── tests ------------------------------> testing units to ensure code works
    ├── __init__.py
    ├── integration
    │   ├── cmor
    │   ├── data_finder.yml
    │   ├── __init__.py
    │   ├── preprocessor
    │   └── test_data_finder.py
    ├── system
    │   ├── config-test.yml
    │   ├── data_simulator.py
    │   ├── esmvaltool_testlib.py
    │   ├── __init__.py
    │   └── test_recipes.py
    └── unit
        ├── cmor
        ├── data_finder
        ├── __init__.py
        ├── preprocessor
        ├── test_lint.py
        └── test_naming.py

36 directories, 78 files



**Where to begin** to understand the code? 

By inspecting the setup.py the console entry is done by *esmvaltool._main:run*, so we will explore this
file to understand how is the step by step from the point of view of the code.

Forgetting the exceptions we get some arguments that are going to main()
```
def run():
    args = get_args()
    conf = main(args)
```

And removing the exceptions and logs parts of the code it basically:
```
def main(args):
    """Define the `esmvaltool` program."""
    recipe = args.recipe
    if not os.path.exists(recipe):
        installed_recipe = os.path.join(
            os.path.dirname(__file__), 'recipes', recipe)
        if os.path.exists(installed_recipe):
            recipe = installed_recipe
    recipe = os.path.abspath(os.path.expandvars(os.path.expanduser(recipe)))

    config_file = os.path.abspath(os.path.expandvars(os.path.expanduser(args.config_file)))

    # Read user config file
    if not os.path.exists(config_file):
        print("ERROR: config file {} does not exist".format(config_file))

    recipe_name = os.path.splitext(os.path.basename(recipe))[0]
    cfg = read_config_user_file(config_file, recipe_name)

    ### Until here we have set up the recipe and read the configuration file
    ### taking in account that recipe.

    ### Now we begin to process things (logs not mentioned)

    # 1. Create run dir
    if os.path.exists(cfg['run_dir']):
        print("ERROR: run_dir {} already exists, aborting to "
              "prevent data loss".format(cfg['output_dir']))
    os.makedirs(cfg['run_dir'])

    # 2. Process information related with downloading if necessary
    cfg['synda_download'] = args.synda_download
    for limit in ('max_datasets', 'max_years'):
        value = getattr(args, limit)
        if value is not None:
            if value < 1:
                raise ValueError("--{} should be larger than 0.".format(
                    limit.replace('_', '-')))
            cfg[limit] = value

    resource_log = os.path.join(cfg['run_dir'], 'resource_usage.txt')

    # 3. Now it begins the processing of the recipe

    with resource_usage_logger(pid=os.getpid(), filename=resource_log):
        process_recipe(recipe_file=recipe, config_user=cfg)

    return cfg

```

Nothing heavy has been already done, but now we have to explore the process_recipe function
but it basically read the recipe file and create a python object (whose class is named recipe)
that has a method named run so at some point we have to understand **recipe.run()**

The **recipe object** is something relatively complex with several methods with decorators that
require carefull reading. My understanding indice that this object has the following components:

```
        """Parse a recipe file into an object."""
        self._cfg = config_user
        self._recipe_file = os.path.basename(recipe_file)
        self._preprocessors = raw_recipe.get('preprocessors', {})
        if 'default' not in self._preprocessors:
            self._preprocessors['default'] = {}
        self._support_ncl = self._need_ncl(raw_recipe['diagnostics'])
        self.diagnostics = self._initialize_diagnostics(
            raw_recipe['diagnostics'], raw_recipe.get('datasets', []))
        self.tasks = self.initialize_tasks() 
```

So now we understand better what a recipe, or at least what components it has:

1. user configuration [ymal file?]
2. recipe file        [yaml file?]
3. preprocessors      [a recipe will have a number of pre-processing steps defined here]
4. support_ncl        [legacy dependence in ncl or to use ncl?]
5. diagnostics        [derived information using diagnostics on raw-recipe and datasets]
6. tasks              [probably all the other aspects are coded in a number of tasks]


