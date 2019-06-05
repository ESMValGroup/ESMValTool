********************************************
How to prepare and run your first diagnostic
********************************************

Instructions for personal diagnostic
====================================

Anyone can run a personal diagnostic, no matter where the location of it;
there is no need to install esmvaltool in developer mode nor is it to
git push or for that matter, do any git operations; the example recipe

.. code-block:: bash
    esmvaltool/recipes/recipe_my_personal_diagnostic.yml

shows the use of running a personal diagnostic; the example

.. code-block:: bash
    esmvaltool/diag_scripts/examples/my_little_diagnostic.py

and any of its alterations may be used as training wheels for the future ESMValTool
diagnostic developer. The purpose of this example is to familiarize the user with
the framework of ESMValTool without the constraints of installing and running the
tool as developer. 

Functionality
=============

`my_little_diagnostic` (or whatever the user will call their diagnostic) makes full use
of ESMValTool's preprocessor output (both phyisical files and run variables); this output
comes in form of a nested dictionary, or config dictionary, see an example below;
it also makes full use of the ability to call any of the preprocessor's functions,
note that relative imports of modules from the esmvaltool package are allowed and
work without altering the $PYTHONPATH.

The user may parse this dictionary so that they execute a number of operations on the
preprocessed data; for example the `my_little_diagnostic.plot_time_series` grabs the
preprocessed data output, computes global area averages for each model, then plots
a time-series for each model. Different manipulation functionalities for grouping,
sorting etc of the data in the config dictionary are available,
please consult ESMValTool User Manual.


Writing a basic recipe
======================
The user will need to write a basic recipe to be able to run their own personal diagnostic.
An example of such a recipe is found in `esmvaltool/recipes/recipe_my_personal_diagnostic.yml`.
For general guidelines with regards to ESMValTool recipes please consult the User Guide;
the specific parameters needed by a recipe that runs a personal diagnostic are:

.. code-block:: bash
    scripts:
      my_diagnostic:
        script: /path/to/your/my_little_diagnostic.py        

i.e. the full path to the personal diagnostic that the user needs to run.

Example of config dictionary
============================
```
{'input_files':
['/group_workspaces/jasmin2/cmip6_prep/esmvaltool_users/valeriu/MyDIAG/recipe_my_personal_diagnostic_20181001_112918/preproc/simple_pp_ta/metadata.yml'],
'log_level': 'info',
'max_data_filesize': 100,
'myDiag': 'my_personal_diagnostic_example',
'myDiagPlace': '/group_workspaces/jasmin2/cmip6_prep/esmvaltool_users/valeriu',
'output_file_type': 'pdf',
'plot_dir': '/group_workspaces/jasmin2/cmip6_prep/esmvaltool_users/valeriu/MyDIAG/recipe_my_personal_diagnostic_20181001_112918/plots/simple/my_diagnostic', 'recipe': 'recipe_my_personal_diagnostic.yml',
'run_dir': '/group_workspaces/jasmin2/cmip6_prep/esmvaltool_users/valeriu/MyDIAG/recipe_my_personal_diagnostic_20181001_112918/run/simple/my_diagnostic',
'script': 'my_diagnostic',
'title': 'My First Diagnostic',
'version': '2.0a1',
'work_dir': '/group_workspaces/jasmin2/cmip6_prep/esmvaltool_users/valeriu/MyDIAG/recipe_my_personal_diagnostic_20181001_112918/work/simple/my_diagnostic',
'write_netcdf': True,
'write_plots': True,
'input_data': {'/group_workspaces/jasmin2/cmip6_prep/esmvaltool_users/valeriu/MyDIAG/recipe_my_personal_diagnostic_20181001_112918/preproc/simple_pp_ta/CMIP5_MPI-ESM-LR_Amon_historical_r1i1p1_T3M_ta_2000-2002.nc':
    {'cmor_table': 'CMIP5',
     'dataset': 'MPI-ESM-LR',
     'diagnostic': 'simple',
     'end_year': 2002,
     'ensemble': 'r1i1p1',
     'exp': 'historical',
     'field': 'T3M',
     'filename': '/group_workspaces/jasmin2/cmip6_prep/esmvaltool_users/valeriu/MyDIAG/recipe_my_personal_diagnostic_20181001_112918/preproc/simple_pp_ta/CMIP5_MPI-ESM-LR_Amon_historical_r1i1p1_T3M_ta_2000-2002.nc',
     'fx_files': {'areacello': '/badc/cmip5/data/cmip5/output1/MPI-M/MPI-ESM-LR/historical/fx/ocean/fx/r0i0p0/latest/areacello/areacello_fx_MPI-ESM-LR_historical_r0i0p0.nc', 'sftlf': '/badc/cmip5/data/cmip5/output1/MPI-M/MPI-ESM-LR/historical/fx/atmos/fx/r0i0p0/latest/sftlf/sftlf_fx_MPI-ESM-LR_historical_r0i0p0.nc', 'sftof': '/badc/cmip5/data/cmip5/output1/MPI-M/MPI-ESM-LR/historical/fx/ocean/fx/r0i0p0/latest/sftof/sftof_fx_MPI-ESM-LR_historical_r0i0p0.nc'},
     'long_name': 'Air Temperature',
     'mip': 'Amon',
     'preprocessor': 'pp',
     'project': 'CMIP5',
     'short_name': 'ta',
     'standard_name': 'air_temperature', 
     'start_year': 2000, 
     'units': 'K'
    } -- end of input_data member value (key: preprocessed file)
  } -- end of input_data dictionary
} -- end of config dictionary
```
