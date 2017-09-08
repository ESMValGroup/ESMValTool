9 August 2017
-------------
First working prototype of the backend integrated with yaml parser and its functionality. The backend runs an example diagnostic end-to-end. To see its workings, download the orchestrator, the data needed:

pr_Amon_MPI-ESM-LR_historical_r1i1p1_185001-200512.nc
ta_Amon_MPI-ESM-LR_historical_r1i1p1_200001-200512.nc

put the data in a test_data directory here, and then run the two examples:

python main.py -n nml/namelist_MyVar.yml -c config.ini
python main.py -n nml/namelist_perfmetrics.yml -c config.ini

config.ini is the configuration file that holds all the GLOBAL variables,
parsed by main.py and integrated in the analysis. You can also do multiple
preprocess instances, initiated from config.ini via preprocess_id parameter
