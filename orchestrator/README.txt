9 August 2017
-------------
First working prototype of the backend integrated with yaml parser and its functionality. The backend runs an example diagnostic end-to-end. To see its workings, download the orchestrator, the data needed:

pr_Amon_MPI-ESM-LR_historical_r1i1p1_185001-200512.nc
ta_Amon_MPI-ESM-LR_historical_r1i1p1_200001-200512.nc

put the data in a test_data directory here, and then run python main.py nml/namelist_MyVar.yml. Tere is a whole lot of functionality still to be added in the bckend but the barebones (input file parsing, cmor reformatting, diagnostic running via NCL and outputting files) now works fine.
