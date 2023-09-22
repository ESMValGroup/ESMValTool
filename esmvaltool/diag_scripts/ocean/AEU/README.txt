Instruction to reproduce the AEU analysis:
# NB: 1. need to change some paths because I moved stuff after all was done
#     2. need a couple of modified esmvaltool scripts:
#        /home/users/gig/git/ESMValTool/esmvaltool/utils/find_AEU_gig.py
#        /home/users/gig/git/ESMValTool/esmvaltool/diag_scripts/ocean/diagnostic_aeu_gig.py
#        for good measure I copy them in /home/users/gig/CRACAB/ 
#        since they were never committed anywhere


# 1. find all the files with 'uo'
python ESMValTool/esmvaltool/utils/find_AEU_gig.py

# 2. find all models that have both hist and ssp and skim ("screma") those that don't
python screma_files.py

# 3. append one file to the recipe template
bash append_one_file_to_all_recipes.sh

# 4. run all recipes
#    (need to change RECIPEDIR=ENSEMBLE_RUN_monthly_all/RECIPES)
#    this execute each recipe in the /RECIPES/ directory
#    this will run diagnostic_aeu_gig.py which saves intermediate steps to .nc files
#    so one can do the means after
bash main_job_array.sh

# Now post process the intermediate files and plot the results
# all of the plot scripts are independent from one another... I think
# the final plotting scripts should be in
# cd CRACAB/FINAL_PLOT/ 

# 5. plot AEU timeseries
#    there's a bunch of versions of the script, this is the last one I used 
python pproc_and_plot_AEU_0623.py

# 6. plot depth profiles
python plot_AEU_vertical_profile.py

# 7. plot seasonality
python pproc_and_plot_AEU_seasonality.py

# 8. plot transect current state vs data
python plot_AEU_transect_v_data.py

# 9. plot AEU transect ssh - hist
python plot_AEU_transects_delta.py

cd /home/users/gig/CRACAB/FINAL_PLOT

# 10. join plots to have the composite 
python join_plots.py

