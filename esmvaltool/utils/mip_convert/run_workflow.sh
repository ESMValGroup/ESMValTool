# eample how to run the workflow
# first set up all relevant mip_convert rose suites: WORKS
mip_convert_setup -c config-user.yml -r recipe_mip_convert.yml -m setup
# then run the mip convert suites (somehow, WIP, replaced with sleep for now)
sleep 5
# check the suites have finished (somehow, WIP, replaced with sleep for now)
sleep 5
# then finally simlink the CMORized netcdf files: WORKS
mip_convert_setup -c config-user.yml -r recipe_mip_convert.yml -m postproc
