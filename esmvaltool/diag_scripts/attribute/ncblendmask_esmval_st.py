import numpy as np
import iris
import esmvalcore.preprocessor as eprep
# Following two functions for blending and masking modified from Cowtan 2015.
# Calculate blended temperatures using general methods
# Source of original Cowtan code:
# http://www-users.york.ac.uk//~kdc3/papers/robust2015/methods.html
# Usage:
#  python ncblendmask.py <mode> tas.nc tos.nc sic.nc sftlf.nc obs.nc dec_warming obs_dec_warming ann_warming gmst_comp_warming diag_name obs ensobs ensobs_diag ensobs_dec_warming
#  <mode> is one of xxx, mxx, xax, max, xxf, mxf, xaf, maf
# max means use masking of anomalies, with time-varying sea ice. See Cowtan website for more details.
# tas.nc. tos.nc, sic.nc, sftlf.nc and obs.nc are names of NetCDF files containing tas, tos, siconc, sftlf from the simulation.
# obs.nc is the name of the observations NetCDF file (the median if using an ensemble obs dataset).
# dec_warming is 2010-2019 warming in GSAT from the model.
# obs_dec_warming is 2010-2019 warming in GMST from the obs.
# ann_warming is a timeseries of annual mean GSAT from the model.
# gmst_comp_warming is 2010-2019 warming globally-complete GMST from the model.
# diag_name is an input diagnostic name e.g. gmst05, hemi10, where the last two digits are the averaging period.
# obs indicates which obs dataset is being used had5/had4 etc.
# ensobs is the partial filename of ensemble obs dataset, if ensemble data is used, otherwise empty string.
# ensobs_diag is the diabnostic requested for each of the ensemble members of an ensemble obs dataset.
# ensobs_dec_warming is 2010-2019 warming in GMST for each for each of the ensemble members of obs dataset.
#Outputs
# diag is the requested diagnostic (e.g. gmst05) for the model.
# obs_diag is the requested dignostic (e.g. gmst05) for the obs.

## Nathan Gillett - Adapted from ncblendmask-nc4.py from Cowtan 201

def ncblendmask_esmval(tas_file,obs_cb,dec_warming,obs_dec_warming,ann_warming,gmst_comp_warming,year_block, ens_cblst='',ensobs_diag=[],ensobs_dec_warming=[],warming_years=[2010,2019],warming_base=[1850,1900]):
# MAIN PROGRAM

# m = mask
# a = blend anomalies
# f = fix ice
# (use x for none)

    # read tas.nc
    tas_cb = iris.load_cube(tas_file)

    if tas_cb.shape != obs_cb.shape:
      #  assuming the difference is in region and not time
        tas_cb = eprep.extract_region(tas_cb, start_longitude=obs_cb.coord('longitude').points.min(), 
                                      end_longitude=obs_cb.coord('longitude').points.max(), 
                                      start_latitude=obs_cb.coord('latitude').points.min(), 
                                      end_latitude=obs_cb.coord('latitude').points.max())      

    cvgmsk = obs_cb.data.mask

    masked_tas_cb = tas_cb 
    masked_tas_cb.data.mask = masked_tas_cb.data.mask | cvgmsk

    # calculate diagnostic
    diag=calc_diag(masked_tas_cb, year_block) #Diagnostic for attribution analysis.
    dec_warming.append(calc_dec_warming(tas_cb,warming_years,warming_base)) #Diagnose SAT warming with global coverage for attributable trends.
    obs_dec_warming.append(calc_dec_warming(obs_cb,warming_years,warming_base))
    
    if ann_warming!=0:
      ann_warming.append(calc_diag(tas_cb,year_block)) #Calculate ann warming.
    if gmst_comp_warming!=0:
      gmst_comp_warming.append(calc_diag(tas_cb,year_block))
    obs_diag=calc_diag(obs_cb,year_block)

    #Repeat obs diagnostics for each member of ensemble observational dataset if ensobs is set.
    #Assume missing data mask is the same as for main obs dataset.
    if ens_cblst != '':
      #Assume 100 member ensemble observations dataset.
      for ens_cb in ens_cblst:
        ensobs_dec_warming.append(calc_dec_warming(ens_cb,warming_years,warming_base))
        ensobs_diag.append(calc_diag(ens_cb, year_block))
    return (diag,obs_diag)


def calc_diag(cube, year_block):
   
    area_av_cb = eprep.area_statistics(cube, 'mean')
    area_av_arr = area_av_cb.data
    reshaped_data = area_av_arr.reshape(-1, 12*year_block)
    block_arr = reshaped_data.mean(axis=1)
    diag = block_arr - block_arr.mean()
   
    return diag


def calc_dec_warming(cube,warming_years,warming_base):
    area_av_cb = eprep.area_statistics(cube, 'mean')
    ann_av_cb = eprep.annual_statistics(area_av_cb)
    warming_cb = eprep.extract_time(ann_av_cb, start_year=warming_years[0],
                                    start_month=1, start_day=1,
                                    end_year=warming_years[1], end_month=12, end_day=31)
    warming_base_cb = eprep.extract_time(ann_av_cb, start_year=warming_base[0],
                                    start_month=1, start_day=1,
                                    end_year=warming_base[1], end_month=12, end_day=31)
    
    warming = np.mean(warming_cb.data) - np.mean(warming_base_cb.data)

    return (warming)