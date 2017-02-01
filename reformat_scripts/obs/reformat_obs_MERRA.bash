#!/bin/bash -eu

###############################################################################
## REFORMAT SCRIPT FOR THE MERRA REANALYSIS DATA SET
###############################################################################
##
## Tier
##    Tier 2: other freely-available dataset.
##
## Source
##    http://gmao.gsfc.nasa.gov/merra/
##
## Last access
##    20150703
##
## Download and processing instructions
##    Followin the links for the MERRA products and chose the "on-the-fly"
##    data subsettter. Make sure to choose NetCDF as download format. This
##    script is written for monthly data input files.
##
## Caveats
##    Only variable 'pr' implemented so far
##
## Modification history
##    20??????-?????????: Written by Grigory Nikulin, SMHI
##    20150703-A_eval_ma: Rewritten sligthly to fit with ESMValTool
##
################################################################################

function log
{
 echo "[ $(date -u '+%Y-%m-%d  %H:%M') ]: " $*
}

info()
{
    log "*II* $*"
}

warning()
{
    log "*WW* $*"
}

error()
{
    log "*EE* $*" 1>&2
    exit 1
}

# VARIABLES
var_in=(pr)

# SIMULATIONS
run_in=(MERRA)
gcm_name=MERRA
experiment_id=reanalysis
gcm_version_id=r1i1p1

#---- FIRST and LAST YEARS ---
fy_in=(1979)   # first year
ly_in=(1984)   # last year

# --- INPUT / OUTPUT ---
path_in="/nobackup/rossby17/sm_maeva/tmp/tmp/in/"
path_out="/nobackup/rossby17/sm_maeva/tmp/tmp/out/"

# --- META INFO ---
ref_time='1950-01-01,00:00'  # reference time
mm_s=(01 02 03 04 05 06 07 08 09 10 11 12) # months
cdo='cdo -s'
freq=mon
table_out=Amon

# --- Number of variables, runs etc..  ---
num_var=${#var_in[@]}
num_run=${#run_in[@]}
num_fy=${#fy_in[@]}
num_ly=${#ly_in[@]}

# --- check if fy_in and ly_in have the same dimension ---
if [ $num_fy -ne $num_ly ]; then
    error "... fy_in and ly doesn't have the same dimension .... TERMINATED"
    exit 1
fi


# --------------------
# ---  MAIN BLOCK  ---
# --------------------
beg_time="$(date +%s)"


for ((run=0;run<=num_run-1;run++)); do  # simulation loop
    for ((var=0;var<=num_var-1;var++)); do   # variable loop

        case ${var_in[$var]} in
           pr) merra_name=prectot
               mask_in='.prod.assim.tavgM_2d_flx_Nx.'
               st_name='precipitation_flux '
               long_name='Precipitation'
               units='kg m-2 s-1'
              ;;

        esac

        # --- YEAR LOOP ---
        for ((nn=0;nn<num_fy;nn++)); do

            fy=${fy_in[$nn]}
            ly=${ly_in[$nn]}

            year_p=$(($ly - $fy + 1))  # number of years

            # --- number of output files ----
            sy=$year_p
            file_p=1



            # --- OUTPUT FILE LOOP ---
            for ((ff=1;ff<=file_p;ff++)); do

                fy_out=$((fy+(ff-1)*sy)) #  first year in output file
                ly_out=$((fy_out+sy-1))  #  last year in output file

                # --- output and temporary files ---
                file_out=$path_out${var_in[$var]}_$table_out'_'$gcm_name'_'$experiment_id'_'$gcm_version_id'_'$fy_out'-'$ly_out'.nc'
                file_tmp0=$path_out'tmp0_'${var_in[$var]}_$table_out'_'$gcm_name'_'$experiment_id'_'$gcm_version_id'_'$fy_out'-'$ly_out'.nc'

                # --- remove output and temporary files ---
                [ -f $file_out ] && rm $file_out
                [ -f $file_tmp0 ] && rm $file_tmp0

                # --- post-processing info ---
                info "... VARIABLE NETCDF      ' '... ${var_in[$var]} '|' $st_name '|' $long_name '|' $units"
                info "... SIMULATION' ... ' ${run_in[$run]}':' $gcm_name '|' $experiment_id '|' $gcm_version_id '|' $fy_out'-'$ly_out '|'"
                info "... OUTPUT FILE ... $file_out"

                for ((yy=fy_out;yy<=ly_out;yy++)); do  # year loop

                    for mm in 0 1 2 3 4 5 6 7 8 9 10 11 ; do  # month loop

                        # --- input file ---
                        if [ $yy -ge 1979 -a $yy -le 1992 ]; then
                            version='MERRA100'
                        fi

                        if [ $yy -ge 1993 -a $yy -le 2001 ]; then
                            version='MERRA200'
                        fi

                        if [ $yy -ge 2001 -a $yy -le 2011 ]; then
                            version='MERRA300'
                        fi

                        if [ $yy -eq 2010 ]; then
                            if  [ $mm -eq 4 -o $mm -eq 5 -o $mm -eq 6 -o $mm -eq 7 ]; then
                                version='MERRA301'
                            fi
                        fi

                        file_in=$path_in$version$mask_in$yy${mm_s[$mm]}'.SUB.nc'

                        if [ -f $file_in ]; then
                            info "... now processing ... $file_in"
                            $cdo cat -selname,$merra_name $file_in $file_tmp0
                        else
                            warning "... not exists $file_in"
                        fi
                    done  # month loop
                done # year loop


                if [ -f "$file_tmp0" ]; then

                    info "... FINALIZING OUTPUT FILE ... $file_out"

                    # FIX TIME
                    $cdo -setmissval,1.e20 -setreftime,$ref_time -shifttime,14day -setcalendar,standard $file_tmp0 $file_out

                    ncatted -a long_name,time,c,c,"time" \
                            -a bounds,time,c,c,"time_bnds" \
                            -a axis,time,c,c,"T" -h $file_out

                    # HORIZONTAL COORDINATES
                    ncrename -d longitude,lon -d latitude,lat -h $file_out
                    ncrename -v longitude,lon_old -v latitude,lat -h $file_out


                    # -180 to 180 --> 0 to 360
                    ncap2 -O -h -s ''${var_in[$var]}'='$merra_name'; '${var_in[$var]}'(:,:,0:269)='$merra_name'(:,:,270:539); '${var_in[$var]}'(:,:,270:539)='$merra_name'(:,:,0:269);' $file_out $file_out
                    ncap2 -O -h -s 'lon=lon_old; lon(0:269)=lon_old(270:539); lon(270:539)=lon_old(0:269);  where(lon < 0) lon=lon+360;' $file_out $file_out
                    ncks -O -x -v lon_old,$merra_name -h $file_out $file_out


                    # VARIABLE ATTRIBUTES
                    ncatted -a grid_name,${var_in[$var]},d,, \
                            -a comments,${var_in[$var]},d,, \
                            -a level_description,${var_in[$var]},d,, \
                            -a time_statistic,${var_in[$var]},d,, \
                            -a standard_name,${var_in[$var]},c,c,"$st_name" \
                            -a long_name,${var_in[$var]},m,c,"$long_name" \
                            -a units,${var_in[$var]},c,c,"$units" \
                            -a missing_value,${var_in[$var]},c,f,1.e+20 \
                            -a cell_methods,${var_in[$var]},c,c,"mean" -h $file_out


                    # ADDING TIME BOUNDS
                    # time bounds files
                    file_dw_tmp=$path_out'time_dw_tmp.nc'
                    file_up=$path_out'time_up.nc'
                    file_dw=$path_out'time_dw.nc'
                    file_bnds=$path_out'time_bnds.nc'
                    file_bnds_tmp=$path_out'tmp_time_bnds.nc'

                    [ -f $file_dw_tmp ] && rm $file_dw_tmp       # remove temporary lower time bound file if exist
                    [ -f $file_dw ] && rm $file_dw               # remove lower time bound file if exist
                    [ -f $file_up ] && rm $file_up               # remove upper time bound file if exist
                    [ -f $file_bnds_tmp ] && rm $file_bnds_tmp   # remove temporary time bound file if exist
                    [ -f $file_bnds ] && rm $file_bnds           # remove time bound file if exist

                    # fixing time bounds
                    ncks -O -v time,${var_in[$var]} -d lat,0 -d lon,0 -h $file_out $file_dw_tmp     # down time limit temporary

                    ncatted -a coordinates,${var_in[$var]},d,, -h $file_dw_tmp
                    ncatted -a grid_mapping,${var_in[$var]},d,, -h $file_dw_tmp
                    ncatted -a bounds,time,d,, -h $file_dw_tmp

                    $cdo -shifttime,-14days $file_dw_tmp $file_dw                                           # lower time limit
                    $cdo -shifttime,1month $file_dw $file_up                                                # upper time limit

                    ncrename -v time,time_dw -h $file_dw
                    ncrename -v time,time_up -h $file_up

                    ncks -v time_dw -h $file_dw $file_bnds_tmp
                    ncks -A -v time_up -h $file_up $file_bnds_tmp

                    ncap2 -h -O -s 'defdim("bnds",2); time_bnds[$time, bnds]=0.; time_bnds(:, 0)=time_dw; time_bnds(:, 1)=time_up;' -h $file_bnds_tmp $file_bnds
                    ncks -A -v time_bnds -h $file_bnds $file_out
                    ncatted -a bounds,time,c,c,"time_bnds" -h $file_out

                    # ---- remove temporary file (time bounds) --------
                    [ -f $file_dw_tmp ] && rm $file_dw_tmp           # remove temporary time bounds file if exist
                    [ -f $file_dw ] && rm $file_dw                   # remove temporary lower time bound file if exist
                    [ -f $file_up ] && rm $file_up                   # remove temporary upper time bound file if exist
                    [ -f $file_bnds ] && rm $file_bnds               # remove time bound file if exist
                    [ -f $file_bnds_tmp ] && rm $file_bnds_tmp       # remove temporary time bound file if exist
                    ncatted -a NCO,global,d,, -h $file_out


                    # -------------------------------------
                    # ---  DATES in FILE NAME: YYYYMMDD ---
                    # -------------------------------------
                    num_time=`$cdo ntime $file_out`
                    f_date=`$cdo showdate -seltimestep,1 $file_out`; f_date=${f_date// /}; f_date=${f_date//-/}
                    f_time=`$cdo showtime -seltimestep,1 $file_out`; f_time=${f_time:1:2}
                    l_date=`$cdo showdate -seltimestep,$num_time $file_out`; l_date=${l_date// /}; l_date=${l_date//-/}
                    l_time=`$cdo showtime -seltimestep,$num_time $file_out`; l_time=${l_time:1:2}
                    file_cmip=$path_out${var_in[$var]}_$table_out'_'$gcm_name'_'$experiment_id'_'$gcm_version_id'_'${f_date:0:6}'-'${l_date:0:6}'.nc'
                    mv $file_out $file_cmip

                    # --- GLOBAL ATTRIBUTES ----
                    ncatted -a Conventions,global,m,c,"CF-1.4" -h $file_cmip
                    ncatted -a source,global,c,c,"MERRA Reanalysis" -h $file_cmip
                    ncatted -a institution,global,c,c,"GSFC" -h $file_cmip
                    ncatted -a references,global,c,c,"http://gmao.gsfc.nasa.gov/merra/" -h $file_cmip
                    ncatted -a calendar,global,d,, -h $file_cmip
                    ncatted -a history,global,d,, -h $file_cmip
                    ncatted -a CDI,global,d,, -h $file_cmip
                    ncatted -a CDO,global,d,, -h $file_cmip
                    ncatted -a NCO,global,d,, -h $file_cmip


                    [ -f $file_tmp0 ] && rm $file_tmp0

                else
                    warning "... TEMPORARY FILE does not exist $file_tmp0"
                fi # if file_out exists
            done  # output file loop
        done  # year loop
    done  # variable loop
done  # simulation loop

end_time="$(date +%s)"
info "ELAPSED TOTAL TIME: "$(expr $end_time - $beg_time)" sec"
info "END"
