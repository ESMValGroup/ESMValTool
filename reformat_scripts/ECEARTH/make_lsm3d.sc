#!/bin/csh
###############################################################################
## CREATE A LAND-SEA MASK FROM SEA SALT FIELD, FOR SURFACE OR MULTILEVEL
###############################################################################
##
## Description
##    This script creates a land-sea mask netcdf file from an EC-Earth/NEMO
##    netcdf file. Where there is no sea salt, there is land, otherwise there is
##    sea. When variable $lsm is set to lsm a surface lsm is created from field
##    sosaline, if $lsm is set to lsm3d a multilevel field is create from field
##    vosaline.
## Notes
##    nco (NetCDF operators) must be installed on your system
##
## Modification history
##    20141204-A_vanu_be: written.
##
###############################################################################
#-------------------------------------------------------------------------------
# User definitions
#-------------------------------------------------------------------------------
# make surface (lsm) or multilevel (lsm3d) land sea mask
set lsm     = lsm # lsm | lsm3d

# set paths
set wdir    = $klad/EMBRACE/ESMVal/ORCA_files
set infile  = historic/2000/Output/ORCA1_MM_20000101_20001231_grid_T.nc
set outfile = ${lsm}_fx_ORCA1_historic_r0i0p0.nc

#-------------------------------------------------------------------------------
# No more user input required below
#-------------------------------------------------------------------------------
# goto working directory and remove any existing output files
cd $wdir
rm -f $outfile
rm -f tmp*_$outfile

# Define which fields to use and the properties
if ($lsm == lsm) then
  set infield = sosaline
  set copyflds = "$lsm,nav_lon,nav_lat"
  set associate = "nav_lat nav_lon"
  set axis = "YX"
else if ($lsm == lsm3d) then
  set infield = vosaline
  set copyflds = "$lsm,nav_lon,nav_lat,deptht"
  set associate = "deptht nav_lat nav_lon"
  set axis = "ZYX"
else
  echo "Unknown lsm type $lsm"
  exit(1)
endif

# create lsm
echo "Creating $lsm from $infile"
ncap2 -s "$lsm=int($infield==0.0)" $infile tmp1_$outfile    || echo "ncap2 failed" && exit(1)
ncwa -a time_counter tmp1_$outfile tmp2_$outfile            || echo "ncwa failed" && exit(1)
ncks -v $copyflds tmp2_$outfile $outfile                    || echo "ncks failed" && exit(1)

# adapt attributes
echo "Setting attributes"
ncatted -O -a associate,$lsm,o,c,"$associate" $outfile      || echo "ncatted  1 failed" && exit(1)
ncatted -O -a axis,$lsm,o,c,$axis $outfile                  || echo "ncatted  2 failed" && exit(1)
ncatted -O -a units,$lsm,o,c,1 $outfile                     || echo "ncatted  3 failed" && exit(1)
ncatted -O -a standard_name,$lsm,o,c,land_sea_mask $outfile || echo "ncatted  4 failed" && exit(1)
ncatted -O -a short_name,$lsm,o,c,land_sea_mask $outfile    || echo "ncatted  5 failed" && exit(1)
ncatted -O -a long_name,$lsm,o,c,"Land sea mask" $outfile   || echo "ncatted  6 failed" && exit(1)
ncatted -O -a valid_min,$lsm,d,, $outfile                   || echo "ncatted  7 failed" && exit(1)
ncatted -O -a valid_max,$lsm,d,, $outfile                   || echo "ncatted  8 failed" && exit(1)
ncatted -O -a missing_value,$lsm,d,, $outfile               || echo "ncatted  9 failed" && exit(1)
ncatted -O -a online_operation,$lsm,d,, $outfile            || echo "ncatted 10 failed" && exit(1)
ncatted -O -a interval_write,$lsm,d,, $outfile              || echo "ncatted 11 failed" && exit(1)
ncatted -O -a interval_operation,$lsm,d,, $outfile          || echo "ncatted 12 failed" && exit(1)

# clean up
rm tmp*_$outfile

# done
echo "Created $outfile"
exit(0)
