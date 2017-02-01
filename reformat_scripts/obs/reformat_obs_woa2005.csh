#!/usr/bin/env csh -eu
###############################################################################
## REFORMAT SCRIPT FOR THE woa2005 OBSERVATIONAL DATA
###############################################################################
##
## Tier
##    Tier 2: other freely-available datasets.
##
## Source
##    Download at
##    http://www.eps.mcgill.ca/~dbianchi/Sites/personal_page/oxygen.html
##
## Download and processing instructions
##    Select "Monthly WOA05 with linear correction based on in situ GLODAP
##    oxygen" (o2_woa05_linear_monthly.nc)
##    Run this script (requires requires NCO, http://nco.sourceforge.net/ and
##    CDO, http://https://code.zmaw.de/projects/cdo/)
##
## Caveats
##
## Modification history
##    20151111-A_laue_ax: written.
##
###############################################################################

set inpath="${ESMValTool_RAWOBSPATH}/Tier2/woa2005"
set outpath="${ESMValTool_OBSPATH}/Tier2/woa2005"

if (! -d $outpath) then
    mkdir -p $outpath
endif

set infile=$inpath/o2_woa05_linear_monthly.nc
set outfile=$outpath/O2_monthly_woa2005_bianchi_reg_2005-2005.nc

if (-e $infile) then
    echo 'input file = '$infile
else
    echo 'error: input file not found - '$infile
    exit
endif

cp $infile tmp.nc
ncatted -a calendar,TIME,c,c,noleap tmp.nc
ncatted -a units,TIME,m,c,"days since 2005-01-01 00:00:00" tmp.nc $outfile
ncrename -v O2_LINEAR,O2 $outfile
ncap2 -s O2="O2*0.001" $outfile O2_monthly_woa2005_bianchi_reg_2005-2005_ncap.nc
cdo chunit,"mmol/m3","mol m-3" O2_monthly_woa2005_bianchi_reg_2005-2005_ncap.nc O2_monthly_woa2005_bianchi_reg_2005-2005_chunit.nc
mv O2_monthly_woa2005_bianchi_reg_2005-2005_chunit.nc $outfile

ncatted -O -a title,global,a,c,"World Ocean Atlas (WOA) 2005 dissolved oxygen concentration data with corrections applied as described in Bianchi et al. (2012)" $outfile
ncatted -O -a source,global,a,c,"http://www.eps.mcgill.ca/~dbianchi/Sites/personal_page/oxygen.html" $outfile
ncatted -O -a tier,global,a,c,"2" $outfile
ncatted -O -a period,global,a,c,"2005" $outfile
ncatted -O -a reference,global,a,c,"Garcia, H. E., R. A. Locarnini, T. P. Boyer, and J. I. Antonov, 2006. World Ocean Atlas 2005, Volume 3: Dissolved Oxygen, Apparent Oxygen Utilization, and Oxygen Saturation. S. Levitus, Ed. NOAA Atlas NESDIS 63, U.S. Government Printing Office, Washington, D.C., 342 pp.; Bianchi, D., Dunne, J. P., Sarmiento, J. L., and Galbraith, E. D.: Data-based estimates of suboxia, denitrification, and N2O production in the ocean and their sensitivities to dissolved O2, Global Biogeochem. Cy., 26, GB2009, doi:10.1029/2011GB004209, 2012." $outfile

if (-e $outfile) then
    rm -f tmp.nc
    rm -f O2_monthly_woa2005_bianchi_reg_2005-2005_ncap.nc
    echo 'created '$outfile
else
    echo 'error: no output written'
endif

