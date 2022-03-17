; CLOUDS_CALIPSO_POSTPROC
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; ############################################################################
; Description
;
; Modification history
;   20211210-lauer_axel: written.
;
; ############################################################################

begin

  inpath1 = "/home/b/b380103/DATA/ESMValTool2/DEVOBS/Tier2/CALIPSO-ICECLOUD/"
  inpath2 = "/home/b/b380103/workesm/esmvaltool_output/" \
            + "recipe_clouds_clim_20220125_151649/preproc/" \
            + "clouds_clim_cmip6/lwp/"

  diag = "clouds_calipso_postproc.ncl"

  refnames = (/"ESACCI-CLOUD", "CLARA-AVHRR", "CERES-EBAF", \
               "ESACCI-WATERVAPOUR", "MODIS", "PATMOS-x", "ISCCP-FH", \
               "MAC-LWP", "SSMI-MERIS", "CLOUDSAT-L2"/)

  fname1 = inpath1 + "OBS_CALIPSO-ICECLOUD_sat_1-00_Amon_cli_200701-201512.nc"
  infile1 = addfile(fname1, "r")

  ; read and average data

  cli = infile1->cli
  plev = infile1->plev

print(plev)

  dp = plev
  np = dimsizes(dp) - 1
  do i = 1, np - 1
    dp(i) = 0.5 * (plev(i - 1) - plev(i + 1))
  end do

  dp(0) = 101325. - 0.5 * (plev(0) + plev(1))
  dp(np) = 0.5 * (plev(np - 1) + plev(np))

print(dp)

  dp4d = conform(cli, dp, 1)
  clivi = dim_sum_n(cli * dp4d, 1)
  copy_VarCoords(cli(:, 0, :, :), clivi)

  outfile = "postproc_calipso.nc"
  system("rm " + outfile)
  debugfile = addfile(outfile, "c")
  debugfile->clivi = clivi
end
