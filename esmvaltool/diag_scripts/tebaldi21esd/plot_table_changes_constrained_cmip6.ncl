; #############################################################################
;  Table of warming level crossing for constrained projections relative to 1986-2005 from cmip6 SSPs
;  Author: Debeire Kevin (DLR, Germany)
; #############################################################################
; Caveats
;    to be completed
;
; Modification history
;    to be completed
; #############################################################################
;diag_scripts="/mnt/lustre02/work/bd1083/b309165/esmvaltool/ESMValTool/esmvaltool/diag_scripts"
load "$diag_scripts/../interface_scripts/interface.ncl"
;
;load "$diag_scripts/shared/latlon.ncl"
;
;load "$diag_scripts/shared/statistics.ncl"
;
;load "$diag_scripts/shared/plot/aux_plotting.ncl"
;load "$diag_scripts/shared/plot/style.ncl"

begin
  file_type = "png"
  ; Output netcdf directory
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)
end

begin
  ; #############################################
  ; # Get parameters from ./variable_defs/*.ncl #
  ; # passed via the 'info' attribute           #
  ; #############################################
  ssp585_data= asciiread("$diag_scripts/tebaldi21esd/constrained_gsat_ssp585.txt",(/2,4/),"float")
  ssp126_data= asciiread("$diag_scripts/tebaldi21esd/constrained_gsat_ssp126.txt",(/2,4/),"float")
  ssp245_data= asciiread("$diag_scripts/tebaldi21esd/constrained_gsat_ssp245.txt",(/2,4/),"float")
  ssp585_2041_2060 = ssp585_data(0,:)
  ssp585_2081_2100 = ssp585_data(1,:)
  ssp245_2041_2060 = ssp245_data(0,:)
  ssp245_2081_2100 = ssp245_data(1,:)
  ssp126_2041_2060 = ssp126_data(0,:)
  ssp126_2081_2100 = ssp126_data(1,:)

  mean_ = (/ (/dim_avg(ssp126_2041_2060(2)), dim_avg(ssp245_2041_2060(2)), dim_avg(ssp585_2041_2060(2))/),\
            (/dim_avg(ssp126_2081_2100(2)), dim_avg(ssp245_2081_2100(2)), dim_avg(ssp585_2081_2100(2))/)/)
  five_perc_ = (/ (/dim_avg(ssp126_2041_2060(1)), dim_avg(ssp245_2041_2060(1)), dim_avg(ssp585_2041_2060(1))/),\
            (/dim_avg(ssp126_2081_2100(1)), dim_avg(ssp245_2081_2100(1)), dim_avg(ssp585_2081_2100(1))/)/)
  ninetyfive_perc_ = (/ (/dim_avg(ssp126_2041_2060(3)), dim_avg(ssp245_2041_2060(3)), dim_avg(ssp585_2041_2060(3))/),\
            (/dim_avg(ssp126_2081_2100(3)), dim_avg(ssp245_2081_2100(3)), dim_avg(ssp585_2081_2100(3))/)/)

  period_list = (/"2041-2060","2081-2100"/)

  n_PER = dimsizes(period_list)
  n_SCEN= 3
  ; ###########################################
  ; *********** PLOTTING ************
  ; Create plot variables
  outfile = "TAS_cmip6_constrained_table_changes"
  wks = gsn_open_wks("png", plot_dir+"/"+outfile)
  res        = True
  res2 = True
  res3 = True
  res4 = True
  res@gsLineThicknessF = 2.
  res2@gsLineThicknessF = 2.
  res3@gsLineThicknessF = 2.
  ncr  = (/n_PER,n_SCEN/)
  ncr2  = (/n_PER+1,1/)
  ncr3= (/1,n_SCEN/)
  text= new((/n_SCEN+1,n_SCEN/),string)
  text2= new((/n_SCEN+1,1/),string)
  text3 = (/"SSP1-2.6","SSP2-4.5","SSP5-8.5"/)
  do i = 0 , n_PER-1
    text2(i+1,0) = period_list(i)
    do j = 0, n_SCEN -1
        text(i,j) = sprintf("%3.2f",mean_(i,j))+"~C~["+sprintf("%3.2f",five_perc_(i,j))+","+sprintf("%3.2f",ninetyfive_perc_(i,j))+"]"

    end do
  end do
  res2@txJust = "CenterCenter"
  res@txJust = "CenterCenter"
  res3@txJust = "CenterCenter"

  res@txFontHeightF  = 0.013
  res2@txFontHeightF= 0.02
  res3@txFontHeightF= 0.02

  res2@gsFillColor = "gray"
  res3@gsFillColor = "gray"

  x = (/0.30,0.95/)
  y = (/0.45,0.65/)
  x2 = (/0.05,0.30/)
  y2 = (/0.45,0.75/)
  x3 = (/0.30,0.95/)
  y3 = (/0.65,0.75/)
  ;Title
  x4 = (/0.05,0.95/)
  y4 = (/0.77,0.79/)
  ncr4 = (/1,1/)
  res4@gsLineColor   = "transparent"
  res4@txFontHeightF = 0.02
  text4 = "Table of TAS change relative to 1986-2005, CMIP6 constrained projections (~S~o~N~C)"

  gsn_table(wks,ncr,x,y,text,res)
  gsn_table(wks,ncr2,x2,y2,text2,res2)
  gsn_table(wks,ncr3,x3,y3,text3,res3)
  gsn_table(wks,ncr4,x4,y4,text4,res4)
  draw(wks)
  frame(wks)
end
