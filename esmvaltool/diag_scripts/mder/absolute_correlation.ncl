; #############################################################################
; mder/absolute_correlation.ncl
; #############################################################################
; MULTIPLE DIAGNOSTIC ENSEMBLE REGRESSION
; Author: Sabrina Wenzel (DLR, Germany)
; ESMVal project PI: Veronika Eyring (DLR, Germany)
; #############################################################################
;
; Description:
;     This routine calculates absolute correlation coefficients s between all
;     diagnostics, calculated for a present-day period, called in the namelist
;     and a future variable of interest.
;
;     Calls one plotscript:
;     Plots the absolute correlation coefficient between the variable of
;     interest and present-day diagnostic.
;
; Required recipe attributes:
;     p_time : time range for projections.
;     p_step : averaging time range for target period.
;     scale_time : a year to scale the timeseries with.
;     time_oper : select for opperation in time_opperation function.
;     time_opt : select for option in time_opperation function.
;     calc_typ :
;
; Optional recipe attributes:
;
; Caveats:
;     Modifies dataset_info@dataset and dataset_info@ensemble.
;
; Modification history:
;   * 20180726_A-schl_ma: ported to v2.0
;   * 201307??_A-wenz_sa: written
;
; #############################################################################

load "interface_scripts/interface.ncl"

load "./diag_scripts/shared/statistics.ncl"
load "./diag_scripts/shared/set_operators.ncl"
load "./diag_scripts/shared/latlon.ncl"

load "./diag_scripts/shared/plot/mder.ncl"

; load "./diag_scripts/shared/meta_data.ncl"

load "./diag_scripts/shared/mder.ncl"


begin

  enter_msg(diag_script, "")
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + variables(0) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  print(dataset_info)
  print(variable_info)
  print(config_user_info)
  print(input_file_info)
  print(diag_script_info)

  ; Variable
  var0 = variables(0)
  field_type0 = field_types(0)

  ; Extract observations from model-list
  IDX_OBS = ind(dataset_info@project .eq. "OBS")
  DIM_OBS = dimsizes(IDX_OBS)

  ; Extract historical simulations from datasets
  IDX_HDAT = ind(input_file_info@start_year .ne. \
                 max(tofloat(input_file_info@start_year)) .and. \
                 dataset_info@project .ne. "OBS")
  DIM_HDAT = dimsizes(IDX_HDAT)

  ; Extract individual RCP 4.5 simulations from datasets (= unique models)
  avgens = False
  if (isatt(diag_script_info, "average_ens")) then
    avgens = diag_script_info@average_ens
  end if
  if (avgens .eq. "True") then
    idx_datasets = UNIQ(dataset_info@dataset(DIM_OBS + DIM_HDAT:)) + \
                   DIM_OBS + DIM_HDAT
    avgens = True
  else
    idx_datasets = ispan(DIM_OBS + DIM_HDAT, \
                         dimsizes(dataset_info@dataset) - 1, 1)
    idx_hist = UNIQ(dataset_info@dataset(DIM_OBS:DIM_OBS + DIM_HDAT))
    if (dimsizes(idx_hist) .gt. 1) then
      all_datasets = unique_labels_min((/"dataset", "ensemble"/))
      dataset_info@dataset := all_datasets
      delete(all_datasets)
    end if
    delete(idx_hist)
    avgens = False
  end if

  ; Get dataset variables
  DIM_DAT = dimsizes(dataset_info@dataset(idx_datasets))
  dataset_info@ensemble := dataset_info@ensemble(idx_datasets)
  obs_name = dataset_info@dataset(IDX_OBS)

end

begin

  ; ---------------------------------------------------------------------------
  ; Read recipe and config data
  ; ---------------------------------------------------------------------------

  ; Read latitude range and level from input files
  f = addfile(input_file_info@filename(0), "r")
  if (isfilevar(f, "plev") .and. isfilevar(f, "lat")) then
    plev = read_data(0, "plev", field_type0)
    lat = read_data(0, "lat", field_type0)
    if (dimsizes(plev) .eq. 1) then
      lev = plev(0)
    else
      error_msg("f", diag_script, "", "ambiguous pressure levels given, " + \
                "only one supported")
    end if
    if (dimsizes(lat) .gt. 1) then
      lat_range = (/min(lat), max(lat)/)
    else
      lat_range = read_data(0, "lat_bnds", field_type0)
    end if
  else
    if (isatt(f, "plev") .and. \
        isatt(f, "lat_range_0") .and. \
        isatt(f, "lat_range_1")) then
      lev = f@plev
      lat_range = (/f@lat_range_0, f@lat_range_1/)
    else
      error_msg("f", diag_script, "", "cannot read plev and latrange")
    end if
  end if
  log_info("Retrieved pressure level " + lev + " and latitude range " + \
           lat_range(0) + " to " + lat_range(1) + " from preprocessed data")

  ; Global attributes
  if (isatt(diag_script_info, "p_time"))then
    p_time = diag_script_info@p_time
  else
    error_msg("f", diag_script, "", "requires attribute p_time in recipe")
  end if
  if (isatt(diag_script_info, "p_step"))then
    p_step = diag_script_info@p_step
  else
    error_msg("f", diag_script, "", "requires attribute p_step in recipe")
  end if
  if (isatt(diag_script_info, "scal_time"))then
    scal_time = diag_script_info@scal_time
  else
    error_msg("f", diag_script, "", "requires attribute scal_time in recipe")
  end if
  if (isatt(diag_script_info, "time_oper"))then
    time_oper = diag_script_info@time_oper
  else
    error_msg("f", diag_script, "", "requires attribute time_oper in recipe")
  end if
  if (isatt(diag_script_info, "time_opt"))then
    time_opt = diag_script_info@time_opt
  else
    error_msg("f", diag_script, "", "requires attribute time_opt in recipe")
  end if
  if (isatt(diag_script_info, "calc_typ"))then
    calc_type = diag_script_info@calc_typ
  else
    error_msg("f", diag_script, "", "requires attribute calc_typ in recipe")
  end if

  ; Output dir
  diag_script_base = basename(diag_script)
  plot_dir = config_user_info@plot_dir
  plot_num = integertochar(ispan(97, 122, 1))

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Output arrays
  xx = new((/dimsizes(p_time), DIM_DAT/), float)

;   ;;===================================================================
;   ;;MAIN
;   ;;Call the variable of interest and merge to one timeline yy
;   do imod = 0,DIM_DAT-1
;     i = idx_datasets(imod)
;
;     ;;Call projection simulations for regression analysis and time series,
;     ;;scale time series and avg over 10-year period (see also cfg_mder.ncl)
;     data = read_data(i, var0, field_type0)
;     tmp = extract_data(i, data, -1, 0, 0)
;     idxH = ind(dataset_info@dataset(i).eq.dataset_info@dataset(:DIM_OBS+DIM_HDAT-1))
;     idxF = ind(dataset_info@dataset(i).eq.dataset_info@dataset(DIM_OBS+DIM_HDAT:))+DIM_OBS+DIM_HDAT
;
;     if (dimsizes(idxH).ge.1) then
;       ;;Call historical simulations for time series
;       tmph = extract_data(idxH(0), data, -1, 0, 0)
;
;       ;;Average over ensemble members if more than one exists (see namlist)
;       if (getenv("ESMValTool_average_ens").eq."True") then
;         if (dimsizes(idxF).gt.1) then
;           if (dimsizes(idxF).eq.dimsizes(idxH)) then
;             tmpAF = new((/dimsizes(idxF),dimsizes(tmp&time)/),typeof(tmp))
;             copy_VarCoords(tmp,tmpAF(0,:))
;             tmpAF(0,:) = tmp
;             delete(tmp)
;
;             tmpAH = new((/dimsizes(idxH),dimsizes(tmph&time)/),typeof(tmph))
;             copy_VarCoords(tmph,tmpAH(0,:))
;             tmpAH(0,:) = tmph
;             delete(tmph)
;
;             ;;Loop over ensemble members, read the data and save in arry
;             do iens = 1,dimsizes(idxF)-1
;               tmpF = extract_data(idxF(iens), data, -1, 0, 0)
;               tmpAF(iens,:) = (/tmpF/)
;
;               tmpH = extract_data(idxH(iens), data, -1, 0, 0)
;               tmpAH(iens,:) = (/tmpH/)
;               delete([/tmpF,tmpH/])
;             end do
;
;             ;;Average over ensemble members
;             tmp = tofloat(dim_avg_n_Wrap(tmpAF, 0))
;             tmph = tofloat(dim_avg_n_Wrap(tmpAH,0))
;             copy_VarCoords(tmpAF(0,:),tmp)
;             copy_VarCoords(tmpAH(0,:),tmph)
;             delete([/tmpAF,tmpAH/])
;           else
;             print("fatal: inconsistend ensemble member size for model " \
;                   + dataset_info@dataset(i) + " in " + diag_script)
;             status_exit(1)
;           end if
;         end if
;       end if
;       tmpa = array_append_record(tmph,tmp,0)
;       tmpy = tmpa
;       delete([/tmp,tmpa,tmph,idxH,idxF/])
;     end if
;
;     do pt = 0, dimsizes(p_time)-1
;
;         ;;Average over var time series
;         tmps  = time_operations(tmpy,p_time(pt), (/p_time(pt) + p_step/), \
;                                 time_oper, time_opt, True)
;         tmpj  = time_operations(tmps,p_time(pt), (/p_time(pt) + p_step/), \
;                                 "average", "yearly", True)
;
;         ;;Calculate the scaling factor between scal_time(0) and scal_time(1)
;         if (scal_time(1).gt.toint(models@end_year(i-DIM_HDAT))) then
;             scal_time(1) = toint(models@end_year(i-DIM_HDAT))
;         end if
;
;         tmpxx = time_operations(tmpy,scal_time(0),scal_time(1),"extract",time_opt,True)
;         tmpyy = time_operations(tmpxx,scal_time(0),scal_time(1),"average","yearly",True)
;         tmppd = tofloat(avg(tmpyy))
;
;         if (calc_type.eq."trend") then
;             xx(pt, imod) = tofloat(avg(tmpj))-tmppd
;         else if (calc_type.eq."pos") then
;             xx(pt, imod) = tofloat(avg(tmpj))
;         end if
;         end if
;
;         variable_info@units = "deg"
;
;         delete([/tmppd,tmps,tmpj/])
;         delete([/tmpyy,tmpxx/])
;
;     end do
;
;     delete(tmpy)
;     scal_time(1) = diag_script_info@scal_time(1)
;
;   end do ;end imod
;
;   ;;Save unique model names
;   tmpm = dataset_info@dataset
;   delete(dataset_info@dataset)
;   dataset_info@dataset  = tmpm(idx_datasets)
;   delete(tmpm)
;
;   ;;Call diagnostics
;   ;; Read diagnostic list and add files together
;   xml_name_char = tochar(getenv("ESMValTool_xml_name"))
;   temp_dir = get_output_dir(getenv("ESMValTool_wrk_dir"), "temp") + \
;             tostring(xml_name_char(:dimsizes(xml_name_char) - 6)) + \
;             "_diagnostic.nc"
;
;   if (isfilepresent(temp_dir)) then
;
;     ;; Read diagnostic list
;     temp_file = addfile(temp_dir, "r")
;     temp_list = temp_file->temp_list
;     file_list = tostring(temp_list)
;     file_list@_FillValue = default_fillvalue(typeof(file_list))
;
;     delete(temp_list)
;
;     ;; Remove duplicates
;     do ii = 0, dimsizes(file_list) - 1
;       curr = file_list(ii)
;       idx = ind(curr.eq.file_list)
;       if (dimsizes(idx).gt.1) then
;         file_list(idx) = file_list@_FillValue
;         file_list(ii) = curr
;       end if
;       delete(curr)
;       delete(idx)
;     end do
;     tmp = file_list(ind(.not.ismissing(file_list)))
;     delete(file_list)
;     file_list = tmp
;     delete(tmp)
;
;     ;;Create diagnostic arreys
;     dim_DIAG = dimsizes(file_list)
;     diag_all = new((/dim_DIAG,DIM_DAT/),float)
;     wdiag = new((/dim_DIAG/),string)           ;;array with diagnostic names
;     obs0  = new((/dim_DIAG/),float)
;
;     ;; Loop over files in the list, append data
;     do ii = 0, dim_DIAG - 1
;         data_temp = ncdf_read(file_list(ii), "diagn")
;
;         do imod = 0,DIM_DAT-1
;             info_output(dataset_info@dataset(imod),verbosity, 2)
;             do dmod = 0,dimsizes(data_temp&models)-1
;                 if (dataset_info@dataset(imod).eq.data_temp&models(dmod)) then
;                     diag_all(ii,imod) = data_temp(0,dmod,0) ;;collect diagnostics in one array
;                     if (isatt(data_temp,"wdiag")) then
;                         wdiag(ii) = data_temp@wdiag
;                     end if
;                 end if
;             end do
;         end do
;         delete(data_temp)
;     end do
;
;     diag_all&diagnostics = wdiag
;   else
;     print("fatal: cannot find diagnostic list in " + diag_script)
;     status_exit(1)
;   end if
;
;   ;;===================================================================
;   ;;Correlation calculations
;
;   corr     = new((/dimsizes(p_time),dim_DIAG/),float)
;   diag_sign = new((/dimsizes(p_time),dim_DIAG/),string)
;
;   do pt = 0, dimsizes(p_time)-1
;
;       corr_sig = new((/2,dim_DIAG/),float)
;       r_sig0  = new((/dim_DIAG/),float)
;       r_sig1  = new((/dim_DIAG/),float)
;       r_wdiag  = new((/dim_DIAG/),string)
;
;       do id = 0,dim_DIAG-1
;           hlp   = diag_all(id,:)
;           good  = ind(hlp.ne. 999.)
;           ngood = dimsizes(good)
;           corr(pt,id) = escorc(hlp(good),xx(pt,good))
;           rc = regline(hlp(good),xx(pt,good))
;           df = rc@nptxy-2
;           z = 0.5*log((1.+abs(corr(pt,id)))/(1.-abs(corr(pt,id))))
;           corr_sig(0,id) = tanh(z-1.96/(sqrt(ngood-3.)))
;           corr_sig(1,id) = tanh(z+1.96/(sqrt(ngood-3.)))
;
;           t0    = abs(corr(pt,id))/sqrt((1-corr(pt,id)^2)/(ngood-2))
;           sign0 = 1.-betainc(df/(df+rc@tval^2), df/2.0, 0.5)
;
;           info_output(wdiag(id) + " corr: " + corr(pt,id) + "  models: " \
;           + ngood + "  stddv: " + sign0(0), verbosity, 1)
;           delete([/hlp,good,ngood,rc,df/])
;
;       end do
;
;       ;;Sort abs(coerr. coef.) to degreacing order in array
;       acorr = abs(corr(pt,:))
;       qsort(acorr)
;       acorr = acorr(::-1)
;
;       do id = 0,dim_DIAG-1
;           i_diag = ind(abs(corr(pt,:)).eq.acorr(id))
;           r_wdiag(id) = wdiag(i_diag)
;           r_sig0(id) = corr_sig(0,i_diag)
;           r_sig1(id) = corr_sig(1,i_diag)
;           diag_sign(pt,id) = where(corr(pt,i_diag).gt.0.,"red","blue")
;           delete(i_diag)
;       end do
;
;       delete([/wdiag, corr_sig/])
;       wdiag = r_wdiag
;       corr_sig = new((/2,dim_DIAG/),float)
;       corr_sig(0,:) = r_sig0
;       corr_sig(1,:) = r_sig1
;       delete([/r_wdiag, r_sig0, r_sig1, id/])
;
;       print(wdiag+": "+acorr)
;
;       acorr!0 = "diag"
;       acorr&diag = wdiag
;
;       ;;===================================================================
;       ;;Plot
;       ;;Call plot scripts
;       plot_file  = diag_script_base + "_" + var0 + "-" + calc_type + "_" + \
;                    models@experiment(DIM_OBS+DIM_HDAT) + "_" + p_time(pt) + \
;                    "-" + (/p_time(pt) + p_step/)
;       wks = gsn_open_wks(file_type, plot_dir + plot_file)
;
;       MainStg = plot_num(pt) + ") " + p_time(pt) + "-" + (/p_time(pt) + p_step/)
;       YStg = "Absolute Correlation Coefficient"
;       plot = error_bar_plot(acorr, corr_sig, wdiag, \
;       diag_sign(pt,:), MainStg, YStg, wks, "draw")
;
;       ;; add meta data to plot (for reporting)
;       climofiles = new(DIM_DAT, string)
;       do ii = 0, DIM_DAT - 1
;           climofiles(ii) = interface_get_inpaths(ii) + "/" \
;           + interface_get_infile(variables(0), field_types(0), ii)
;       end do
;
;       alltags = array_append_record(tags, (/"ST_mder","ST_spectrum", "ST_corr", "ST_range", \
;       "PT_errorbar", "PT_scatter", "DM_et", \
;       "DM_global", "DM_polar", "DM_sh", "DM_trop"/), 0)
;       allmodelnames = dataset_info@dataset
;       caption = "Absolute values of the correlation coefficient between " + variables(0) + \
;       " and pressent-day diagnostics, similar to Wenzel et al., JClim, 2016, fig. 2"
;       id = diag_script + "_" + variables(0)
;       contrib_authors = (/"A_wenz_sa"/)
;
;       ESMValMD(plot_dir+plot_file + "." + file_type, alltags, caption, id, variables(0), \
;       allmodelnames, climofiles, diag_script, contrib_authors)
;
;       delete([/alltags, caption, id, climofiles, allmodelnames, wks, plot, MainStg, acorr/])
;
;       print("gv "+plot_dir+plot_file+".ps")
;
;   end do

  leave_msg(diag_script, "")

end
