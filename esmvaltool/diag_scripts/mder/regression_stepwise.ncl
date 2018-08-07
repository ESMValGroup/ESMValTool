; #############################################################################
; mder/regression_stepwise.ncl
; #############################################################################
; MULTIPLE DIAGNOSTIC ENSEMBLE REGRESSION
; Author: Sabrina Wenzel (DLR, Germany)
; ESMVal project PI: Veronika Eyring (DLR, Germany)
; #############################################################################
;
; Description:
;     This routine computes weights for a set of models listed in the
;     recipe recipes/recipe_wenzel16jclim.xml following the multiple diagnostic
;     ensemble regression (MDER) described by Karpechko et al. 2013.
;     The script calls a variable of interest from different models (which
;     shall be weighted) and additionally several diagnostics which are
;     specified in the recipe for computig the weights.
;     Calls three plotscripts:
;     + Plots the absolute correlation coefficient between the variable of
;       interest and present-day diagnostic.
;     + Scatter plot of the climatological mean of the variable of interest
;       and the result of the MDER.
;     + Timeseries of the variable of interest including the multi-model
;       mean, the weighted and unweighted model ensemble.
;
; Required recipe attributes:
;     lat : latitude for area averaging
;     lon : longitude for area averaging
;     p_time : time range for projections
;     scale_time : a year to scale the timeseries with
;     smooth : (optional) smoothes time period with 1-2-1 filter
;     iter : (only if 'smooth' eq True) iteration for smoothing
;     area_oper : select for operation in area_operation function
;     time_oper : select for operation in time_operation function
;     time_opt : select for option in time_operation function
;     units : units of target variable
;     wdiag : array with diagnostic names
;     pvalue : prediction interval
;
; Optional recipe attributes:
;     loglin : perform log-linear regression (default: False)
;
; Caveats:
;
; Modification history
;   * 20180806_A-schl_ma: ported to v2.0
;   * 201808??_A-senf_da: added option to plot log-linear regression
;   * 201303??_A-wenz_sa: written
;
; #############################################################################

load "interface_scripts/interface.ncl"

load "./diag_scripts/shared/statistics.ncl"
load "./diag_scripts/shared/set_operators.ncl"
load "./diag_scripts/shared/latlon.ncl"

load "./diag_scripts/shared/plot/mder.ncl"
load "./diag_scripts/shared/plot/aux_plotting.ncl"

; load "./diag_scripts/shared/meta_data.ncl"

load "./diag_scripts/shared/mder.ncl"
; load "./diag_scripts/aux/SeaIce/seaice.ncl"






;; Read data
begin
  verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
  info_output("<<<<<<<< Entering " + diag_script, \
                verbosity, 2)

  i_MOD = dimsizes(models@name)
  dim_VAR = dimsizes(variables)
  field_type0 = field_types(0)
  var0 = variables(0)
  target_var = var0
  if (var0.eq."sic") then
    target_var = "SIE"  ; sea ice extent
  end if

  ;; Check field type: need 1-d field
  supported_vars = (/"sic"/)  ; List of non-1-d variables that are supported
  if (.not.isStrSubset(field_type0,"0")) then
    if (.not.any(var0.eq.supported_vars)) then
      error_msg("f", diag_script, "", "Variable " + var0 + "not supported yet.")
    end if
  end if

  ;; Check optional diag_script_info attributes
  if .not.(isatt(diag_script_info,"units")) then
    diag_script_info@units = variable_info@units
  end if

  ;;Extract observations from model-list
  i_OBS = str_match_ind_ic(models@project,"OBS")  ; Case insensitive
  dim_OBS = dimsizes(i_OBS)

  ;;Extract historical simulations from model-list
  i_hMOD_tmp = ind(models@start_year.ne.max(tofloat(models@start_year)))
  i_hMOD     = i_hMOD_tmp(ind(i_hMOD_tmp.ne.i_OBS))
  dim_hMOD   = dimsizes(i_hMOD)
  delete(i_hMOD_tmp)

  ;;Extract future simulations and unique model names from model-list
  if (getenv("ESMValTool_average_ens").eq."True") then
    idx_model = UNIQ(models@name(dim_OBS+dim_hMOD:))+dim_OBS+dim_hMOD
  else
    idx_model = ispan(dim_OBS+dim_hMOD,dimsizes(models@name)-1,1)
    imodel = UNIQ(models@name(dim_OBS:dim_OBS+dim_hMOD))
    if (dimsizes(imodel).eq.1) then
      modelnames = unique_labels_min((/"name","ensemble"/))
      delete(models@name)
      models@name = modelnames
      delete(modelnames)
    end if
    delete(imodel)
  end if
  dim_MOD = dimsizes(models@name(idx_model))
  tmpm = models@name
  tmpe = models@ensemble
  delete(models@ensemble)
  models@ensemble = tmpe(idx_model)
  obs_name        = tmpm(i_OBS)
  delete([/tmpm,tmpe,i_hMOD/])

  ;;read data
  do idx_mod = 0, i_MOD - 1
    data = read_data(idx_mod, var0, field_type0)
  end do
end

;; Load var_info
begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

;; Main
begin
  ;; File type for ploting
  file_type = getenv("ESMValTool_output_file_type")
  if (ismissing(file_type)) then
    file_type = "PS"
  end if

  ;;Global attributes
  if (isatt(diag_script_info,"lat"))then
    latrange = diag_script_info@lat
  end if
  if (isatt(diag_script_info,"lon"))then
    lonrange = diag_script_info@lon
  end if
  if (isatt(diag_script_info,"lev"))then
    level = diag_script_info@lev
  end if
  if (isatt(diag_script_info,"p_time"))then
    p_time = diag_script_info@p_time
    p_step = diag_script_info@pstep;p_time(1) - p_time(0)-1
  end if
  if (isatt(diag_script_info,"scal_time"))then
    scal_time = diag_script_info@scal_time
  end if

  smooth = False
  if (isatt(diag_script_info,"smooth"))then
    smooth = diag_script_info@smooth
  end if
  if (isatt(diag_script_info,"iter"))then
    iter = diag_script_info@iter
  end if
  if (isatt(diag_script_info,"area_oper"))then
    area_oper = diag_script_info@area_oper
  end if
  if (isatt(diag_script_info,"time_oper"))then
    time_oper = diag_script_info@time_oper
  end if
  if (isatt(diag_script_info,"time_opt"))then
    time_opt = diag_script_info@time_opt
  end if
  if (isatt(diag_script_info,"cross_validation_mode"))then
    cross_validation_mode = diag_script_info@cross_validation_mode
  end if
  if (isatt(diag_script_info,"typ"))then
    calc_type = diag_script_info@typ
  end if
  if (isatt(diag_script_info,"loglin"))then
    loglin = diag_script_info@loglin
  else
    loglin = False
  end if

  ;;Output dir
  plot_dir = getenv("ESMValTool_plot_dir")
  climo_dir = getenv("ESMValTool_climo_dir")
  diag_script_base = basename(diag_script)
  diri_out = get_output_dir(plot_dir, diag_script_base)
  plot_num = integertochar(ispan(97,122,1))

  ;;Create file
  plot_file  = models@project(dim_OBS) + "_" + target_var + \
               "-" + calc_type + "_" + models@experiment(0+dim_OBS+dim_hMOD) + \
               "_" + (p_step+1) + "ystep"

  ;;Create variable arrays
  nyMax = max(toint(models@end_year))
  nyMin = min(toint(models@start_year(dim_OBS:)))
  nMax  = nyMax - nyMin + 1
  xx    = new((/dim_MOD,dimsizes(p_time)/),float)
  xx_val = new((/dim_MOD,dimsizes(p_time),2/),float)
  yy     = new((/dim_MOD+3,nMax/),float) ;;incl. 2 placeholders uMMM and wMMM
  yyw    = new((/2,nMax/),float)
  yywt   = new((/2,nMax/),float)
  yy_obs = new((/dim_OBS,nMax/),float)
  weights_t = new((/1,dim_MOD,dimsizes(p_time)/),float)
  cfint_t   = new((/dim_MOD,dimsizes(p_time)/),float)
  coefYest  = new((/dim_MOD,dimsizes(p_time)/),float)
  Yweight   = new((/dimsizes(p_time)/),float)
  MCORR     = new((/dimsizes(p_time)/),float)
  ;;===========================================================================
  ;;MAIN
  ;;Call the variable of interest and save in one array yy
  do imod = 0,dim_MOD-1
    i = idx_model(imod)
    info_output(models@name(i),verbosity, 3)

    ;;Call projection simulations for regression analysis and time series,
    ;;scale time series and avg over 10-year period (see also cfg_mder.ncl)
    data = read_data(i, var0, field_type0)
    tmp = extract_data(i, data, -1, 0, 0)
    idxH = ind(models@name(i).eq.models@name(:dim_OBS+dim_hMOD-1))
    idxF = ind(models@name(i).eq.models@name(dim_OBS+dim_hMOD:)) + \
           dim_OBS + dim_hMOD

    ;; Calculate derived var (i.e., make it 1-d)
    if (var0.eq."sic") then  ; Calculate sea ice extent
      if (isatt(diag_script_info,"region")) then
        region = diag_script_info@region
      else
        error_msg("w", diag_script, "", "No region for sea ice extent " + \
                  "selected in cfg file, selecting Arctic.")
        region = "Arctic"
      end if
      tmp := sea_ice_area(tmp, i, data, region, "all", "extent", True)
    end if

    if (dimsizes(idxH).ge.1) then

      ;;Call historical simulations for time series
      tmph = extract_data(idxH(0), data, -1, 0, 0)
      if (var0.eq."sic") then  ; Calculate sea ice extent
        tmph := sea_ice_area(tmph, idxH(0), data, region, "all", "extent", True)
      end if

      ;;Average over ensemble members if more than one exists (see namlist)
      if (getenv("ESMValTool_average_ens").eq."True") then
        if (dimsizes(idxF).gt.1) then
          if (dimsizes(idxF).eq.dimsizes(idxH)) then
            tmpAF = new((/dimsizes(idxF),dimsizes(tmp&time)/),typeof(tmp))
            copy_VarCoords(tmp,tmpAF(0,:))
            tmpAF(0,:) = tmp
            delete(tmp)

            tmpAH = new((/dimsizes(idxH),dimsizes(tmph&time)/),typeof(tmph))
            copy_VarCoords(tmph,tmpAH(0,:))
            tmpAH(0,:) = tmph
            delete(tmph)

            ;;Loop over ensemble members, read the data and save in arry
            do iens = 1,dimsizes(idxF)-1
              tmpF = extract_data(idxF(iens), data, -1, 0, 0)
              if (var0.eq."sic") then  ; Calculate sea ice extent
                tmpF := sea_ice_area(tmpF, idxF(iens), data, region, "all", \
                                     "extent", True)
              end if
              tmpAF(iens,:) = (/tmpF/)

              tmpH = extract_data(idxH(iens), data, -1, 0, 0)
              if (var0.eq."sic") then  ; Calculate sea ice extent
                tmpH := sea_ice_area(tmpH, idxH(iens), data, region, "all", \
                                     "extent", True)
              end if
              tmpAH(iens,:) = (/tmpH/)
              delete([/tmpF,tmpH/])
            end do

            ;;Average over ensemble members
            tmp = tofloat(dim_avg_n_Wrap(tmpAF, 0))
            tmph = tofloat(dim_avg_n_Wrap(tmpAH,0))
            copy_VarCoords(tmpAF(0,:),tmp)
            copy_VarCoords(tmpAH(0,:),tmph)
            delete([/tmpAF,tmpAH/])
          else
            print("fatal: inconsistent ensemble member size for model " + \
                  models@name(i) + " in " + diag_script)
            status_exit(1)
          end if
        end if
      end if

      ;;Append furure array to historical array
      tmpa = array_append_record(tmph,tmp,0)
      tmpy = tmpa
      delete([/tmp,tmpa,tmph,idxH,idxF/])
    end if

    ;;Calculate average over starting period for trend calculation
    tmpxx = time_operations(tmpy, scal_time(0), scal_time(1), time_oper, \
                            time_opt,True)
    tmpyy = time_operations(tmpxx, scal_time(0), scal_time(1), "average", \
                            "yearly",True)
    tmppd = tofloat(avg(tmpyy))

    ;;Loop over different target periods
    do pt = 0,dimsizes(p_time)-1

      ;;Average over var time series
      tmps  = time_operations(tmpy, p_time(pt), p_time(pt)+p_step, time_oper, \
                              time_opt,True)
      tmpj  = time_operations(tmps, p_time(pt), p_time(pt)+p_step, "average", \
                              "yearly",True)
      if (calc_type.eq."trend") then
        xx(imod,pt) = tofloat(avg(tmpj)-tmppd)
        xx_val(imod,pt,0) = tofloat(avg(tmpj)-tmppd)
        xx_val(imod,pt,1) = tofloat(stddev(tmpj)/sqrt(dimsizes(tmpj))-tmppd)
        variable_info@units = diag_script_info@units + "/dec"
      else if (calc_type.eq."pos") then
        xx(imod,pt) = tofloat(avg(tmpj))
        xx_val(imod,pt,0) = tofloat(avg(tmpj))
        xx_val(imod,pt,1) = tofloat(stddev(tmpj)/sqrt(dimsizes(tmpj)))
        variable_info@units = diag_script_info@units
      end if
      end if
      delete([/tmps,tmpj/])
    end do

    ;;Extract season for var time series
    tmpst = time_operations(tmpy, -1, -1, time_oper, time_opt, True)
    tmpsc = time_operations(tmpst, -1, -1, "average", "yearly", True)

    ;;Smoothing var time series with 1-2-1 filter
    if (smooth.eq.True) then
      tmpsm = filter121(tmpsc,iter)
    else
      tmpsm = tmpsc
    end if
    if (calc_type.eq."trend") then
      yy(imod+3,0:dimsizes(tmpsm)-1) = tofloat(tmpsm-tmppd)
    else if (calc_type.eq."pos") then
      yy(imod+3,0:dimsizes(tmpsm)-1) = tofloat(tmpsm)
    end if
    end if

    delete([/tmppd,tmpyy,tmpxx/])
    delete([/tmpy,tmpsm,tmpsc,tmpst/])
    ;delete([/aveX,varX,sX,iflag/])
  end do

  yy!1 = "time"
  yy&time = ispan(nyMin,nyMax,1);tmpsc&year
  yy&time@calendar = "standard"

  ;;Save unique model names
  tmpm = models@name
  delete(models@name)
  models@name  = tmpm(idx_model)
  delete(tmpm)

  ;;Call OBS of the variable of interest and save in one array yy_obs
  do iobs = 0, dim_OBS-1
    info_output(obs_name(iobs),verbosity, 3)

    ;;Call projection simulations for regression analysis and time series,
    ;;scale time series and avg over 10-year period (see also cfg_mder.ncl)
    tmpy = extract_data(iobs, data, -1, 0, 0)

    ;; Calculate derived var (i.e., make it 1-d)
    if (var0.eq."sic") then  ; Calculate sea ice extent
      SIE = sea_ice_area(tmpy, iobs, data, region, "all", "extent", True)
      tmpy := SIE
      delete(SIE)
    end if

    ;;Scale var time series
    tmpoy = time_operations(tmpy,scal_time(0),-1,"average","yearly",True)
    tmppd = tofloat(avg(tmpoy))
    tmpyy = time_operations(tmpy,-1,-1,time_oper,time_opt,True)
    tmpsm = time_operations(tmpyy,-1,-1,"average","yearly",True)
    ymin = min(toint(models@start_year(iobs)))-nyMin;-1

    if (calc_type.eq."trend") then
      yy_obs(iobs,ymin:ymin+dimsizes(tmpsm)-1) = tofloat(tmpsm-tmppd)
    else if (calc_type.eq."pos") then
      yy_obs(iobs,ymin:ymin+dimsizes(tmpsm)-1) = tofloat(tmpsm)
    end if
    end if

    ;yy_obs(iobs,:dimsizes(tmpsm)-1) = tofloat(tmpsm)

    delete([/tmpy,tmpyy,tmpsm,ymin,tmppd,tmpoy/])
  end do
  yy_obs!1 = "time"
  yy_obs&time = yy&time
  yy_obs&time@calendar = "standard"
  yy_obs@name = obs_name
  delete(data)

  ;;Call diagnostics listed in diag_script_info@fname_diag
  ;; Read diagnostic list and add files together
  xml_name_char = tochar(getenv("ESMValTool_xml_name"))
  temp_dir = get_output_dir(getenv("ESMValTool_wrk_dir"), "temp") + \
             tostring(xml_name_char(:dimsizes(xml_name_char) - 6)) + \
             "_diagnostic.nc"

  script_DIAG = getenv("ESMValTool_diag_script")
  dim_SCRIPT  = dimsizes(script_DIAG)

  if (isfilepresent(temp_dir)) then

    ;; Read diagnostic list
    temp_file = addfile(temp_dir, "r")
    temp_list = temp_file->temp_list
    file_list = tostring(temp_list)
    file_list@_FillValue = default_fillvalue(typeof(file_list))

    if (dim_SCRIPT.eq.1.or.script_DIAG(dim_SCRIPT).eq.diag_script) then
      delete(temp_list)
    end if

    ;; Remove duplicates
    do ii = 0, dimsizes(file_list) - 1
      curr = file_list(ii)
      idx = ind(curr.eq.file_list)
      if (dimsizes(idx).gt.1) then
        file_list(idx) = file_list@_FillValue
        file_list(ii) = curr
      end if
      delete([/curr,idx/])
    end do
    tmp = file_list(ind(.not.ismissing(file_list)))
    delete(file_list)
    file_list = tmp
    delete(tmp)

    ;;Create diagnostic arrays
    dim_DIAG = dimsizes(file_list)
    diag_a   = new((/dim_DIAG,dim_MOD/),float)
    diag_a_std = new((/dim_DIAG,dim_MOD/),float)
    wdiag    = new((/dim_DIAG/),string)           ;;array with diagnostic names
    obs0     = new((/dim_DIAG/),float)
    obs0_std = new((/dim_DIAG/),float)

    ;; Loop over files in the list, append data
    do ii = 0, dim_DIAG - 1
      data_temp = ncdf_read(file_list(ii), "diag")

      do imod = 0,dim_MOD-1
        do dmod = 0,dimsizes(data_temp&models)-1
          if (models@name(imod).eq.data_temp&models(dmod)) then
            diag_a(ii,imod) = data_temp(0,dmod,0)     ;;collect diagnostics in one array
            diag_a_std(ii,imod) = data_temp(0,dmod,1) ;;collect diagnostics in one array
            if (isatt(data_temp,"wdiag")) then
              wdiag(ii) = data_temp@wdiag
            end if
          end if
        end do
      end do
      if (isatt(data_temp,"ref_model"))then
        if (any(data_temp&models.eq.data_temp@ref_model)) then
            idx_ref = ind(data_temp&models.eq.data_temp@ref_model)
            obs0(ii) = data_temp(0,idx_ref,0)         ;;collect obs in seperate array
            obs0_std(ii) = data_temp(0,idx_ref,1)     ;;collect obs in seperate array
            delete(idx_ref)
        end if
      end if
      delete(data_temp)
      info_output(wdiag(ii)+": OBS +/- 95% conv.: "+obs0(ii)+" +/- "+(cdft_t(0.05/2., 26)*obs0_std(ii))+\
                  " MODELS: "+(avg(diag_a(ii,:)))+" +/- "+(stddev(diag_a(ii,:))),verbosity,1)

      ;;Only account for good diagnostics where no model is missing!
      if (.not.any(ismissing(diag_a(ii,:)))) then
          if (isvar("good_diag_ind").eq.False) then
              good_diag_ind = ii
          else
              tmp_gd = good_diag_ind
              delete(good_diag_ind)
              good_diag_ind = new((/dimsizes(tmp_gd)+1/),integer)
              good_diag_ind(:dimsizes(tmp_gd)-1) = tmp_gd
              good_diag_ind(dimsizes(tmp_gd)) = ii
              delete(tmp_gd)
          end if
      else
          misMOD = ind(ismissing(diag_a(ii,:)))
          error_msg("w", diag_script, "", \
                    "MDER was not succesfull!!! Diagnostic input for " +\
                    " some models are missing:")
          print(models@name(misMOD))
          error_msg("f", diag_script, "", \
                    "Delete files in work/temp/* and try again!")
      end if
    end do  ; ii

    tmp_wdiag = wdiag(good_diag_ind)
    delete(wdiag)
    diag_all = diag_a(good_diag_ind,:)
    diag_all_std = diag_a_std(good_diag_ind,:)
    wdiag = tmp_wdiag
    delete(tmp_wdiag)
    diag_all&diagnostics = wdiag;(good_diag_ind)
    diag_all_std&diagnostics = wdiag;(good_diag_ind)
    dim_DIAG = dimsizes(good_diag_ind)
    dim_DIAG_all = dimsizes(diag_a(:,0))
    delete(good_diag_ind)
  else
    print("fatal: cannot find diagnostic list in " + diag_script)
    status_exit(1)
;      print(diag_a)
  end if

  ;;---------------------------------------------------------------------------
  ;;Loop over different target periods
  do pt = 0,dimsizes(p_time)-1
    if (any(ismissing(xx(:,pt)))) then
      continue
    end if

    ;;Regression model
    ;;The cross-validation cycle
    if (cross_validation_mode.eq.True) then
      do im = 0,dim_MOD-1
        good   = ind(models@name .ne. models@name(im))
        coef   = regress_stepwise(diag_all(:,good),xx(good,pt))
        const  = coef@const
        terms  = toint(coef@terms)
        nterms = coef@nterms
        YFIT   = coef@YFIT
        mcorr  = coef@MCORRELATION
        corr   = coef@correlation

        info_output("Cross-validation mode. The model: "+models@name(im),verbosity, 2)
        info_output("Constant: "+const,verbosity, 2)
        if (verbosity.ge.2) then
          if (nterms .gt. 0) then
            do it = 0, dimsizes(coef)-1
              print("info:  "+wdiag(terms(it))+" , "+coef(it)+" , "+corr(it))
            end do
          end if
        end if
        info_output("Multiple correlation: "+mcorr,verbosity,2)
        info_output("",verbosity, 2)

        delete([/const,terms,nterms,YFIT,mcorr,corr,coef,good/])
      end do
    end if

    ;;Building the final regression
    if (all(isvar((/"coef","terms"/))).eq.True) then
      delete([/coef,terms/])
    end if
    coef = regress_stepwise(diag_all,xx(:,pt))

    const  = coef@const
    terms  = toint(coef@terms)
    nterms = coef@nterms
    YFIT   = coef@YFIT
    mcorr  = coef@MCORRELATION
    corr   = coef@correlation

    info_output("The final model ("+scal_time(0)+"-"+(p_time(pt)+p_step)+"):",verbosity, 0)
    info_output("Constant:"+const,verbosity, 0)
    if (verbosity.ge.0) then
      if (nterms.gt.0) then
        do it = 0,dimsizes(coef)-1
          print("info:   "+wdiag(terms(it))+" , "+coef(it)+" , "+corr(it))
        end do
      end if
    end if
    info_output("Multiple correlation: "+mcorr,verbosity, 0)
    info_output("",verbosity, 0)

    ;;Model weighting
    x_diag = transpose(diag_all(terms,:))
    good_diag = ind(.not.ismissing(x_diag(:,0)))
    x = x_diag(good_diag,:)

    ;tmp_MOD = dim_MOD
    tmp_MOD = dimsizes(x(:,0))

    ;;Create Design Matrix D
    if (isvar("D").eq.True) then
      delete(D)
    end if
    D = new((/tmp_MOD,dimsizes(terms)+1/),typeof(diag_all))
    D = 1
    do id = 0,dimsizes(x(:,0))-1
      D(id,1:) = x(id,:)
    end do

    xobs = new((/1,dimsizes(terms)+1/),typeof(diag_all),getFillValue(diag_all))
    xobs = 1.
    xobs(0,1:) = obs0(terms)

    df = tmp_MOD-dimsizes(xobs(0,:))

    I = new((/tmp_MOD,tmp_MOD/),typeof(diag_all))
    I = 0.
    do imod = 0,tmp_MOD-1
      I(imod,imod) = 1.
    end do

    y = xx(good_diag,pt)

    m  = new((/tmp_MOD,1/),float)  ;;column vector ;REPLICATE(1.,nmods)
    m  = 1.
    n  = new((/tmp_MOD,1/),float)  ;;column vector ;REPLICATE(1./nmods,nmods)
    n  = (/1./tmp_MOD/)


    weights = (transpose(n)+(xobs(:,1:)-transpose(n)#x)#inverse_matrix(transpose(x)#x-transpose(x)#m#transpose(n)#x)#(transpose(x)-transpose(x)#m#transpose(n)))
print(weights)
    cfint  = cdft_t(0.05/2., df)^2*(1.+xobs#inverse_matrix(transpose(D)#D)#transpose(xobs))
    SSE1   = y#(I-D#inverse_matrix(transpose(D)#D)#transpose(D))
    sigmae = (1./df)*SSE1
    tmp_cf = cfint(0,0)*sigmae
    delete(cfint)

    ;;Weighted mean and std
    y_weight = (weights#y)
    cf = sqrt(tmp_cf#y)

    uM = avg(y)
    uM_std = stddev(y)/sqrt(dim_MOD)

    info_output("uMMM     = "+uM+" +/- "+uM_std,verbosity,0)
    info_output("y_weight = "+y_weight,verbosity,0)
    info_output("cf       = "+cf,verbosity,0)
    info_output(" ",verbosity,1)

    do imod = 0,dim_MOD-1
      do gdiad = 0,dimsizes(good_diag)-1
        if (imod.eq.good_diag(gdiad))
          weights_t(:,imod,pt) = weights(:,gdiad)
          cfint_t(imod,pt)     = tmp_cf(gdiad)
        end if
      end do
    end do

    info_output("Model weights: ",verbosity,0)
    do imod = 0,dim_MOD-1
      info_output(weights_t(:,imod,pt)+" : "+models@name(imod),verbosity,0)
    end do
    info_output(" ",verbosity,0)

    ;;Save Y-estimates and multiple correllation coeficients
    coefYest(:,pt) = coef@Yest
    MCORR(pt) = mcorr
    Yweight(pt) = y_weight

    ;;Save coefficients and terms of regression model in array
    if (isvar("mderrm")) then
      idx = max((/dimsizes(coef),dimsizes(mderrm(0,1:))/))
      tmprm = mderrm
      tmptr = termsrm
      tmpt  = termst
      delete([/mderrm,termsrm,termst/])

      mderrm  = new((/pt+1,idx+1/),typeof(tmprm))
      termsrm = new((/pt+1,idx+1/),typeof(tmptr))
      termst  = new((/pt+1,idx+1/),typeof(tmpt))

      mderrm(:pt-1,:dimsizes(tmprm(0,:))-1) = tmprm
      termsrm(:pt-1,:dimsizes(tmptr(0,:))-1) = tmptr
      termst(:pt-1,:dimsizes(tmpt(0,:))-1) = tmpt
      mderrm(pt,0) = const
      mderrm(pt,1:dimsizes(coef)) = coef
      termsrm(pt,1:dimsizes(wdiag(terms))) = wdiag(terms)
      termst(pt,1:dimsizes(terms)) = terms
      delete([/tmprm,tmptr,tmpt,idx/])
    else
      mderrm  = new((/pt+1,dimsizes(coef)+1/),float)
      termsrm = new((/pt+1,dimsizes(coef)+1/),string)
      termst  = new((/pt+1,dimsizes(coef)+1/),typeof(terms))

      mderrm(pt,0)  = const
      mderrm(pt,1:) = coef
      termsrm(pt,1:) = wdiag(terms)
      termst(pt,1:) = terms
    end if

    delete([/y,corr,mcorr/])
    delete([/x,x_diag,good_diag/])
    delete([/xobs,cf,tmp_cf,weights,uM,uM_std/]);,cfint,t_diff,tmp_t/])
  end do  ;;END Loop over different target periods

  ;;Calculate diverse multy model means
  tmpm = models@name
  delete(models@name)
  models@name = new(dim_MOD+3,string)
  models@name(0)  = "uMMM"
  models@name(1)  = "MDER time dependent"
  models@name(2)  = "MDER"
  models@name(3:) = tmpm
  delete(tmpm)

  do pt = 0,dimsizes(p_time)-1
    ;ptt = 0
    tmp_yy  = new(dimsizes(yy),typeof(yy))
    tmp_yywt = new(dimsizes(yywt),typeof(yywt))
    do tt = 0,nMax-2
      yy(0,tt) = avg(yy(3:,tt))
      if (tt.ge.p_time(pt)-nyMin.and.tt.lt.(p_time(pt)+p_step+1-nyMin)) then
          tmp_yy(1,tt)  = sum(yy(3:,tt)*weights_t(0,:,pt))
          tmp_yywt(0,tt) = tmp_yy(1,tt)+sqrt(abs(sum(yy(3:,tt)*cfint_t(:,pt))))
          tmp_yywt(1,tt) = tmp_yy(1,tt)-sqrt(abs(sum(yy(3:,tt)*cfint_t(:,pt))))
      end if
      if (tt.eq.(p_time(pt)+p_step-nyMin)) then
          yy(1,tt-(p_step+1)/2)  = avg(tmp_yy(1,tt-p_step:tt))
          yywt(0,tt-(p_step+1)/2) = avg(tmp_yywt(0,tt-p_step:tt))
          yywt(1,tt-(p_step+1)/2) = avg(tmp_yywt(1,tt-p_step:tt))
      end if
      yy(2,tt) = sum(yy(3:,tt)*weights_t(0,:,pt))
      yyw(0,tt) = yy(2,tt)+sqrt(abs(sum(yy(3:,tt)*cfint_t(:,pt))))
      yyw(1,tt) = yy(2,tt)-sqrt(abs(sum(yy(3:,tt)*cfint_t(:,pt))))
    end do
  end do

  tmpm = models@name
  delete(models@name)
  models@name = tmpm(3:)
  delete(tmpm)

  ;;---------------------------------------------------------------------------
  ;;This is for plotting the regression line
;  nfit     = fspan(0,49,50)
  diag_reg = new((/dimsizes(p_time),dim_MOD+2/),float)

  nfit     = fspan(0, dimsizes(diag_reg(0,2:))-1,dimsizes(diag_reg(0,2:)))

  diag_std = new((/dimsizes(p_time),dim_MOD+2/),float)
  wy       = new((/dimsizes(p_time),dim_MOD+2/),float)
  wy_std   = new((/dimsizes(p_time),dim_MOD+2/),float)
  yfit_reg = new((/dimsizes(p_time),dimsizes(nfit)/),float)
  yfit_std = new((/dimsizes(p_time),2,dimsizes(nfit)/),float)

  do ppt = 0, dimsizes(p_time)-1

    xgood    = ind(.not.ismissing(termst(ppt,:)))
    hlp_tmp  = new((/dimsizes(mderrm(ppt,xgood)),dimsizes(nfit)/),float)
    hlp      = new((/1,dimsizes(mderrm(ppt,xgood))+1/),float)
    diag_reg(ppt,:) = mderrm(ppt,0);const
    diag_std(ppt,:) = 0.
    yfit_reg(ppt,:) = mderrm(ppt,0);const
    hlp             = 1.

    wy(ppt,0) = avg(xx(:,ppt))
    wy(ppt,1) = Yweight(ppt)
    wy(ppt,2:) = xx(:,ppt)
    wy_std(ppt,2:) = xx_val(:,ppt,1);cfint_t(:,ppt)

    do ii = 1,dimsizes(mderrm(ppt,:))-1
      if (.not.ismissing(mderrm(ppt,ii))) then
        diag_reg(ppt,1)  = diag_reg(ppt,1) + mderrm(ppt,ii) * \
                                             obs0(termst(ppt,ii))
        diag_std(ppt,1)  = diag_std(ppt,1) + (mderrm(ppt,ii) * \
                                             (obs0_std(termst(ppt,ii))))^2
        diag_reg(ppt,2:) = coefYest(:,ppt)
                           ;diag_reg(2:)+coef(ii)*diag_all(terms(ii),:)
        diag_std(ppt,2:) = diag_std(ppt,2:) + (mderrm(ppt,ii) * \
                                          (diag_all_std(termst(ppt,ii),:)))^2

        minx = min(diag_all(termst(ppt,ii),:)) - 1.1 * \
              (max(diag_all(termst(ppt,ii),:)) - \
               min(diag_all(termst(ppt,ii),:)))
        maxx = max(diag_all(termst(ppt,ii),:)) + 1.1 * \
              (max(diag_all(termst(ppt,ii),:)) - \
               min(diag_all(termst(ppt,ii),:)))
        diag_var = minx+(maxx-minx)*nfit/dimsizes(nfit)

        hlp_tmp(ii-1,:) = diag_var
        delete([/diag_var,minx,maxx/])
      end if
    end do

    diag_std = sqrt(diag_std)

    x_diag = transpose(diag_all(termst(ppt,xgood),:))
    good_diag = ind(.not.ismissing(x_diag(:,0)))
    x = x_diag(good_diag,:)

    ;;Create Design Matrix D and Einheits Matrix I for pt = 1
    if (isvar("D").eq.True) then
      delete(D)
    end if
    D = new((/dimsizes(x(:,0)),dimsizes(termst(ppt,xgood))+1/),typeof(diag_all))
    D = 1
    do id = 0,dimsizes(x(:,0))-1
      D(id,1:) = x(id,:)
    end do

    df = tmp_MOD-dimsizes(obs0(termst(ppt,xgood)))+1

    I = new((/dimsizes(x(:,0)),dimsizes(x(:,0))/),typeof(diag_all))
    I = 0.
    do imod = 0,tmp_MOD-1
      I(imod,imod) = 1.
    end do

    SSE=xx(:,ppt)#(I-D#inverse_matrix(transpose(D)#D)#transpose(D))#transpose(xx(:,ppt))
    sigmae2 = sqrt((1./df)*SSE)

    ;;Calculate 95% confidence interval bounds
    if (.not.loglin) then
      miny = min(xx(:,ppt)) - 1.1 * (max(xx(:,ppt)) - min(xx(:,ppt)))
      maxy = max(xx(:,ppt)) + 1.1 * (max(xx(:,ppt)) - min(xx(:,ppt)))
      yfit_reg(ppt,:) = miny+(maxy-miny)*nfit/dimsizes(nfit)
      delete([/miny, maxy/])
      do ti = 0,dimsizes(nfit)-1
        hlp(0,1:) = (/hlp_tmp(:,ti)/)
        yfit_std(ppt,0,ti) = yfit_reg(ppt,ti) + cdft_t(0.05/2., df) * sigmae2 * \
          sqrt((1.+hlp#inverse_matrix(transpose(D)#D)#transpose(hlp)))
        yfit_std(ppt,1,ti) = yfit_reg(ppt,ti) - cdft_t(0.05/2., df) * sigmae2 * \
          sqrt((1.+hlp#inverse_matrix(transpose(D)#D)#transpose(hlp)))
      end do
    else  ; log-linear regression
      isort_tmp = dim_pqsort_n(diag_reg(ppt,2:), 1, 0)
      isort     = new(dimsizes(isort_tmp)+2, integer)
      isort(0)  = 0
      isort(1)  = 1
      isort(2:) = isort_tmp + 2
      diag_reg(ppt,:) = diag_reg(ppt,isort)
      diag_std(ppt,:) = diag_std(ppt,isort)
      wy(ppt,:)       = wy(ppt,isort)
      wy_std(ppt,:)   = wy_std(ppt,isort)
      model_numbers   = ispan( 0, dimsizes(diag_reg(ppt,:))-1, 1)
      diag_reg@model_numbers = model_numbers(isort)  ; memorize order of models
      rc = regline_stats(diag_reg(ppt,2:), log(wy(ppt,2:)))
      yfit_reg(ppt,:)   = exp(rc@b(0) + rc@b(1) * diag_reg(ppt,2:))
      yfit_std(ppt,:,:) = (/ exp( rc@y95(0) + rc@b95(0) * diag_reg(ppt,2:) ), \
                             exp( rc@y95(1) + rc@b95(1) * diag_reg(ppt,2:) ) /)
      wy(ppt,1)         = exp(rc@b(0) + rc@b(1) * diag_reg(ppt,1))  ; yfit(obs)
      delete([/isort_tmp, isort, model_numbers, rc/])
    end if

    ;; Clean up
    delete([/xgood,hlp_tmp,hlp,x_diag,good_diag,x,D,df,SSE,sigmae2/])

  end do  ;ppt

  ;;---------------------------------------------------------------------------
  ;;This is for plotting the squared differences between ensemble mean
  ;;and the total change in pseudoreality
  grade_all = new((/dim_DIAG_all,dim_MOD-1,dim_MOD/),float)
              ;;1. #grades, 2. #pseudorealities
  xx_pre    = new((/dim_MOD,dimsizes(p_time)/),float)

  grade_trans  = dim_avg_n_Wrap(grade_all(:,:,:),0)
  do im = 0,dim_MOD-1
    good = ind(models@name(im).ne.models@name)
    do pt = 0, dimsizes(p_time)-1
      xx_pre(im,pt)   = avg(xx(good,pt))
    end do
    delete(good)
  end do

  ;Calculating statistics
  x25  = round(.25*dim_MOD,3)-1 ;;-1 to account for NCL indexing starting at 0
  x75  = round(.75*dim_MOD,3)-1

  presq     = new((/dimsizes(p_time),dim_MOD/),float)
  premdersq = new((/dimsizes(p_time),dim_MOD/),float)
  BSS = new((/dimsizes(p_time)/),float)

  do pt = 0, dimsizes(p_time)-1
    presq(pt,:)     = (xx(:,pt)-xx_pre(:,pt))^2
    premdersq(pt,:) = (xx(:,pt)-coefYest(:,pt))^2

    qsort(presq(pt,:))
    qsort(premdersq(pt,:))

    tmpRDS = (/(/presq(pt,x25),presq(pt,x25),dim_median(presq(pt,:)),presq(pt,x75),presq(pt,x75)/),\
                       (/premdersq(pt,x25),premdersq(pt,x25),dim_median(premdersq(pt,:))\
                        ,premdersq(pt,x75),premdersq(pt,x75)/)/)
    tmpAVG = (/avg(presq(pt,:)),avg(premdersq(pt,:))/)
    tmpPRE = (/presq(pt,:),premdersq(pt,:)/)

    if (all(isvar((/"preRDS","preAVG","prePRE"/)))) then
      aarRDS = preRDS
      aarAVG = preAVG
      aarPRE = prePRE
      delete([/preRDS,preAVG,prePRE/])
      preRDS = array_append_record(aarRDS,sqrt(tmpRDS),0)
      preAVG = array_append_record(aarAVG,sqrt(tmpAVG),0)
      prePRE = array_append_record(aarPRE,sqrt(tmpPRE),0)
      delete([/aarRDS,aarAVG,aarPRE/])
    else
      preRDS = sqrt(tmpRDS)
      preAVG = sqrt(tmpAVG)
      prePRE = sqrt(tmpPRE)
    end if
    delete([/tmpRDS,tmpAVG,tmpPRE/])

    info_output("RMSE ("+p_time(pt)+"-"+(p_time(pt)+p_step)+"), uMMM: "+sqrt(avg(presq(pt,:)))+", MDER: "+sqrt(avg(premdersq(pt,:))),verbosity, 0)
  end do

  tmpm = models@name
  delete(models@name)
  models@name = new(dim_MOD+3,string)
  models@name(0)  = "uMMM"
  models@name(1)  = "MDER time dependent"
  models@name(2)  = "MDER"
  models@name(3:) = tmpm
  delete(tmpm)
  dim_MOD = dim_MOD+2

  ;;===========================================================================
  ;;Write regression model equations to ASCII table
  RegStg = new((/dimsizes(p_time)/),string)
  PStg = new((/dimsizes(p_time)/),string)
  do pt = 0, dimsizes(p_time)-1
    RegStg(pt) = sprintf("%4.2f",mderrm(pt,0))
    do ii = 1,dimsizes(mderrm(pt,:))-1
      if (.not.ismissing(mderrm(pt,ii))) then
        if (mderrm(pt,ii).gt.0.) then
          RegStg(pt) = RegStg(pt)+" +"+sprintf("%4.2f",mderrm(pt,ii))+" x "+termsrm(pt,ii)
        else
          RegStg(pt) = RegStg(pt)+" "+sprintf("%4.2f",mderrm(pt,ii))+" x "+termsrm(pt,ii)
        end if
      end if
    end do
    p1_time = p_time(pt)+p_step-2000;p_time(0)
    PStg(pt) = p_time(pt)+"-"+sprinti("%0.2hi",p1_time)
  end do

  header = (/"Target Period             Model Equation  R^2",\
             "                          "+models@experiment(0+dim_OBS+dim_hMOD)+"         "/)
  hlist = [/header/]
  alist = [/PStg,RegStg,MCORR/]

  write_table(diri_out+plot_file+".txt", "w", hlist, "%s")
  write_table(diri_out+plot_file+".txt", "a", alist, "%s %s %f")
  delete([/header,hlist,alist/])

  asciiwrite(climo_dir + "/" + target_var + "_" + \
             models@experiment(0+dim_OBS+dim_hMOD) + "_" +(p_step+1) + \
             "ystep_diagnostics.txt", termsrm)
  ;;--------------------------------------------------------------------------
  ;;Write weights to NetCDF
  model_list = models@name(3:)
  do pt = 0, dimsizes(p_time)-1
    weights_out = weights_t(0,:,pt)
    weights_filename = diri_out + plot_file + "_weights_" + \
                       p_time(pt) + "-" + (p_time(pt)+p_step)
    weights_out!0           = "model"
    weights_out&model       = model_list
    weights_out@diag_script = diag_script
    weights_out@var         = "weight"
    weights_out@description = "Weights from MDER"
    if (target_var.eq."SIE") then
      weights_out@description = "Weights for September Arctic sea ice extent"
    end if
    weights_file            = ncdf_write(weights_out, weights_filename+".nc")
    delete([/weights_out, weights_filename, weights_file/])
  end do
  ;;--------------------------------------------------------------------------
  ;;PLOT

  ;; add meta data to plot (for reporting)
  climofiles = new(dim_MOD, string)
  do ii = 0, dim_MOD - 1
      climofiles(ii) = interface_get_inpaths(ii) + "/" \
                       + interface_get_infile(variables(0), field_types(0), ii)
  end do
  allmodelnames = models@name
  zz = diag_script + "_" + variables(0)
  contrib_authors = (/"A_wenz_sa"/)

  ;; Line plot showing the time series including MDER weights
  wks     = get_wks(0, diag_script, plot_file + "_FIG1")
  MainStg = ""
  YStg    = target_var+" ["+variable_info@units+"]"
  plot    = weight_lin_plot(yy,yywt,yyw,yy_obs,MainStg,YStg,wks,"draw")

  alltags = array_append_record(tags, (/"ST_mder", "ST_corr", "ST_range", \
                                        "ST_stddev", "PT_times", "DM_sh"/), 0)
  caption = "Time series of " + variables(0) + " for " + models@ensemble(0) + \
            " scenario, similar to Wenzel et al., JClim, 2016, fig. 5"

  plot_fname = wks@fullname
  ESMValMD(plot_fname, alltags, caption, zz, variables(0), \
           allmodelnames, climofiles, diag_script, contrib_authors)
  delete([/alltags, caption, wks, plot_fname/])

  do ppt=0, dimsizes(p_time)-1
    wks     = get_wks(0, diag_script, plot_file + "_FIG2" + plot_num(ppt))
    MainStg = plot_num(ppt)+") "+p_time(ppt)+"-"+(p_time(ppt)+p_step)
    XStg    = RegStg(ppt)+" ["+variable_info@units+"]"
    YStg    = target_var+" ["+variable_info@units+"]"
print(diag_reg(ppt,:))
print(diag_std(ppt,:))
print(wy(ppt,:))
print(wy_std(ppt,:))
print(yfit_reg(ppt,:))
print(yfit_std(ppt,:,:))
    plot2   = mder_scatter_plot(diag_reg(ppt,:), diag_std(ppt,:), \
                                wy(ppt,:), wy_std(ppt,:), \
                                yfit_reg(ppt,:), yfit_std(ppt,:,:), \
                                MainStg, XStg, YStg, wks, "draw")
    alltags = array_append_record(tags, (/"ST_mder", "ST_corr", "ST_range", \
                                          "ST_stddev", "ST_corr", "PT_scatter",\
                                          "T_EC", "DM_sh"/), 0)
    caption = "Scatterplot showing the correlation between future " + \
              variables(0) + " and the quantity resulting from MDER analysis" +\
              " on the x-axes, similar to Wenzel et al., JClim, 2016, fig. 3"
    plot_fname = wks@fullname
    ESMValMD(plot_fname, alltags, caption, zz, variables(0), \
             allmodelnames, climofiles, diag_script, contrib_authors)
    delete([/alltags, caption, wks, plot_fname/])
  end do  ; ppt
exit

  YStg   = "RMSE projection error ["+variable_info@units+"]"
  XBname = new(dimsizes(PStg)*2,string)

  if (dimsizes(XBname).eq.2) then
      XBname = (/"uMMM","MDER"/)
  else if (dimsizes(XBname).eq.4) then
      XBname(0::2) = (/"uMMM~C~"+PStg/)
      XBname(1::2) = (/"MDER~C~"+PStg/)
  else
      XBname(0::2) = (/" "/)
      XBname(1::2) = PStg
  end if
  end if

  ;; Boxplot of RMSE of different weigting methods inkl pseudoreality check
  wks     = get_wks(0, diag_script, plot_file + "_FIG3")
  MainStg = ""
  Xtitle  = PStg(0)
  do pt = 1, dimsizes(p_time)-1
    Xtitle = Xtitle+"  "+PStg(pt)
  end do
  plot3   = sqrt_error_plot(preRDS,preAVG,prePRE,\
                           MainStg,XBname,Xtitle,YStg,wks,"draw")

  alltags = array_append_record(tags, (/"ST_mder", "ST_rmsd", "ST_stddev", \
                                "ST_mean", "PT_bar", "PT_errorbar", "DM_sh"/),0)
  caption = "RMSE difference between the ensemble mean future climatological" +\
            " mean " + variables(0) + "and its equivalent in pseudoreality " +\
            "under " + models@ensemble(0) + " scenario, similar to Wenzel " +\
            "et al., JClim, 2016, fig. 4"
  plot_fname = wks@fullname
  ESMValMD(plot_fname, alltags, caption, zz, variables(0), \
           allmodelnames, climofiles, diag_script, contrib_authors)
  delete([/alltags, caption, zz, climofiles, allmodelnames, plot_fname/])

  ;; Remove temporary diagnostic list
  system("rm -f " + temp_dir)

end
