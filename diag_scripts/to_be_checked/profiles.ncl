;*********************************************************************
; profiles.ncl
;*********************************************************************
; Plot vertical profiles from model and raobs data.
; Written by Mark Stevens Feb. 2001
; last update: 07/17/01
;*********************************************************************
; The original NCL script, profiles.ncl, was written to make...
; vertical profile plot at 56 raobs stations around the globe (set 12)
; for the AMWG diagnostics package.
; http://www.cgd.ucar.edu/cms/diagnostics/
;
; The NCL script described above has been modified for use with the CCMVal
; diagnostic software package.
; (PIs: Andrew Gettelman & Veronika Eyring, co-I: Hisako Shiona)
;
; This script takes input from the var_att.ncl files
; info attributes required:
; profile_stnid         filepath name to station_id file
; profile_raobs         filepath name to raobs file
; profile_allstn        (optional) set True for plotting all 56 stations
;*********************************************************************
load "./plot_type/functions_contrib.ncl"
load "./plot_type/uniform_coordinates.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
load "./interface_data/ncl.interface"
;*********************************************************************

begin

   print("++++++++++++++++++++++++++++++++++++++++++")
   print("plot SET12 - profiles (var: "+variable+")")
   print("++++++++++++++++++++++++++++++++++++++++++")

   field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
   plot_fields = (/"T3M","C3M"/)
   ind_f = ind(plot_fields .eq. field_number)
   if (all(ismissing(ind_f))) then
      print("fatal:can't use field, "+field_number+" for this plot_type")
      status_exit(1)
   end if
   delete(plot_fields)
   delete(ind_f)

   wkdir = plot_dir+"/profiles/"
   compare = "MODEL"
   file_type = "PS"

   dim_MOD = dimsizes(models@name)

   ;******************************************************************
   ;check the input parameters
  
   if (isatt(info,"profile_allstn")) then
      plotall = info@profile_allstn
   else
      plotall = False
   end if

   if (isatt(info,"profile_stnid")) then
      stnid_file = info@profile_stnid
      data_line = asciiread(stnid_file,-1,"string")
      cdata = stringtochar(data_line)
      ; read station name and id numbers for plotting
      station_names = charactertostring(cdata(:,39:))
      stnid2plot = stringtointeger(charactertostring(cdata(:,37:37)))
      station_ids = ind(stnid2plot.eq.0)
      delete(data_line)
      delete(cdata)
      delete(stnid2plot)
   else
      print("fatal:plot requires variable attribute 'profile_stnid'")
      status_exit(1)
   end if

   if plotall then       ; do all stations
      delete (station_ids)
      station_ids = ispan (0,55,1)       ; all station ids
   end if
   nstations = dimsizes(station_ids)

   if (isatt(info,"profile_raobs")) then
      raobs_file = info@profile_raobs
      if (isfilepresent_esmval(raobs_file))
	 fin = addfile(raobs_file, "r")   ; open raobs data file
      else
         print("fatal:File Not Found: " + raobs_file)
         status_exit(1)
      end if
      fin_vars = getfilevarnames(fin)
      idx_var = ind(fin_vars .eq. variable)
      if (all(ismissing(idx_var))) then
         print("fatal:Variable Not Found: " + variable)
         status_exit(1)
      else
         vout = fin->$fin_vars(idx_var)$
         dimNames = getvardims(vout)
         idx_pt = ind(dimNames .eq. "plev_trop")
         if (all(ismissing(idx_pt))) then   ; surface to TOA
            p0_raobs = fin->plev_atmo           ; raobs pressure levels
            tropOnly = False
         else
            p0_raobs = fin->plev_trop           ; raobs pressures in troposphere only
            tropOnly = True
         end if
         nlev0 = dimsizes(p0_raobs)
         slat = fin->slat                    ; raobs station latitudes
         slon = fin->slon                    ; raobs station longitudes

         ; reverse plev
         vout = vout(:,:,::-1)
         p0_raobs = p0_raobs(::-1)

         delete(dimNames)
         delete(idx_pt)
         delete(fin)
      end if
      delete(fin_vars)
      delete(idx_var)
   else
      print("fatal:plot requires variable attribute 'profile_raobs'")
      status_exit(1)
   end if

   ;******************************************************************
   ; constants
   Cp = 1.00464e3   ; J/(K*kg)
   g  = 9.80616     ; m/s^2
   L  = 2.501e6     ; J/kg
   p0 = 1000.       ; mb required by vinth2p
   TropMinP = 225.0   ; mb

   idx_mon = (/0,3,6,9/)
   month_name = (/"JAN","APR","JUL","OCT"/)

   ;******************************************************************
   dim_TOTAL = dim_MOD+1
   lgLabels = new((/dim_TOTAL/),string)
   colors = new((/dim_TOTAL/),string)
   dashs = new((/dim_TOTAL/),"integer")
   LineThick = new((/dim_TOTAL/),integer)
   LineMode = new((/dim_TOTAL/),string)
   Markers = new((/dim_TOTAL/),integer)

   lgLabels(0) = "raobs"
   lgLabels(1:) = (/models@name/)
   color_mod = CCMVal2_ColorTable(models@name)
   colors(0) = "black"
   colors(1:) = (/color_mod/)
   dash_mod = CCMVal2_DashTable(models@name)
   dashs(0) = 0
   dashs(1:) = (/dash_mod/)
   LineThick(:) = 2
   LineMode(0) = "Markers"
   LineMode(1:) = "Lines"
   Markers(0) = 16
   Markers(1:) = 0

   ;------------------------------------------------------------------
   ; common plot resources
   res                   = True
   res@trYReverse        = True
   res@trYMaxF           = 1025.0
   res@xyLineColors      = colors
   res@xyLineThicknesses = LineThick
   res@xyMarkLineModes   = LineMode
   res@xyDashPatterns    = dashs
   res@xyMarkers         = Markers
   res@xyMarkerColor     = "maroon"
   res@xyMarkerSizeF     = 0.01
   res@pmLegendWidthF    = 0.15
   if (dim_TOTAL.le.4) then
      res@pmLegendHeightF = 0.08
      res@lgLabelFontHeightF = 0.016
   else
      if (dim_TOTAL.le.8) then
         res@pmLegendHeightF = 0.15
         res@lgLabelFontHeightF = 0.012
      else
         res@pmLegendHeightF = 0.30
         res@lgLabelFontHeightF = 0.008
      end if
   end if
   res@lgPerimOn         = True
   res@txFontHeightF     = 0.022
   res@xyExplicitLegendLabels = lgLabels

   ;==================================================================
   ; setup data holders for plots
   dimXY = new((/dim_TOTAL/),"integer")
   dimXY(0) = nlev0
   do imod=0,dim_MOD-1
      A0 = extract_data(models,data,imod,-1,0,0)
      uniform_coordinates(A0,field_number)
      plvl = A0&plev
      if tropOnly then
        itop = ind_nearest_coord(TropMinP,plvl,0) ; nearest index of top level
        nlvl = itop + 1                         ; number of pressures to plot
      else
        nlvl = dimsizes(plvl)
      end if
      dimXY(imod+1) = nlvl
      delete(A0)
      delete(plvl)
   end do

   nMax = max(dimXY)
   Xdata = new((/dim_TOTAL,nMax,nstations,dimsizes(idx_mon)/),"float")
   Ydata = new((/dim_TOTAL,nMax/),"float")
   Ydata!0 = "line"
   Ydata!1 = "lev"

   ; get data into arrays (X) & (Y)
   do i=0,dim_MOD

      ; assign pressure level data (Y)
      if (i .eq. 0) then
         Ydata(0,0:dimXY(0)-1) = (/p0_raobs/)   ; raobs
      else
         imod = i-1
         A0 = extract_data(models,data,imod,-1,0,0)
         uniform_coordinates(A0,field_number)
         if isStrSubset(field_number,"T") then
            A = clmMonTLLL(A0)
         else
            A = A0
         end if
         plvl = A0&plev
         if tropOnly then
            Ydata(i,0:dimXY(i)-1) = (/plvl(:dimXY(i)-1)/)
            tmp = A(:,:dimXY(i)-1,:,:)
            delete(A)
            A = tmp
            delete(tmp)
         else
            Ydata(i,0:dimXY(i)-1) = (/plvl/)
         end if
         delete(A0)
         delete(plvl)
      end if

      ; assign data (X)
      do ns = 0, nstations-1         ;   loop over stations

         sid = station_ids(ns)        ; station id number
         latpt = slat(sid)            ; station latitude
         lonpt = slon(sid)            ; station longitude

         if (sid .eq. 3) then         ; Western Europe
            latpt = 51.626
         end if
         if (sid .eq. 26) then        ; San Francisco
            lonpt = 236.250
         end if
         if (sid .eq. 31) then        ; Sydney Australia
            lonpt = 151.875
         end if
         if (sid .eq. 32) then        ; Christchurch NZ
            lonpt = 174.375
         end if
         if (sid .eq. 34) then        ; Miami FL
            latpt = 26.511
         end if
         if (sid .eq. 47) then        ; Gibraltor
            latpt = 37.673
         end if

         if (latpt .ge. -87.8638) then ; model can't do South Pole
            do m = 0, 3                 ; loop over months
              if (i .eq. 0) then
                 Xdata(0,0:dimXY(0)-1,ns,m) = (/vout(sid,m,:)/)
              else
                 Xdata(i,0:dimXY(i)-1,ns,m) = (/A(idx_mon(m),:,{latpt},{lonpt})/)
              end if
            end do ;m
         end if
      end do ;ns

      if (i .ne. 0) then
         delete(A)
      end if
   end do  ; i

   delete(vout)
   delete(p0_raobs)

   ;==================================================================
   ; plot each stations
   do ns = 0, nstations-1         ;   loop over stations

      sid = station_ids(ns)        ; station id number
      latpt = slat(sid)            ; station latitude
      lonpt = slon(sid)            ; station longitude

      if (sid .eq. 3) then         ; Western Europe
         latpt = 51.626
      end if
      if (sid .eq. 26) then        ; San Francisco
         lonpt = 236.250
      end if
      if (sid .eq. 31) then        ; Sydney Australia
         lonpt = 151.875
      end if
      if (sid .eq. 32) then        ; Christchurch NZ
         lonpt = 174.375
      end if
      if (sid .eq. 34) then        ; Miami FL
         latpt = 26.511
      end if
      if (sid .eq. 47) then        ; Gibraltor
         latpt = 37.673
      end if

      if (latpt .ge. -87.8638) then ; model can't do South Pole

         print ("plot(" + (ns+1) + " of " + nstations + "): " + station_names(sid))
         print("=========================================")
         wks = gsn_open_wks(file_type,wkdir+"profiles_"+station_names(sid)+"_"+variable)
         plot = new(4,"graphic")
         res@tiXAxisString = info@long_name + " (" + info@units + ")"
         do m = 0, 3                 ; loop over months
            XdataF = Xdata(:,:,ns,m)
            XdataF!0 = "line"
            XdataF!1 = "lev"
            ;----------------------------------------------------------
            ; convert to g/kg
            if (variable.eq."Q") then
               XdataF = (/XdataF*1000./)        ; raobs in g/kg
               res@tiXAxisString = info@long_name + " (g/kg)"
            end if

            ;**********************************************************
            if (month_name(m).eq."JAN" .or. month_name(m).eq."JUL") then
               res@tiYAxisString = "Pressure (mb)"
            else
               delete(res@tiYAxisString)
            end if
            res@gsnLeftString = month_name(m)

            if (month_name(m).eq."APR") then
               res@pmLegendDisplayMode    = "Always"
               res@pmLegendSide           = "Right"
	       if (variable.eq."H") then
                  res@pmLegendOrthogonalPosF    = -.45
                  if (dim_TOTAL.le.4) then
                     res@pmLegendParallelPosF   = .10
                  else
                     if (dim_TOTAL.le.8) then
                        res@pmLegendParallelPosF = .15
                     else
                        res@pmLegendParallelPosF = .30
                     end if
                  end if
               end if
               if (variable.eq."T" .or. variable.eq."ta") then
                  res@pmLegendOrthogonalPosF    = -1.04
                  if (dim_TOTAL.le.4) then
                     res@pmLegendParallelPosF   = .10
                  else
                     if (dim_TOTAL.le.8) then
                        res@pmLegendParallelPosF = .15
                     else
                        res@pmLegendParallelPosF = .30
                     end if
                  end if
               end if
               if (variable.eq."Q") then

                  res@pmLegendOrthogonalPosF    = -.45
                  if (dim_TOTAL.le.4) then
                     res@pmLegendParallelPosF   = .90
                  else
                     if (dim_TOTAL.le.8) then
                        res@pmLegendParallelPosF = .80
                     else
                        res@pmLegendParallelPosF = .70
                     end if
                  end if
               end if
            else
               res@pmLegendDisplayMode    = "Never"
            end if

            res@gsnDraw = False
            res@gsnFrame = False
            plot(m) = gsn_csm_xy(wks,XdataF,Ydata,res)  ; create the plot

            delete(res@gsnLeftString)
            delete(XdataF)
         end do           ; loop over months

         ;**********************************************************
         ; add lat long info to plot
         txres               = True
         txres@txFontHeightF = 0.017

         if (latpt.lt.0) then
            prntlat = fabs(latpt)+ " S"
         else
            prntlat = latpt + " N"
         end if
         if (lonpt.gt.180) then
            prntlon = 360-lonpt + " W"
         else
            prntlon = lonpt + " E"
         end if
         title = station_names(sid)+"  latitude= "+prntlat+" longitude= "+prntlon
         gsn_text_ndc(wks,title,.50,.98,txres)

         pan = True
         pan@gsnFrame = False
         pan@gsnMaximize = True
         pan@gsnPanelTop = 0.96
         gsn_panel (wks,plot,(/2,2/),pan)
         frame(wks)

         delete(title)
         delete(res@tiXAxisString)
         delete(txres)
         delete(plot)
         delete(pan)
      end if               ; check for south pole
   end do               ; loop over station

   delete(res)
   delete(station_ids)
   delete(station_names)
   delete(lgLabels)
   delete(colors)
   delete(dashs)
   delete(LineThick)
   delete(LineMode)
   delete(Markers)
   delete(Xdata)
   delete(Ydata)

end
