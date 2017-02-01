;*********************************************************************
; table_one.ncl
;*********************************************************************
; NCL script
; tables.ncl
; Mark Stevens, Sept 2001
;*********************************************************************
; The original NCL script, plot_surfaces_cons.ncl, was written to make...
; tables of global and regional means, differences and RMSE (set 1)
; for the AMWG diagnostics package.
; http://www.cgd.ucar.edu/cms/diagnostics/
;
; The NCL script described above has been modified for use with the CCMVal
; diagnostic software package.
; (PIs: Andrew Gettelman & Veronika Eyring, co-I: Hisako Shiona)
;
; The switch 'table_one_type' stands for the output option. If present
; and set to "DIFF", then differences with the other models are calculated.
; If present and set to "RMSE", then weighted area root-mean-square-
; differences (rmse) with the other models are calculated.
; In both cases, requires 'refModel' which stands for the reference model.
; If omitted, mean of all models are calculated.
;
; This script creates a file, table_one.nc, which contains the results
; and then the file gets called by another script, make_table_one.ncl,
; to make summary table(s).
;
; This script takes input from the var_att.ncl files
; info attributes required:
; table_one_type	(optional) output option ("MEAN", "DIFF", or
;			"RMSE") Note: if not set, use default "MEAN"
; refModel		reference model name (used with table_one_type)
; seasons		(optional) array of for averaging. (DJF, MAM,
;			JJA, SON or ANN for annual mean)
;			Note: if not set, use default (/"DJF","JJA","ANN"/)
; regions		(optional) array of regions.
;			"GLBL": GLOBAL
;			"TROP": TROPICS (20N-20S)
;			"NEXT": NORTHERN EXTRATROPICS (90N-20N)
;			"SEXT": SOUTHERN EXTRATROPICS (90S-20S
;			Note: if not set, use default, all four regions
;			(/"GLBL","TROP","NEXT","SEXT"/)
;*********************************************************************
load "./plot_type/calc_seasonal.ncl"
load "./plot_type/uniform_coordinates.ncl"
load "./interface_data/ncl.interface"
;*********************************************************************

begin

   print("++++++++++++++++++++++++++++++++++++++++++")
   print("plot SET1 - table_one (variable: "+variable+")")
   print("++++++++++++++++++++++++++++++++++++++++++")

   field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
   plot_fields = (/"T2Ms","T2Ds","C2Ms"/)
   ind_f = ind(plot_fields .eq. field_number)
   if (all(ismissing(ind_f))) then
      print("fatal:can't use field, "+field_number+" for this plot_type")
      status_exit(1)
   end if
   delete(plot_fields)
   delete(ind_f)

   wkdir = plot_dir+"/table_one/"
   compare = "MODEL"

   ; constants required for varaible "AAM"
   re = 6.3712e6         ; radius of earth
   pi = 3.14159265

   ; extract settings from materfile if exists
   masterfile = wkdir+"table_one.nc"
   if (isfilepresent_esmval(masterfile)) then
      fin = addfile(masterfile,"r")
      ds = dimsizes(fin->VARNAMES)
      nvars = ds(0)
      outTYPE = fin@out_type
      indexREF =  stringtoint(fin@indexREF)
      noMODEL =  stringtoint(fin@numMODELS)
      noSEASON = stringtoint(fin@numSEASON)
      labelSEA = new(noSEASON,"string")
      do i=0,noSEASON-1
         vName = "SEASON" + sprinti("%0.2i", i)
	 labelSEA(i) = fin@$vName$
      end do
      noREGION = stringtoint(fin@numREGION)
      labelREG = new(noREGION,"string")
      do i=0,noREGION-1
         vName = "REGION" + sprinti("%0.2i", i)
	 labelREG(i) = fin@$vName$
      end do
      delete(fin)
   else
      nvars = 0
   end if

   labelMOD = models@name
   ;labelMOD = models@name + " " + models@case_name
   dim_MOD = dimsizes(labelMOD)
   if (nvars .gt. 0) then
      if (dim_MOD .ne. noMODEL) then
         print("fatal:plot_type(table_one): NUMBER OF MODELS is different between variables")
         status_exit(1)
      end if
      delete(noMODEL)
   end if

   ;******************************************************************
   ;check the input parameters and compare them with masterfile settings

   if (isatt(info,"table_one_type")) then
      out_type = info@table_one_type
      if ((out_type .eq. "DIFF") .or. (out_type .eq. "RMSE")) then
         if (isatt(info,"refModel")) then
            idx_ref = ind(models@name .eq. info@refModel)
            if (all(ismissing(idx_ref))) then
               print("fatal:refModel isn't included in model list (namelist)")
               status_exit(1)
            end if
         else
            print("fatal:plot requires variable attribute 'refModel'")
            status_exit(1)
         end if
      end if
   else
      out_type = "MEAN"
      idx_ref = -1
   end if
   if (nvars .gt. 0) then
      if ((out_type .ne. outTYPE) .or. (idx_ref .ne. indexREF)) then
         print("fatal:plot_type(table_one): OUTPUT_TYPE is different between variables")
         status_exit(1)
      end if
      delete(outTYPE)
      delete(indexREF)
   end if

   if (isatt(info,"seasons")) then
      seasons = info@seasons
   else
      seasons = (/"ANN","DJF","JJA"/)
   end if
   ns = dimsizes(seasons)
   if (nvars .gt. 0) then
      if (ns .ne. noSEASON) then
        print("fatal:plot_type(table_one): SEASON DIMSIZE is different between variables")
        status_exit(1)
      else
         do i=0,noSEASON-1
            if (labelSEA(i) .ne. seasons(i)) then
               print("fatal:plot_type(table_one): SEASON NAMES is different between variables")
               status_exit(1)
            end if
         end do
      end if
      delete(noSEASON)
      delete(labelSEA)
   end if

   if (isatt(info,"regions")) then
      regions = info@regions
   else
      regions = (/"GLBL","TROP","NEXT","SEXT"/)
   end if
   nreg = dimsizes(regions)
   if (nvars .gt. 0) then
      if (nreg .ne. noREGION) then
        print("fatal:plot_type(table_one): REGION DIMSIZE is different between variables")
        status_exit(1)
      else
         do i=0,noREGION-1
            if (labelREG(i) .ne. regions(i)) then
               print("fatal:plot_type(table_one): REGION NAMES is different between variables")
               status_exit(1)
            end if
         end do
      end if
      delete(noREGION)
      delete(labelREG)
   end if

   ;******************************************************************
   nrec = 1
   if (idx_ref .ne. -1) then
      ; (b) difference/rmse of each model from 'reference'
      nmod = dim_MOD-1
   else
      ; (a) full field
      nmod = dim_MOD
   end if
   stats  = new((/nrec,ns,nreg,nmod/),"float")

   ;============================================================
   ; step1: compute the means, diff and rmse of the two fields.
   ;============================================================

   ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ;loop through seasons
   do s=0,dimsizes(seasons)-1
      season = seasons(s)

      if (idx_ref .ne. -1) then
         ; out_type = "DIFF"/"RMSE"

         ; CASE 1
         A0 = extract_data(models,data,idx_ref,-1,0,0)
         uniform_coordinates(A0,field_number)
         ; convert to T2Ms
         if (field_number .eq. "T2Ds") then
            tmp = day_to_month(A0,"T2Ds")
            delete(A0)
            A0 = tmp
            delete(tmp)
         end if

         lat1 = A0&lat
         nlat1 = dimsizes(lat1)
         lon1 = A0&lon
         nlon1 = dimsizes(lon1)
         if (lon1(0).lt.0.0) then
            lon1 = lon1 + 180.0
            A0&lon = (/lon1/)
         end if

         ;gw1 = inptr1->gw
         if (mod(nlat1,2).eq.1) then
            gw1_tmp = latGauWgt(nlat1-1,"lat","gaussian weights","none")  ; no gw in file
            gw1 = new(nlat1,float)
            gw1(0:(nlat1/2)) = gw1_tmp(0:(nlat1/2))
            gw1(((nlat1/2)+1):(nlat1-1)) = gw1_tmp((nlat1/2):(nlat1-2))
            delete(gw1_tmp)
         else
            gw1 = latGauWgt(nlat1,"lat","gaussian weights","none")  ; no gw in file
         end if
         gw1!0 = "lat"
         gw1&lat = lat1

         A = calc_SEASON(A0,season)
         A&lat@units = "degrees_north"
         A&lon@units = "degrees_east"
         delete(A0)

         if (variable.eq."ICEAREA") then
            ;if (isfilevar(inptr1,"slat")) then   ; finite volume lat/lon grid
            SLAT1 = True
         end if

         jmod = 0
         do imod=0,dim_MOD-1
	    if (imod .ne. idx_ref) then
               ; CASE 2
               case2 = models@name(imod) + " " + models@case_name(imod)
               B0 = extract_data(models,data,imod,-1,0,0)
               uniform_coordinates(B0,field_number)
               ; convert to T2Ms
               if (field_number .eq. "T2Ds") then
                  tmp = day_to_month(B0,"T2Ds")
                  delete(B0)
                  B0 = tmp
                  delete(tmp)
               end if

               lat2 = B0&lat
               nlat2 = dimsizes(lat2)
               lon2 = B0&lon
               nlon2 = dimsizes(lon2)
               if (lon2(0).lt.0.0) then
                  lon2 = lon2 + 180.0
                  B0&lon = (/lon2/)
               end if

               ;gw2 = inptr2->gw
               if (mod(nlat2,2).eq.1) then
                  gw2_tmp = latGauWgt(nlat2-1,"lat","gaussian weights","none")  ; no gw in file
                  gw2 = new(nlat2,float)
                  gw2(0:(nlat2/2)) = gw2_tmp(0:(nlat2/2))
                  gw2(((nlat2/2)+1):(nlat2-1)) = gw2_tmp((nlat2/2):(nlat2-2))
                  delete(gw2_tmp)
               else
                  gw2 = latGauWgt(nlat2,"lat","gaussian weights","none")  ; no gw in file
               end if
               gw2!0 = "lat"
               gw2&lat = lat2

               B = calc_SEASON(B0,season)
               B&lat@units = "degrees_north"
               B&lon@units = "degrees_east"
               delete(B0)

               if (variable.eq."ICEAREA") then
                  ;if (isfilevar(inptr2,"slat")) then   ; finite volume lat/lon grid
                  SLAT2 = True
               end if

               if (out_type .eq. "DIFF") then

                  ;regions = (/"GLBL","TROP","NEXT","SEXT"/)
                  ; GLOBAL
                  idx_reg = ind(regions .eq. "GLBL")
                  if (.not.ismissing(idx_reg)) then
                     if (variable.eq."ICEAREA") then
                        if (SLAT1) then   ; finite volume lat/lon grid
                           WgtAveA = iceareaFIX(A,0)+iceareaFIX(A,1)     ; north+south
                        else
                           WgtAveA = iceareaGAU(A,0)+iceareaGAU(A,1)     ; north+south
                        end if
                        if (SLAT2) then
                           WgtAveB = iceareaFIX(B,0)+iceareaFIX(B,1)
                        else
                           WgtAveB = iceareaGAU(B,0)+iceareaGAU(B,1)
                        end if
                     else
                        WgtAveA = wgt_areaave(A,gw1,1.0,0)            ; case1
                        WgtAveB = wgt_areaave(B,gw2,1.0,0)            ; case2
                     end if
                     WgtAveC = WgtAveA - WgtAveB
                     if (variable.eq."AAM") then
                        WgtAveC = WgtAveC*4.*pi*re^2
                     end if
                     stats(0,s,idx_reg,jmod) = (/WgtAveC/)
                  end if

                  ; TROPICAL
                  idx_reg = ind(regions .eq. "TROP")
                  if (.not.ismissing(idx_reg)) then
                     if ((variable.eq."ICEAREA") .or. (variable.eq."SNOWHLND") .or. (variable.eq."AAM")) then
                        WgtAveA = 0.0
                        WgtAveB = 0.0
                     else
                        WgtAveA = wgt_areaave(A({-20:20},:),gw1({-20:20}),1.0,0)
                        WgtAveB = wgt_areaave(B({-20:20},:),gw2({-20:20}),1.0,0)
                      end if
                     WgtAveC = WgtAveA - WgtAveB
                     stats(0,s,idx_reg,jmod) = (/WgtAveC/)
                  end if

                  ; NORTH EXTRATROPICAL
                  idx_reg = ind(regions .eq. "NEXT")
                  if (.not.ismissing(idx_reg)) then
                     if (variable.eq."ICEAREA") then
                        if (SLAT1) then   ; finite volume lat/lon grid
                           WgtAveA = iceareaFIX(A,1)     ; north
                        else
                           WgtAveA = iceareaGAU(A,1)     ; north
                        end if
                        if (SLAT2) then
                           WgtAveB = iceareaFIX(B,1)
                        else
                           WgtAveB = iceareaGAU(B,1)
                        end if
                     else
                        if (variable.eq."AAM") then
                           WgtAveA = 0.0
                           WgtAveB = 0.0
                        else
                           WgtAveA = wgt_areaave(A({20:90},:),gw1({20:90}),1.0,0)
                           WgtAveB = wgt_areaave(B({20:90},:),gw2({20:90}),1.0,0)
                        end if
                     end if
                     WgtAveC = WgtAveA - WgtAveB
                     stats(0,s,idx_reg,jmod) = (/WgtAveC/)
                  end if

                  ; SOUTH EXTRATROPICAL
                  idx_reg = ind(regions .eq. "SEXT")
                  if (.not.ismissing(idx_reg)) then
                     if (variable.eq."ICEAREA") then
                        if (SLAT1) then   ; finite volume lat/lon grid
                           WgtAveA = iceareaFIX(A,1)     ; south
                        else
                           WgtAveA = iceareaGAU(A,1)     ; south
                        end if
                        if (SLAT2) then
                           WgtAveB = iceareaFIX(B,1)
                        else
                           WgtAveB = iceareaGAU(B,1)
                        end if
                     else
                        if (variable.eq."AAM") then
                           WgtAveA = 0.0
                           WgtAveB = 0.0
                        else
                           WgtAveA = wgt_areaave(A({-90:-20},:),gw1({-90:-20}),1.0,0)
                           WgtAveB = wgt_areaave(B({-90:-20},:),gw2({-90:-20}),1.0,0)
                        end if
                     end if
                     WgtAveC = WgtAveA - WgtAveB
                     stats(0,s,idx_reg,jmod) = (/WgtAveC/)
                  end if
               else
                  ; out_type = "RMSE"
                  if (nlat1 .ge. nlat2) then           ; interp A to nlat2
                     tmpA = B
                     tmpA = linint2(lon1,lat1,A,False,lon2,lat2,0)
                     tmpB = B
                     gwC = gw2
                  else
                     tmpB = A
                     tmpB = linint2(lon2,lat2,B,False,lon1,lat1,0)
                     tmpA = A
                     gwC = gw1
                  end if

                  ; GLOBAL
                  idx_reg = ind(regions .eq. "GLBL")
                  if (.not.ismissing(idx_reg)) then
                     stats(0,s,idx_reg,jmod) = wgt_arearmse (tmpA,tmpB,gwC,1.0,0)
                  end if

                  ; TROPICAL
                  idx_reg = ind(regions .eq. "TROP")
                  if (.not.ismissing(idx_reg)) then
                     stats(0,s,idx_reg,jmod) = wgt_arearmse (tmpA({-20:20},:),tmpB({-20:20},:),gwC({-20:20}),1.0,0)
                  end if

                  ; NORTH EXTRATROPICAL
                  idx_reg = ind(regions .eq. "NEXT")
                  if (.not.ismissing(idx_reg)) then
                     stats(0,s,idx_reg,jmod) = wgt_arearmse (tmpA({20:90},:),tmpB({20:90},:),gwC({20:90}),1.0,0)
                  end if

                  ; SOUTH EXTRATROPICAL
                  idx_reg = ind(regions .eq. "SEXT")
                  if (.not.ismissing(idx_reg)) then
                     stats(0,s,idx_reg,jmod) = wgt_arearmse (tmpA({-90:-20},:),tmpB({-90:-20},:),gwC({-90:-20}),1.0,0)
                  end if
                  delete(tmpA)
                  delete(tmpB)
                  delete(gwC)

               end if

               delete(lat2)
               delete(nlat2)
               delete(lon2)
               delete(nlon2)
               delete(gw2)
               delete(B)

               jmod = jmod + 1
            end if
         end do  ; imod

         delete(lat1)
         delete(nlat1)
         delete(lon1)
         delete(nlon1)
         delete(gw1)
         delete(A)

      else
         ; out_type = "MEAN"
         do imod=0,dim_MOD-1
            ; CASE 1
            A0 = extract_data(models,data,imod,-1,0,0)
            uniform_coordinates(A0,field_number)
            ; convert to T2Ms
            if (field_number .eq. "T2Ds") then
               tmp = day_to_month(A0,"T2Ds")
               delete(A0)
               A0 = tmp
               delete(tmp)
            end if

            lat1 = A0&lat
            nlat1 = dimsizes(lat1)

            ;gw1 = inptr1->gw
            if (mod(nlat1,2).eq.1) then
               gw1_tmp = latGauWgt(nlat1-1,"lat","gaussian weights","none")  ; no gw in file
               gw1 = new(nlat1,float)
               gw1(0:(nlat1/2)) = gw1_tmp(0:(nlat1/2))
               gw1(((nlat1/2)+1):(nlat1-1)) = gw1_tmp((nlat1/2):(nlat1-2))
               delete(gw1_tmp)
            else
               gw1 = latGauWgt(nlat1,"lat","gaussian weights","none")  ; no gw in file
            end if
            gw1!0 = "lat"
            gw1&lat = lat1

            A = calc_SEASON(A0,season)
            A&lat@units = "degrees_north"
            A&lon@units = "degrees_east"
            delete(A0)

            if (variable.eq."ICEAREA") then
               ;if (isfilevar(inptr1,"slat")) then   ; finite volume lat/lon grid
               SLAT1 = True
            end if

            ; GLOBAL
            idx_reg = ind(regions .eq. "GLBL")
            if (.not.ismissing(idx_reg)) then
               if (variable.eq."ICEAREA") then
                  if (SLAT1) then   ; finite volume lat/lon grid
                     WgtAveA = iceareaFIX(A,0)+iceareaFIX(A,1)     ; north+south
                  else
                     WgtAveA = iceareaGAU(A,0)+iceareaGAU(A,1)     ; north+south
                  end if
               else
                  WgtAveA = wgt_areaave(A,gw1,1.0,0)            ; case1
               end if
               if (variable.eq."AAM") then
                  WgtAveA = WgtAveA*4.*pi*re^2
               end if
               stats(0,s,idx_reg,imod) = (/WgtAveA/)
            end if

            ; TROPICAL
            idx_reg = ind(regions .eq. "TROP")
            if (.not.ismissing(idx_reg)) then
               if ((variable.eq."ICEAREA") .or. (variable.eq."SNOWHLND") .or. (variable.eq."AAM")) then
                  WgtAveA = 0.0
               else
                  WgtAveA = wgt_areaave(A({-20:20},:),gw1({-20:20}),1.0,0)
               end if
               stats(0,s,idx_reg,imod) = (/WgtAveA/)
            end if

            ; NORTH EXTRATROPICAL
            idx_reg = ind(regions .eq. "NEXT")
            if (.not.ismissing(idx_reg)) then
               if (variable.eq."ICEAREA") then
                  if (SLAT1) then   ; finite volume lat/lon grid
                     WgtAveA = iceareaFIX(A,1)     ; north
                  else
                     WgtAveA = iceareaGAU(A,1)     ; north
                  end if
               else
                  if (variable.eq."AAM") then
                     WgtAveA = 0.0
                  else
                     WgtAveA = wgt_areaave(A({20:90},:),gw1({20:90}),1.0,0)
                  end if
               end if
               stats(0,s,idx_reg,imod) = (/WgtAveA/)
            end if

            ; SOUTH EXTRATROPICAL
            idx_reg = ind(regions .eq. "SEXT")
            if (.not.ismissing(idx_reg)) then
               if (variable.eq."ICEAREA") then
                  if (SLAT1) then   ; finite volume lat/lon grid
                     WgtAveA = iceareaFIX(A,1)     ; south
                  else
                     WgtAveA = iceareaGAU(A,1)     ; south
                  end if
               else
                  if (variable.eq."AAM") then
                     WgtAveA = 0.0
                  else
                     WgtAveA = wgt_areaave(A({-90:-20},:),gw1({-90:-20}),1.0,0)
                  end if
               end if
               stats(0,s,idx_reg,imod) = (/WgtAveA/)
            end if

            delete(lat1)
            delete(nlat1)
            delete(gw1)
            delete(A)
         end do  ; imod

      end if

   ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   end do ;s
   ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   ;============================================================
   ; step2: create/open output file and write/append the results
   ;============================================================

   rec = nvars
   rec!0 = "record"

   nlen = 32
   varNames = new((/1,nlen/),"character")
   varC = stringtochar(variable)
   len_varC = dimsizes(varC)-1
   varNames(0,0:len_varC) = varC(0:len_varC)
   varNames!0 = "record"
   varNames!1 = "len"
   varNames&record = rec

   var_names3D = new(nreg,"string")
   varvar_types3D = new(nreg,"string")
   varvar_types3D(:) = "double"

   ;var_names3D = (/"GLOBAL", "TROPICS", "NEXTROP", "SEXTROP"/)
   ; GLOBAL
   idx_reg = ind(regions .eq. "GLBL")
   if (.not.ismissing(idx_reg)) then
      var_names3D(idx_reg) = "GLOBAL"
   end if

   ; TROPICAL
   idx_reg = ind(regions .eq. "TROP")
   if (.not.ismissing(idx_reg)) then
      var_names3D(idx_reg) = "TROPICS"
   end if

   ; NORTH EXTRATROPICAL
   idx_reg = ind(regions .eq. "NEXT")
   if (.not.ismissing(idx_reg)) then
      var_names3D(idx_reg) = "NEXTROP"
   end if

   ; SOUTH EXTRATROPICAL
   idx_reg = ind(regions .eq. "SEXT")
   if (.not.ismissing(idx_reg)) then
      var_names3D(idx_reg) = "SEXTROP"
   end if

   varfile = wkdir+"table_one_"+variable+".nc"
   fout = addfile(varfile,"c")

   ; predefine coordinate variables
   dimNames = (/"record","seasons","t1models","len"/)
   dimSizes = (/-1,ns,nmod,nlen/)
   dimUnlim = (/True,False,False,False/)
   filedimdef(fout,dimNames,dimSizes,dimUnlim)

   ; predefine variable names, type, and dimensions
   filevardef(fout,"record",(/"integer"/),(/"record"/))
   filevardef(fout,"t1models",(/"integer"/),(/"t1models"/))
   filevardef(fout,"seasons",(/"integer"/),(/"seasons"/))
   filevardef(fout,var_names3D,varvar_types3D,(/"record", "seasons", "t1models"/))
   filevardef(fout,"VARNAMES",typeof(varNames),getvardims(varNames))

   ; define global attributes.
   globalAtt = True
   globalAtt@out_type = out_type
   globalAtt@indexREF = sprinti("%0.2i", idx_ref)
   globalAtt@numMODELS = sprinti("%0.2i", dim_MOD)
   do imod=0,dim_MOD-1
      vName = "MODEL" + sprinti("%0.2i", imod)
      globalAtt@$vName$ = labelMOD(imod)
   end do
   globalAtt@numSEASON = sprinti("%0.2i", ns)
   do i=0,ns-1
      vName = "SEASON" + sprinti("%0.2i", i)
      globalAtt@$vName$ = seasons(i)
   end do
   globalAtt@numREGION = sprinti("%0.2i", nreg)
   do i=0,nreg-1
      vName = "REGION" + sprinti("%0.2i", i)
      globalAtt@$vName$ = regions(i)
   end do
   fileattdef(fout,globalAtt)

   ; output data values only
   fout->record   = (/rec/)
   fout->seasons  = ispan(0,ns-1,1)
   fout->t1models = ispan(0,nmod-1,1)
   fout->VARNAMES = (/varNames/)

   ; GLOBAL
   idx_reg = ind(regions .eq. "GLBL")
   if (.not.ismissing(idx_reg)) then
      fout->GLOBAL   = (/stats(:,:,idx_reg,:)/)
   end if

   idx_reg = ind(regions .eq. "TROP")
   if (.not.ismissing(idx_reg)) then
      fout->TROPICS  = (/stats(:,:,idx_reg,:)/)
   end if

   idx_reg = ind(regions .eq. "NEXT")
   if (.not.ismissing(idx_reg)) then
      fout->NEXTROP  = (/stats(:,:,idx_reg,:)/)
   end if

   idx_reg = ind(regions .eq. "SEXT")
   if (.not.ismissing(idx_reg)) then
      fout->SEXTROP  = (/stats(:,:,idx_reg,:)/)
   end if

   if (nvars .eq. 0) then
      move_cmd = "mv "+ varfile + " " + masterfile
      system(move_cmd)
      print ("move_cmd="+move_cmd)
   else
      ncecat_cmd = "ncrcat -O " + masterfile + " " + varfile + " " + masterfile
      system(ncecat_cmd)
      print ("ncecat_cmd="+ncecat_cmd)
      remove_cmd = "rm "+ varfile
      system(remove_cmd)
      print ("remove_cmd="+remove_cmd)
   end if

   delete(seasons)
   delete(regions)
   delete(labelMOD)

end
