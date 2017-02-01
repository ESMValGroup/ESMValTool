;*********************************************************************
; make_table_one.ncl
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
; This script makes tables of seasonal & regional means, differences,
; or RMSE depends on the file created by script table_one.ncl
;*********************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "./plot_type/taylor_metrics_table.ncl"
;*********************************************************************

begin
   print("")
   print("++++++++++++++++++++++++++++++++++++++++++")
   print("plot SET1 - make_table_one")
   print("++++++++++++++++++++++++++++++++++++++++++")
   print("")

   plot_dir = getenv("ESMValTool_plot_dir")
   wkdir = plot_dir+"/table_one/"
   file_type = "PS"

   ;============================================================
   masterfile = wkdir+"table_one.nc"
   fin = addfile(masterfile,"r")
   vars = fin->VARNAMES
   varStr = charactertostring(vars)
   nvars = dimsizes(varStr)

   noSEASON = stringtoint(fin@numSEASON)
   SEASONS = new(noSEASON,"string")
   do i=0,noSEASON-1
      vName = "SEASON" + sprinti("%0.2i", i)
      SEASONS(i) = fin@$vName$
   end do
   noREGION = stringtoint(fin@numREGION)
   regions = new(noREGION,"string")
   do i=0,noREGION-1
      vName = "REGION" + sprinti("%0.2i", i)
      labelREG = fin@$vName$
      ;regions = (/"GLOBAL","TROPICS","NORTH_EXT","SOUTH_EXT"/)
      if (labelREG .eq. "GLBL") then
         regions(i) = "GLOBAL"
      end if
      if (labelREG .eq. "TROP") then
         regions(i) = "TROPICS"
      end if
      if (labelREG .eq. "NEXT") then
         regions(i) = "NORTH_EXT"
      end if
      if (labelREG .eq. "SEXT") then
         regions(i) = "SOUTH_EXT"
      end if
   end do
   idx_ref = stringtoint(fin@indexREF)
   out_type = fin@out_type
   noMODEL = stringtoint(fin@numMODELS)
   if (idx_ref .eq. -1) then
      labels1 = new(noMODEL,"string")
      do imod=0,noMODEL-1
         vName = "MODEL" + sprinti("%0.2i", imod)
         labels1(imod) = fin@$vName$
      end do
      labels2 = labels1
   else
      labels1 = new(noMODEL-1,"string")
      labels2 = new(noMODEL-1,"string")
      vName = "MODEL" + sprinti("%0.2i", idx_ref)
      refMODEL = fin@$vName$

      jmod = 0
      do imod=0,noMODEL-1
         vName = "MODEL" + sprinti("%0.2i", imod)
	 if (imod .ne. idx_ref) then
	    labels1(jmod) = fin@$vName$
            if (out_type .eq. "DIFF") then
               labels2(jmod) = refMODEL + "-" + fin@$vName$
            else
               labels2(jmod) = refMODEL + " vs " + fin@$vName$
            end if
            jmod = jmod + 1
	 end if
      end do
   end if
   nLabel = dimsizes(labels1)

   ;============================================================

   do r=0,dimsizes(regions)-1

      region =  regions(r)

      opt = True
      opt@pltType = file_type
      opt@tableTitle = region + "(" + out_type + ")"
      opt@color0 = "white"
      opt@color1 = "white"

      if (region .eq. "GLOBAL") then
         table_data = fin->GLOBAL
      end if
      if (region .eq. "TROPICS") then
         table_data = fin->TROPICS
      end if
      if (region .eq. "NORTH_EXT") then
         table_data = fin->NEXTROP
      end if
      if (region .eq. "SOUTH_EXT") then
         table_data = fin->SEXTROP
      end if

      do s=0,dimsizes(SEASONS)-1

         season = (/SEASONS(s)/)

         tabVal = table_data(seasons|s:s,t1models|:,record|:)
         table_name = wkdir + "table_one_" + season + "_" + region + "_" + out_type
         if (out_type .eq. "MEAN") then
            case_name = season
         else
            if (out_type .eq. "DIFF") then
               case_name = season + " ["+ refMODEL + " - modelX]"
            else
               case_name = season + " ["+ refMODEL + " vs modelX]"
            end if
         end if
         taylor_metrics_table(table_name, varStr, (/case_name/), labels1, tabVal, opt)

         delete(tabVal)
      end do  ; s

      if (dimsizes(SEASONS) .gt. 1) then
         tabValAll = table_data(t1models|:,seasons|:,record|:)
         summary_table = wkdir + "table_one_" + region + "_" + out_type
         taylor_metrics_table(summary_table, varStr, labels2, SEASONS, tabValAll, opt)
         delete(tabValAll)
      end if

      delete(opt)
      delete(table_data)

   end do   ; r

   delete(regions)
   delete(SEASONS)
   delete(fin)
   delete(vars)
   delete(varStr)
   delete(nvars)
   delete(labels1)
   delete(labels2)

   system("rm "+masterfile)

end
