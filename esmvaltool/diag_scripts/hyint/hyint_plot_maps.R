######################################################
#---------Maps plotting routine for HyInt------------#
#-------------E. Arnone (September 2017)-------------#
######################################################

# DECLARING THE FUNCTION: EXECUTION IS AT THE BOTTOM OF THE SCRIPT

hyint.plot.maps<-function(work_dir,plot_dir,ref_dir,ref_idx,season) {

# load settings
source('esmvaltool/diag_scripts/aux/hyint/hyint_parameters.r')
for (myname in names(settings)) { temp=get(myname,settings); assign(myname,temp) }

#source('interface_data/r.interface')
# diag_base <- "HyInt"
 
# setting up path and parameters
dataset_ref=models_name[ref_idx]
model_exp_ref=models_experiment[ref_idx]
model_ens_ref=models_ensemble[ref_idx]
year1_ref=models_start_year[ref_idx]
year2_ref=models_end_year[ref_idx]
years_ref <- year1_ref:year2_ref

# check path to reference dataset
#if (ref_dir!=work_dir) {
#  ref_dir=file.path(ref_dir,diag_base,paste0(year1_ref,"_",year2_ref),season)
#} else {
#  ref_dir=file.path(work_dir,dataset_ref,paste0(year1_ref,"_",year2_ref),season)
#}

# Define fields to be used
if (selfields[1]!=F) {
  field_names=field_names[selfields,drop=F]
  levels_m=levels_m[selfields,,drop=F]
  tlevels_m=tlevels_m[selfields,,drop=F]
  title_unit_m=title_unit_m[selfields,,drop=F] 
}
nfields=length(field_names)

# Define quantity (exp, ref, exp-ref) to be plotted depending on plot_type 
nquantity=c(1,3,3,1) # 1=exp_only, 3=exp/ref/exp-ref

# Define regions to be used
nregions=length(selregions)
if (nregions > dim(regions)[1]) { stop(paste(diag_base,": requesting regions outside list"))}


# ------- loading reference data ----------
# load topography if needed
if (maskSeaLand) {
  topo_file_idx=paste0(topography_file,ref_idx,'.nc')
  gridfile<-getfilename.indices(work_dir,diag_base,ref_idx,grid=T)
  if (!file.exists(topo_file_idx)) { createLandSeaMask(regrid=gridfile,regridded_topo=topo_file_idx,topo_only=T) }
  relevation=ncdf.opener(topo_file_idx,"topo","lon","lat",rotate="no")
}
if (highreselevation) { highresel=get.elevation(elev_range=c(highreselevation,9000)) }

# produce desert areas map if required from reference file (mean annual precipitation <0.5 mm, Giorgi et al. 2014)
if (removedesert) {
  ref_filename<-getfilename.indices(ref_dir,diag_base,ref_idx,season)
  pry<-ncdf.opener(ref_filename,"pry","lon","lat",rotate="no")
  retdes=which(pry<0.5)
  pry[retdes]=NA
  retdes2D=apply(pry*0,c(1,2),sum)+1 # create mask with NAs for deserts and 1's for not-desert
  retdes3D=replicate(dim(pry)[length(dim(pry))],retdes2D) # replicate for number of years
}

# open reference field
ref_filename<-getfilename.indices(ref_dir,diag_base,ref_idx,season)
print(paste("Reading reference ",ref_filename))
for (field in field_names) {
  field_ref=ncdf.opener(ref_filename,field,"lon","lat",rotate="no")
  if (removedesert) { field_ref<-field_ref*retdes3D }
  if (maskSeaLand) {field_ref<-apply.elevation.mask(field_ref,relevation,seaLandElevation)}
  if (rmultiyear_mean) { # if requested calculate multiyear average and store at time=1
    field_ref[,,1]<-apply(field_ref,c(1,2),mean,na.rm=T)
  }
  assign(paste(field,"_ref",sep=""),field_ref)
}

# Loop over models
for (model_idx in c(1:(length(models_name)))) {

  # setting up path and parameters
  exp    <- models_name[model_idx]
  model_exp <- models_experiment[model_idx]
  model_ens <- models_ensemble[model_idx]
  year1  <- models_start_year[model_idx]
  year2  <- models_end_year[model_idx]

  # set main paths
  work_dir_exp=work_dir #file.path(work_dir,exp,paste0(year1,"_",year2),season)
  plot_dir_exp=plot_dir # file.path(plot_dir,exp,paste0(year1,"_",year2),season)
  dir.create(plot_dir_exp,recursive=T)

  # Years to be considered based on namelist and cfg_file
  years <- year1:year2
  if (ryearplot[1] == "ALL") {
    years <- year1:year2
  } else if (ryearplot[1] == "FIRST") {
    years <- year1
  } else {
    years <- years[match(ryearplot,years)] ; years <- years[!is.na(years)] 
  }
  nyears <- length(years)

  #-----------------Loading data-----------------------#

  # open experiment field
  for (field in field_names) {
    filename<-getfilename.indices(work_dir_exp,diag_base,model_idx,season)
    print(paste("Reading experiment ",filename))
    field_exp=ncdf.opener(filename,field,"lon","lat",rotate="no")
    if (removedesert) { field_exp<-field_exp*retdes3D }
    if (maskSeaLand) { field_exp<-apply.elevation.mask(field_exp,relevation,seaLandElevation)}
    if (rmultiyear_mean) { # if requested calculate multiyear average and store it at time=1
      field_exp[,,1]<-apply(field_exp,c(1,2),mean,na.rm=T)
    } 
    if (highreselevation_only) { field_exp[] <- NA }
    assign(paste(field,"_exp",sep=""),field_exp)        
  }

 
  #---------------Multiyear mean-----#
  if (rmultiyear_mean) { nyears <- 1 }

  #-----------------Producing figures------------------------#

  print(paste0(diag_base,": starting figures"))

  if (boxregion!=0) { nregions = 1 } # boxregion will plot region boxes over a global map of the selected field

  # LOOP over selected regions 
  for (iselregion in 1:nregions) {
    iregion=selregions[iselregion]
    print(paste("region: ",region_names[iregion]))
  
    # Startup graphics for multiple years in one figure
    if (plot_type == 4) {
      field_label=paste(field_names,collapse="-") 
      figname=getfilename.figure(plot_dir_exp,field_label,year1,year2,model_idx,season,
              "multiyear",region_codes[iregion],label,"map",output_file_type) 
      graphics.startup(figname,output_file_type,diag_script_cfg)
      par(mfrow=c(nyears,nfields),cex.main=1.3,cex.axis=1.2,cex.lab=1.2,mar=c(2,2,2,2),oma=c(1,1,1,1))
    }
    # LOOP over years defined in namelist and cfg_file
    for (iyear in c(1:nyears)) {
      if (ryearplot_ref[1] == "EXP") {iyear_ref=iyear} else {iyear_ref<-match(ryearplot_ref,years_ref)}
      time_label=years[iyear]; time_label_ref=years[iyear_ref] ; time_label_fig=time_label
      if (rmultiyear_mean) { time_label=paste(year1,year2,sep="-") ; time_label_ref=paste(year1_ref,year2_ref,sep="-") ; time_label_fig="myearmean"}      
      print(paste0(diag_base, ": plotting data for  ",region_names[iregion]," ",time_label))
    
      #standard properties
      info_exp=paste(exp,time_label)#,season)
      info_ref=paste(dataset_ref,time_label_ref)#,season)

      # Startup graphics for multiple fields/quantities in one figure 
      if (plot_type == 3) {
        field_label=paste(field_names,collapse="-")
        figname=getfilename.figure(plot_dir_exp,field_label,year1,year2,model_idx,
                season,time_label_fig,region_codes[iregion],label,"map",output_file_type)
        graphics.startup(figname,output_file_type,diag_script_cfg)
        par(mfrow=c(nfields,3),cex.main=1.3,cex.axis=1.2,cex.lab=1.2,mar=c(2,2,2,2),oma=c(1,1,1,1))
      }
      # LOOP over fields
      for (field in field_names) {
       # print(paste0("working on ",field))
        ifield<-which(field == field_names)
        if (anyNA(title_unit_m[ifield,1:3])) {title_unit_m[ifield,1:3]<-field;title_unit_m[ifield,4]<-""}

        # get fields
        field_ref=get(paste(field,"_ref",sep=""))
        field_exp=get(paste(field,"_exp",sep=""))
 
        # MAPS: select required year (if requested, multiyear average is stored at iyear=1) 
        field_ref<-field_ref[,,iyear]
        field_exp<-field_exp[,,iyear_ref] 
        tmp.field<-field_exp
 
        # define quantity-dependent properties (exp, ref, exp-ref)
        tmp.colorbar<-c(F,T,T)
        if (plot_type == 1) { tmp.colorbar<-T }
        tmp.palette<-palette_giorgi2011
        if (is.na(levels_m[ifield,1])|is.na(levels_m[ifield,2])) { 
          print("No value for range: assigning min and max")
          tmp.levels<-seq(min(field_ref,na.rm=T),max(field_ref,na.rm=T),len=nlev)
        } else { tmp.levels<-seq(levels_m[ifield,1],levels_m[ifield,2],len=nlev) }
        if (highreselevation_only) {title_unit_m[ifield,1]<-"Elevation" }
        tmp.titles<-paste0(title_unit_m[ifield,1],": ",region_names[iregion]," ",c(info_exp,info_ref,"Difference"))
        if (plot_type == 4) { tmp.titles<-paste(title_unit_m[ifield,1],time_label) }

        # Startup graphics for individual fields and multi quantities in each figure
        if (plot_type == 2) {
          figname=getfilename.figure(plot_dir_exp,field,year1,year2,model_idx,
                  season,time_label_fig,region_codes[iregion],label,"comp_map",output_file_type)
          graphics.startup(figname,output_file_type,diag_script_cfg)
          par(mfrow=c(3,1),cex.main=2,cex.axis=1.5,cex.lab=1.5,mar=c(5,5,4,8),oma=c(1,1,1,1))
        }

        # --- MAPS ----
        # LOOP over quantity (exp,ref,exp-ref difference) to be plotted 
        for (iquantity in c(1:nquantity[plot_type])) {
          if (iquantity == 2) { tmp.field<-field_ref }
          if (iquantity == 3) { 
            tmp.palette<-palette2
            tmp.field<- field_exp-field_ref             
            if (is.na(levels_m[ifield,3])|is.na(levels_m[ifield,4])) { 
              tmp.levels<-seq(min(field_ref,na.rm=T),max(field_ref,na.rm=T),len=nlev) 
            } else { tmp.levels<-seq(levels_m[ifield,3],levels_m[ifield,4],len=nlev) }
          }
          # Startup graphics for individual field in each figure
          if (plot_type == 1) {
            figname=getfilename.figure(plot_dir_exp,field,year1,year2,model_idx,
                    season,time_label_fig,region_codes[iregion],label,"map",output_file_type)
            graphics.startup(figname,output_file_type,diag_script_cfg)
            par(cex.main=2,cex.axis=1.5,cex.lab=1.5,mar=c(5,5,4,8),oma=c(1,1,1,1))
          }

          # set active panel
          if (plot_type == 3) { par(mfg=c(ifield,iquantity,nfields,3)) }
          if (plot_type == 4) { par(mfg=c(iyear,ifield,nyears,nfields)) } 
          # contours
          filled.contour3(ics,ipsilon,tmp.field,xlab="Longitude",ylab="Latitude",
               main=tmp.titles[iquantity],levels=tmp.levels,color.palette=tmp.palette,
               xlim=c(regions[iregion,1],regions[iregion,2]), ylim=c(regions[iregion,3],regions[iregion,4]), axes=F)
          # continents
          continents_col="white"; if (map_continents<=0) {continents_col="gray30"}
          map("world",regions=".",interior=map_continents_regions,exact=F,boundary=T,add=T,col=continents_col,lwd=abs(map_continents))
          # grid points
        
          if (oplot_grid) { 
            # build up grid if needed 
            ics2=replicate(length(ipsilon),ics)
            ipsilon2=t(replicate(length(ics),ipsilon))
            print(str(ics2))
            print(str(ipsilon2))
            print(str(tmp.field))
            points(ics2,ipsilon2,pch=1,col="grey40",cex=oplot_grid) 
          }        
 
          # add highres elevation contours
          if (highreselevation) {
            palette(terrain.colors(10))
            contour(highresel$lon_el,highresel$lat_el,highresel$elevation,
                 levels=seq(500,5000,length.out=10),col=1:10,add=T)
          }
          # boxes
          box(col="grey60") 
          if (boxregion!=0) {
            box_col="white"; if (boxregion<=0) {box_col="grey30"}
            for (ireg in 2:length(selregions)) {  
              iselreg=selregions[ireg] 
              rect(regions[iselreg,1],regions[iselreg,3],regions[iselreg,2],
                   regions[iselreg,4],border=box_col,lwd=abs(boxregion))
              text(regions[iselreg,1],regions[iselreg,3],paste0("         ",
                   region_codes[iselreg]),col=box_col,pos=3,offset=0.5)
            }
          }
          # axis
          if (plot_type <= 2) {
            axis(1,col="grey40") 
            axis(2,col="grey40") 
          } else if (plot_type == 3) {
            if (iquantity == 1) { axis(2,col="grey40") }
            if (ifield == length(field_names)) { axis(1,col="grey40") }
          } else if (plot_type == 4) {
            if (iyear == nyears) { axis(1,col="grey40") }
            if (field == "int_norm") { axis(2,col="grey40") }
          }
          #colorbar
          if ((tmp.colorbar[iquantity]) & add_colorbar) { 
            image.scale3(volcano,levels=tmp.levels,color.palette=tmp.palette,colorbar.label=paste(title_unit_m[ifield,1],title_unit_m[ifield,4]),
                          cex.colorbar=1.0,cex.label=1.0,colorbar.width=1,line.label=legend_distance,line.colorbar=1.0)
          }
        } # close loop over quantity
        if (plot_type <= 2) { graphics.close(figname) }
   
      } # close loop over field 
      if (plot_type == 3) { graphics.close(figname) }
    } # close loop over years
    if (plot_type == 4) { graphics.close(figname) }
  } # close loop over regions
} # close loop over models
} # close function

