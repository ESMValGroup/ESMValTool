######################################################
#--------Trend plotting routine for HyInt------------#
#-------------E. Arnone (September 2017)-------------#
######################################################

hyint.plot.trends<-function(work_dir,plot_dir,ref_dir,ref_idx,season) {

# load settings
source('esmvaltool/diag_scripts/aux/hyint/hyint_parameters.r')
for (myname in names(settings)) { temp=get(myname,settings); assign(myname,temp) }



#source('interface_data/r.interface')
var_type = c("tseries","tseries-sd","trend","trend-stat")
# diag_base <- "HyInt"
 
# Number of models
nmodels=length(models_name)

# Define regions to be used
nregions=length(selregions)
if ((plot_type == 13)|(plot_type == 15)) { nregions <- 1 } # if plotting multiple models use only one region at the time (first of list) 

# Define fields to be used (note that the routine is optimized for 6 fields in 3x2 panels per multi-panel figures)
if (selfields[1]!=F) {
  field_names=field_names[selfields,drop=F]
  levels_m=levels_m[selfields,,drop=F]
  tlevels_m=tlevels_m[selfields,,drop=F]
  title_unit_m=title_unit_m[selfields,,drop=F] 
}

# Update number of panels and columns if selfields has one element only 
if (length(selfields) == 1) { 
	npancol=1
	npanrow=1
}

# Define array to store plotting limits for each panel of multi-panel figures 
plot_limits=array(NaN,c(4,length(field_names)))

# Load parameters for reference dataset
dataset_ref=models_name[ref_idx]
model_exp_ref=models_experiment[ref_idx]
model_ens_ref=models_ensemble[ref_idx]
year1_ref=models_start_year[ref_idx]
year2_ref=models_end_year[ref_idx]
plot_dir_ref=plot_dir #file.path(plot_dir,model_exp_ref,paste0(year1_ref,"_",year2_ref),season)
dir.create(plot_dir_ref,recursive=T)

# Handle label tag when overplotting data from tseries files with different labels in plot_type 14,15,16
label_figname=label[1]
if (length(label)>1 & plot_type>=10 ) { label_figname=paste0(label[1],"-plus") }

# Startup graphics for multi-model timeseries
if ((plot_type == 13)|(plot_type == 15)) {
  field_label=paste(field_names,collapse="-")
  tseries_trend_tag="timeseries"
  if (plot_type == 15) { tseries_trend_tag="trend_summary" }
  figname=getfilename.figure(plot_dir_ref,field_label,year1_ref,year2_ref,ref_idx,season
          ,"",region_codes[selregions[1]],label_figname,tseries_trend_tag,output_file_type,multimodel=T)
  graphics.startup(figname,output_file_type,diag_script_cfg)
  par(mfrow=c(npanrow,npancol),cex.main=1.3,cex.axis=1.2,cex.lab=1.2,mar=c(5,5,5,5),oma=c(1,1,1,1))
}

# Loop over models
for (model_idx in 1:nmodels) {
  # setting up path and parameters
  exp    <- models_name[model_idx]
  model_exp <- models_experiment[model_idx]
  model_ens <- models_ensemble[model_idx]
  year1  <- models_start_year[model_idx]
  year2  <- models_end_year[model_idx]

  # set main paths
  work_dir_exp=work_dir # file.path(work_dir,exp,paste0(year1,"_",year2),season)
  plot_dir_exp=plot_dir # file.path(plot_dir,exp,paste0(year1,"_",year2),season)
  dir.create(plot_dir_exp,recursive=T)

#  # check path to reference dataset
   if (!exists("ref_dir")) { ref_dir = work_dir}
   if (file.exists(ref_dir) == F) { stop(paste(diag_base,": reference file path not found")) }
#  if (ref_dir!=work_dir) {
#    ref_dir=file.path(ref_dir,diag_base)
#  } else {
#    ref_dir=file.path(work_dir,dataset_ref,paste0(year1_ref,"_",year2_ref),season)
#  }

  # Years to be considered based on namelist and cfg_file
  years <- year1:year2
  years_ref <- year1_ref:year2_ref
  if (ryearplot[1] == "ALL") {
    years <- year1:year2
  } else if (ryearplot[1] == "FIRST") {
    years <- year1
  } else {
    years <- years[match(ryearplot,years)] ; years <- years[!is.na(years)] 
  }
  nyears <- length(years)
  if (plot_type >= 14) { add_trend<-F } # do not plot trend line for plot 14 or 15

  # Startup graphics for multi-region timeseries
  if (plot_type == 12) {
    field_label=paste(field_names,collapse="-")
    figname=getfilename.figure(plot_dir_exp,field_label,year1,year2,model_idx,season,"","regions",label_figname,"timeseries",output_file_type) 
    graphics.startup(figname,output_file_type,diag_script_cfg)
    par(mfrow=c(npanrow,npancol),cex.main=1.3,cex.axis=1.2,cex.lab=1.2,mar=c(5,5,5,5),oma=c(1,1,1,1))
  }
  # Startup graphics for bar plot of trend coefficients 
  if (plot_type == 14) {
    field_label=paste(field_names,collapse="-")
    figname=getfilename.figure(plot_dir_exp,field_label,year1,year2,model_idx,season,"","regions",label_figname,"trend_summary",output_file_type)
    graphics.startup(figname,output_file_type,diag_script_cfg)
    par(mfrow=c(npanrow,npancol),cex.main=1.3,cex.axis=1.2,cex.lab=1.2,mar=c(8,8,2,2),oma=c(1,1,1,1))
  }

 if (model_idx == 1) {store_label=label} 
  #------------ Loop over label when plotting more timeseries files in the same panel ----------
  for (ilabel in 1:length(store_label)) {
    label = store_label[ilabel]
  #-----------------Loading data-----------------------#

  # open timeseries and trends for exp and ref
  infile<-getfilename.trends(work_dir_exp,label,model_idx,season)
  print(paste("HyInt_trends: reading file ",infile))
  field_long_names<-array(NaN,length(field_names))
  field_units<-array(NaN,length(field_names))

  for (var in field_names) {
    ivar<-which(field_names==var)
    for (stype in var_type[1:2]) { 
      svar=paste0(var,"_",stype)
      rfield=ncdf.opener(infile,svar,"region",timedimname,rotate="no")
      assign(svar,rfield) # assign field data to field name       
      nc<-nc_open(infile)
      dlname <- ncatt_get(nc,svar,"long_name")
      dunits <- ncatt_get(nc,svar,"units")
      field_long_names[ivar]<-dlname$value      
      field_units[ivar]<-dunits$value
      nc_close(nc)
    } 
    for (stype in var_type[3:4]) {
      svar=paste0(var,"_",stype)
      rfield=ncdf.opener(infile,svar,"region","coefficients",rotate="no")
      assign(svar,rfield) # assign field data to field name        
    } 
  }

  # store size of time and region arrays
  time<-ncdf.opener(infile,timedimname,timedimname,rotate="no") + 1950
  regions<-ncdf.opener(infile,"regions","region","boundaries",rotate="no")

  #-----------------Producing figures------------------------#

  print(paste0(diag_base,": starting figures"))
  
  # LOOP over fields
  for (field in field_names) {
    ifield<-which(field == field_names)
    if (anyNA(title_unit_m[ifield,1:3])) { 
      title_unit_m[ifield,1]<-field_names[ifield]
      title_unit_m[ifield,2:3]<-field_long_names[ifield]
      title_unit_m[ifield,4]<-field_units[ifield]
    }

    # TIMESERIES: get timeseries and trends     
    tfield_exp<-get(paste0(field,"_",var_type[1]))
    tfield_exp_sd<-get(paste0(field,"_",var_type[2]))
    trend_exp<-get(paste0(field,"_",var_type[3]))
    trend_exp_stat<-get(paste0(field,"_",var_type[4]))

    if (length(dim(tfield_exp))<2) { # reshape data to matrix if regions has only one element
      tfield_exp<-array(tfield_exp,c(1,length(tfield_exp)))
      tfield_exp_sd<-array(tfield_exp_sd,c(1,length(tfield_exp_sd)))
      trend_exp<-array(trend_exp,c(1,length(trend_exp)))
      trend_exp_stat<-array(trend_exp_stat,c(1,length(trend_exp_stat)))
    }
    if (is.na(levels_m[ifield,1])|is.na(levels_m[ifield,2])) {
      print("No value for range: assigning min and max")
      tmp.levels<-seq(min(tfield_exp,na.rm=T),max(tfield_exp,na.rm=T),len=nlev)
    } else { tmp.levels<-seq(levels_m[ifield,1],levels_m[ifield,2],len=nlev) }
    # setup time array 
    rettimes = which(!is.na(time))
    if (trend_years[1] != F) { # apply trend to limited time interval if required
      rettimes = which((time >= trend_years[1]) & time <= trend_years[2])          
      if (length(trend_years) == 4) { # apply trend also to second time interval if required
        rettimes2 = which((time >= trend_years[3]) & time <= trend_years[4])
      } 
    }
    xlim=c(min(time),max(time))
    if (trend_years_only&(trend_years[1] != F)) { xlim=trend_years[1:2] }
   
    # Startup graphics for one timeseries in one figure 
    if (plot_type == 11) {
      figname=getfilename.figure(plot_dir_exp,field,year1,year2,model_idx,season,
              "",region_codes[iregion],label_figname,"timeseries_single",output_file_type)
      graphics.startup(figname,output_file_type,diag_script_cfg)
      par(cex.main=1.3,cex.axis=1.2,cex.lab=1.2,mar=c(4,4,2,2),oma=c(1,1,1,1))
    } 

    # Actual plotting
    if ((plot_type == 11)|(plot_type == 12)|(plot_type == 13)) {

      # set active panel
      par_row=(ifield-1)%/%npancol+1
      par_col=(ifield-1)%%npancol+1
      par(mfg=c(par_row,par_col,npanrow,npancol)) 

      # Base plot
      if (!(plot_type == 13 & model_idx >1)&ilabel==1) {
        ylab=paste0(title_unit_m[ifield,1])
        if (title_unit_m[ifield,4]!="") { ylab=paste0(ylab,"(",title_unit_m[ifield,4],")")}  
        # plot(time,tfield_exp[1,],ylim=c(tmp.levels[1],tmp.levels[length(tmp.levels)]),xlim=xlim,
        plot(time,type="n",ylim=c(tmp.levels[1],tmp.levels[length(tmp.levels)]),xlim=xlim,
             xlab="Year",ylab=ylab,main=title_unit_m[ifield,3])               
        # store panel plot limits
        plot_limits[,ifield]=par("usr")
      }

      # Update plot limits in case panel has changed
      par(usr=plot_limits[,ifield])

      # LOOP over regions to plot timeseries 
      if (add_trend_sd_shade) {
        for (ireg in 1:nregions) { 
          iselreg=selregions[ireg] 
          shade_area=c(tfield_exp[ireg,]+tfield_exp_sd[ireg,],rev(tfield_exp[ireg,]-tfield_exp_sd[ireg,]))
          shade_area[shade_area < tmp.levels[1]] <- tmp.levels[1]            
          polygon(c(time,rev(time)),shade_area,col="grey95", border=NA)
        }
      }
      for (ireg in 1:nregions) { 
        iselreg=selregions[ireg]  
        col_ts=ireg
	if (length(label) > 1) {col_ts=c("dodgerblue4","darkseagreen4","goldenrod4","coral4","grey","mediumorchid1","black")[ilabel]} 
        if (plot_type == 13) { col_ts=model_idx } 
        if (add_trend_sd) {
          lines(time,tfield_exp[ireg,]+tfield_exp_sd[ireg,],lty=3,col=col_ts)
          lines(time,tfield_exp[ireg,]-tfield_exp_sd[ireg,],lty=3,col=col_ts)       
        }
        if (add_tseries_lines) { lines(time,tfield_exp[ireg,],col=col_ts) }          
        points(time,tfield_exp[ireg,],col=col_ts)
        if (add_trend) {
          lines(time[rettimes],trend_exp[ireg,1]+trend_exp[ireg,2]*time[rettimes],col=col_ts,lwd=2)  
          if (length(trend_years) == 4) { # apply trend also to second time interval if required
            lines(time[rettimes2],trend_exp[ireg,3]+trend_exp[ireg,4]*time[rettimes2],col=col_ts,lwd=2) 
          } 
        }
      }
      if (abs(add_legend)&((plot_type == 11)|(plot_type == 12))&(ifield==1)) {
        pos_legend=c(plot_limits[1,ifield]+(plot_limits[2,ifield]-plot_limits[1,ifield])*xy_legend[1],
                     plot_limits[3,ifield]+(plot_limits[4,ifield]-plot_limits[3,ifield])*xy_legend[2])
        ncol=1       
        #    text((xlim[1]+(xlim[2]-xlim[1])*ireg/nregions),
        #         tmp.levels[1],region_codes[iselreg],col=col_ts,offset=0.5)
        if (add_legend < 0) { ncol=nregions }
        if (add_legend > 1) { ncol=add_legend }
        legend(pos_legend[1],pos_legend[2],region_codes[selregions],text.col=(1:nregions),ncol=ncol)
      }          
      if (plot_type == 11) { graphics.close(figname) }
    }
    if ((plot_type == 14)|(plot_type == 15)) { # plot trend coefficients for different regions, one panel per field
      #ylim=c(min(trend_exp[,2]-trend_exp_stat[,2]),max(trend_exp[,2]+trend_exp_stat[,2])) 

      if (anyNA(tlevels_m[ifield,])) {
        print("No value for range: assigning min and max")
        ylim<-c(min(trend_exp[,2]-trend_exp_stat[,2],na.rm=T),max(trend_exp[,2]+trend_exp_stat[,2],na.rm=T))
      } else { ylim=tlevels_m[ifield,] }
  
      if (trend_years[1]!=F) { xlim=trend_years[1:2] }
      ylab=paste0("Avg trend")
      # change y scale to % and 1/100 years          
      if (scalepercent & (field != "hyint")) {
        trend_exp=trend_exp*100   # trend coefficients 
        trend_exp_stat[,2]=trend_exp_stat[,2]*100 # standard error
        ylab=paste0(ylab," (%)")
        ylim=ylim*100
      }
      if (scale100years) {
        trend_exp=trend_exp*100   # trend coefficients 
        trend_exp_stat[,2]=trend_exp_stat[,2]*100 # standard error
        ylab=paste0(ylab," (1/100 years)")
        ylim=ylim*100
      }
      nx=nregions
      xlab="Regions"
      xlabels=region_codes[selregions]
      if (plot_type==15) { 
        nx=nmodels 
        xlab=""#"Models"
        xlabels=models_name
      }
      xregions=1:nx # hereafter xregions is the x which also holds models for plot_type 15

      # Actual plotting
      # set active panel
      par_row=(ifield-1)%/%npancol+1
      par_col=(ifield-1)%%npancol+1
      par(mfg=c(par_row,par_col,npanrow,npancol))

      # Base plot
      if (!(plot_type == 15 & model_idx >1)&ilabel==1) {
        #plot(xregions,trend_exp[,2],type="n",pch=22,axes=F,xlab=xlab,ylab=ylab,
        plot(xregions,xregions,type="n",pch=22,axes=F,xlab=xlab,ylab=ylab,
            ylim=ylim, main=(paste0(title_unit_m[ifield,1]," trend (",xlim[1],"-",xlim[2],")")))          
        box()         
        # store panel plot limits
        plot_limits[,ifield]=par("usr")
      }

      # Update plot limits in case panel has changed
      par(usr=plot_limits[,ifield])
      for (ireg in 1:nregions) {
        iregion=selregions[ireg]
        ixregion = iregion
        if (plot_type==15) { ixregion = model_idx }
        # add errorbar (standard error) 
	if (!anyNA(trend_exp_stat[iregion,])) {
          arrows(xregions[ixregion], trend_exp[iregion,2]-trend_exp_stat[iregion,2], xregions[ixregion], trend_exp[iregion,2]+trend_exp_stat[iregion,2], 
                 length=0.05, angle=90, code=3)
          points(xregions[ixregion], trend_exp[iregion,2], pch=22, col="grey40", bg="white",cex=2)
          # add filled points for significant (95% level)
          col90=c("dodgerblue3","darkseagreen3","goldenrod3","coral3","grey","mediumorchid1","black")
          col95=c("dodgerblue4","darkseagreen4","goldenrod4","coral4","grey","mediumorchid1","black")
          if (trend_exp_stat[iregion,4]<=0.1) { points(xregions[ixregion], trend_exp[iregion,2], pch=22, col=col90[ilabel], bg=col90[ilabel],cex=2)}
          points(xregions[ixregion], trend_exp[iregion,2], pch=22, col=col95[ilabel], bg=col95[ilabel],cex=2)
          if (trend_exp_stat[iregion,4]<=0.05) { points(xregions[ixregion], trend_exp[iregion,2], pch=22, col=col95[ilabel], bg=col95[ilabel],cex=2)}
        } else {
          print(paste("MISSING VALUES in index ",field,", region ",region_codes[iregion]))
          print(trend_exp_stat[iregion,])
        }
        # retsig90=which(trend_exp_stat[,4]<0.1)
        # if (!is.na(retsig90[1])) { points(xregions[retsig90], trend_exp[retsig90,2], pch=22, col="grey70", bg="grey70",cex=2) }
        # retsig95=which(trend_exp_stat[,4]<0.05)
        # if (!is.na(retsig95[1])) { points(xregions[retsig95], trend_exp[retsig95,2], pch=22, col="dodgerblue3", bg="dodgerblue3",cex=2) }

      }
      box()         
      if (!((plot_type == 15)&(model_idx>1))){
        if (add_zeroline&(ylim[1]!=0)) { lines(c(-1,nx+1),c(0,0),lty=2,lwd=1.5,col="grey40") }  
        las = 1
        cex.axis = 1 
        if (plot_type == 15) { 
          las = 2 
          cex.axis = 0.8
        }
        axis(1,labels=xlabels,at=xregions,las=las,cex.axis=cex.axis)
        axis(2)                  
      }
    } 
  } # close loop over field 

  } # close loop over label
  if ((plot_type == 12)|(plot_type == 14)) { graphics.close(figname) }
} # close loop over model

# Legend for plot_type 13
if (abs(add_legend)&(plot_type==13)) {   
  ncol=1
  if (add_legend > 1) { ncol=add_legend }
  if (add_legend < 0) { ncol=nmodels }
#  for (ifield in 1:length(field_names)) {  
    ifield=1
    # set active panel
    par_row=(ifield-1)%/%npancol+1
    par_col=(ifield-1)%%npancol+1
    par(mfg=c(par_row,par_col,npanrow,npancol), usr=plot_limits[,ifield])
    pos_legend=c(plot_limits[1,ifield]+(plot_limits[2,ifield]-plot_limits[1,ifield])*xy_legend[1],
                 plot_limits[3,ifield]+(plot_limits[4,ifield]-plot_limits[3,ifield])*xy_legend[2])
              #    text((xlim[1]+(xlim[2]-xlim[1])*model_idx/nmodels),tmp.levels[1],
              #         paste(exp,model_exp),col=col_ts,offset=0.5)
    legend_label=""
    if (tag_legend[1]) legend_label=models_name
    if (tag_legend[2]) legend_label=paste(legend_label,models_experiments,sep=" ")
    if (tag_legend[3]) legend_label=paste(legend_label,models_ensemble,sep=" ")
    legend(pos_legend[1],pos_legend[2],legend=legend_label,text.col=(1:nmodels),ncol=ncol,cex=0.9)
print(legend_label)
print("legend_label")
#  }
} 
if ((plot_type == 13)|(plot_type == 15)) { graphics.close(figname) }
} # close function
