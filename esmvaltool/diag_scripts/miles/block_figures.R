######################################################
#------Blocking routines plotting for MiLES----------#
#-------------P. Davini (May 2017)-------------------#
######################################################

#DECLARING THE FUNCTION: EXECUTION IS AT THE BOTTOM OF THE SCRIPT
miles.block.figures<-function(dataset,expid,ens,year1,year2,
			      dataset_ref,expid_ref,ens_ref,year1_ref,year2_ref,
			      season,FIGDIR,FILESDIR,REFDIR,CFGSCRIPT,PROGDIR) {

#source function scripts
source(file.path(PROGDIR,"script/basis_functions.R"))

#figures configuration files
source(CFGSCRIPT)

#which fieds to load/plot
fieldlist=c("InstBlock","ExtraBlock","Z500","MGI","BI","CN","ACN","BlockEvents","LongBlockEvents","DurationEvents","NumberEvents","TM90")

##########################################################
#-----------------Loading datasets-----------------------#
##########################################################

#open field
for (field in fieldlist) {	

	#use file.builder function
	nomefile=file.builder(FILESDIR,"Block","BlockClim",dataset,expid,ens,year1,year2,season)
	field_exp=ncdf.opener(nomefile,namevar=field,rotate="no")
	assign(paste(field,"_exp",sep=""),field_exp)
}

#open reference field
for (field in fieldlist) {
	
	# check for REFDIR==FILESDIR, i.e. if we are using the climatology provided by MiLES or another dataset MiLES-generated	
	if (REFDIR!=FILESDIR) {
		nomefile_ref=paste0(file.path(REFDIR,"Block"),"/BlockClim_",dataset_ref,"_",year1_ref,"_",year2_ref,"_",season,".nc")
	} else { 
		
		#use file.builder to create the path of the blocking files
		nomefile_ref=file.builder(FILESDIR,"Block","BlockClim",dataset_ref,expid_ref,ens_ref,year1_ref,year2_ref,season)
	}
    
	field_ref=ncdf.opener(nomefile_ref,namevar=field,rotate="no")
    	assign(paste(field,"_ref",sep=""),field_ref)
}

##########################################################
#-----------------Produce figures------------------------#
##########################################################

#standard properties
info_exp=info.builder(dataset,expid,ens,year1,year2,season)
info_ref=info.builder(dataset_ref,expid_ref,ens_ref,year1_ref,year2_ref,season)

#loop on fields
for (field in fieldlist) {
	
	#define field-dependent properties
	fp=field.details(field)
	
	#get fields
    	field_ref=get(paste(field,"_ref",sep=""))
    	field_exp=get(paste(field,"_exp",sep=""))

    	#create figure names with ad-hoc function
    	figname=fig.builder(FIGDIR,"Block",field,dataset,expid,ens,year1,year2,season,output_file_type)
    	print(figname)

	#special treatment for TM90: it is a 1D field!
	if (field=="TM90") {

        	open.plot.device(figname,output_file_type,CFGSCRIPT,special=TRUE)

		#panels option
	        par(cex.main=2,cex.axis=1.5,cex.lab=1.5,mar=c(5,5,4,3),oma=c(0,0,0,0))

		#rotation to simplify the view (90 deg to the west)
		n=(-length(ics)/4)
		ics2=c(tail(ics,n),head(ics,-n)+360)
		field_exp2=c(tail(field_exp,n),head(field_exp,-n))
		field_ref2=c(tail(field_ref,n),head(field_ref,-n))

		#plot properties
		lwdline=4
		tm90cols=fp$color_field
		plot(ics2,field_exp2,type="l",lwd=lwdline,ylim=fp$lev_field,main=fp$title_name,xlab="Longitude",ylab=fp$legend_unit,col=tm90cols[1])
		points(ics2,field_ref2,type="l",lwd=lwdline,lty=1,col=tm90cols[2])
		grid()
        	legend(100,30,legend=c(info_exp,info_ref),lwd=lwdline,lty=c(1,1),col=tm90cols,bg="white",cex=1.5)
	
		#par(new=TRUE)	
		#plot(ics2,field_exp2,type="n",ylim=c(0,90),xlab="",ylab="",axes=F)
		#map("world",regions=".",interior=F,exact=F,boundary=T,add=T,ylim=c(40,80))
		
		dev.off()

		#skip other part of the script
		next()

	}
	
	# Chose output format for figure - by JvH
    	open.plot.device(figname,output_file_type,CFGSCRIPT)

	#plot options	
	par(plotpar)

	#main experiment plot
	im=plot.prepare(ics,ipsilon,field_exp,proj=map_projection,lat_lim=lat_lim)
	filled.contour3(im$x,im$y,im$z,xlab=im$xlab,ylab=im$ylab,main=paste(info_exp),levels=fp$lev_field,color.palette=fp$color_field,xlim=im$xlim,ylim=im$ylim,axes=im$axes)
	mtext(fp$title_name,side=3,line=.5,outer=TRUE,cex=2,font=2)
	proj.addland(proj=map_projection)

	#reference field plot
	im=plot.prepare(ics,ipsilon,field_ref,proj=map_projection,lat_lim=lat_lim)
	filled.contour3(im$x,im$y,im$z,xlab=im$xlab,ylab=im$ylab,main=paste(info_ref),levels=fp$lev_field,color.palette=fp$color_field,xlim=im$xlim,ylim=im$ylim,axes=im$axes)
	proj.addland(proj=map_projection)
	image.scale3(volcano,levels=fp$lev_field,color.palette=fp$color_field,colorbar.label=fp$legend_unit,cex.colorbar=imgscl_colorbar,cex.label=imgscl_label,colorbar.width=1*af,line.label=fp$legend_distance)

	#delta field plot
	im=plot.prepare(ics,ipsilon,field_exp-field_ref,proj=map_projection,lat_lim=lat_lim)
	filled.contour3(im$x,im$y,im$z,xlab=im$xlab,ylab=im$ylab,main=paste("Difference"),levels=fp$lev_diff,color.palette=fp$color_diff,xlim=im$xlim,ylim=im$ylim,axes=im$axes)
	proj.addland(proj=map_projection)
	image.scale3(volcano,levels=fp$lev_diff,color.palette=fp$color_diff,colorbar.label=fp$legend_unit,cex.colorbar=imgscl_colorbar,cex.label=imgscl_label,colorbar.width=1*af,line.label=fp$legend_distance)

	dev.off()
	}


}

#blank line
cat("\n\n\n")

# REAL EXECUTION OF THE SCRIPT 
# read command line
args <- commandArgs(TRUE)

# number of required arguments from command line
name_args=c("dataset","expid","ens","year1","year2","dataset_ref","expid_ref","ens_ref","year1_ref","year2_ref","season","FIGDIR","FILESDIR","REFDIR","CFGSCRIPT","PROGDIR")
req_args=length(name_args)

# print error message if uncorrect number of command 
if (length(args)!=0) {
    if (length(args)!=req_args) {
    	print(paste("Not enough or too many arguments received: please specify the following",req_args,"arguments:"))
		print(name_args)
    } else {
# when the number of arguments is ok run the function()
	for (k in 1:req_args) {assign(name_args[k],args[k])}
	miles.block.figures(dataset,expid,ens,year1,year2,
			    dataset_ref,expid_ref,ens_ref,year1_ref,year2_ref,
			    season,FIGDIR,FILESDIR,REFDIR,CFGSCRIPT,PROGDIR) 
    }
}


