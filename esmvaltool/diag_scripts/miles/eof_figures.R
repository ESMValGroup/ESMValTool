######################################################
#--------Routines for EOFs plotting for MiLES--------#
#-------------P. Davini (May 2017)-------------------#
######################################################

#DECLARING THE FUNCTION: EXECUTION IS AT THE BOTTOM OF THE SCRIPT

miles.eof.figures<-function(exp,year1,year2,dataset_ref,year1_ref,year2_ref,season,FIGDIR,FILESDIR,REFDIR,CFGSCRIPT,tele)
{

#R configuration file 
source(CFGSCRIPT)

#correct folder to experiment dependent
EOFDIR=file.path(FILESDIR,exp,"EOFs",tele,paste0(year1,"_",year2),season)
FIGDIREOF=file.path(FIGDIR,exp,"EOFs",tele,paste0(year1,"_",year2),season)
dir.create(FIGDIREOF,recursive=T)

#check path for reference dataset
#if (dataset_ref=="ERAINTERIM" & year1_ref=="1979" & year2_ref=="2014")
if (REFDIR!=FILESDIR) 
 	{REFDIR=file.path(REFDIR,"EOFs",tele,season) } else {REFDIR=file.path(FILESDIR,dataset_ref,"EOFs",tele,paste0(year1_ref,"_",year2_ref),season)}


#EOFs to plot (depends on how many computed by CDO!)
neofs=4

##########################################################
#-----------------Loading datasets-----------------------#
##########################################################

#loading anomalies and variances of experiment
nomefile=paste0(EOFDIR,"/Z500_monthly_anomalies_",exp,"_",year1,"_",year2,"_",season,".nc")
anomalies_exp=ncdf.opener(nomefile,"zg","lon","lat",rotate="full")
nomefile=paste0(EOFDIR,"/",tele,"_Z500_eigenvalues_",exp,"_",year1,"_",year2,"_",season,".nc")
variance=ncdf.opener(nomefile,"zg")
variance_exp=round(variance[1:neofs]/sum(variance)*100,1)

#loading reference field
nomefile=paste0(REFDIR,"/Z500_monthly_anomalies_",dataset_ref,"_",year1_ref,"_",year2_ref,"_",season,".nc")
anomalies_ref=ncdf.opener(nomefile,"zg","lon","lat",rotate="full")
nomefile=paste0(REFDIR,"/",tele,"_Z500_eigenvalues_",dataset_ref,"_",year1_ref,"_",year2_ref,"_",season,".nc")
variance=ncdf.opener(nomefile,"zg")
variance_ref=round(variance[1:neofs]/sum(variance)*100,1)


##########################################################
#-----------------Produce figures------------------------#
##########################################################

#plot properties
info_exp=paste(exp,year1,"-",year2,season)
info_ref=paste(dataset_ref,year1_ref,"-",year2_ref,season)
lat_lim=c(20,90)
lev_field=seq(-150,150,20)
lev_diff=seq(-95,95,10)
nlev_field=length(lev_field)-1
nlev_diff=length(lev_diff)-1

#loop on number of EOFs
for (neof in 1:neofs) {

	#loading PCs of experiment and normalize
	nomefile=paste(EOFDIR,"/",tele,"_monthly_timeseries_",exp,"_",year1,"_",year2,"_",season,"_0000",neof-1,".nc",sep="")
	timeseries_exp0=ncdf.opener(nomefile,"zg")
	timeseries_exp=standardize(timeseries_exp0)

	#linear regression on Z500 anomalies for experiment (faster function)
	#linear_exp=apply(anomalies_exp,c(1,2),function(linreg) lm(linreg ~ timeseries_exp,na.action=na.exclude)$coef[2])
	linear_exp=apply(anomalies_exp,c(1,2),function(linreg) lin.fit(as.matrix(timeseries_exp,ncol=1),linreg)$coefficients)
	
	#loading PC of reference and normalize
	nomefile=paste(REFDIR,"/",tele,"_monthly_timeseries_",dataset_ref,"_",year1_ref,"_",year2_ref,"_",season,"_0000",neof-1,".nc",sep="")
	timeseries_ref0=ncdf.opener(nomefile,"zg")
	timeseries_ref=standardize(timeseries_ref0)

	#linear regression on Z500 anomalies for reference (faster lm.fit function)
	#linear_ref=apply(anomalies_ref,c(1,2),function(linreg) lm(linreg ~ timeseries_ref,na.action=na.exclude)$coef[2])
	linear_ref=apply(anomalies_ref,c(1,2),function(linreg) lin.fit(as.matrix(timeseries_ref,ncol=1),linreg)$coefficients)

	#check and flip signs (to be in agreement with reference field)
	if (cor(c(linear_ref),c(linear_exp))<0) {linear_exp=(-linear_exp)}
	
	#-----plotting-------#
	
	#plot properties
	if (tele=="NAO") {region="North Atlantic"}
	if (tele=="AO") {region="Northern Hemisphere"}
	title_name=paste0(region," EOF",neof)

	#final plot production
	figname=paste0(FIGDIREOF,"/EOF",neof,"_",exp,"_",year1,"_",year2,"_",season,".",output_file_type)
	print(figname)

	# Chose output format for figure - by JvH
        if (tolower(output_file_type) == "png") {
           png(filename = figname, width=png_width, height=png_height)
        } else if (tolower(output_file_type) == "pdf") {
            pdf(file=figname,width=pdf_width,height=pdf_height,onefile=T)
        } else if (tolower(output_file_type) == "eps") {
            setEPS(width=pdf_width,height=pdf_height,onefile=T,paper="special")
            postscript(figname)
        }

	#plot properties
	par(mfrow=c(3,1),cex.main=2,cex.axis=1.5,cex.lab=1.5,mar=c(5,5,4,8),oma=c(1,1,1,1))

	filled.contour3(ics,ipsilon,linear_exp,xlab="Longitude",ylab="Latitude",main=paste(title_name,info_exp),levels=lev_field,color.palette=palette3,ylim=lat_lim)
	map("world",regions=".",interior=F,exact=F,boundary=T,add=T)
	text(120,85,paste("Variance Explained: ",variance_exp[neof],"%",sep=""),cex=2)

	filled.contour3(ics,ipsilon,linear_ref,xlab="Longitude",ylab="Latitude",main=paste(title_name,info_ref),levels=lev_field,color.palette=palette3,ylim=lat_lim)
	map("world",regions=".",interior=F,exact=F,boundary=T,add=T)
	image.scale3(volcano,levels=lev_field,color.palette=palette0,colorbar.label="m",cex.colorbar=1.2,cex.label=1.5,colorbar.width=1,line.label=3)
	text(120,85 ,paste("Variance Explained: ",variance_ref[neof],"%",sep=""),cex=2)

	#delta field plot
	filled.contour3(ics,ipsilon,linear_exp-linear_ref,xlab="Longitude",ylab="Latitude",main=paste(title_name,"Difference"),levels=lev_diff,color.palette=palette2,ylim=lat_lim)
	map("world",regions=".",interior=F,exact=F,boundary=T,add=T)
	image.scale3(volcano,levels=lev_diff,color.palette=palette2,colorbar.label="m",cex.colorbar=1.2,cex.label=1.5,colorbar.width=1,line.label=3)
	
	dev.off()
	}

}

# REAL EXECUTION OF THE SCRIPT 
# read command line
args <- commandArgs(TRUE)

# number of required arguments from command line
name_args=c("exp","year1","year2","dataset_ref","year1_ref","year2_ref","season","FIGDIR","FILESDIR","REFDIR","CFGSCRIPT","PROGDIR","tele")
req_args=length(name_args)

# print error message if uncorrect number of command 
if (length(args)!=0) {
    if (length(args)!=req_args) {
        print(paste("Not enough or too many arguments received: please specify the following",req_args,"arguments:"))
        print(paste("If running from bash, please specify the",req_args,"arguments here below:"))
        print(name_args)
     } else {
# when the number of arguments is ok run the function()
        for (k in 1:req_args) {assign(name_args[k],args[k])}
        source(paste0(PROGDIR,"/script/basis_functions.R"))
        miles.eof.figures(exp,year1,year2,dataset_ref,year1_ref,year2_ref,season,FIGDIR,FILESDIR,REFDIR,CFGSCRIPT,tele)
     }
}


