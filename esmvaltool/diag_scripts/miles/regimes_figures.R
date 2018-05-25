######################################################
#------Regimes routines figures for MiLES------------#
#-------------P. Davini (May 2017)-------------------#
######################################################

#DECLARING THE FUNCTION: EXECUTION IS AT THE BOTTOM OF THE SCRIPT

miles.regimes.figures<-function(exp,year1,year2,dataset_ref,year1_ref,year2_ref,season,FIGDIR,FILESDIR,REFDIR,CFGSCRIPT,nclusters)
{

if (nclusters!=4 | season!="DJF") {stop("Beta version: unsupported season and/or number of clusters")}

#R configuration file 
source(CFGSCRIPT)

#correct folder to experiment dependent
REGIMESDIR=file.path(FILESDIR,exp,"Regimes",paste0(year1,"_",year2),season)
FIGDIRREGIMES=file.path(FIGDIR,exp,"Regimes",paste0(year1,"_",year2),season)
dir.create(FIGDIRREGIMES,recursive=T)

#check path for reference dataset
#if (dataset_ref=="ERAINTERIM" & year1_ref=="1979" & year2_ref=="2014")
if (REFDIR!=FILESDIR)
        {REFDIR=file.path(REFDIR,"Regimes") } else {REFDIR=file.path(FILESDIR,dataset_ref,"Regimes",paste0(year1_ref,"_",year2_ref),season)}

##########################################################
#-----------------Loading datasets-----------------------#
##########################################################

#loading anomalies and variances of experiment
nomefile=paste0(REGIMESDIR,"/RegimesPattern_",exp,"_",year1,"_",year2,"_",season,".nc")
frequencies_exp=ncdf.opener(nomefile,"Frequencies")
regimes_exp=ncdf.opener(nomefile,"Regimes","Lon","Lat",rotate="no")

#loading reference field
nomefile=paste0(REFDIR,"/RegimesPattern_",dataset_ref,"_",year1_ref,"_",year2_ref,"_",season,".nc")
frequencies_ref=ncdf.opener(nomefile,"Frequencies")
regimes_ref=ncdf.opener(nomefile,"Regimes","Lon","Lat",rotate="no")

#try to assign the 4 standard regimes names to the dataset using the distance between 
#the location of the maximum/minimum of the pattern and 4 "standard" locations
#when something is wrong (i.e. multiple assignments) general "Regime X" names are set
#It is not perfect, it is just aimed at simplying the plots 
for (ii in c("ref","exp"))
{
	compose=get(paste0("regimes_",ii))
	names=paste("Regimes",1:nclusters)
        position=rbind(c(-50,60),c(-40,50),c(0,60),c(-15,60))
        rownames(position)<-c("NAO-","Atlantic Ridge","Scandinavian Blocking","NAO+")
        for (i in 1:nclusters)
                {
                MM=which(compose[,,i]==max(compose[,,i],na.rm=T),arr.ind=T)
                mm=which(compose[,,i]==min(compose[,,i],na.rm=T),arr.ind=T)
                if (max(compose[,,i],na.rm=T)>abs(min(compose[,,i],na.rm=T)))
                        {
                        distMM=dist(rbind(c(ics[MM[1]],ipsilon[MM[2]]),position))
                        } else {
                        distMM=dist(rbind(c(ics[mm[1]],ipsilon[mm[2]]),position))
                        }
                #print(distMM)
                names[i]=rownames(position)[which.min(distMM[1:nclusters])]
		if (i>1 & any(names[i]==names[1:max(c(1,i-1))])) {names[i]=paste("Regime",i)	}
		#print(names[i])
                }

assign(paste0("names_",ii),names)
}

#plot properties
lev_field=seq(-250,250,20)
lev_diff=seq(-150,150,20)
lat_lim=c(20,90)
info_exp=paste(exp,year1,"-",year2,season)
info_ref=paste(dataset_ref,year1_ref,"-",year2_ref,season)

kk0=1
# loop on regimes
for (name in names_ref)
{
	#-----plotting-------#

	#a bit complicated but it is used to compare similar regimes even if they not
	#equale percentage of occurrence.
	ii=which(name==names_exp)
	jj=which(name==names_ref)
	if (length(ii)==0) {ii=which(setdiff(names_exp,names_ref)[kk0]==names_exp); kk0=kk0+1}
	print(name)

        #final plot production
        figname=paste(FIGDIRREGIMES,"/Regime",ii,"_",exp,"_",year1,"_",year2,"_",season,".",output_file_type,sep="")
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

        filled.contour3(ics,ipsilon,regimes_exp[,,ii],xlab="Longitude",ylab="Latitude",main=paste(names_exp[ii],info_exp),levels=lev_field,color.palette=palette3,ylim=lat_lim)
        map("world",regions=".",interior=F,exact=F,boundary=T,add=T)
        text(120,85,paste("Frequencies: ",round(frequencies_exp[ii],2),"%",sep=""),cex=2)

        filled.contour3(ics,ipsilon,regimes_ref[,,jj],xlab="Longitude",ylab="Latitude",main=paste(names_ref[jj],info_ref),levels=lev_field,color.palette=palette3,ylim=lat_lim)
        map("world",regions=".",interior=F,exact=F,boundary=T,add=T)
        image.scale3(volcano,levels=lev_field,color.palette=palette0,colorbar.label="m",cex.colorbar=1.2,cex.label=1.5,colorbar.width=1,line.label=3)
        text(120,85 ,paste("Frequencies: ",round(frequencies_ref[jj],2),"%",sep=""),cex=2)

        #delta field plot
        filled.contour3(ics,ipsilon,regimes_exp[,,ii]-regimes_ref[,,jj],xlab="Longitude",ylab="Latitude",main=paste(name,"Difference"),levels=lev_diff,color.palette=palette2,ylim=lat_lim)
        map("world",regions=".",interior=F,exact=F,boundary=T,add=T)
        image.scale3(volcano,levels=lev_diff,color.palette=palette2,colorbar.label="m",cex.colorbar=1.2,cex.label=1.5,colorbar.width=1,line.label=3)

        dev.off()
        }

}


# REAL EXECUTION OF THE SCRIPT 
# read command line
args <- commandArgs(TRUE)

# number of required arguments from command line
name_args=c("exp","year1","year2","dataset_ref","year1_ref","year2_ref","season","FIGDIR","FILESDIR","REFDIR","CFGSCRIPT","PROGDIR","nclusters")
req_args=length(name_args)

# print error message if uncorrect number of command 
if (length(args)!=0) {
    if (length(args)!=req_args) {
        print(paste("Not enough or too many arguments received: please specify the following",req_args,"arguments:"))
        print(name_args)
    } else {
# when the number of arguments is ok run the function()
        for (k in 1:req_args) {assign(name_args[k],args[k])}
        source(paste0(PROGDIR,"/script/basis_functions.R"))
        miles.regimes.figures(exp,year1,year2,dataset_ref,year1_ref,year2_ref,season,FIGDIR,FILESDIR,REFDIR,CFGSCRIPT,nclusters)
    }
}



