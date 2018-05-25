######################################################
#------Regimes routines computation for MiLES--------#
#-------------P. Davini (May 2017)-------------------#
######################################################

miles.regimes.fast<-function(exp,year1,year2,season,z500filename,FILESDIR,nclusters=nclusters)
{

#t0
t0<-proc.time()

#region boundaries for North Atlantic 
if (nclusters!=4 | season!="DJF") {
	stop("Beta version: unsupported season and/or number of clusters")
}

#test function to smooth seasonal cycle: it does not work fine yet, keep it false
smoothing=F
xlim=c(-80,40)
ylim=c(30,87.5)


#setting up main variables
REGIMESDIR=file.path(FILESDIR,exp,"Regimes",paste0(year1,"_",year2),season)
dir.create(REGIMESDIR,recursive=T)

#setting up time domain
years=year1:year2
timeseason=season2timeseason(season)

#if (smoothing) {
#	timeseason0=timeseason
#	if (season=="DJF") {
#		timeseason=sort(c(timeseason,11,3)) 
#	} else {
#	timeseason=sort(c(timeseason[1]-1,timeseason,timeseason[length(timeseason)]+1))
#	}
#}

#new file opening
nomefile=z500filename
fieldlist=ncdf.opener.time(nomefile,"zg",tmonths=timeseason,tyears=years,rotate="full")

#time array
datas=fieldlist$time
etime=power.date.new(datas)

#declare variable
Z500=fieldlist$field

print("Compute anomalies based on daily mean")
Z500cycle=apply(Z500,c(1,2),ave,etime$month,etime$day)

if (!smoothing) {
	Z500anom=Z500-aperm(Z500cycle,c(2,3,1)) 
}

#if (smoothing) {
#	print("running mean")
#	rundays=5
#	runZ500cycle=apply(Z500cycle,c(2,3),filter,rep(1/rundays,rundays),sides=2)
#	Z500anom0=Z500-aperm(runZ500cycle,c(2,3,1))
#	whichdays=which(as.numeric(format(datas,"%m")) %in% timeseason0)
#	Z500anom=Z500anom0[,,whichdays]
#	etime=power.date.new(datas[whichdays])
#	}

#compute weather regimes
weather_regimes=regimes(ics,ipsilon,Z500anom,ncluster=nclusters,ntime=1000,neof=4,xlim,ylim,alg="Hartigan-Wong")

t1=proc.time()-t0
print(t1)

##########################################################
#------------------------Save to NetCDF------------------#
##########################################################

#saving output to netcdf files
print("saving NetCDF climatologies...")
savefile1=paste(REGIMESDIR,"/RegimesPattern_",exp,"_",year1,"_",year2,"_",season,".nc",sep="")

# dimensions definition
TIME=paste("days since ",year1,"-",timeseason[1],"-01 00:00:00",sep="")
LEVEL=50000
fulltime=as.numeric(etime$data)-as.numeric(etime$data)[1]
x <- ncdim_def( "Lon", "degrees", ics)
x0 <- ncdim_def( "Lon0", "degrees", 0)
y <- ncdim_def( "Lat", "degrees", ipsilon)
y0 <- ncdim_def( "Lat0", "degrees", 0)
z <- ncdim_def( "Lev", "Pa", LEVEL)
cl <- ncdim_def( "Time0", TIME, 1:nclusters)
t <- ncdim_def( "Time", TIME, fulltime,unlim=T)

unit="m"; longvar="Weather Regimes Pattern"
pattern_ncdf=ncvar_def("Regimes",unit,list(x,y,z,cl),-999,longname=longvar,prec="single",compression=1)
unit=paste0("0-",nclusters); longvar="Weather Regimes Cluster Undex"
cluster_ncdf=ncvar_def("Indices",unit,list(x0,y0,z,t),-999,longname=longvar,prec="single",compression=1)
unit="%"; longvar="Weather Regimes Frequencies"
frequencies_ncdf=ncvar_def("Frequencies",unit,list(cl),-999,longname=longvar,prec="single",compression=1)

ncfile1 <- nc_create(savefile1,list(pattern_ncdf,cluster_ncdf,frequencies_ncdf))
ncvar_put(ncfile1, "Regimes", weather_regimes$regimes, start = c(1, 1, 1, 1),  count = c(-1,-1,-1,-1))
ncvar_put(ncfile1, "Indices", weather_regimes$cluster, start = c(1, 1, 1, 1),  count = c(-1,-1,-1,-1))
ncvar_put(ncfile1, "Frequencies", weather_regimes$frequencies, start = c(1),  count = c(-1))
nc_close(ncfile1)

}

# REAL EXECUTION OF THE SCRIPT 
# read command line
args <- commandArgs(TRUE)

# number of required arguments from command line
name_args=c("exp","year1","year2","season","z500filename","FILESDIR","PROGDIR","nclusters")
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
        miles.regimes.fast(exp,year1,year2,season,z500filename,FILESDIR,nclusters)
    }
}


