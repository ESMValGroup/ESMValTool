#basis functions


##########################################################
#------------------------Packages------------------------#
##########################################################

#loadin packages
library("maps")
library("ncdf4")
library("PCICt")

#check if fast linear fit is operative (after R 3.1): 3x faster than lm.fit, 36x faster than lm
if (exists(".lm.fit")) {lin.fit=.lm.fit} else {lin.fit=lm.fit}

##########################################################
#-----------------Basic functions------------------------#
##########################################################

#normalize a time series
standardize<-function(timeseries)
{
	out=(timeseries-mean(timeseries,na.rm=T))/sd(timeseries,na.rm=T)
	return(out)
}


#detect ics ipsilon lat-lon
whicher<-function(axis,number)
{
	out=which.min(abs(axis-number))
	return(out)
}

#produce a 2d matrix of area weight
area.weight<-function(ics,ipsilon,root=T)
{
field=array(NA,dim=c(length(ics),length(ipsilon)))
if (root==T)
{
        for (j in 1:length(ipsilon))
        {field[,j]=sqrt(cos(pi/180*ipsilon[j]))}
}

if (root==F)
{
        for (j in 1:length(ipsilon))
        {field[,j]=cos(pi/180*ipsilon[j])}
}
return(field)

}

##########################################################
#--------------Time Based functions----------------------#
##########################################################

# to convert season charname to months number
season2timeseason<-function(season)
{
	if (season=="ALL")  {timeseason=1:12}
	if (season=="JJA")  {timeseason=6:8}
	if (season=="DJF")  {timeseason=c(1,2,12)}
	if (season=="MAM")  {timeseason=3:5}
	if (season=="SON")  {timeseason=9:11}
	if (!exists("timeseason")) {stop("wrong season selected!")}
	return(timeseason)
}

#leap year treu/false function
is.leapyear=function(year)
{
	return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
}

power.date.new<-function(datas)
{
whichdays=as.numeric(format(datas,"%m"))
#create a "season" for continuous time, used by persistance tracking
seas=whichdays*1; ss=1
for (i in 1:(length(whichdays)-1))
       {
       if (diff(whichdays)[i]>1)  {ss=ss+1}
       seas[i+1]=ss
       }

etime=list(day=as.numeric(format(datas,"%d")),month=as.numeric(format(datas,"%m")),year=as.numeric(format(datas,"%Y")),data=datas,season=seas)
print("Time Array Built")
print(paste("Length:",length(seas)))
print(paste("From",datas[1],"to",datas[length(seas)]))
return(etime)
}

power.date<-function(season,ANNO1,ANNO2)
{
#evalute the number of days that will analyze in order
#to create arrays of the needed dimensions

	#create continous calendar
	p1<-as.Date(paste0(ANNO1,"-01-01"))
	p2<-as.Date(paste0(ANNO2,"-12-31"))
	datas=seq(p1,p2,by="day")

	#select only days correspondeing to the needed season
	timeseason=season2timeseason(season)
	month=as.numeric(format(datas,"%m"))
	whichdays=which(month %in% timeseason)

	#create a "season" for continuous time, used by persistance tracking
	seas=whichdays*1; ss=1
	for (i in 1:(length(whichdays)-1)) 
		{
		if (diff(whichdays)[i]>1)  {ss=ss+1}
		seas[i+1]=ss
		}

	#produce a final timeseries of dates
	datas=datas[whichdays]
	dataline=list(day=as.numeric(format(datas,"%d")),month=as.numeric(format(datas,"%m")),year=as.numeric(format(datas,"%Y")),season=seas,data=datas)
	print("Time Array Built")
	print(paste("Length:",length(seas),"days for",season,"season"))
	print(paste("From",datas[1],"to",datas[length(seas)]))

	return(dataline)
}

power.date.no.leap<-function(season,ANNO1,ANNO2)
{
	#apply to power.date object to clean out elements for leap years
	e=power.date(season,ANNO1,ANNO2)
	leap.days=which(e$month==2 & e$day==29)
	dataline.leap=list(day=e$day[-leap.days],month=e$month[-leap.days],year=e$year[-leap.days],season=e$season[-leap.days],data=e$data[-leap.days])
	print("FIXED FOR NO LEAP CALENDAR: Time Array Built")
	print(paste("Length:",length(dataline.leap$season),"days for",season,"season"))
	print(paste("From",dataline.leap$data[1],"to",dataline.leap$data[length(dataline.leap$season)]))
	return(dataline.leap)
}

power.date.30day<-function(season,ANNO1,ANNO2)
{
	#apply to power.date object to clean out elements for leap years
	nmonths=length(season2timeseason(season))
	nyears=as.numeric(ANNO2)-as.numeric(ANNO1)+1
	dd=rep(seq(1,30),nmonths*nyears)
	mm=rep(rep(season2timeseason(season),each=30),nyears)
	#create a "season" for continuous time, used by persistance tracking
	seas=mm*0+1; ss=1
	for (i in 1:(length(mm)-1))
	        {
	        if (diff(mm)[i]>1)  {ss=ss+1}
	        seas[i+1]=ss
	        }
	dataline.30day=list(day=dd,month=mm,season=seas)
	print("SIMPLIFIED CALENDAR FOR 30-day CALENDAR: Time Array Built")
	print(paste("Length:",length(dataline.30day$season),"days for",season,"season"))
	return(dataline.30day)
}

##########################################################
#--------------NetCDF loading function-------------------#
##########################################################

#function to open ncdf files (much more refined, with CDO-based interpolation)
ncdf.opener<-function(namefile,namevar=NULL,namelon="lon",namelat="lat",rotate="full",interp2grid=F,grid="r144x73",remap_method="remapcon2")
{
#function to open netcdf files. It uses ncdf4 library. support only 1D (t), 2D (x,y) or 3D (x,y,t) data in any netcdf format.
#automatically rotate matrix to place greenwich at the center (flag "rotate") and flip the latitudes in order to have increasing
#if require (flag "interp2grid") additional interpolation with CDO can be used. "grid" can be used to specify the grid name
require(ncdf4)

if (rotate=="full") {rot=T; move1=move2=1/2} #180 degrees rotation of longitude
if (rotate=="half") {rot=T; move1=1/4; move2=3/4} #90 degree rotation (useful for TM90)
if (rotate=="no") {rot=F} #keep as it is, breaking at Greemwich

#interpolation made with CDO: second order conservative remapping
if (interp2grid)
        {
        print(paste("Remapping with CDO on",grid,"grid"))
        filename=basename(normalizePath(namefile))
	filedir=dirname(normalizePath(namefile))
	cdo=Sys.which("cdo")
        tempfile=paste0(file.path(filedir,paste0("tempfile_",filename)))
        #system(paste0(cdo," ",remap_method,",",grid," ",namefile," ",tempfile))
	system2(cdo,args=c(paste0(remap_method,",",grid),namefile,tempfile))
        namefile=tempfile 
        }

#define rotate function (faster than with apply)
rotation<-function(line) {
vettore=line; dims=length(dim(vettore))
if (dims==1) #for longitudes
{ll=length(line); line[(ll*move1):ll]=vettore[1:(ll*move2+1)]; line[1:(ll*move1-1)]=vettore[(ll*move2+2):ll]-360}
if (dims==2) #for x,y data
{ll=length(line[,1]); line[(ll*move1):ll,]=vettore[1:(ll*move2+1),]; line[1:(ll*move1-1),]=vettore[(ll*move2+2):ll,]}
if (dims==3) #for x,y,t data
{ll=length(line[,1,1]); line[(ll*move1):ll,,]=vettore[1:(ll*move2+1),,]; line[1:(ll*move1-1),,]=vettore[(ll*move2+2):ll,,]}
return(line)    }

#define flip function ('cos rev/apply is not working)
flipper<-function(field) {
dims=length(dim(field))
if (dims==2) {ll=length(daily[1,]); field=field[,ll:1]} #for x,y data
if (dims==3) {ll=length(daily[1,,1]); field=field[,ll:1,]} #for x,y,t data
return(field) }

#opening file: getting variable (if namevar is given, that variable is extracted)
a=nc_open(namefile)
if (is.null(namevar)) {daily=ncvar_get(a)} else {daily=ncvar_get(a,namevar)}

#check for dimensions (presence or not of time dimension)
dimensions=length(dim(daily))

#if dimensions are multiple, get longitude, latitude
#if needed, rotate and flip the array
if (dimensions>1)
{
        #read attributes
        ics=ncvar_get(a,namelon); ipsilon=ncvar_get(a,namelat)
        
	#longitute rotation around Greenwich
        if (rot)     {ics=rotation(ics); daily=rotation(daily) }
        if (ipsilon[2]<ipsilon[1] & length(ipsilon)>1)
                if (length(ics)>1)
                {ipsilon=sort(ipsilon); daily=flipper(daily) }
                else
                {ipsilon=sort(ipsilon); daily=flipper.zonal(daily) }

        #exporting variables to the main program
        assign("ics",ics, envir = .GlobalEnv)
        assign("ipsilon",ipsilon, envir = .GlobalEnv)

} 

if (dimensions>3)
{stop("This file is more than 3D file")}

#close connection
nc_close(a)

#remove interpolated file
if (interp2grid) {system2("rm",tempfile)}

#showing array properties
#print(dim(daily))

return(daily)
}


#function to open ncdf files (much more refined, with CDO-based interpolation)
ncdf.opener.time<-function(namefile,namevar=NULL,namelon=NULL,namelat=NULL,tmonths=NULL,tyears=NULL,rotate="full",interp2grid=F,grid="r144x73",remap_method="remapcon2")
{
#function to open netcdf files. It uses ncdf4 library
#time selection of month and years needed
#automatically rotate matrix to place greenwich at the center (flag "rotate") and flip the latitudes in order to have increasing
#if require (flag "interp2grid") additional interpolation with CDO can be used. "grid" can be used to specify the grid name
require(ncdf4)
require(PCICt)

if (is.null(tyears) | is.null(tmonths)) {stop("Please specify both months and years to load")}

if (rotate=="full") {rot=T; move1=move2=1/2} #180 degrees rotation of longitude
if (rotate=="half") {rot=T; move1=1/4; move2=3/4} #90 degree rotation (useful for TM90)
if (rotate=="no") {rot=F} #keep as it is, breaking at Greemwich

#interpolation made with CDO: second order conservative remapping
if (interp2grid)
        {
        print(paste("Remapping with CDO on",grid,"grid"))
        filename=basename(normalizePath(namefile))
        filedir=dirname(normalizePath(namefile))
        cdo=Sys.which("cdo")
        tempfile=paste0(file.path(filedir,paste0("tempfile_",filename)))
        system2(cdo,args=c(paste0(remap_method,",",grid),namefile,tempfile))
        namefile=tempfile
        }

#define rotate function (faster than with apply)
rotation<-function(line) {
vettore=line; dims=length(dim(vettore))
if (dims==1) #for longitudes
{ll=length(line); line[(ll*move1):ll]=vettore[1:(ll*move2+1)]; line[1:(ll*move1-1)]=vettore[(ll*move2+2):ll]-360}
if (dims==2) #for x,y data
{ll=length(line[,1]); line[(ll*move1):ll,]=vettore[1:(ll*move2+1),]; line[1:(ll*move1-1),]=vettore[(ll*move2+2):ll,]}
if (dims==3) #for x,y,t data
{ll=length(line[,1,1]); line[(ll*move1):ll,,]=vettore[1:(ll*move2+1),,]; line[1:(ll*move1-1),,]=vettore[(ll*move2+2):ll,,]}
return(line)    }

#define flip function ('cos rev/apply is not working)
flipper<-function(field) {
dims=length(dim(field))
if (dims==2) {ll=length(field[1,]); field=field[,ll:1]} #for x,y data
if (dims==3) {ll=length(field[1,,1]); field=field[,ll:1,]} #for x,y,t data
return(field) }


#opening file: getting variable (if namevar is given, that variable is extracted)
print(paste("opening file:",namefile))
a=nc_open(namefile)

#load axis
naxis=names(a$dim)[1:min(c(4,length(a$dim)))]
for (axis in naxis) {print(axis); assign(axis,ncvar_get(a,axis))}

print("selecting years and months")
#extracting time (BETA)
#origin=strsplit(ncatt_get(a,"time","units")$value," ")[[1]][3]


#if (substring(origin, 1, 1)=="%") {
#	print("qui")
#	timeline=strptime(time,format="%Y%m%d") #time with format
#} else {
#	origin.pcict <- as.PCICt(origin, cal)
#	timeline=origin.pcict + (time * 86400)	#time with origin
#}
#select time needed

#based on preprocessing of CDO time format: get calendar type and use PCICt package for irregular data
cal=ncatt_get(a,"time","calendar")$value
timeline=as.PCICt(as.character(time),format="%Y%m%d",cal=cal)

# break if the calendar has not been recognized
if (any(is.na(timeline))) {
	stop("Calendar from NetCDF is unsupported or not present. Stopping!!!")
}
#break if the data requested is not there
lastday=as.PCICt(paste0(max(tyears),"-",max(tmonths),"-28"),cal="standard",format="%Y-%m-%d")
firstday=as.PCICt(paste0(min(tyears),"-",min(tmonths),"-01"),cal="standard",format="%Y-%m-%d")
if (max(timeline)<lastday | min(timeline)>min(timeline)) {
	stop("You requested a time interval that is not present in the NetCDF")
}

#time selection and variable loading
print("loading full field...")
#if no name provided load the only variable available
if (is.null(namevar)) {namevar=names(a$var)}
field=ncvar_get(a,namevar)
print(str(field))

# select data we need
select=which(as.numeric(format(timeline,"%Y")) %in% tyears & as.numeric(format(timeline,"%m")) %in% tmonths)
field=field[,,select]
time=timeline[select]

print(paste("This is a",cal,"calendar"))
print(paste(length(time),"days selected from",time[1],"to",time[length(time)]))

#check for dimensions (presence or not of time dimension)
dimensions=length(dim(field))

#if dimensions are multiple, get longitude, latitude
#if needed, rotate and flip the array
if (dimensions>1)
{
	#assign ics and ipsilon 
	if (is.null(namelon)) {
		xlist=c("lon","Lon","longitude","Longitude")
		if (any(xlist %in% naxis))  {
			  ics=get(names(a$dim[which(naxis %in% xlist)]))} else {stop("No lon found")}
		} else {
		ics=ncvar_get(a,namelon)
		}
	if (is.null(namelat)) {
		ylist=c("lat","Lat","latitude","Latitude")
		if (any(ylist %in% naxis))  {
			ipsilon=get(names(a$dim[which(naxis %in% ylist)]))} else {stop("No lat found")}
		} else {
		ipsilon=ncvar_get(a,namelat)
		}
		
	print("flipping and rotating")
        #longitute rotation around Greenwich
        if (rot)     {ics=rotation(ics); field=rotation(field) }
        if (ipsilon[2]<ipsilon[1] & length(ipsilon)>1 )
                if (length(ics)>1)
                {ipsilon=sort(ipsilon); field=flipper(field) }

        #exporting variables to the main program
        assign("ics",ics, envir = .GlobalEnv)
        assign("ipsilon",ipsilon, envir = .GlobalEnv)
	assign(names(a$dim[which(naxis %in% xlist)]),ics)
	assign(names(a$dim[which(naxis %in% ylist)]),ipsilon)

}


if (dimensions>3)
{stop("This file is more than 3D file")}

#close connection
nc_close(a)

#remove interpolated file
if (interp2grid) {system2("rm",tempfile)}

#showing array properties
print(paste(dim(field)))
print(paste("From",time[1],"to",time[length(time)]))

return(mget(c("field",naxis)))
}


##########################################################
#--------------Plotting functions------------------------#
##########################################################


#extensive filled.contour function 
filled.contour3 <-
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), extend=TRUE, plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar, ...) 
{
  # modification by Ian Taylor of the filled.contour function
  # to remove the key and facilitate overplotting with contour()
  # further modified by Carey McGilliard and Bridget Ferris
  # to allow multiple plots on one page
  # modification to allow plot outside boundaries

  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")

 #trim extremes for nicer plots
 if (extend) {
	#extendvalue=10^8
	#levels=c(-extendvalue,levels,extendvalue)
	#col=c(col[1],col,col[length(col)])
	z[z<min(levels)]=min(levels)
	z[z>max(levels)]=max(levels)
	}

 plot.new()
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
  .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                          col = col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1,...)
      Axis(y, side = 2,...)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}

image.scale3 <- function(z,levels,color.palette=heat.colors,colorbar.label="image.scale",extend=T,
line.label=2,line.colorbar=0,cex.label=1,cex.colorbar=1,colorbar.width=1,...){

 #save properties from main plotting region
 old.par <- par(no.readonly = TRUE)
 mfg.save <- par()$mfg
 old.fig=par()$fig

 #defining plotting region with proper scaling
 #print(old.fig)
 xscal=(old.fig[2]-old.fig[1]); yscal=(old.fig[4]-old.fig[3]); lw=colorbar.width; lp=line.colorbar/100
 new.fig=c(old.fig[2]-0.07*xscal*lw-lp,old.fig[2]-0.03*xscal-lp,old.fig[3]+0.1*yscal,old.fig[4]-0.1*yscal)
 #print(new.fig)

 if (missing(levels)) { levels=seq(min(z),max(z),,12)}
 #fixing color palette
 col=color.palette(length(levels)-1)

 #starting plot
 par(mar=c(1,1,1,1),fig=new.fig,new=TRUE)

 #creating polygons for legend 
 poly <- vector(mode="list", length(col))
 for(i in seq(poly))
  {poly[[i]] <- c(levels[i], levels[i+1], levels[i+1], levels[i])}
  
 xlim<-c(0,1)
 if (extend) {
	longer=1.5
	dl=diff(levels)[1]*longer
	ylim<-c(min(levels)-dl,max(levels)+dl)
	} else {
	ylim<-range(levels)
 }
 plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)
 for(i in seq(poly))
   {polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)}

if (extend) {
  polygon(c(0,1,1/2), c(levels[1], levels[1], levels[1]-dl), 
                col = col[1],border=NA) 
  polygon(c(0,1,1/2), c(levels[length(levels)], levels[length(levels)], levels[length(levels)]+dl), 
                col = col[length(col)],border=NA)
  polygon(c(0,0,1/2,1,1,1/2),c(levels[1],levels[length(levels)],levels[length(levels)]+dl, levels[length(levels)],levels[1],
		levels[1]-dl),border="black",lwd=2)
  ylim0=range(levels)
  prettyspecial=pretty(ylim0); prettyspecial=prettyspecial[prettyspecial<=max(ylim0) & prettyspecial>=min(ylim0)]
  axis(4,las=1,cex.axis=cex.colorbar,at=prettyspecial,labels=prettyspecial,...)
} else {
  box()
  axis(4,las=1,cex.axis=cex.colorbar,...)
}

 #box, axis and leged
 mtext(colorbar.label,line=line.label,side=4,cex=cex.label,...)

 #resetting properties for starting a new plot (mfrow style)
 par(old.par)
 par(mfg = mfg.save, new = FALSE)
 invisible()

}

##########################################################
#------------Blocking Tracking Functions-----------------#
##########################################################

#time persistence (used for longitude filter too)
time.persistence<-function(timeseries,persistence=5)
{
rr=rle(timeseries)
rr$values[which(rr$values==1 & rr$length<persistence)]=0
nn=rep(rr$values,rr$length)
return(nn)
}


#blocking 5 days tracking
blocking.persistence<-function(field,persistence=5,time.array)
{

#function for persistence
#pers<-function(timeseries,persistence,time.array)
#{
#        xx=NULL
#        for (s in min(time.array$season):max(time.array$season))
#                {
#                yy=timeseries[which(time.array$season==s)]
#                nn=time.persistence(yy,persistence)
#                xx=append(xx,nn)
#                }
#        return(xx)
#}

#function for persistence
pers2<-function(timeseries,persistence,time.array)
{     
dd=min(time.array$season):max(time.array$season) 
nn=sapply(dd, function(x) {time.persistence(timeseries[which(time.array$season==x)],persistence)})
xx=c(unlist(nn))
return(xx)
}

# check for etime
if (length(time.array$month)!=length(field[1,1,])) { stop("Wrong time array! Exiting...") }

print("Time filtering...")
newfield=apply(field,c(1,2),function(x) pers2(x,persistence=5,time.array))
newfield=aperm(newfield,c(2,3,1))
print("Mean field...")
meanfield=apply(newfield,c(1,2),mean,na.rm=T)*100


print("Events detection...")
maxdim=max(apply(newfield,c(1,2),function(x) length(rle(x)$length[which(rle(x)$values==1)])))
events=apply(newfield,c(1,2),function(x) c(rle(x)$lengths[which(rle(x)$values==1)],rep(NA,maxdim-length(rle(x)$length[which(rle(x)$values==1)]))))
events=aperm(events,c(2,3,1))
print("Mean Duration...")
duration=apply(events,c(1,2),mean,na.rm=T)
print("Number of Events...")
nevents=apply(events,c(1,2),function(x) length(x[!is.na(x)]))

out=list(track=newfield,percentage=meanfield,duration=duration,events=events,nevents=nevents)
return(out)
}


#large scale extension with further implementation
largescale.extension.if<-function(ics,ipsilon,field)
{
print("Large Scale Extension based on fixed angle")
fimin=30 #southern latitude to be analyzed
fimax=75 #northern latitude to be analyzed
deltaics=diff(ics)[1]
deltaips=diff(ipsilon)[1]
passo=round(5/deltaics)  #horizontal movemenent
vertical=round(2.5/deltaips) #vertical movement
#time=1:length(field[1,1,]) #elements of the length of the dataset
time=which(apply(field,3,max)!=0) #elements length of the dataset (removing no blocked days)

print(paste("Box dimension:",passo*2*deltaics,"° lon x ",vertical*2*deltaips,"° lat"))

short<-function(ics,ipsilon,field,passo,vertical) {
	control=field
	range=which.min(abs(ipsilon-fimin)):which.min(abs(ipsilon-fimax)) #check range for latitude excursion
	#range=range[(1+vertical):(length(range)-vertical)] #reduce range considering border effect
	
	new=rbind(field,field,field) #bind domain for cross-date line
	for (i in 1:length(ics))
		{
		ii=i+length(ics)
		if (!all(new[(ii-passo):(ii+passo),]==0)) #check to speed up
			{
			for (j in range)
				{
				control[i,j]=mean(new[(ii-passo):(ii+passo),(j-vertical):(j+vertical)],na.rm=T)
				}
			}
		}
	control[control>0]=1
	return(control)
}


for (t in time)
{
	if (any(t==round(seq(0,length(field[1,1,]),,11))))
        	{print(paste("--->",round(t/length(field[1,1,])*100),"%"))}
			{field[,,t]=short(ics,ipsilon,field[,,t],passo,vertical)}
}
return(field)
}

#large scale extension
largescale.extension2<-function(ics,ipsilon,field)
{
print("Large Scale Extension based on fixed angle")
deltaics=(ics[20]-ics[19])
deltaips=(ipsilon[20]-ipsilon[19])
passo=round(5/(ics[20]-ics[19]))
vertical=round(2.5/(ipsilon[20]-ipsilon[19]))

print(paste("Box dimension:",passo*2*deltaics,"° lon x ",vertical*2*deltaips,"° lat"))

short<-function(ics,ipsilon,field,passo,vertical)
{
out=field
startipsilon=which.min(abs(ipsilon-30))
estension=round((75-30)/(ipsilon[20]-ipsilon[19]))
new=rbind(field,field,field)
for (i in 1:length(ics))
{ii=i+length(ics)
for (j in startipsilon:(startipsilon+estension))
{
control=mean(new[(ii-passo):(ii+passo),(j-vertical):(j+vertical)],na.rm=T)
if (control>0)
{out[i,j]=1}
}
}
return(out)
}

for (t in 1:length(field[1,1,]))
{
if (any(t==round(seq(0,length(field[1,1,]),,11))))
        {print(paste("--->",round(t/length(field[1,1,])*100),"%"))}
if (all(!is.na(field[,,t])))
{field[,,t]=short(ics,ipsilon,field[,,t],passo,vertical)}
}
return(field)
}

#Longitude filter for minimum extension
longitude.filter<-function(ics,ipsilon,field)
{
print("Longitude filter based on fixed angle")
out=field
deltaics=(ics[20]-ics[19])
startipsilon=which.min(abs(ipsilon-30))
estension=round((75-30)/(ipsilon[20]-ipsilon[19]))
passo=round(15/(ics[20]-ics[19]))

print(paste("Continous longitude contrain",passo*deltaics,"° lon"))

for (t in 1:length(field[1,1,]))
{
if (any(t==round(seq(0,length(field[1,1,]),,11))))
        {print(paste("--->",round(t/length(field[1,1,])*100),"%"))}

new=rbind(field[,,t],field[,,t],field[,,t])
for (j in startipsilon:((startipsilon+estension)))
{
new[,j]=time.persistence(new[,j],persistence=passo)
}
field[,,t]=new[length(ics)+(1:length(ics)),]
}
return(field)
}


##########################################################
#------------EOFs and regims functions-------------------#
##########################################################

eofs<-function(lon,lat,field,neof=4,xlim,ylim,method="SVD",do_standardize=F,do_regression=F)
{
# R tool for computing EOFs based on Singular Value Decomposition ("SVD", default)
# or with the eigenvectors of the covariance matrix ("covariance", slower) 
# If requested, computes linear regressions and standardizes the PCs
# If you want to use the regressions, remember to standardize the PCs
# Take as input a 3D anomaly field.
# Requires "personal" functions area.weight, whicher and standardize

#area weighting, based on the root of cosine
print("Area Weighting...")
ww=area.weight(lon,lat,root=T)
wwfield=sweep(field,c(1,2),ww,"*")

#selection of the box
box=wwfield[whicher(lon,xlim[1]):whicher(lon,xlim[2]),whicher(lat,ylim[1]):whicher(lat,ylim[2]),]
slon=lon[whicher(lon,xlim[1]):whicher(lon,xlim[2])]
slat=lat[whicher(lat,ylim[1]):whicher(lat,ylim[2])]

#transform 3D field in a matrix
new_box=array(box,dim=c(dim(box)[1]*dim(box)[2],dim(box)[3]))

#calling SVD
if (method=="SVD")
{       
        print("Calling SVD...")
        SVD=svd(new_box,nu=neof,nv=neof)
        
        #extracting EOFs (loading pattern), expansions coefficient and variance explained
        pattern=array(SVD$u,dim=c(dim(box)[1],dim(box)[2],neof))
        coefficient=SVD$v
        variance=(SVD$d[1:neof])^2/sum((SVD$d)^2)
        if (do_standardize)
                { coefficient=apply(coefficient,c(2),standardize) }
                else
                { coefficient=sweep(coefficient,c(2),sqrt(variance),"*") }
}

#calling covariance matrix
if (method=="covariance")
{       
        print("Calling eigenvectors of the covariance matrix...")
        covma=cov(t(new_box))
        eig=eigen(covma)
        coef=(t(new_box)%*%eig$vector)[,1:neof]
        pattern=array(eig$vectors,dim=c(dim(box)[1],dim(box)[2],dim(box)[3]))[,,1:neof]
        variance=eig$values[1:neof]/sum(eig$values)
           if (do_standardize)
                { coefficient=apply(coef,c(2),standardize) }
                else
                { coefficient=coef }
}

#linear regressions on anomalies
regression=NULL
if (do_regression)
{       
        print("Linear Regressions (it can takes a while)... ")
        regression=array(NA,dim=c(length(lon),length(lat),neof))
        for (i in 1:neof) {regression[,,i]=apply(field,c(1,2),function(x) coef(lm(x ~ coefficient[,i]))[2])}
}

#preparing output
print("Finalize...")
pattern=list(x=slon,y=slat,z=pattern)
out=list(pattern=pattern,coeff=coefficient,variance=variance,regression=regression)
return(out)
}


regimes<-function(lon,lat,field,ncluster=4,ntime=1000,neof=10,xlim,ylim,alg="Hartigan-Wong")
{
# R tool to compute cluster analysis based on k-means.
# Requires "personal" function eofs
# Take as input a 3D anomaly field

#Reduce the phase space with EOFs: use SVD and do not standardize PCs
print("Launching EOFs...")
t0=proc.time()
reducedspace=eofs(lon,lat,field,neof=neof,xlim=xlim,ylim=ylim,method="SVD",do_regression=F,do_standardize=F)
t1=proc.time()-t0
#print(t1)

#extract the principal components
PC=reducedspace$coeff
print(str(PC))

#k-means computation repeat for ntime to find best solution. 
print("Computing k-means...")
t0=proc.time()
print(str(ncluster))
regimes=kmeans(PC,as.numeric(ncluster),nstart=ntime,iter.max=1000,algorithm=alg)
t1=proc.time()-t0
#print(t1)

#Extract regimes frequencyr and timeseries of occupation
cluster=regimes$cluster
frequencies=regimes$size/dim(field)[3]*100
print(frequencies[order(frequencies,decreasing=T)])
#print(regimes$tot.withinss)

print("Creating Composites...")
compose=aperm(apply(field,c(1,2),by,cluster,mean),c(2,3,1))

#sorting from the more frequent to the less frequent
kk=order(frequencies,decreasing=T)
cluster=cluster+10
for (ss in 1:ncluster) {cluster[cluster==(ss+10)]=which(kk==ss)}

#prepare output
print("Finalize...")
out=list(cluster=cluster,frequencies=frequencies[kk],regimes=compose[,,kk],tot.withinss=regimes$tot.withinss)
return(out)
}

