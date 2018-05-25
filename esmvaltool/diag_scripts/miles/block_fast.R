######################################################
#-----Blocking routines computation for MiLES--------#
#-------------P. Davini (Oct 2014)-------------------#
######################################################
miles.block.fast<-function(exp,year1,year2,season,z500filename,FILESDIR)
{

#t0
t0<-proc.time()

#setting up main variables
BLOCKDIR=file.path(FILESDIR,exp,"Block",paste0(year1,"_",year2),season)
dir.create(BLOCKDIR,recursive=T)

#setting up time domain
years=year1:year2
timeseason=season2timeseason(season)

#new file opening
nomefile=z500filename
fieldlist=ncdf.opener.time(nomefile,"zg",tmonths=timeseason,tyears=years,rotate="full")
print(str(fieldlist))

#time array
datas=fieldlist$time
#etime=list(day=as.numeric(format(datas,"%d")),month=as.numeric(format(datas,"%m")),year=as.numeric(format(datas,"%Y")),data=datas)
etime=power.date.new(datas)
totdays=length(datas)

#declare variable
Z500=fieldlist$field


##########################################################
#--------------Tibaldi and Molteni 1990------------------#
##########################################################

print("Tibaldi and Molteni (1990) index...")
# TM90: parametres for blocking detection (beta)
tm90_fi0=60 #central_lat
tm90_fiN=tm90_fi0+20; tm90_fiS=tm90_fi0-20 #south and north lath
tm90_central=whicher(ipsilon,tm90_fi0)
tm90_south=whicher(ipsilon,tm90_fiS)
tm90_north=whicher(ipsilon,tm90_fiN)
tm90_range=c(-2:2)  #5 degrees to the north, 5 to the south (larger than TM90 or D'Andrea et al 1998)

#TM90: beta version, the amazing power of R vectorization!
#6 lines to get the climatology
tm90_ghgn=(Z500[,tm90_north+tm90_range,]-Z500[,tm90_central+tm90_range,])/(tm90_fiN-tm90_fi0)
tm90_ghgs=(Z500[,tm90_central+tm90_range,]-Z500[,tm90_south+tm90_range,])/(tm90_fi0-tm90_fiS)
tm90_check=(tm90_ghgs>0 & tm90_ghgn<(-10)) # TM90 conditions
tm90_check[tm90_check==T]=1; tm90_check[tm90_check==F]=0
totTM90=apply(tm90_check,c(1,3),max,na.rm=T)
TM90=apply(totTM90,1,mean)*100
print("Done!")

##########################################################
#--------------Davini et al. 2012------------------------#
##########################################################

# decleare main variables to be computed (considerable speed up!)
totrwb=totmeridional=totBI=Z500*NA
totblocked=totblocked2=Z500*0

# Davini et al. 2012: parameters to be set for blocking detection
fi0=30                          #lowest latitude to be analyzed
jump=15                         #distance on which compute gradients
step0=round(jump/diff(ipsilon)[1])   #number of grid points to be used
central=which.min(abs(ipsilon-fi0))             #lowest starting latitude
north=central+step0                             #lowest north latitude
south=central-step0                                     #lowest sourth latitude
maxsouth=central-2*step0
fiN=ipsilon[north]
fiS=ipsilon[south]
range=round((90-fi0-jump)/diff(ipsilon)[1])  #escursion to the north for computing blocking

print("--------------------------------------------------")
print("Davini et al. (2012) index and diagnostics...")
print(c("distance for gradients:",step0*diff(ics)[1]))
print(paste("range of latitudes ",fi0,"-",90-step0*diff(ics)[1]," N",sep=""))

##########################################################
#--------------Istantaneous Blocking---------------------#
##########################################################

#----COMPUTING BLOCKING INDICES-----
for (t in 1:totdays)
		{
		if (any(t==round(seq(0,totdays,,11))))
                {print(paste("--->",round(t/totdays*100),"%"))}

		#multidim extension
                new_field=rbind(Z500[,,t],Z500[,,t],Z500[,,t])

                for (delta in 0:range)  # computing blocking for different latitudes
                {
                ghgn=(Z500[,north+delta,t]-Z500[,central+delta,t])/(fiN-fi0)
                ghgs=(Z500[,central+delta,t]-Z500[,south+delta,t])/(fi0-fiS)
                gh2gs=(Z500[,south+delta,t]-Z500[,maxsouth+delta,t])/(fi0-fiS)
                check1=which(ghgs>0 & ghgn<(-10))
                check2=which(ghgs>0 & ghgn<(-10) & gh2gs<(-5))  #supplementary condition

		if (length(check2)>0) {
			totblocked2[check2,central+delta,t]=1
		}

                if (length(check1)>0)
                {
                # 1-MATRIX FOR INSTANTANEOUS BLOCKING
                totblocked[check1,central+delta,t]=1
		

                # 2-PART ON COMPUTATION OF ROSSBY WAVEBREAKING
                r=check1+length(ics)
                rwb_jump=jump/2
                steprwb=round(rwb_jump/(ics[20]-ics[19]))
                rwb_west=new_field[(r-steprwb),south+delta+steprwb]
                rwb_east=new_field[(r+steprwb),south+delta+steprwb]
                fullgh=(rwb_west-rwb_east)

                totrwb[check1[fullgh<0],central+delta,t]=(-10)     # gradient decreasing: cyclonic RWB     
                totrwb[check1[fullgh>0],central+delta,t]=10        # gradient increasing: anticyclonic RWB                 

                # 4-part about adapted version of blocking intensity by Wiedenmann et al. (2002)
                step=round(60/(ics[2]-ics[1]))
                ii=check1+length(ics)
                zu=zd=NULL
                for (ll in ii)
                        {
                        zu=c(zu,min(new_field[(ll-step):ll,central+delta]))
                        zd=c(zd,min(new_field[ll:(ll+step),central+delta]))
                        }
                mz=Z500[check1,central+delta,t]
                rc=0.5*((zu+mz)/2+(zd+mz)/2)
                totBI[check1,central+delta,t]=100*(mz/rc-1)

                #5 - part about meridional gradient index
                totmeridional[check1,central+delta,t]=ghgs[check1]


                }}}

                print(paste("Total # of days:",t))
                print("-------------------------")

##########################################################
#--------------------Mean Values-------------------------#
##########################################################

#compute mean values (use rowMeans that is faster when there are no NA values)
frequency=rowMeans(totblocked,dims=2)*100               #frequency of Instantaneous Blocking days
frequency2=rowMeans(totblocked2,dims=2)*100               #frequency of Instantaneous Blocking days with GHGS2
Z500mean=rowMeans(Z500,dims=2)                          #Z500 mean value
BI=apply(totBI,c(1,2),mean,na.rm=T)                     #Blocking Intensity Index as Wiedenmann et al. (2002)
MGI=apply(totmeridional,c(1,2),mean,na.rm=T)   		 #Value of meridional gradient inversion

#anticyclonic and cyclonic averages RWB
CN=apply(totrwb,c(1,2),function(x) sum(x[x==(-10)],na.rm=T))/(totdays)*(-10)
ACN=apply(totrwb,c(1,2),function(x) sum(x[x==(10)],na.rm=T))/(totdays)*(10)

t1=proc.time()-t0
print(t1)

print("Instantaneous blocking and diagnostics done!")

##########################################################
#--------------------Time filtering----------------------#
##########################################################

#spatial filtering on fixed longitude distance
spatial=longitude.filter(ics,ipsilon,totblocked)
#CUT=apply(spatial,c(1,2),sum,na.rm=T)/ndays*100

#large scale extension on 10x5 box
large=largescale.extension.if(ics,ipsilon,spatial)
#LARGE=apply(large,c(1,2),sum,na.rm=T)/ndays*100

#5 days persistence filter
block=blocking.persistence(large,persistence=5,time.array=etime)

tf=proc.time()-t1
print(tf)


##########################################################
#------------------------Save to NetCDF------------------#
##########################################################

#saving output to netcdf files
print("saving NetCDF climatologies...")
savefile1=paste0(BLOCKDIR,"/BlockClim_",exp,"_",year1,"_",year2,"_",season,".nc")
savefile2=paste0(BLOCKDIR,"/BlockFull_",exp,"_",year1,"_",year2,"_",season,".nc")

#which fieds to plot/save
fieldlist=c("TM90","InstBlock","ExtraBlock","Z500","MGI","BI","CN","ACN","BlockEvents","DurationEvents","NumberEvents")
full_fieldlist=c("TM90","InstBlock","ExtraBlock","Z500","MGI","BI","CN","ACN","BlockEvents")

# dimensions definition
TIME=paste("days since ",year1,"-",timeseason[1],"-01 00:00:00",sep="")
LEVEL=50000
fulltime=as.numeric(etime$data)-as.numeric(etime$data)[1]
x <- ncdim_def( "Lon", "degrees", ics)
y <- ncdim_def( "Lat", "degrees", ipsilon)
z <- ncdim_def( "Lev", "Pa", LEVEL)
t1 <- ncdim_def( "Time", TIME, 1,unlim=T)
t2 <- ncdim_def( "Time", TIME, fulltime,unlim=T)

for (var in fieldlist)
{
        #name of the var
	if (var=="TM90")
                {longvar="TibaldiMolteni 1990 Instantaneous Blocking frequency"; unit="%"; field=TM90; full_field=totTM90}
        if (var=="InstBlock")
                {longvar="Instantaneous Blocking frequency"; unit="%"; field=frequency; full_field=totblocked}
        if (var=="ExtraBlock")
                {longvar="Instantaneous Blocking frequency (GHGS2)"; unit="%"; field=frequency2; full_field=totblocked2}
        if (var=="Z500")
                {longvar="Geopotential Height"; unit="m"; field=Z500mean; full_field=Z500}
        if (var=="BI")
                {longvar="BI index"; unit=""; field=BI; full_field=totBI}
        if (var=="MGI")
                {longvar="MGI index"; unit=""; field=MGI; full_field=totmeridional}
        if (var=="ACN")
                {longvar="Anticyclonic RWB frequency"; unit="%"; field=ACN; full_field=totrwb/10; full_field[full_field==(-1)]=NA}
        if (var=="CN")
                {longvar="Cyclonic RWB frequency"; unit="%"; field=CN; full_field=totrwb/10; full_field[full_field==(1)]=NA}
        if (var=="BlockEvents")
                {longvar="Blocking Events frequency"; unit="%"; field=block$percentage; full_field=block$track}
        if (var=="DurationEvents")
                {longvar="Blocking Events duration"; unit="days"; field=block$duration}
        if (var=="NumberEvents")
                {longvar="Blocking Events number"; unit=""; field=block$nevents}

	#fix eventual NaN	
	field[is.nan(field)]=NA

        #variable definitions
        if (var=="TM90") {
		var_ncdf=ncvar_def(var,unit,list(x,t=t1),-999,longname=longvar,prec="single",compression=1)
		full_var_ncdf=ncvar_def(var,unit,list(x,t=t2),-999,longname=longvar,prec="single",compression=1)
	} else {
		var_ncdf=ncvar_def(var,unit,list(x,y,z,t=t1),-999,longname=longvar,prec="single",compression=1)
		full_var_ncdf=ncvar_def(var,unit,list(x,y,z,t=t2),-999,longname=longvar,prec="single",compression=1)
	}
	
        assign(paste0("var",var),var_ncdf)
	assign(paste0("full_var",var),full_var_ncdf)
        assign(paste0("field",var),field)
	assign(paste0("full_field",var),full_field)
}

#Climatologies Netcdf file creation
print(savefile1)
namelist1=paste0("var",fieldlist)
nclist1 <- mget(namelist1)
ncfile1 <- nc_create(savefile1,nclist1)
for (var in fieldlist)
{
        # put variables into the ncdf file
	#ncvar_put(ncfile1, fieldlist[which(var==fieldlist)], get(paste0("field",var)), start = c(1, 1, 1, 1),  count = c(-1,-1,-1,-1))
	ndims=get(paste0("var",var))$ndims
	ncvar_put(ncfile1, var, get(paste0("field",var)), start = rep(1,ndims),  count = rep(-1,ndims))
}
nc_close(ncfile1)

#Fullfield Netcdf file creation
print(savefile2)
namelist2=paste0("full_var",full_fieldlist)
nclist2 <- mget(namelist2)
ncfile2 <- nc_create(savefile2,nclist2)
for (var in full_fieldlist)
{
        # put variables into the ncdf file
        #ncvar_put(ncfile2, full_fieldlist[which(var==full_fieldlist)], get(paste0("full_field",var)), start = c(1, 1, 1, 1),  count = c(-1,-1,-1,-1))
	ndims=get(paste0("full_var",var))$ndims
	ncvar_put(ncfile2, var, get(paste0("full_field",var)), start = rep(1,ndims),  count = rep(-1,ndims))

}
nc_close(ncfile2)

}

# REAL EXECUTION OF THE SCRIPT 
# read command line
args <- commandArgs(TRUE)

# number of required arguments from command line
name_args=c("exp","year1","year2","season","z500filename","FILESDIR","PROGDIR")
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
        miles.block.fast(exp,year1,year2,season,z500filename,FILESDIR)
    }
}
