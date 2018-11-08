library(yaml)
library(ncdf4)
library(SPEI)
library(RColorBrewer)

leap_year <- function(year) {
  return(ifelse((year %%4 == 0 & year %%100 != 0) | 
                 year %%400 == 0, TRUE, FALSE))
}

getnc <- function(yml, m, lat=FALSE) {
  id <- nc_open(yml[m][[1]]$filename, readunlim=FALSE)
  if(lat){
     v <- ncvar_get(id, 'lat')
  }else{
   v <- ncvar_get(id, yml[m][[1]]$short_name)
   #print(dim(v))
   if(yml[m][[1]]$short_name == "tas") v <- v - 273.15
   if(yml[m][[1]]$short_name == "pr"){
     time <- ncvar_get(id,"time")
     tcal <- ncatt_get(id,"time",attname="calendar")
     tunits <- ncatt_get(id,"time",attname="units")
     tustr <- strsplit(tunits$value, " ")
     stdate <- as.Date(time[1],origin=unlist(tustr)[3])
     nddate <- as.Date(time[length(time)],origin=unlist(tustr)[3])
     if(tcal$value=="365_day"){
       # Correct for missing leap years in nddate
       diff <- as.numeric(nddate-stdate, units="days")
       dcorr <- floor((diff/365 - diff/365.25)*365.25)
       nddate <- nddate + dcorr
     }
     if(tcal$value=="360_day"){
       v <- v * 30 * 24 * 3600.
     }else{
       cnt <- 1
       monarr <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
       date <- stdate
       while(date <= nddate){
         year <- as.numeric(substr(date,1,4))
	 lpyear <- leap_year(year)
         month <- as.numeric(substr(date,6,7))
         mdays = monarr[month]
	 pdays <- mdays
         if(tcal$value != "365_day" & month==2 & lpyear == TRUE){
           mdays <- 29
	   pdays <- 29
         }else{
           pdays <- 29
         }
         v[,,cnt] <- v[,,cnt] * mdays * 24 * 3600.
	 date <- date + pdays
       }
     }
   }
  }
  nc_close(id)
  return(v)
}
ncwrite <- function(yml, m, data, wdir){
  fnam <- strsplit(yml[mod][[1]]$filename,"/")[[1]]
  pcs <- strsplit(fnam[length(fnam)],"_")[[1]]
  pcs[which(pcs==yml[m][[1]]$short_name)] = "spei"
  onam <- paste0(wdir,"/",paste(pcs,collapse="_"))
  system(paste0("cp ", yml[mod][[1]]$filename, " ",onam))
  idw <- nc_open(onam, write=TRUE)
    ncvar_put(idw, yml[m][[1]]$short_name, data)
  nc_close(idw)
}
dothornthwaite <- function(v, lat){
 print("Estimating PET wit Thorntwaite method.")
 dpet <- v*NA
 d <- dim(dpet)
 for(i in 1:d[2]){
  tmp <- v[,i,]
  tmp2 <- thornthwaite(t(tmp), rep(lat[i],d[1]), na.rm = TRUE)
  d2 <- dim(tmp2)
  tmp2 <- as.numeric(tmp2)
  dim(tmp2) <- d2
  dpet[,i,] <- t(tmp2) #thornthwaite(t(v), rep(lat[i],d[1]))
 }
 return(dpet)
}

args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
wdir <- params$work_dir
dir.create(wdir, recursive = TRUE)
pdir <- params$plot_dir
dir.create(pdir, recursive = TRUE)
var1_input <- read_yaml(params$input_files[1])
var2_input <- read_yaml(params$input_files[2])
nmods <- length(names(var1_input))

histbrks <- c(-99999, -2, -1.5, -1, 1, 1.5, 2, 99999)
histnams <- c("Extremely dry", "Moderately dry", "Dry",
              "Neutral",
              "Wet", "Moderately wet", "Extremely wet")

print(names(var1_input))
for(n in 1:nmods){
   if(var1_input[n][[1]]$cmor_table == "OBS"){
      ref <- getnc(var1_input, n)
      refmsk <- apply(ref,c(1,2),FUN=mean)
      refmsk[refmsk > 10000] = NA
      refmsk[!is.na(refmsk)]=1
      nref <- n
   }
}

histarr <- array(NA, c(nmods, length(histnams)))
for(mod in 1:nmods){
   lat <- getnc(var1_input, mod, TRUE)
   v1 <- getnc(var1_input, mod, FALSE)
   v2 <- getnc(var2_input, mod, FALSE)
   if(var1_input[1][[1]]$short_name == "pr") prtas <- TRUE else prtas <- FALSE
   if(prtas){
     #print(prtas)
     pet <- dothornthwaite(v2, lat)
     pme <- v1 - pet
   }else{
     #print(prtas)
     pet <- dothornthwaite(v1, lat)
     pme <- v2 - pet
   }
   print(var1_input[mod][[1]]$cmor_table)
   d <- dim(pme)
   pme_spei <- pme*NA
   for(i in 1:d[1]){
    tmp <- v1[i,,]
    pme_spei[i,,] <- t(spei(t(tmp), 1, na.rm=TRUE)$fitted)
   }
   print(summary(as.numeric(pme_spei)))
   pme_spei[is.infinite(pme_spei)] = NA
   print(summary(as.numeric(pme_spei)))
   pme_spei[pme_spei > 10000] = NA
   print(summary(as.numeric(pme_spei)))

   ncwrite(var1_input, mod, pme_spei, wdir)
   for(t in 1:d[3]){
      tmp <- pme_spei[,,t]
      tmp[is.na(refmsk)] = NA
      pme_spei[,,t] <- tmp
   }#t
   pme_spei[is.infinite(pme_spei)] = NA
   pme_spei[pme_spei > 10000] = NA
   # Should weight against latitude!
   print(summary(as.numeric(pme_spei)))
   print(histbrks)
   h <- hist(pme_spei, breaks=histbrks,
                         plot=FALSE)$density
   histarr[mod,] <- h/sum(h)
}#mod
save(histarr, file=paste0(params$work_dir,
                          "/histarr.rsav"))

bhistarr <- array(NA,c(nmods-1,7))
marr <- c(1:nmods)[c(1:nmods) != nref]
cnt <- 1
for(m in marr){
  bhistarr[cnt,] <- histarr[m,]-histarr[nref,]
  cnt <- cnt + 1
}
parr <- c(nref,marr)

mnam <- c(1:nmods)*NA
for(m in 1:nmods){mnam[m] <- var1_input[m][[1]]$dataset}

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                           rownames(qual_col_pals)))
cols=c("black",sample(col_vector, nmods-1))

png(paste0(params$plot_dir,"/histplot.png"),
    width=1000, height=500)
 par(mfrow=c(2,1),oma=c(3,3,3,13),mar=c(2,1,1,1))
 barplot(histarr[parr,], beside=1, names.arg=histnams,
         col=cols, xaxs="i")
 box()
 mtext("Probability", side=2, line=2.1)
 barplot(bhistarr, beside=1, names.arg=histnams,
         col=cols[2:17], xaxs="i")
 box()
 mtext("Absolute difference", side=2, line=2.1)
 mtext("Standardized precipitation-evapotranspiration index",
        outer=TRUE, cex=2, font=2)
 par(fig=c(0.8,.95,0.1,0.9),new=T,oma=c(1,1,1,1)*0,mar=c(0,0,0,0))
 legend("topright", mnam[parr], fill=cols)
dev.off()



