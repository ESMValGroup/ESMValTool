library(yaml)
library(ncdf4)
library(SPEI)
getnc <- function(yml, m, lat=FALSE) {
  id <- nc_open(yml[m][[1]]$filename, readunlim=FALSE)
  if(lat){
     v <- ncvar_get(id, 'lat')
  }else{
   v <- ncvar_get(id, yml[m][[1]]$short_name)
   if(yml[m][[1]]$short_name == "tas") v <- v - 273.15
   if(yml[m][[1]]$short_name == "pr"){
     #timeatt <- ncatt_get(id, 'time')
     #timeunit <- timeatt$units
     #if(timeatt$calendar == 360_day){
     v <- v*30*24*3600
     print("Warning: assuming 30 days per month!")
     #}else if(timeatt$calendar == 365_day){
     #  units(time) <- make_unit(timeunit)
     #}else{
     #}
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
print("dothorn...")
 dpet <- v*NA
 d <- dim(dpet)
 for(i in 1:d[2]){
  tmp <- v[,i,]
  tmp2 <- thornthwaite(t(tmp), rep(lat[i],d[1]))
  d2 <- dim(tmp2)
  tmp2 <- as.numeric(tmp2)
  dim(tmp2) <- d2
  dpet[,i,] <- t(tmp2) #thornthwaite(t(v), rep(lat[i],d[1]))
 }
 return(dpet)
}

#Parsing input file paths and creating output dirs
args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
wdir <- params$work_dir
dir.create(wdir, recursive = TRUE)
var1_input <- read_yaml(params$input_files[1])
var2_input <- read_yaml(params$input_files[2])
nmods <- length(names(var1_input))
if(length(names(var2_input)) != nmods){
  stop(paste0("Number of models not equal for both variables! Aborting!"))
}

for(mod in 1:nmods){
   lat <- getnc(var1_input, mod, TRUE)
   v1 <- getnc(var1_input, mod, FALSE)
   v2 <- getnc(var2_input, mod, FALSE)
   if(var1_input[1][[1]]$short_name == "pr") prtas <- TRUE else prtas <- FALSE
   if(prtas){
   print(prtas)
     pet <- dothornthwaite(v2, lat)
     pme <- v1 - pet
   }else{
   print(prtas)
     pet <- dothornthwaite(v1, lat)
     pme <- v2 - pet
   }
   d <- dim(pme)
   v1_spei <- pme*NA
   for(i in 1:d[1]){
    print(i)
    tmp <- pme[i,,]
    v1_spei[i,,] <- t(spei(t(tmp), 1, na.rm=TRUE)$fitted)
   }
   #if(dim(d)>2){
   #  dim(pme) <- c(d[1]*d[2],d[3])
   #}else{
   #  dim(pme) <- c(d[1]*d[2])
   #}
   #spei <- spei(pme, 1)
   ncwrite(var1_input, mod, v1_spei, wdir)
}


