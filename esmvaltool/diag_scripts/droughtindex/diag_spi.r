library(yaml)
library(ncdf4)
library(SPEI)
getnc <- function(yml, m) {
  id <- nc_open(yml[m][[1]]$filename, readunlim=FALSE)
  v <- ncvar_get(id, yml[m][[1]]$short_name)
  nc_close(id)
  return(v)
}
ncwrite <- function(yml, m, data, wdir){
  fnam <- strsplit(yml[mod][[1]]$filename,"/")[[1]]
  pcs <- strsplit(fnam[length(fnam)],"_")[[1]]
  pcs[which(pcs==yml[m][[1]]$short_name)] = "spi"
  onam <- paste0(wdir,"/",paste(pcs,collapse="_"))
  system(paste0("cp ", yml[mod][[1]]$filename, " ",onam))
  idw <- nc_open(onam, write=TRUE)
    ncvar_put(idw, yml[m][[1]]$short_name, data)
  nc_close(idw)
}

args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
wdir <- params$work_dir
dir.create(wdir, recursive = TRUE)
var1_input <- read_yaml(params$input_files[1])
nmods <- length(names(var1_input))
for(mod in 1:nmods){
   v1 <- getnc(var1_input, mod)
   d <- dim(v1)
   v1_spi <- v1*NA
   for(i in 1:d[1]){
    print(i)
    tmp <- v1[i,,]
    v1_spi[i,,] <- t(spi(t(tmp), 1, na.rm=TRUE)$fitted)
   }
   ncwrite(var1_input, mod, v1_spi, wdir)
}
