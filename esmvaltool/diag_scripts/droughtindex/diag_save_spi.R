library(yaml)
library(ncdf4)
library(SPEI)
library(RColorBrewer) # nolint

getnc <- function(yml, m, lat = FALSE) {
  id <- nc_open(yml[m][[1]]$filename, readunlim = FALSE)
  if (lat){
    v <- ncvar_get(id, "lat")
  }else{
    v <- ncvar_get(id, yml[m][[1]]$short_name)
  }
  nc_close(id)
  return(v)
}

ncwritespi <- function(yml, m, data, wdir){
  fnam <- strsplit(yml[m][[1]]$filename, "/")[[1]]
  pcs <- strsplit(fnam[length(fnam)], "_")[[1]]
  pcs[which(pcs == yml[m][[1]]$short_name)] <- "spi"
  onam <- paste0(wdir, "/", paste(pcs, collapse = "_"))
  ncid_in <- nc_open(yml[m][[1]]$filename)
  # var <- ncid_in$var[[yml[m][[1]]$short_name]]
  xdim <- ncid_in$dim[["lon"]]
  ydim <- ncid_in$dim[["lat"]]
  tdim <- ncid_in$dim[["time"]]
  allatt <- ncatt_get(ncid_in, "pr")
  fillvalue <- ncatt_get(ncid_in,"pr","_FillValue")
  globat <- ncatt_get(ncid_in, 0)
  fillfloat <- 1.e+20
  as.single(fillfloat)
  var_spi <- ncvar_def("spi", "1", list(xdim, ydim, tdim), fillfloat)
  idw <- nc_create(onam, var_spi)
  ncvar_put(idw, "spi", data)
  cntatt <- 1
  for (thisattname in names(globat)){
    ncatt_put(idw, 0, thisattname, globat[[cntatt]])
    cntatt <- cntatt + 1
  }
  nc_close(idw)
  nc_close(ncid_in)
  return(onam)
}

whfcn <- function(x, ilow, ihigh){
  return(length(which(x >= ilow & x < ihigh)))
}

args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
metadata <- read_yaml(params$input_files)
modfile <- names(metadata)
wdir <- params$work_dir
dir.create(wdir, recursive = TRUE)
rundir <- params$run_dir
pdir <- params$plot_dir
dir.create(pdir, recursive = TRUE)
var1_input <- read_yaml(params$input_files[1])

nmods <- length(names(var1_input))

fillfloat <- 1.e+20
as.single(fillfloat)

# setup provenance file and list
provenance_file <- paste0(rundir, "/", "diagnostic_provenance.yml")
provenance <- list()


refnam <- var1_input[1][[1]]$reference_dataset
n <- 1
while (n <= nmods){
  if (var1_input[n][[1]]$dataset == refnam) break
  n <- n + 1
}
nref <- n
lat <- getnc(var1_input, nref, lat = TRUE)
if (max(lat) > 90){
  print(paste0("Latitude must be [-90,90]: min=",
  min(lat), " max=", max(lat)))
  stop("Aborting!")
}
ref <- getnc(var1_input, nref, lat = FALSE)
refmsk <- apply(ref, c(1, 2), FUN = mean, na.rm = TRUE)
refmsk[refmsk > 10000] <- fillfloat
refmsk[!is.na(refmsk)] <- 1

xprov <- list(
  ancestors = list(""),
  authors = list("weigel_katja"),
  references = list("mckee93proc"),
  projects = list("eval4cmip"),
  caption = "",
  statistics = list("other"),
  realms = list("atmos"),
  themes = list("phys"),
  domains = list("global")
)

for (mod in 1:nmods){
   v1 <- getnc(var1_input, mod)
   d <- dim(v1)
   v1_spi <- array(fillfloat, dim=d)
   for (i in 1:d[1]){
     wh <- which(!is.na(refmsk[i,]))
     if (length(wh) > 0){
       tmp <- v1[i,wh,]
       v1_spi[i,wh,] <- t(spi(t(tmp), params$smooth_month, na.rm = TRUE,
                        distribution = params$distribution)$fitted)
     }
   }
   v1_spi[is.infinite(v1_spi)] <- fillfloat
   v1_spi[is.na(v1_spi)] <- fillfloat
   v1_spi[v1_spi > 10000] <- fillfloat
   filename <- ncwritespi(var1_input, mod, v1_spi, wdir)
   xprov$caption <- "SPI index per grid point."
   xprov$ancestors <- list(modfile[mod])
   provenance[[filename]] <- xprov
   print("provenance[[filename]] kwnew")
   print(provenance[[filename]])
}

print("provenance_file kwnew")
print(provenance_file)
print("provenance kwnew")
print(provenance)
write_yaml(provenance, provenance_file)
