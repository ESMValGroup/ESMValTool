library(yaml)
library(ncdf4)
library(SPEI)
library(RColorBrewer) # nolint

getnc <- function(yml, m, lat = FALSE) {
  id <- nc_open(yml[m][[1]]$filename, readunlim = FALSE)
  if (lat) {
    v <- ncvar_get(id, "lat")
  } else {
    v <- ncvar_get(id, yml[m][[1]]$short_name)
  }
  nc_close(id)
  return(v)
}

ncwritenew <- function(yml, m, hist, wdir, bins) {
  fnam <- strsplit(yml[m][[1]]$filename, "/")[[1]]
  pcs <- strsplit(fnam[length(fnam)], "_")[[1]]
  pcs[which(pcs == yml[m][[1]]$short_name)] <- "spi"
  onam <- paste(pcs, collapse = "_")
  onam <- paste0(wdir, "/", strsplit(onam, ".nc"), "_hist.nc")
  ncid_in <- nc_open(yml[m][[1]]$filename)
  var <- ncid_in$var[[yml[m][[1]]$short_name]]
  xdim <- ncid_in$dim[["lon"]]
  ydim <- ncid_in$dim[["lat"]]
  hdim <- ncdim_def("bins", "level", bins[1:(length(bins) - 1)])
  hdim2 <- ncdim_def("binsup", "level", bins[2:length(bins)])
  var_hist <-
    ncvar_def("hist", "counts", list(xdim, ydim, hdim), NA)
  idw <- nc_create(onam, var_hist)
  ncvar_put(idw, "hist", hist)
  nc_close(idw)
  return(onam)
}

whfcn <- function(x, ilow, ihigh) {
  return(length(which(x >= ilow & x < ihigh)))
}

args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
metadata <- read_yaml(params$input_files)
modfile <- names(metadata)
wdir <- params$work_dir
rundir <- params$run_dir
dir.create(wdir, recursive = TRUE)
pdir <- params$plot_dir
dir.create(pdir, recursive = TRUE)
var1_input <- read_yaml(params$input_files[1])
nmods <- length(names(var1_input))

# setup provenance file and list
provenance_file <- paste0(rundir, "/", "diagnostic_provenance.yml")
provenance <- list()

histbrks <- c(-99999, -2, -1.5, -1, 1, 1.5, 2, 99999)
histnams <- c(
  "Extremely dry",
  "Moderately dry",
  "Dry",
  "Neutral",
  "Wet",
  "Moderately wet",
  "Extremely wet"
)
refnam <- var1_input[1][[1]]$reference_dataset
n <- 1
while (n <= nmods) {
  if (var1_input[n][[1]]$dataset == refnam) {
    break
  }
  n <- n + 1
}
nref <- n
lat <- getnc(var1_input, nref, lat = TRUE)
if (max(lat) > 90) {
  print(paste0(
    "Latitude must be [-90,90]: min=",
    min(lat), " max=", max(lat)
  ))
  stop("Aborting!")
}
ref <- getnc(var1_input, nref, lat = FALSE)
refmsk <- apply(ref, c(1, 2), FUN = mean, na.rm = TRUE)
refmsk[refmsk > 10000] <- NA
refmsk[!is.na(refmsk)] <- 1

xprov <- list(
  ancestors = list(""),
  authors = list("berg_peter"),
  references = list("mckee93proc"),
  projects = list("c3s-magic"),
  caption = "",
  statistics = list("other"),
  realms = list("atmos"),
  themes = list("phys"),
  domains = list("global")
)

histarr <- array(NA, c(nmods, length(histnams)))
for (mod in 1:nmods) {
  v1 <- getnc(var1_input, mod)
  print(var1_input[mod][[1]]$cmor_table)
  d <- dim(v1)
  v1_spi <- v1 * NA
  for (i in 1:d[1]) {
    wh <- which(!is.na(refmsk[i, ]))
    if (length(wh) > 0) {
      tmp <- v1[i, wh, ]
      v1_spi[i, wh, ] <- t(spi(t(tmp), 1,
        na.rm = TRUE,
        distribution = "PearsonIII"
      )$fitted)
    }
  }
  v1_spi[is.infinite(v1_spi)] <- NA
  v1_spi[v1_spi > 10000] <- NA
  hist_spi <- array(NA, c(d[1], d[2], length(histbrks) - 1))
  for (nnh in 1:(length(histbrks) - 1)) {
    hist_spi[, , nnh] <- apply(v1_spi,
      c(1, 2),
      FUN = whfcn,
      ilow = histbrks[nnh],
      ihigh = histbrks[nnh + 1]
    )
  }
  filename <- ncwritenew(var1_input, mod, hist_spi, wdir, histbrks)
  # Set provenance for output files
  xprov$caption <- "Histogram of SPI index per grid point."
  xprov$ancestors <- modfile[mod]
  provenance[[filename]] <- xprov
  # Weight against latitude
  h <- seq_along(histnams) * 0
  for (j in 1:d[2]) {
    h <- h + hist(v1_spi[j, , ],
      breaks = histbrks,
      plot = FALSE
    )$counts * cos(lat[j] * pi / 180.)
  }
  histarr[mod, ] <- h / sum(h, na.rm = TRUE)
}
filehist <- paste0(params$work_dir, "/", "histarr.rsav")
save(histarr, file = filehist)
plot_file <- paste0(params$plot_dir, "/", "histplot.png")
xprov$caption <- "Global latitude-weighted histogram of SPI index."
xprov$ancestors <- modfile
provenance[[plot_file]] <- xprov
provenance[[filehist]] <- xprov
write_yaml(provenance, provenance_file)

bhistarr <- array(NA, c(nmods - 1, 7))
marr <- c(1:nmods)[c(1:nmods) != nref]
cnt <- 1
for (m in marr) {
  bhistarr[cnt, ] <- histarr[m, ] - histarr[nref, ]
  cnt <- cnt + 1
}
parr <- c(nref, marr)

mnam <- c(1:nmods) * NA
for (m in 1:nmods) {
  mnam[m] <- var1_input[m][[1]]$dataset
}

qual_col_pals <-
  brewer.pal.info[brewer.pal.info$category == "qual", ] # nolint
col_vector <-
  unlist(mapply(
    brewer.pal, qual_col_pals$maxcolors, # nolint
    rownames(qual_col_pals)
  ))
cols <- c("black", sample(col_vector, nmods - 1))

png(plot_file, width = 1000, height = 500)
par(
  mfrow = c(2, 1),
  oma = c(3, 3, 3, 13),
  mar = c(2, 1, 1, 1)
)
barplot(
  histarr[parr, ],
  beside = 1,
  names.arg = histnams,
  col = cols,
  xaxs = "i"
)
box()
mtext("Probability", side = 2, line = 2.1)
barplot(
  bhistarr,
  beside = 1,
  names.arg = histnams,
  col = cols[2:nmods],
  xaxs = "i"
)
box()
mtext("Absolute difference", side = 2, line = 2.1)
mtext(
  "Standardized precipitation index",
  outer = TRUE,
  cex = 2,
  font = 2
)
par(
  fig = c(0.8, .95, 0.1, 0.9),
  new = T,
  oma = c(0, 0, 0, 0),
  mar = c(0, 0, 0, 0)
)
legend("topright", mnam[parr], fill = cols)
dev.off()
