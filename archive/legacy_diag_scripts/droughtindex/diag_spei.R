library(yaml)
library(ncdf4)
library(SPEI)
library(RColorBrewer) # nolint

leap_year <- function(year) {
  return(ifelse((year %% 4 == 0 & year %% 100 != 0) |
    year %% 400 == 0, TRUE, FALSE))
}

getnc <- function(yml, m, lat = FALSE) {
  id <- nc_open(yml[m][[1]]$filename, readunlim = FALSE)
  if (lat) {
    v <- ncvar_get(id, "lat")
  } else {
    v <- ncvar_get(id, yml[m][[1]]$short_name)
    if (yml[m][[1]]$short_name == "tas") {
      v <- v - 273.15
    }
    if (yml[m][[1]]$short_name == "pr") {
      time <- ncvar_get(id, "time")
      tcal <- ncatt_get(id, "time", attname = "calendar")
      tunits <- ncatt_get(id, "time", attname = "units")
      tustr <- strsplit(tunits$value, " ")
      stdate <- as.Date(time[1], origin = unlist(tustr)[3])
      nddate <-
        as.Date(time[length(time)], origin = unlist(tustr)[3])
      if (tcal$value == "365_day") {
        # Correct for missing leap years in nddate
        diff <- as.numeric(nddate - stdate, units = "days")
        dcorr <- floor((diff / 365 - diff / 365.25) * 365.25)
        nddate <- nddate + dcorr
      }
      if (tcal$value == "360_day") {
        v <- v * 30 * 24 * 3600.
      } else {
        cnt <- 1
        monarr <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
        date <- stdate
        while (date <= nddate) {
          year <- as.numeric(substr(date, 1, 4))
          lpyear <- leap_year(year)
          month <- as.numeric(substr(date, 6, 7))
          mdays <- monarr[month]
          pdays <- mdays
          if (month == 2 & lpyear == TRUE) {
            pdays <- 29
            if (tcal$value != "365_day") {
              mdays <- 29
            } else {
              mdays <- 28
            }
          }
          v[, , cnt] <- v[, , cnt] * mdays * 24 * 3600.
          date <- date + pdays
          cnt <- cnt + 1
        }
      }
    }
  }
  nc_close(id)
  return(v)
}

ncwritenew <- function(yml, m, hist, wdir, bins) {
  fnam <- strsplit(yml[m][[1]]$filename, "/")[[1]]
  pcs <- strsplit(fnam[length(fnam)], "_")[[1]]
  pcs[which(pcs == yml[m][[1]]$short_name)] <- "spei"
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

dothornthwaite <- function(v, lat) {
  print("Estimating PET with Thornthwaite method.")
  dpet <- v * NA
  d <- dim(dpet)
  for (i in 1:d[2]) {
    tmp <- v[, i, ]
    tmp2 <- thornthwaite(t(tmp), rep(lat[i], d[1]), na.rm = TRUE)
    d2 <- dim(tmp2)
    tmp2 <- as.numeric(tmp2)
    dim(tmp2) <- d2
    dpet[, i, ] <- t(tmp2)
  }
  return(dpet)
}

args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
metadata1 <- read_yaml(params$input_files[1])
metadata2 <- read_yaml(params$input_files[2])
modfile1 <- names(metadata1)
modfile2 <- names(metadata2)
wdir <- params$work_dir
rundir <- params$run_dir
dir.create(wdir, recursive = TRUE)
pdir <- params$plot_dir
dir.create(pdir, recursive = TRUE)
var1_input <- read_yaml(params$input_files[1])
var2_input <- read_yaml(params$input_files[2])
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
  references = list("vicente10jclim"),
  projects = list("c3s-magic"),
  caption = "",
  statistics = list("other"),
  realms = list("atmos"),
  themes = list("phys"),
  domains = list("global")
)

histarr <- array(NA, c(nmods, length(histnams)))
for (mod in 1:nmods) {
  lat <- getnc(var1_input, mod, TRUE)
  v1 <- getnc(var1_input, mod, FALSE)
  v2 <- getnc(var2_input, mod, FALSE)
  if (var1_input[1][[1]]$short_name == "pr") {
    prtas <- TRUE
  } else {
    prtas <- FALSE
  }
  if (prtas) {
    pet <- dothornthwaite(v2, lat)
    pme <- v1 - pet
  } else {
    pet <- dothornthwaite(v1, lat)
    pme <- v2 - pet
  }
  print(var1_input[mod][[1]]$cmor_table)
  d <- dim(pme)
  pme_spei <- pme * NA
  for (i in 1:d[1]) {
    wh <- which(!is.na(refmsk[i, ]))
    if (length(wh) > 0) {
      tmp <- pme[i, wh, ]
      pme_spei[i, wh, ] <- t(spei(t(tmp), 1, na.rm = TRUE)$fitted)
    }
  }
  pme_spei[is.infinite(pme_spei)] <- NA
  pme_spei[pme_spei > 10000] <- NA
  hist_spei <- array(NA, c(d[1], d[2], length(histbrks) - 1))
  for (nnh in 1:(length(histbrks) - 1)) {
    hist_spei[, , nnh] <- apply(
      pme_spei,
      c(1, 2),
      FUN = whfcn,
      ilow = histbrks[nnh],
      ihigh = histbrks[nnh + 1]
    )
  }
  filename <-
    ncwritenew(var1_input, mod, hist_spei, wdir, histbrks)
  # Set provenance for output files
  xprov$caption <- "Histogram of SPEI index per grid point."
  xprov$ancestors <- list(modfile1[mod], modfile2[mod])
  provenance[[filename]] <- xprov
  for (t in 1:d[3]) {
    tmp <- pme_spei[, , t]
    tmp[is.na(refmsk)] <- NA
    pme_spei[, , t] <- tmp
  }
  pme_spei[is.infinite(pme_spei)] <- NA
  pme_spei[pme_spei > 10000] <- NA
  # Weight against latitude
  h <- seq_along(histnams) * 0
  for (j in 1:d[2]) {
    h <- h + hist(pme_spei[j, , ],
      breaks = histbrks,
      plot = FALSE
    )$counts * cos(lat[j] * pi / 180.)
  }
  histarr[mod, ] <- h / sum(h, na.rm = TRUE)
}
filehist <- paste0(params$work_dir, "/", "histarr.rsav")
save(histarr, file = filehist)
plot_file <- paste0(params$plot_dir, "/", "histplot.png")
xprov$caption <- "Global latitude-weighted histogram of SPEI index."
xprov$ancestors <- c(modfile1, modfile2)
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
  "Standardized precipitation-evapotranspiration index",
  outer = TRUE,
  cex = 2,
  font = 2
)
par(
  fig = c(0.8, .95, 0.1, 0.9),
  new = T,
  oma = c(1, 1, 1, 1) * 0,
  mar = c(0, 0, 0, 0)
)
legend("topright", mnam[parr], fill = cols)
dev.off()
