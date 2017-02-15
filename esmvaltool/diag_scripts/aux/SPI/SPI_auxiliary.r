##
## Standardized Precipitation Index helper functions
##
fit.gamma <- function (v) 
{
    ## Strategy: Try gamma fit with fitdistr, if doesn't work, use bordi
    ##           et al, 2001, anal geofisica
    alpha <- beta <- q <- gamma <- NA
    if (sum(!is.na(v)) >= 20) { ## Require at least 20 usable values
        v <- v[!is.na(v)]
        q <- sum(v == 0)/length(v)
        v <- v[v > 0]
        m <- mean(v)

        suppressPackageStartupMessages(require("MASS"))
        gamma <- try(fitdistr(v, "gamma"), silent = TRUE)
        {
            if (length(gamma) == 5) {
                alpha <- gamma$est[1]
                beta <- gamma$est[2]
            }
            else {
                A <- log(m) - sum(log(v))/length(v)
                alpha <- 1/(4 * A) * (1 + sqrt(1 + 4 * A/3))
                beta <- m/alpha
            }
        }
    }
    ## Return
    list(alpha = alpha, beta = beta, q = q, l = (gamma))
}

transform.to.spi <- function (v, gamma) 
{
    ## Cummulative function accounting for 0 precip
    Hv <- gamma$q + (1 - gamma$q) * pgamma(v, gamma$alpha, gamma$beta)
    qnorm(Hv)
}

calc.spi <- function (my.p, ref.inds, timescale) 
{
    info_output("<<<<<<<< Entering calc.spi()", verbosity, 8)
    {
        if (length(my.p) > 1) {
            ## Aggregate my.p over timescale
            agg.p <- my.p
            if (timescale > 1) {
                for (j in timescale:ncol(my.p)) { 
                    agg.p[, j] <- apply(my.p[, (j - timescale + 1):j], 1, sum)
                }
                agg.p[, 1:(timescale - 1)] <- NA
            }
            ## Fit gamma and spi per month
            gammas <- list()
            spis <- matrix(nrow = nrow(agg.p), ncol = ncol(agg.p))
            for (m in 1:12) {
                ## Gamma
                gamma.inds <- seq(m + ref.inds[1] - 1, max(ref.inds), by = 12)
                ## Remove the gamma.inds for which the agg.p-columns are set to NA
                gamma.inds <- setdiff(gamma.inds, 1:(timescale - 1))
                gammas[[m]] <- apply(agg.p[, gamma.inds], 1, fit.gamma)


                ## SPI
                spi.inds <- seq(m, ncol(agg.p), by = 12)
                spis[, spi.inds] <- t(sapply(1:nrow(agg.p), function(i) transform.to.spi(agg.p[i, 
                  spi.inds], gammas[[m]][[i]])))

                info_output(paste0("m=", m, " of 12"), verbosity, 9)
            }
            info_output(">>>>>>>> Leaving calc.spi()", verbosity, 8)
            spis
        } else { 
            info_output(">>>>>>>> Leaving calc.spi()", verbosity, 8)
            NA
       }
    }
}

##
## Function: Read data and calculate spi
##

## Input: Filename of monthly precip netcdf file, first and last time
##        step of this file yyyy/mm/dd, timescale in months, first and
##        last year of the reference period.

generate.spi <- function (fullpath_filename, filename,
        start_year, end_year, timescale, 
        begin.ref.year, end.ref.year, 
        model_idx) 
{
    info_output(filename, verbosity, 1)
    info_output(fullpath_filename, verbosity, 1)
    info_output("<<<<<<<< Entering generate.spi()", verbosity, 6)
    suppressPackageStartupMessages(require(ncdf4))
    suppressPackageStartupMessages(require(chron))
    precip.nc <- nc_open(fullpath_filename)

    ## Try to figure out lon/lat/time naming
    pot.lon.names <- c("lon", "Lon", "Longitude", "longitude", "X", "x")
    pot.lat.names <- c("lat", "Lat", "Latitude", "latitude", "Y", "y")
    pot.time.names <- c("time", "Time", "date", "Date", "T", "t")

    lon.name <- pot.lon.names[pot.lon.names %in% names(precip.nc$dim)]
    lat.name <- pot.lat.names[pot.lat.names %in% names(precip.nc$dim)]
    time.name <- pot.time.names[pot.lat.names %in% names(precip.nc$dim)]

    if (!(length(lon.name) > 0 && length(lat.name) > 0 && length(time.name) > 0
        && is.character(lon.name) > 0 && is.character(lat.name) > 0
        && is.character(time.name) > 0)) {

        stop(paste("fatal:Cannot guess dim names of", fullpath_filename))
    }

    lon <- ncvar_get(precip.nc, lon.name)
    lat <- ncvar_get(precip.nc, lat.name)
    time <- ncvar_get(precip.nc, time.name)

    ## Define time axis
    start_date <- as.Date(paste0(start_year, "-01-15"))
    end_date <- as.Date(paste0(end_year, "-12-15"))
    dates <- seq(start_date, end_date, by = "1 month")
    if (length(dates) != length(time)) 
        stop("fatal:Specified begin/end date and length of time dimension mismatch.")

    ## Check if reference period is covered by data 
    if (begin.ref.year < start_year || end.ref.year > end_year)
        stop("fatal:Specified reference period outside temporal coverage.")
    
    ## Reference period
    begin.ref <- as.Date(paste(begin.ref.year, "01/01", sep = "/"))
    end.ref <- as.Date(paste(end.ref.year, "12/31", sep = "/"))
    ref.inds <- which(findInterval(dates, c(begin.ref, end.ref)) == 1)

    ## Calculate spi
    my.p <- matrix(ncvar_get(precip.nc), ncol = length(dates))

    curr_model_str <- models_name[model_idx]
    output_str <- paste0("Now calculating spi for model ", curr_model_str)
    info_output(output_str, verbosity, 1)

    spi <- calc.spi(my.p, ref.inds, timescale)
    spi[!is.finite(spi)] <- NA

    ## Output to ncdf
    description <- "_values_wrt_"
    out.filename <- processed_spi_info(timescale, description,
                                       begin.ref.year, end.ref.year,
                                       filename)
    x <- ncdim_def("Lon", "degreesE", lon)
    y <- ncdim_def("Lat", "degreesN", lat)
    t <- ncdim_def("Time", "days since 1970-01-01", as.numeric(dates), unlim = TRUE)

    spi.var <- ncvar_def("SPI", "unitless", list(x, y, t), 1e+30)
    info_output(out.filename, verbosity, 1)
    out.filename <- file.path(climo_dir, out.filename)
    info_output(out.filename, verbosity, 1)
    ncid <- nc_create(out.filename, list(spi.var))
    ncvar_put(ncid, spi.var, spi, start = c(1, 1, 1), count = c(-1, -1, -1))
    nc_close(ncid)

    info_output(">>>>>>>> Leaving generate.spi()", verbosity, 6)
    ## Return spi, dates, lon, lat
    list(spi = spi, dates = dates, lon = lon, lat = lat, timescale = timescale, 
        begin.ref.year = begin.ref.year, end.ref.year = end.ref.year, 
        filename = filename)
}

##
## Function: Evaluate spi
##

## Input: A list like the one returned by generate.spi()

evaluate.spi <- function (spi.info) 
{
    info_output("<<<<<<<< Entering evaluate.spi()", verbosity, 6)
    if (!is.list(spi.info) || length(spi.info) != 8) 
        stop("fatal:Wrong input to evaluate spi.")

    spi <- spi.info$spi
    dates <- spi.info$dates
    lon <- spi.info$lon
    lat <- spi.info$lat
    timescale <- spi.info$timescale
    begin.ref.year <- spi.info$begin.ref.year
    end.ref.year <- spi.info$end.ref.year
    filename <- spi.info$filename

    if (length(lon) * length(lat) != nrow(spi) 
            || begin.ref.year >= end.ref.year 
            || length(dates) != ncol(spi)) {
        stop("fatal:Wrong input to evaluate spi.")
    }

    suppressPackageStartupMessages(require(ncdf4))
    suppressPackageStartupMessages(require(chron))

    ## Determine seasonal and annual time series
    my.seasons <- quarters(dates)
    my.seasons <- c(my.seasons[2:length(my.seasons)], NA)
    my.years <- years(dates)
    season.inds <- list(ann = 1:length(dates), 
        mam = which(my.seasons == "Q2"), jja = which(my.seasons == "Q3"),
        son = which(my.seasons == "Q4"), djf = which(my.seasons == "Q1"))

    ## Find grid points with not just missing values
    good.gps <- which(is.finite(apply(spi, 1, max, na.rm = TRUE)))

    ## Aggregate
    spi.agg <- list()
    for (s in names(season.inds)) {
        info_output(paste("aggregating for", s), verbosity, 4)
        spi.agg[[s]] <- matrix(nrow = nrow(spi), ncol = length(levels(my.years)))
        ii <- season.inds[[s]]
        for (i in good.gps) {
            if (s != "djf") 
                spi.agg[[s]][i, ] <- tapply(spi[i, ii], my.years[ii], 
                  mean)
            else {
                djf.years <- c((my.years[2:(length(my.years))]), 
                  NA)
                spi.agg[[s]][i, ] <- tapply(spi[i, ii], djf.years[ii], 
                  mean)
            }
        }
    }

    ## Calculate  linear trends for annual and seasonal time scale
    trends <- p.values <- list()
    for (s in names(season.inds)) {
        info_output(paste("calculating linear trends for", s), verbosity, 4)
        trends[[s]] <- p.values[[s]] <- rep(NA, nrow(spi.agg[[s]]))

        ## Find grid points with at least 70% coverage
        good.gps <- which(apply(!is.na(spi.agg[[s]]), 1, sum)/ncol(spi.agg[[s]]) > 0.7)
        for (i in good.gps) {
            lm <- lm(spi.agg[[s]][i, ] ~ as.numeric(levels(my.years)))
            ## Decadal change
            trends[[s]][i] <- coef(lm)[[2]] * 10
            ## p.value
            p.values[[s]][i] <- summary(lm)$coef[, "Pr(>|t|)"][[2]]
        }
    }

    ## Save to ncdf
    
    ## Create "time steps", indicating the beginning of the season and
    ## using last year as year.
    last.year <- max(levels(my.years))
    out.time <- c(paste(last.year,"01-01",sep="-"),
                  paste(last.year,"03-01",sep="-"),paste(last.year,"06-01",sep="-"),
                  paste(last.year,"09-01",sep="-"),paste(last.year,"12-01",sep="-"))
    out.time <- as.Date(out.time)

    description <- "_trends-n-pvals_wrt_"
    out.filename <- processed_spi_info(timescale, description,
                                       begin.ref.year, end.ref.year,
                                       filename)

    x <- ncdim_def("Lon", "degreesE", lon)
    y <- ncdim_def("Lat", "degreesN", lat)
    t <- ncdim_def("Time", "days since 1970-01-01", as.numeric(out.time), unlim = TRUE)

    trend.out <- as.vector(unlist(trends))
    pval.out <- as.vector(unlist(p.values))

    trend.var <- ncvar_def("Decadal_Trend", "unitless", list(x, y, t), 1e+30)
    pval.var <- ncvar_def("p_Value", "unitless", list(x, y, t), 1e+30)
    out.filename <- file.path(climo_dir, out.filename)
    ncid <- nc_create(out.filename, list(trend.var, pval.var))

    ncvar_put(ncid, trend.var, trend.out, start = c(1, 1, 1), count = c(-1, -1, -1))
    ncvar_put(ncid, pval.var, pval.out, start = c(1, 1, 1), count = c(-1, -1, -1))
    nc_close(ncid)

    info_output(">>>>>>>> Leaving evaluate.spi()", verbosity, 6)
    ## Return -- nothing
}

##
## Function: Define filenamed for SPI processed data
##
processed_spi_info <- function(timescale,
                               description,
                               begin.ref.year,
                               end.ref.year,
                               filename) {
    info_output("<<<<<<<< Entering processed_spi_info()", verbosity, 8)
    info_output(">>>>>>>> Leaving processed_spi_info()", verbosity, 8)
    return_processed_info <- paste("SPI", timescale, description, begin.ref.year, 
        "-", end.ref.year, "_", filename, sep = "")
}
