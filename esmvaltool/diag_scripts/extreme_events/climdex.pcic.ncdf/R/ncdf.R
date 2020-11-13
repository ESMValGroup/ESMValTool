# nolint start
#' climdex.pcic.ncdf, a package to calculate Climdex indices from NetCDF files.
#'
#' This package implements code to facilitate computation of Climdex indices
#' from NetCDF input files.
#'
#' The Climdex climate extremes indices have historically been calculated using
#' Fortran code. This has a number of problems:\itemize{
#' \item{Difficult to test}
#' \item{Difficult to modify (for instance, to add NetCDF file I/O)}
#' \item{Difficult to parallelize}
#' }
#' The \code{climdex.pcic} package provides an easy interface to efficient
#' computation of Climdex indices. This package is complementary to it, providing
#' easy access to functions to compute indices in parallel, using NetCDF files as
#' input and output. It implements chunked processing of input files to keep memory
#' usage reasonable; it implements parallel computation using the \code{snow}
#' library; and it includes a test suite to verify correctness of the implementation.
#' Furthermore, the package has a modular design, allowing for easy extension to
#' allow for adaptation to changing or custom requirements.
#'
#' Users of this package should pay particular attention to the
#' \code{\link{create.indices.from.files}} function, which computes Climdex indices
#' given NetCDF input files; and \code{\link{create.thresholds.from.file}}, which
#' computes thresholds for use with threshold-based indices given NetCDF input files.
#' Many of the other functions exposed by the package are intended to provide for
#' extensibility, but are unlikely to be routinely used by users of this package.
#'
#' @name climdex.pcic.ncdf
#' @aliases climdex.pcic.ncdf-package
#' @docType package
#' @seealso \code{\link{create.indices.from.files}}, \code{\link{create.thresholds.from.file}}
#' @references \url{http://etccdi.pacificclimate.org/list_27_indices.shtml}
#'
#' Karl, T.R., N. Nicholls, and A. Ghazi, 1999: CLIVAR/GCOS/WMO workshop on
#' indices and indicators for climate extremes: Workshop summary. Climatic
#' Change, 42, 3-7.
#'
#' Peterson, T.C., and Coauthors: Report on the Activities of the Working Group
#' on Climate Change Detection and Related Rapporteurs 1998-2001. WMO, Rep.
#' WCDMP-47, WMO-TD 1071, Geneve, Switzerland, 143pp.
#'
#' Zhang, X., 2005: Avoiding inhomogeneity in percentile-based indices of
#' temperature extremes. Journal of Climate 18.11 (2005):1641-.
#' @keywords climate ts
#' @importClassesFrom climdex.pcic climdexInput
#' @import snow PCICt
NULL

## Parallel lapply across 'x', running remote.func, and filtering with local.filter.func .
## Processing is incremental, not batch, to improve parallel throughput and reduce memory consumption.
parLapplyLBFiltered <- function(cl, x, remote.func, ..., local.filter.func=NULL) {
  snow::checkCluster(cl)
  cluster.size <- length(cl)
  num.tasks <- length(x)
  if(num.tasks == 0)
    return(list())
  if(cluster.size == 0)
    stop("Impossible happened; cluster size = 0")

  data.to.return <- vector("list", num.tasks)

  submit.job <- function(cluster.id, task.id) {
    snow::sendCall(cl[[cluster.id]], remote.func, args=c(x[task.id], list(...)), tag=task.id)
  }

  ## Fire off jobs, filling in the cur.task table as we go.
  for(i in 1:min(cluster.size, num.tasks))
    submit.job(i, i)

  next.task <- min(cluster.size, num.tasks)

  ## Stalk and feed jobs
  for(i in 1:num.tasks) {
    d <- snow::recvOneResult(cl)
    next.task <- next.task + 1

    ## Feed the finished node another task if we have one.
    if(next.task <= num.tasks)
      submit.job(d$node, next.task)

    if(!is.null(local.filter.func))
      data.to.return[d$tag] <- list(local.filter.func(d$value, x[[d$tag]]))
    else
      data.to.return[d$tag] <- list(d$value)

    rm(d)
  }

  ## Return data when complete
  return(data.to.return)
}

put.history.att <- function(f, v, definemode=FALSE) {
  history.string <- paste("Created by climdex.pcic", packageVersion("climdex.pcic"), "on", date())
  ncdf4::ncatt_put(f, v, "history", history.string, definemode=definemode)
  invisible(0)
}

put.ETCCDI.atts <- function(f, freq, orig.title, author.data, definemode=FALSE) {
  if("institution" %in% names(author.data))
    ncdf4::ncatt_put(f, 0, "ETCCDI_institution", author.data$institution, definemode=definemode)
  if("institution_id" %in% names(author.data))
    ncdf4::ncatt_put(f, 0, "ETCCDI_institution_id", author.data$institution_id, definemode=definemode)
  if("indices_archive" %in% names(author.data))
    ncdf4::ncatt_put(f, 0, "ETCCDI_indices_archive", author.data$indices_archive, definemode=definemode)

  ncdf4::ncatt_put(f, 0, "ETCCDI_software", "climdex.pcic", definemode=definemode)
  ncdf4::ncatt_put(f, 0, "ETCCDI_software_version", as.character(packageVersion("climdex.pcic")), definemode=definemode)

  if("contact" %in% names(author.data))
    ncdf4::ncatt_put(f, 0, "contact", author.data$contact, definemode=definemode)
  if("references" %in% names(author.data))
    ncdf4::ncatt_put(f, 0, "references", author.data$references, definemode=definemode)

  ncdf4::ncatt_put(f, 0, "frequency", freq, definemode=definemode)
  ncdf4::ncatt_put(f, 0, "creation_date", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="GMT"), definemode=definemode)
  ncdf4::ncatt_put(f, 0, "title", paste("ETCCDI indices computed on", orig.title), definemode=definemode)
  invisible(0)
}

all.the.same <- function(dat) {
  ifelse(length(dat) == 1, TRUE, all(unlist(lapply(dat, identical, dat[[1]]))))
}

#' Creates a list of CMIP5-compliant filenames reflecting the input data.
#'
#' Creates a list of CMIP5-compliant filenames reflecting the input data.
#'
#' This function takes a split filename (as created by \code{get.split.filename.cmip5}) and a list of variables and creates corresponding filenames for the given variables.
#'
#' @param fn.split A vector containing named components, as created by \code{get.split.filename.cmip5}.
#' @param vars.list A vector containing names of variables, as created by \code{\link{get.climdex.variable.list}}.
#' @return A vector containing filenames corresponding to the variables and filename bits supplied.
#'
#' @examples
#' \dontrun{
#' library(ncdf4.helpers)
#' ## Split out filename bits for use below...
#' fn <- "pr_day_BCCAQ+ANUSPLIN300+MRI-CGCM3_historical+rcp85_r1i1p1_19500101-21001231.nc"
#' fn.split <- get.split.filename.cmip5(fn)
#'
#' ## Create filenames with time data and variable appropriately replaced.
#' filenames <- create.climdex.cmip5.filenames(fn.split, c("rx5dayETCCDI_mon", "tn90pETCCDI_yr"))
#' }
#'
#' @export
create.climdex.cmip5.filenames <- function(fn.split, vars.list) {
  time.res <- c("yr", "mon")[grepl("_mon$", vars.list) + 1]
  time.range <- substr(fn.split[c('tstart', 'tend')], 1, 4)

  paste(paste(vars.list, fn.split['model'], fn.split['emissions'], fn.split['run'], sapply(time.res, function(x) { paste(time.range, switch(x, yr=c("", ""), mon=c("01", "12")), sep="", collapse="-") }), sep="_"), ".nc", sep="")
}

#' Returns a list of Climdex variables given constraints
#'
#' Returns a list of Climdex variables given constraints.
#'
#' This function takes a character vector which specifies what source data is present and a time resolution, and generates a list of names consisting of the variable and the time resolution, separated by an underscore.
#'
#' @param source.data.present A vector of strings naming the data that's present; at least one of (tmin, tmax, prec, tavg).
#' @param time.resolution The time resolutions to compute indices at. See \code{\link{create.indices.from.files}}.
#' @param climdex.vars.subset A character vector of lower-case names of Climdex indices to calculate (eg: tr, fd, rx5day). See \code{\link{create.indices.from.files}}.
#' @return A character vector containing variable names with time resolutions appended.
#'
#' @seealso \code{\link{create.indices.from.files}}
#' @examples
#' ## Get all variables which require tmin and/or tmax, for all time resolutions.
#' var.list1 <- get.climdex.variable.list(c("tmax", "tmin"))
#'
#' ## Get all variables which require prec with an annual time resolution.
#' var.list2 <- get.climdex.variable.list("prec", time.resolution="annual")
#'
#' ## Get the intersection of a set list of vars and available data.
#' sub.vars <- c("su", "id", "tr", "fd", "gsl", "csdi", "wsdi", "r10mm")
#' var.list3 <- get.climdex.variable.list("tmax", climdex.vars.subset=sub.vars)
#'
#' @export
get.climdex.variable.list <- function(source.data.present, time.resolution=c("all", "annual", "monthly"), climdex.vars.subset=NULL) {
  time.res <- match.arg(time.resolution)
  annual.only <- c("fdETCCDI", "suETCCDI", "idETCCDI", "trETCCDI", "gslETCCDI", "wsdiETCCDI", "csdiETCCDI", "sdiiETCCDI", "r10mmETCCDI", "r20mmETCCDI", "r1mmETCCDI", "cddETCCDI", "cwdETCCDI", "r95pETCCDI", "r99pETCCDI", "prcptotETCCDI", "altcddETCCDI", "altcwdETCCDI", "altcsdiETCCDI", "altwsdiETCCDI")
  vars.by.src.data.reqd <- list(tmax=c("suETCCDI", "idETCCDI", "txxETCCDI", "txnETCCDI", "tx10pETCCDI", "tx90pETCCDI", "wsdiETCCDI", "altwsdiETCCDI"),
                                tmin=c("fdETCCDI", "trETCCDI", "tnxETCCDI", "tnnETCCDI", "tn10pETCCDI", "tn90pETCCDI", "csdiETCCDI", "altcsdiETCCDI"),
                                prec=c("rx1dayETCCDI", "rx5dayETCCDI", "sdiiETCCDI", "r10mmETCCDI", "r20mmETCCDI", "r1mmETCCDI", "cddETCCDI", "cwdETCCDI", "r95pETCCDI", "r99pETCCDI", "prcptotETCCDI", "altcddETCCDI", "altcwdETCCDI"),
                                tavg=c("gslETCCDI", "dtrETCCDI") )

  if(any(!(source.data.present %in% c("tmin", "tmax", "tavg", "prec"))))
    stop("Invalid variable listed in source.data.present.")

  if(all(c("tmax", "tmin") %in% source.data.present) && !("tavg" %in% source.data.present))
    source.data.present <- c(source.data.present, "tavg")

  climdex.vars <- unlist(vars.by.src.data.reqd[source.data.present])
  if(!is.null(climdex.vars.subset))
    climdex.vars <- climdex.vars[climdex.vars %in% paste(climdex.vars.subset, "ETCCDI", sep="")]

  freq.lists <- list(c("mon", "yr"), c("yr"))
  dat <- switch(time.res,
                all=unlist(lapply(climdex.vars, function(x) { paste(x, freq.lists[[(x %in% annual.only) + 1]], sep="_") })),
                annual=paste(climdex.vars, "yr", sep="_"),
                monthly=paste(climdex.vars[!(climdex.vars %in% annual.only)], "mon", sep="_"))

  names(dat) <- NULL

  return(dat)
}

#' Returns a list of Climdex functions, with parameters curried in.
#'
#' Returns a list of Climdex functions, with parameters curried in.
#'
#' This function takes a variable list (as created by \code{\link{get.climdex.variable.list}}) and creates a list of functions corresponding to the specified indices, with parameters such as time resolution curried in. This allows for these functions to be called with just the \code{climdexInput} object as an argument, easing the automation of computing indices.
#'
#' @param vars.list The variable list, as created by \code{\link{get.climdex.variable.list}}.
#' @param fclimdex.compatible Whether to create fclimdex compatible functions.
#' @return A list of functions, named by the variable they compute.
#'
#' @examples
#' ## Get Climdex functions for a variable list with all appropriate params
#' ## curried in, so that all they take is a ClimdexInput object.
#' cdx.funcs <- get.climdex.functions(get.climdex.variable.list(c("tmax", "tmin")))
#'
#' @export
get.climdex.functions <- function(vars.list, fclimdex.compatible=TRUE) {
  func.names <- c("climdex.fd", "climdex.su", "climdex.id", "climdex.tr", "climdex.gsl",
                 "climdex.txx", "climdex.tnx", "climdex.txn", "climdex.tnn", "climdex.tn10p", "climdex.tx10p", "climdex.tn90p", "climdex.tx90p",
                 "climdex.txx", "climdex.tnx", "climdex.txn", "climdex.tnn", "climdex.tn10p", "climdex.tx10p", "climdex.tn90p", "climdex.tx90p",
                 "climdex.wsdi", "climdex.csdi",
                 "climdex.dtr", "climdex.rx1day", "climdex.rx5day",
                 "climdex.dtr", "climdex.rx1day", "climdex.rx5day",
                 "climdex.sdii", "climdex.r10mm", "climdex.r20mm", "climdex.rnnmm", "climdex.cdd", "climdex.cwd", "climdex.r95ptot", "climdex.r99ptot", "climdex.prcptot",
                 "climdex.cdd", "climdex.cwd", "climdex.csdi", "climdex.wsdi")

  el <- list()
  af <- list(freq="annual")
  mf <- list(freq="monthly")
  cwdd.opts <- list(spells.can.span.years=TRUE)
  altcwdd.opts <- list(spells.can.span.years=FALSE)
  wcsdi.opts <- list(spells.can.span.years=FALSE)
  altwcsdi.opts <- list(spells.can.span.years=TRUE)
  rx5day.opts <- list(center.mean.on.last.day=fclimdex.compatible)
  r1mm.opts <- list(threshold=1)
  options <- list(el, el, el, el, el,
                  mf, mf, mf, mf, mf, mf, mf, mf,
                  af, af, af, af, af, af, af, af,
                  wcsdi.opts, wcsdi.opts,
                  mf, mf, c(mf, rx5day.opts),
                  af, af, c(af, rx5day.opts),
                  el, el, el, r1mm.opts, cwdd.opts, cwdd.opts, el, el, el,
                  altcwdd.opts, altcwdd.opts, altwcsdi.opts, altwcsdi.opts)

  func <- lapply(1:length(func.names), function(n) do.call(functional::Curry, c(list(getFromNamespace(func.names[n], 'climdex.pcic')), options[[n]])))
  names(func) <- c("fdETCCDI_yr", "suETCCDI_yr", "idETCCDI_yr", "trETCCDI_yr", "gslETCCDI_yr",
                   "txxETCCDI_mon", "tnxETCCDI_mon", "txnETCCDI_mon", "tnnETCCDI_mon", "tn10pETCCDI_mon", "tx10pETCCDI_mon", "tn90pETCCDI_mon", "tx90pETCCDI_mon",
                   "txxETCCDI_yr", "tnxETCCDI_yr", "txnETCCDI_yr", "tnnETCCDI_yr", "tn10pETCCDI_yr", "tx10pETCCDI_yr", "tn90pETCCDI_yr", "tx90pETCCDI_yr",
                   "wsdiETCCDI_yr", "csdiETCCDI_yr",
                   "dtrETCCDI_mon", "rx1dayETCCDI_mon", "rx5dayETCCDI_mon",
                   "dtrETCCDI_yr", "rx1dayETCCDI_yr", "rx5dayETCCDI_yr",
                   "sdiiETCCDI_yr", "r10mmETCCDI_yr", "r20mmETCCDI_yr", "r1mmETCCDI_yr", "cddETCCDI_yr", "cwdETCCDI_yr", "r95pETCCDI_yr", "r99pETCCDI_yr", "prcptotETCCDI_yr",
                   "altcddETCCDI_yr", "altcwdETCCDI_yr", "altcsdiETCCDI_yr", "altwsdiETCCDI_yr")

  return(func[vars.list])
}

#' Returns metadata for specified Climdex variables
#'
#' Returns metadata for specified Climdex variables.
#'
#' This function returns metadata suitable for use in NetCDF files for the specified variables.
#'
#' @param vars.list The list of variables, as returned by \code{\link{get.climdex.variable.list}}.
#' @param template.filename The filename template to be used when generating filenames.
#' @return A data frame containing the following:
#' \itemize{
#' \item{long.name}{Long names for the variable}
#' \item{var.name}{Variable name for use in the file}
#' \item{units}{Units for the variable}
#' \item{annual}{Whether the variable is annual}
#' \item{base.period.attr}{Whether to include a base period attribute}
#' \item{standard.name}{Standard name to use for the variable}
#' \item{filename}{Filename to be written out}
#' }
#'
#' @examples
#' ## Get metadata (including filenames) for specified variables.
#' fn <- "pr_day_BCCAQ+ANUSPLIN300+MRI-CGCM3_historical+rcp85_r1i1p1_19500101-21001231.nc"
#' var.list2 <- get.climdex.variable.list("prec", time.resolution="annual")
#' md <- get.climdex.variable.metadata(var.list2, fn)
#'
#' @export
get.climdex.variable.metadata <- function(vars.list, template.filename) {
  all.data <- data.frame(long.name=c("Number of Frost Days", "Number of Summer Days", "Number of Icing Days", "Number of Tropical Nights", "Growing Season Length",
                          "Monthly Maximum of Daily Maximum Temperature", "Monthly Maximum of Daily Minimum Temperature",
                          "Monthly Minimum of Daily Maximum Temperature", "Monthly Minimum of Daily Minimum Temperature",
                          "Percentage of Days when Daily Minimum Temperature is Below the 10th Percentile", "Percentage of Days when Daily Maximum Temperature is Below the 10th Percentile",
                          "Percentage of Days when Daily Minimum Temperature is Above the 90th Percentile", "Percentage of Days when Daily Maximum Temperature is Above the 90th Percentile",
                          "Annual Maximum of Daily Maximum Temperature", "Annual Maximum of Daily Minimum Temperature",
                          "Annual Minimum of Daily Maximum Temperature", "Annual Minimum of Daily Minimum Temperature",
                          "Percentage of Days when Daily Minimum Temperature is Below the 10th Percentile", "Percentage of Days when Daily Maximum Temperature is Below the 10th Percentile",
                          "Percentage of Days when Daily Minimum Temperature is Above the 90th Percentile", "Percentage of Days when Daily Maximum Temperature is Above the 90th Percentile",
                          "Warm Spell Duration Index", "Cold Spell Duration Index",
                          "Mean Diurnal Temperature Range", "Monthly Maximum 1-day Precipitation", "Monthly Maximum Consecutive 5-day Precipitation",
                          "Mean Diurnal Temperature Range", "Annual Maximum 1-day Precipitation", "Annual Maximum Consecutive 5-day Precipitation",
                          "Simple Precipitation Intensity Index", "Annual Count of Days with At Least 10mm of Precipitation",
                          "Annual Count of Days with At Least 20mm of Precipitation", "Annual Count of Days with At Least 1mm of Precipitation",
                          "Maximum Number of Consecutive Days with Less Than 1mm of Precipitation", "Maximum Number of Consecutive Days with At Least 1mm of Precipitation",
                          "Annual Total Precipitation when Daily Precipitation Exceeds the 95th Percentile of Wet Day Precipitation",
                          "Annual Total Precipitation when Daily Precipitation Exceeds the 99th Percentile of Wet Day Precipitation", "Annual Total Precipitation in Wet Days",
                          "Maximum Number of Consecutive Days Per Year with Less Than 1mm of Precipitation", "Maximum Number of Consecutive Days Per Year with At Least 1mm of Precipitation",
                          "Cold Spell Duration Index Spanning Years", "Warm Spell Duration Index Spanning Years"),
                        var.name=c("fdETCCDI", "suETCCDI", "idETCCDI", "trETCCDI", "gslETCCDI",
                          "txxETCCDI", "tnxETCCDI", "txnETCCDI", "tnnETCCDI", "tn10pETCCDI", "tx10pETCCDI", "tn90pETCCDI", "tx90pETCCDI",
                          "txxETCCDI", "tnxETCCDI", "txnETCCDI", "tnnETCCDI", "tn10pETCCDI", "tx10pETCCDI", "tn90pETCCDI", "tx90pETCCDI",
                          "wsdiETCCDI", "csdiETCCDI",
                          "dtrETCCDI", "rx1dayETCCDI", "rx5dayETCCDI",
                          "dtrETCCDI", "rx1dayETCCDI", "rx5dayETCCDI",
                          "sdiiETCCDI", "r10mmETCCDI", "r20mmETCCDI", "r1mmETCCDI", "cddETCCDI", "cwdETCCDI", "r95pETCCDI", "r99pETCCDI", "prcptotETCCDI",
                          "altcddETCCDI", "altcwdETCCDI", "altcsdiETCCDI", "altwsdiETCCDI"),
                        units=c("days", "days", "days", "days", "days",
                          "degrees_C", "degrees_C", "degrees_C", "degrees_C", "%", "%", "%", "%",
                          "degrees_C", "degrees_C", "degrees_C", "degrees_C", "%", "%", "%", "%",
                          "days", "days",
                          "degrees_C", "mm", "mm",
                          "degrees_C", "mm", "mm",
                          "mm d-1", "days", "days", "days", "days", "days", "mm", "mm", "mm",
                          "days", "days", "days", "days"),
                        annual=c(T, T, T, T, T,
                          F, F, F, F, F, F, F, F,
                          T, T, T, T, T, T, T, T,
                          T, T,
                          F, F, F,
                          T, T, T,
                          T, T, T, T, T, T, T, T, T,
                          T, T, T, T),
                        base.period.attr=c(F, F, F, F, F,
                          F, F, F, F, T, T, T, T,
                          F, F, F, F, T, T, T, T,
                          T, T,
                          F, F, F,
                          F, F, F,
                          F, F, F, F, F, F, T, T, F,
                          F, F, T, T),
                         row.names=c("fdETCCDI_yr", "suETCCDI_yr", "idETCCDI_yr", "trETCCDI_yr", "gslETCCDI_yr",
                           "txxETCCDI_mon", "tnxETCCDI_mon", "txnETCCDI_mon", "tnnETCCDI_mon", "tn10pETCCDI_mon", "tx10pETCCDI_mon", "tn90pETCCDI_mon", "tx90pETCCDI_mon",
                           "txxETCCDI_yr", "tnxETCCDI_yr", "txnETCCDI_yr", "tnnETCCDI_yr", "tn10pETCCDI_yr", "tx10pETCCDI_yr", "tn90pETCCDI_yr", "tx90pETCCDI_yr",
                           "wsdiETCCDI_yr", "csdiETCCDI_yr",
                           "dtrETCCDI_mon", "rx1dayETCCDI_mon", "rx5dayETCCDI_mon",
                           "dtrETCCDI_yr", "rx1dayETCCDI_yr", "rx5dayETCCDI_yr",
                           "sdiiETCCDI_yr", "r10mmETCCDI_yr", "r20mmETCCDI_yr", "r1mmETCCDI_yr", "cddETCCDI_yr", "cwdETCCDI_yr", "r95pETCCDI_yr", "r99pETCCDI_yr", "prcptotETCCDI_yr",
                           "altcddETCCDI_yr", "altcwdETCCDI_yr", "altcsdiETCCDI_yr", "altwsdiETCCDI_yr"),
                         stringsAsFactors=FALSE)

  standard.name.lookup <- c(fdETCCDI="number_frost_days", suETCCDI="number_summer_days", idETCCDI="number_icing_days", trETCCDI="number_tropical_nights", gslETCCDI="growing_season_length",
                            txxETCCDI="maximum_daily_maximum_temperature", tnxETCCDI="maximum_daily_minimum_temperature", txnETCCDI="minimum_daily_maximum_temperature", tnnETCCDI="minimum_daily_minimum_temperature",
                            tn10pETCCDI="percent_days_when_daily_minimum_temperature_below_10p", tx10pETCCDI="percent_days_when_daily_maximum_temperature_below_10p",
                            tn90pETCCDI="percent_days_when_daily_minimum_temperature_above_90p", tx90pETCCDI="percent_days_when_daily_maximum_temperature_above_90p",
                            wsdiETCCDI="warm_spell_duration_index", csdiETCCDI="cold_spell_duration_index", dtrETCCDI="diurnal_temperature_range",
                            altwsdiETCCDI="warm_spell_duration_index", altcsdiETCCDI="cold_spell_duration_index",
                            rx1dayETCCDI="maximum_1day_precipitation", rx5dayETCCDI="maximum_5day_precipitation", sdiiETCCDI="simple_precipitation_intensity_index",
                            r10mmETCCDI="count_days_more_than_10mm_precipitation", r20mmETCCDI="count_days_more_than_20mm_precipitation", r1mmETCCDI="count_days_more_than_1mm_precipitation",
                            cddETCCDI="maximum_number_consecutive_dry_days", cwdETCCDI="maximum_number_consecutive_wet_days",
                            altcddETCCDI="maximum_number_consecutive_dry_days", altcwdETCCDI="maximum_number_consecutive_wet_days",
                            r95pETCCDI="total_precipitation_exceeding_95th_percentile", r99pETCCDI="total_precipitation_exceeding_99th_percentile", prcptotETCCDI="total_wet_day_precipitation")

  all.data$standard.name <- standard.name.lookup[all.data$var.name]

  all.data$filename <- create.climdex.cmip5.filenames(ncdf4.helpers::get.split.filename.cmip5(template.filename), rownames(all.data))
  return(all.data[vars.list,])
}

get.output.time.data <- function(ts, time.origin.PCICt, time.units, time.dim.name, time.bnds.name, bnds.dim, res=c("year", "month"), origin="1970-01-01") {
  res <- match.arg(res)
  time.bounds <- ncdf4.helpers::nc.make.time.bounds(ts, res)
  time.series <- PCICt::as.PCICt.numeric((unclass(time.bounds[1,]) + unclass(time.bounds[2,])) / 2, cal=attr(time.bounds, "cal"), origin=origin)
  time.bounds.days <- as.numeric(julian(time.bounds, origin=time.origin.PCICt))
  time.days <- as.numeric(julian(time.series, origin=time.origin.PCICt))
  time.dim <- ncdf4::ncdim_def(time.dim.name, units=time.units, vals=time.days, unlim=TRUE, longname='')
  time.bnds.var <- ncdf4::ncvar_def(time.bnds.name, '', list(bnds.dim, time.dim), longname='', prec="double")
  return(list(time.dim=time.dim, time.bnds.var=time.bnds.var, time.bnds.data=time.bounds.days))
}

#' Creates output files for Climdex variables.
#'
#' Creates output files for Climdex variables.
#'
#' This function creates a set of output files for the set of variable parameters passed in \code{cdx.dat}, as created by \code{\link{get.climdex.variable.metadata}}. It copies metadata from input files as appropriate and adds new metadata as required.
#'
#' @param cdx.dat The variable description data, as created by \code{\link{get.climdex.variable.metadata}}.
#' @param f The file(s) being used as input.
#' @param v.f.idx A mapping from variables to files, as created by \code{\link{get.var.file.idx}}.
#' @param variable.name.map A mapping from standardized names (tmax, tmin, prec) to NetCDF variable names.
#' @param ts The associated time data, as created by \code{nc.get.time.series}.
#' @param time.origin The time origin, as specified in the source NetCDF file(s).
#' @param base.range The base range; a vector of two numeric years.
#' @param out.dir The output directory name.
#' @param author.data A vector containing named elements describing the author; see \code{\link{create.indices.from.files}}.
#' @return A list of objects of type \code{ncdf4}.
#'
#' @examples
#' \donttest{
#' ## Establish basic inputs.
#' author.data <- list(institution="Looney Bin", institution_id="LBC")
#' input.files <- c("pr_NAM44_CanRCM4_ERAINT_r1i1p1_1989-2009.nc")
#'
#' ## Prepare derived inputs.
#' f <- lapply(input.files, ncdf4::nc_open)
#' variable.name.map <- c(tmax="tasmax", tmin="tasmin", prec="pr")
#' f.meta <- create.file.metadata(f, variable.name.map)
#' climdex.var.list <- get.climdex.variable.list(names(f.meta$v.f.idx), "all", NULL)
#' cdx.meta <- get.climdex.variable.metadata(climdex.var.list, input.files[1])
#'
#' ## Create output files
#' cdx.ncfile <- create.ncdf.output.files(cdx.meta, f, f.meta$v.f.idx, variable.name.map,
#'                                        f.meta$ts, get.time.origin(f, f.meta$dim.axes),
#'                                        c(1981,1990), "/foo", author.data)
#' }
#'
#' @export
create.ncdf.output.files <- function(cdx.dat, f, v.f.idx, variable.name.map, ts, time.origin, base.range, out.dir, author.data) {
  f.example <- f[[v.f.idx[1]]]
  v.example <- variable.name.map[names(v.f.idx)[1]]
  time.dim.name <- ncdf4.helpers::nc.get.dim.for.axis(f.example, v.example, "T")$name
  old.time.bnds.att <- ncdf4::ncatt_get(f.example, time.dim.name, "bounds")
  time.bnds.name <- if(old.time.bnds.att$hasatt) old.time.bnds.att$value else paste(time.dim.name, "bnds", sep="_")

  ## Create new time dimensions
  time.origin.PCICt <- PCICt::as.PCICt.default(time.origin, cal=attr(ts, "cal"))
  time.units <- paste("days since", time.origin)

  input.bounds <- ncdf4.helpers::nc.get.dim.bounds.var.list(f.example, v.example)
  ## FIXME: I'm not sure how solid the assumption about the location of bnds here is.
  bnds <- if(length(input.bounds) > 0) f.example$var[[input.bounds[1]]]$dim[[1]] else ncdf4::ncdim_def("bnds", "", 1:2, create_dimvar=FALSE)
  time.dat <- list(annual=get.output.time.data(ts, time.origin.PCICt, time.units, time.dim.name, time.bnds.name, bnds, res="year"),
                   monthly=get.output.time.data(ts, time.origin.PCICt, time.units, time.dim.name, time.bnds.name, bnds, res="month"))

  grid.mapping.att <- ncdf4::ncatt_get(f.example, v.example, "grid_mapping")
  vars.to.copy <- c(input.bounds[input.bounds != time.bnds.name], names(ncdf4.helpers::nc.get.coordinate.axes(f.example, v.example)), if(grid.mapping.att$hasatt) grid.mapping.att$value)
  vars.to.clone.atts.for <- c(vars.to.copy, ncdf4.helpers::nc.get.dim.names(f.example, v.example))
  vars.ncvars <- sapply(vars.to.copy, function(x) { f.example$var[[x]] }, simplify=FALSE)
  vars.data <- lapply(vars.ncvars, function(ncvar) { if(length(ncvar$dim) == 0) NULL else ncdf4::ncvar_get(f.example, ncvar) })

  return(lapply(1:length(cdx.dat$var.name), function(x) {
    annual <- cdx.dat$annual[x]
    time.for.file <- time.dat[[c("monthly", "annual")[1 + annual]]]

    ## Establish variables, create file
    nc.var.list <- c(vars.ncvars, list(time.for.file$time.bnds.var, ncdf4::ncvar_def(name=cdx.dat$var.name[x], units=cdx.dat$units[x], dim=c(f.example$var[[v.example]]$dim[1:2], list(time.for.file$time.dim)), missval=1e20, longname=cdx.dat$long.name[x])))
    new.file <- ncdf4::nc_create(paste(out.dir, cdx.dat$filename[x], sep="/"), nc.var.list, force_v4=TRUE)

    ## Copy attributes for all variables plus global attributes
    att.rename <- c("frequency"="input_frequency", "creation_date"="input_creation_date", "title"="input_title", "tracking_id"="input_tracking_id")
    inst.id <- ncdf4::ncatt_get(f.example, 0, "institution_id")
    if(inst.id$hasatt) {
      att.rename.inst <- c("contact"="contact", "references"="references")
      names(att.rename.inst) <- paste(inst.id$value, names(att.rename.inst), sep="_")
      att.rename <- c(att.rename, att.rename.inst)
    }

    ## Copy attributes with renaming and exclusions.
    ncdf4.helpers::nc.copy.atts(f.example, 0, new.file, 0, definemode=TRUE, rename.mapping=att.rename)
    ncdf4.helpers::nc.copy.atts(f.example, v.example, new.file, cdx.dat$var.name[x], definemode=TRUE, exception.list=c("units", "long_name", "standard_name", "base_period", "missing_value", "_FillValue", "add_", "valid_min", "valid_max", "valid_range", "scale_factor", "add_offset", "signedness", "history"))
    for(v in vars.to.clone.atts.for) {
      ncdf4.helpers::nc.copy.atts(f.example, v, new.file, v, definemode=TRUE)
    }
    ncdf4::ncatt_put(new.file, time.dim.name, "units", time.units, definemode=TRUE)

    ## Put additional attributes.
    put.history.att(new.file, cdx.dat$var.name[x], definemode=TRUE)
    put.ETCCDI.atts(new.file, c("mon", "yr")[1 + annual], ncdf4::ncatt_get(f.example, 0, "title")$value, author.data, definemode=TRUE)
    if(cdx.dat$base.period.attr[x])
      ncdf4::ncatt_put(new.file, cdx.dat$var.name[x], "base_period", paste(base.range, collapse="-"), definemode=TRUE)
    ncdf4::nc_enddef(new.file)

    ## Copy data from vars.to.copy and put time bounds.
    ncdf4::ncvar_put(new.file, time.bnds.name, time.for.file$time.bnds.data)
    for(v in vars.to.copy)
      if(!is.null(vars.data[[v]]))
         ncdf4::ncvar_put(new.file, v, vars.data[[v]])

    new.file
  }))
}

## Get dim sizes, with checking to make sure sizes are all the same.
get.dim.size <- function(f, v.f.idx, variable.name.map) {
  dim.size.list <- lapply(names(v.f.idx), function(i) { f[[v.f.idx[i]]]$var[[variable.name.map[i]]]$varsize })
  stopifnot(all.the.same(dim.size.list))
  dim.size.list[[1]]
}

## Get dim axes, with checking to make sure they all have same axes.
get.dim.axes <- function(f, v.f.idx, variable.name.map) {
  dim.axes.list <- lapply(names(v.f.idx), function(i) { ncdf4.helpers::nc.get.dim.axes(f[[v.f.idx[i]]], variable.name.map[i]) })
  stopifnot(all.the.same(dim.axes.list))
  dim.axes.list[[1]]
}

## Get timeseries (as PCICt), with error checking to ensure input files have same TS.
## FIXME: This will need to be revised for fixed time dimensions. Should be identified by axis.
get.ts <- function(f) {
  ts.list <- lapply(lapply(f, ncdf4.helpers::nc.get.time.series), trunc, "days")
  stopifnot(all.the.same(ts.list))
  ts.list[[1]]
}

## Compute all indices for a single grid box
#' Compute Climdex indices using provided data.
#'
#' Compute Climdex indices using provided data.
#'
#' Given the provided data and functions, compute the Climdex indices defined by the functions.
#'
#' @param in.dat The input data to compute indices on.
#' @param cdx.funcs The functions to be applied to the data, as created by \code{\link{get.climdex.functions}}.
#' @param ts The associated time data, as created by \code{nc.get.time.series}.
#' @param base.range The base range; a vector of two numeric years.
#' @param fclimdex.compatible Whether to make the results identical to those of fclimdex; this affects how the data in the base period is padded.
#' @return A list of data for each index.
#'
#' @examples
#' library(climdex.pcic)
#'
#' ## Prepare input data
#' in.dat <- list(tmax=ec.1018935.tmax$MAX_TEMP)
#' cdx.funcs <- get.climdex.functions(get.climdex.variable.list(names(in.dat)))
#' in.dat$northern.hemisphere <- TRUE
#' ts <- as.PCICt(do.call(paste, ec.1018935.tmax[,c("year", "jday")]),
#'                format="%Y %j", cal="gregorian")
#'
#' ## Compute indices
#' res <- compute.climdex.indices(in.dat, cdx.funcs, ts, c(1981, 1990), FALSE)
#'
#' @export
compute.climdex.indices <- function(in.dat, cdx.funcs, ts, base.range, fclimdex.compatible) {
  ci <- climdex.pcic::climdexInput.raw(
                        in.dat$tmax, in.dat$tmin, in.dat$prec,
                        if(is.null(in.dat$tmax)) NULL else ts,
                        if(is.null(in.dat$tmin)) NULL else ts,
                        if(is.null(in.dat$prec)) NULL else ts,
                        tavg=in.dat$tavg, tavg.dates=if(is.null(in.dat$tavg)) NULL else ts,
                        base.range=base.range, northern.hemisphere=in.dat$northern.hemisphere,
                        quantiles=in.dat$quantiles)

  ## NOTE: Names must be stripped here because it increases memory usage on the head by a factor of 8-9x (!)
  return(lapply(cdx.funcs, function(f) { d <- f(ci=ci); names(d) <- NULL; d }))
}

#' Flatten the X and Y dimensions down to a space dimension.
#'
#' Flatten the X and Y dimensions down to a space dimension.
#'
#' This function takes input data, a vector of dimensions to reduce to 1 dimension, and optionally a subset of dimnames to copy. It returns the data with the specified dimensions shrunk down to 1 dimension.
#'
#' @param dat The data to operate on.
#' @param reduce.dims The names or indices of the dimensions to reduce to 1 dimension.
#' @param names.subset Optionally, a subset of dimension names to copy.
#' @return The data with the specified dimensions reduced to 1 dimension.
#'
#' @note The dimensions to reduce must be adjoining dimensions.
#'
#' @examples
#' ## Take example data and flatten the last two dims down to one.
#' dat <- structure(1:8, .Dim=c(2, 2, 2))
#' dat.flat <- flatten.dims(dat, 2:3)
#'
#' @export
flatten.dims <- function(dat, reduce.dims, names.subset) {
  stopifnot(all(diff(reduce.dims) == 1))
  dat.dim <- dim(dat)
  if(!missing(names.subset))
    dat.dimnames <- dimnames(dat)
  before.reduce <- 1:length(dat.dim) < min(reduce.dims)
  after.reduce <- 1:length(dat.dim) > max(reduce.dims)
  new.dims <- c(dat.dim[before.reduce], prod(dat.dim[reduce.dims]), dat.dim[after.reduce])
  dim(dat) <- new.dims
  if(!missing(names.subset))
    dimnames(dat) <- dat.dimnames[names.subset]
  return(dat)
}

## FIXME: Handle time-minor data gracefully.
#' Retrieve and convert data to correct units and dimensions.
#'
#' Retrieve and convert data to correct units and dimensions.
#'
#' This function retrieves NetCDF data for the specified subset from the specified file and variable; converts from \code{src.units} to \code{dest.units}, transposes the data to (T, S) dimensionality, and returns the result.
#'
#' @param f The NetCDF file to retrieve data from; an object of class \code{ncdf4}.
#' @param v The variable to retrieve data from.
#' @param subset The subset to retrieve.
#' @param src.units The source units to convert data from.
#' @param dest.units The destination units to convert to.
#' @param dim.axes The dimension axes to be used.
#' @return The retrieved and massaged data.
#'
#' @examples
#' \donttest{get.data(f, "pr", list(Y=3), "kg m-2 s-1", "kg m-2 s-1", c(X="lon",Y="lat",T="time"))}
#'
#' @export
get.data <- function(f, v, subset, src.units, dest.units, dim.axes) {
  stopifnot(inherits(f, "ncdf4"))
  dat <- if(!missing(src.units) && !missing(dest.units))
    udunits2::ud.convert(ncdf4.helpers::nc.get.var.subset.by.axes(f, v, subset), src.units, dest.units)
  else
    ncdf4.helpers::nc.get.var.subset.by.axes(f, v, subset)

  reduce.dims <- which(dim.axes %in% c("X", "Y", "Z"))
  return(t(flatten.dims(dat, reduce.dims=reduce.dims)))
}

## Produce slab of northern.hemisphere booleans of the same shape as the data.
#' Determine what portions of a subset are within the northern hemisphere.
#'
#' Determine what portions of a subset are within the northern hemisphere.
#'
#' Given a subset, a file, a variable, and a projection, determine what positions are within the northern hemisphere, returning the result as an array of booleans.
#'
#' @param subset The subset to use.
#' @param f The NetCDF file to use; an object of class \code{ncdf4}.
#' @param v The variable in question.
#' @param projection The proj4 string to use; NULL if the data is not in a projected coordinate space.
#' @return An array of booleans corresponding to the subset containing TRUE if the point is within the northern hemisphere, and FALSE otherwise.
#'
#' @examples
#' \donttest{
#' ## Open files, etc.
#' input.files <- c("tasmax_NAM44_CanRCM4_ERAINT_r1i1p1_1989-2009.nc")
#' f <- list(nc_open(input.files))
#' f.v <- lapply(f, ncdf4.helpers::nc.get.variable.list, min.dims=2)
#' bools <- get.northern.hemisphere.booleans(list(X=1:2, Y=1:2), f[[1]], f.v[[1]], NULL)
#' }
#'
#' @export
get.northern.hemisphere.booleans <- function(subset, f, v, projection) {
  y.dim <- ncdf4.helpers::nc.get.dim.for.axis(f, v, "Y")
  x.dim <- ncdf4.helpers::nc.get.dim.for.axis(f, v, "X")
  y.subset.vals <- rep(y.dim$vals[if(is.null(subset$Y)) 1:y.dim$len else subset$Y],
                       each=(if(is.null(subset$X)) x.dim$len else length(subset$X)))
  if(!is.null(projection)) {
    x.subset.vals <- rep(x.dim$vals[if(is.null(subset$X)) 1:x.dim$len else subset$X],
                         (if(is.null(subset$Y)) y.dim$len else length(subset$Y)))
    dat <- proj4::project(list(x=x.subset.vals, y=y.subset.vals), projection, inverse=TRUE, ellps.default=NA)
    return(dat$y >= 0)
  } else
    return(y.subset.vals >= 0)
}

#' Extract a single quantiles object from a set of thresholds.
#'
#' Extract a single quantiles object from a set of thresholds.
#'
#' From a set of thresholds as retrieved from one or more NetCDF files containing thresholds, this function extracts a single point and converts the format to one suitable for passing to \code{climdexInput.raw}.
#'
#' @param thresholds The thresholds, as extracted by \code{\link{get.thresholds.chunk}}.
#' @param idx The index to extract.
#' @return A quantiles object suitable for passing to \code{climdexInput.raw} as the \code{quantiles} argument.
#'
#' @examples
#' \donttest{
#' ## Define mappings and filenames.
#' thresholds.name.map <- c(tx10thresh="tx10thresh", tn10thresh="tn10thresh", tx90thresh="tx90thresh",
#'                          tn90thresh="tn90thresh", r95thresh="r95thresh", r99thresh="r99thresh")
#' thresh.files <- "thresholds.nc"
#'
#' ## Open files, etc.
#' cdx.funcs <- get.climdex.functions(get.climdex.variable.list("tmax"))
#' thresholds.netcdf <- lapply(thresh.files, nc_open)
#' t.f.idx <- get.var.file.idx(thresholds.name.map, lapply(thresholds.netcdf,
#'                             ncdf4.helpers::nc.get.variable.list, min.dims=2))
#'
#' ## Get thresholds chunk.
#' dat <- get.thresholds.chunk(list(Y=1), cdx.funcs, thresholds.netcdf, t.f.idx, thresholds.name.map)
#'
#' ## Get quantiles object for index 2
#' q <- get.quantiles.object(dat, 2)
#' }
#'
#' @export
get.quantiles.object <- function(thresholds, idx) {
  if(is.null(thresholds))
    return(NULL)

  thresh.path.2d <- list(tx10thresh=c("tmax", "outbase", "q10"),
                         tx90thresh=c("tmax", "outbase", "q90"),
                         tn10thresh=c("tmin", "outbase", "q10"),
                         tn90thresh=c("tmin", "outbase", "q90"))
  thresh.path.1d <- list(r95thresh=c("prec", "q95"),
                         r99thresh=c("prec", "q99"))
  result <- list()


  recursive.append <- function(x, l, data) {
    if(length(x) == 0) return(data)
    if(is.null(l)) l <- list()
    return(c(l[!(names(l) %in% x[1])], structure(list(recursive.append(tail(x, n=-1), l[[x[1]]], data)), .Names=x[1])))
  }


  for(threshold.var in names(thresh.path.2d)[names(thresh.path.2d) %in% names(thresholds)])
    result <- recursive.append(thresh.path.2d[[threshold.var]], result, thresholds[[threshold.var]][,idx])

  for(threshold.var in names(thresh.path.1d)[names(thresh.path.1d) %in% names(thresholds)]) {
    thresh.path <- thresh.path.1d[[threshold.var]]
    result[[thresh.path[1]]] <- c(result[[thresh.path[1]]], structure(thresholds[[threshold.var]][idx], .Names=thresh.path[2]))
  }

  return(result)
}

#' Compute Climdex indices for a subset / stripe
#'
#' Compute Climdex indices for a subset / stripe
#'
#' Given a subset, a set of Climdex functions (as created by \code{\link{get.climdex.functions}}), and ancillary data, load and convert data, create a climdexInput object for each point, run all of the functions in \code{cdx.funcs} on that data, and return the result.
#'
#' @param subset The subset to use.
#' @param cdx.funcs The functions to be applied to the data, as created by \code{\link{get.climdex.functions}}.
#' @param ts The associated time data, as created by \code{nc.get.time.series}.
#' @param base.range The base range; a vector of two numeric years.
#' @param dim.axes The dimension axes for the input data.
#' @param v.f.idx A mapping from variables to files, as created by \code{\link{get.var.file.idx}}.
#' @param variable.name.map A mapping from standardized names (tmax, tmin, prec) to NetCDF variable names.
#' @param src.units The source units to convert data from.
#' @param dest.units The destination units to convert to.
#' @param t.f.idx A mapping from threshold variables to threshold files, as created by \code{\link{get.var.file.idx}}.
#' @param thresholds.name.map A mapping from standardized names (tx10thresh, tn90thresh, etc) to NetCDF variable names.
#' @param fclimdex.compatible Whether to make the results identical to those of fclimdex; this affects how the data in the base period is padded.
#' @param projection A proj4 string representing the projection the data is in.
#' @param f A list of objects of type \code{ncdf4}, consisting of the open input files. If missing, will be pulled from the global namespace.
#' @param thresholds.netcdf A list of objects of type \code{ncdf4}, consisting of the open threshold files. If missing, will be pulled from the global namespace.
#'
#' @note This function relies on an object named 'f' and containing the opened NetCDF files being part of the global namespace.
#'
#' @examples
#' \donttest{
#' ## Define mappings and filenames.
#' author.data <- list(institution="Looney Bin", institution_id="LBC")
#' input.files <- c("pr_NAM44_CanRCM4_ERAINT_r1i1p1_1989-2009.nc")
#' variable.name.map <- c(tmax="tasmax", tmin="tasmin", prec="pr")
#'
#' ## Open files, etc.
#' cdx.funcs <- get.climdex.functions(get.climdex.variable.list("tmax"))
#' f <- lapply(input.files, ncdf4::nc_open)
#' f.meta <- create.file.metadata(f, variable.name.map)
#' climdex.var.list <- get.climdex.variable.list(names(f.meta$v.f.idx), "all", NULL)
#' cdx.meta <- get.climdex.variable.metadata(climdex.var.list, input.files[1])
#'
#' ## Compute indices for stripe
#' cdx <- compute.indices.for.stripe(list(Y=1), cdx.funcs, f.meta$ts, c(1981, 1990), f.meta$dim.axes,
#'                            f.meta$v.f.idx, variable.name.map, f.meta$src.units, f.meta$dest.units,
#'                            t.f.idx, NULL, f=f, thresholds.netcdf=NULL)
#' }
#'
#' @export
compute.indices.for.stripe <- function(subset, cdx.funcs, ts, base.range, dim.axes, v.f.idx, variable.name.map, src.units, dest.units, t.f.idx, thresholds.name.map, fclimdex.compatible=TRUE, projection=NULL, f, thresholds.netcdf) {
  f <- if(missing(f)) get("f", .GlobalEnv) else f
  thresholds.netcdf <- if(missing(thresholds.netcdf)) get("thresholds.netcdf", .GlobalEnv) else thresholds.netcdf

  ## Dimension order: Time, Space for each Var in list
  data.list <- sapply(names(v.f.idx), function(x) { gc(); get.data(f[[v.f.idx[x]]], variable.name.map[x], subset, src.units[x], dest.units[x], dim.axes) }, simplify=FALSE)
  gc()

  northern.hemisphere <- get.northern.hemisphere.booleans(subset, f[[v.f.idx[1]]], variable.name.map[names(v.f.idx)[1]], projection)

  thresholds <- if(is.null(thresholds.netcdf)) NULL else get.thresholds.chunk(subset, cdx.funcs, thresholds.netcdf, t.f.idx, thresholds.name.map)
  return(lapply(1:(dim(data.list[[1]])[2]), function(x) {
    dat.list <- sapply(names(data.list), function(name) { data.list[[name]][,x] }, simplify=FALSE)
    ## Fast-path the all-NA case.
    if(all(sapply(dat.list, function(x) { all(is.na(x)) }))) {
      ## We don't need to pad this out to full length; cbind will do that for us.
      return(structure(as.list(rep(NA, length(cdx.funcs))), .Names=names(cdx.funcs)))
    } else {
      indices.input <- c(dat.list, northern.hemisphere=northern.hemisphere[x], list(quantiles=get.quantiles.object(thresholds, x)))
      return(compute.climdex.indices(indices.input, cdx.funcs, ts, base.range, fclimdex.compatible))
  }
  }))
}

#' Retrieve thresholds for a subset
#'
#' Retrieve thresholds for a subset
#'
#' Given a subset, a set of Climdex functions (as created by \code{\link{get.climdex.functions}}), and ancillary data, load the thresholds required for the functions being called and return them.
#'
#' @param subset The subset to use.
#' @param cdx.funcs The functions to be applied to the data, as created by \code{\link{get.climdex.functions}}.
#' @param thresholds.netcdf One or more NetCDF files containing thresholds.
#' @param t.f.idx A mapping from threshold variables to threshold files, as created by \code{\link{get.var.file.idx}}.
#' @param thresholds.name.map A mapping from standardized names (tx10thresh, tn90thresh, etc) to NetCDF variable names.
#'
#' @examples
#' \donttest{
#' ## Define mappings and filenames.
#' thresholds.name.map <- c(tx10thresh="tx10thresh", tn10thresh="tn10thresh", tx90thresh="tx90thresh",
#'                          tn90thresh="tn90thresh", r95thresh="r95thresh", r99thresh="r99thresh")
#' thresh.files <- "thresholds.nc"
#'
#' ## Open files, etc.
#' cdx.funcs <- get.climdex.functions(get.climdex.variable.list("tmax"))
#' thresholds.netcdf <- lapply(thresh.files, nc_open)
#' t.f.idx <- get.var.file.idx(thresholds.name.map, lapply(thresholds.netcdf,
#'                             ncdf4.helpers::nc.get.variable.list, min.dims=2))
#'
#' ## Get thresholds chunk.
#' dat <- get.thresholds.chunk(list(Y=1), cdx.funcs, thresholds.netcdf, t.f.idx, thresholds.name.map)
#' }
#'
#' @export
get.thresholds.chunk <- function(subset, cdx.funcs, thresholds.netcdf, t.f.idx, thresholds.name.map) {
  var.thresh.map <- list(tx10thresh=c("tx10p"), tx90thresh=c("tx90p", "WSDI"), tn10thresh=c("tn10p", "CSDI"), tn90thresh=c("tn90p"), r95thresh=c("r95p"), r99thresh=c("r99p"))

  cdx.names <- names(cdx.funcs)
  thresh.var.needed <- names(var.thresh.map)[sapply(var.thresh.map, function(x) { any(unlist(lapply(x, function(substr) { any(grepl(substr, cdx.names)) }))) })]
  stopifnot(all(thresh.var.needed %in% names(t.f.idx)))
  return(sapply(thresh.var.needed, function(threshold.var) {
    dim.axes <- ncdf4.helpers::nc.get.dim.axes(thresholds.netcdf[[t.f.idx[threshold.var]]], thresholds.name.map[threshold.var]);
    return(get.data(thresholds.netcdf[[t.f.idx[threshold.var]]], thresholds.name.map[threshold.var], subset, dim.axes=dim.axes))
  }, simplify=FALSE))
}

## Write out results for variables computed
#' Write out computed climdex results
#'
#' Write out computed climdex results
#'
#' Given a set of Climdex results, a subset, a set of files, and dimension sizes, write out the data to the appropriate files.
#'
#' @param climdex.results The results to write out.
#' @param chunk.subset The corresponding subset.
#' @param cdx.ncfile The list of NetCDF files to write the results out to.
#' @param dim.size The overall size of the input data.
#' @param cdx.varname The list of NetCDF variable names for the files in \code{cdx.ncfile}.
#'
#' @examples
#' \donttest{
#' ## Define mappings and filenames.
#' author.data <- list(institution="Looney Bin", institution_id="LBC")
#' input.files <- c("pr_NAM44_CanRCM4_ERAINT_r1i1p1_1989-2009.nc")
#' variable.name.map <- c(tmax="tasmax", tmin="tasmin", prec="pr")
#'
#' ## Open files, etc.
#' cdx.funcs <- get.climdex.functions("tmax")
#' f <- lapply(input.files, ncdf4::nc_open)
#' f.meta <- create.file.metadata(f, variable.name.map)
#' climdex.var.list <- get.climdex.variable.list(names(f.meta$v.f.idx), "all", NULL)
#' cdx.meta <- get.climdex.variable.metadata(climdex.var.list, input.files[1])
#'
#' ## Create output files
#' cdx.ncfile <- create.ncdf.output.files(cdx.meta, f, f.meta$v.f.idx, variable.name.map,
#'                                        f.meta$ts, get.time.origin(f, f.meta$dim.axes),
#'                                        c(1981,1990), "/foo", author.data)
#'
#' ## Compute indices for stripe
#' cdx <- compute.indices.for.stripe(list(Y=1), cdx.funcs, f.meta$ts, c(1991, 2000), f.meta$dim.axes,
#'                            f.meta$v.f.idx, variable.name.map, f.meta$src.units, f.meta$dest.units,
#'                            t.f.idx, NULL, f=f, thresholds.netcdf=NULL)
#'
#' ## Write out indices
#' write.climdex.results(cdx, list(Y=1), cdx.ncfile, f.meta$dim.size, cdx.meta$varname)
#' }
#'
#' @export
write.climdex.results <- function(climdex.results, chunk.subset, cdx.ncfile, dim.size, cdx.varname) {
  xy.dims <- dim.size[1:2]
  if(!is.null(chunk.subset$X))
    xy.dims[1] <- length(chunk.subset$X)
  if(!is.null(chunk.subset$Y))
    xy.dims[2] <- length(chunk.subset$Y)

  ## Write out results, variable by variable
  lapply(1:length(cdx.ncfile), function(v) {
    dat <- t(do.call(cbind, lapply(climdex.results, function(cr) { cr[[v]] })))
    t.dim.len <- ncdf4.helpers::nc.get.dim.for.axis(cdx.ncfile[[v]], cdx.varname[v], "T")$len

    ## If data is of length 1, it's an error.
    if(length(dat) == 1)
      stop(dat)

    ## Special case of an entire slab missing values... repeat such that we have full data.
    if(prod(dim(dat)) != prod(c(xy.dims, t.dim.len)))
      dat <- rep(dat, t.dim.len)

    dim(dat) <- c(xy.dims, t.dim.len)
    ncdf4.helpers::nc.put.var.subset.by.axes(cdx.ncfile[[v]], cdx.varname[v], dat, chunk.subset)
  })
  invisible(0)
}

#' Compute Climdex thresholds for a subset / stripe
#'
#' Compute Climdex thresholds for a subset / stripe
#'
#' Given a subset and ancillary data, load and convert data, get the out-of-base quantiles for the data for each point, and return the result.
#'
#' @param subset The subset to use.
#' @param ts The associated time data, as created by \code{nc.get.time.series}.
#' @param base.range The base range; a vector of two numeric years.
#' @param dim.axes The dimension axes for the input data.
#' @param v.f.idx A mapping from variables to files, as created by \code{\link{get.var.file.idx}}.
#' @param variable.name.map A mapping from standardized names (tmax, tmin, prec) to NetCDF variable names.
#' @param src.units The source units to convert data from.
#' @param dest.units The destination units to convert to.
#' @param f A list of objects of type \code{ncdf4}, consisting of the open input files. If missing, will be pulled from the global namespace.
#'
#' @note This function relies on an object named 'f' and containing the opened NetCDF files being part of the global namespace.
#'
#' @examples
#' \donttest{
#' ## Establish basic inputs.
#' author.data <- list(institution="Looney Bin", institution_id="LBC")
#' input.files <- c("pr_NAM44_CanRCM4_ERAINT_r1i1p1_1989-2009.nc")
#'
#' ## Prepare derived inputs.
#' f <- lapply(input.files, ncdf4::nc_open)
#' variable.name.map <- c(tmax="tasmax", tmin="tasmin", prec="pr")
#' f.meta <- create.file.metadata(f, variable.name.map)
#' threshold.dat <- get.thresholds.metadata(names(f.meta$v.f.idx))
#'
#' ## Create output file
#' thresh.file <- create.thresholds.file("thresh.nc", f, f.meta$ts, f.meta$v.f.idx, variable.name.map,
#'                                       c(1991, 2000), f.meta$dim.size, f.meta$dim.axes,
#'                                       threshold.dat, author.data)
#'
#' ## Compute threshold quantiles for stripe
#' q <- get.quantiles.for.stripe(list(Y=1), f.meta$ts, c(1991, 2000), f.meta$dim.axes,
#'                               f.meta$v.f.idx, variable.name.map, f.meta$src.units,
#'                               f.meta$dest.units, f)
#' }
#'
#' @export
get.quantiles.for.stripe <- function(subset, ts, base.range, dim.axes, v.f.idx, variable.name.map, src.units, dest.units, f) {
  f <- if(missing(f)) get("f", .GlobalEnv) else f
  data.list <- sapply(names(v.f.idx), function(x) { gc(); get.data(f[[v.f.idx[x]]], variable.name.map[x], subset, src.units[x], dest.units[x], dim.axes) }, simplify=FALSE)
  gc()

  r <- 1:(dim(data.list[[1]])[2])
  if(!is.null(data.list$tmax)) {
    if(!is.null(data.list$tmin)) {
      if(!is.null(data.list$prec)) {
        return(lapply(r, function(x) climdex.pcic::get.outofbase.quantiles(data.list$tmax[,x], data.list$tmin[,x], data.list$prec[,x], ts, ts, ts, base.range)))
      } else {
        return(lapply(r, function(x) climdex.pcic::get.outofbase.quantiles(data.list$tmax[,x], data.list$tmin[,x], NULL, ts, ts, NULL, base.range)))
      }
    } else {
      if(!is.null(data.list$prec)) {
        return(lapply(r, function(x) climdex.pcic::get.outofbase.quantiles(data.list$tmax[,x], NULL, data.list$prec[,x], ts, NULL, ts, base.range)))
      } else {
        return(lapply(r, function(x) climdex.pcic::get.outofbase.quantiles(data.list$tmax[,x], NULL, NULL, ts, NULL, NULL, base.range)))
      }
    }
  } else {
    if(!is.null(data.list$tmin)) {
      if(!is.null(data.list$prec)) {
        return(lapply(r, function(x) climdex.pcic::get.outofbase.quantiles(NULL, data.list$tmin[,x], data.list$prec[,x], NULL, ts, ts, base.range)))
      } else {
        return(lapply(r, function(x) climdex.pcic::get.outofbase.quantiles(NULL, data.list$tmin[,x], NULL, NULL, ts, NULL, base.range)))
      }
    } else {
      if(!is.null(data.list$prec)) {
        return(lapply(r, function(x) climdex.pcic::get.outofbase.quantiles(NULL, NULL, data.list$prec[,x], NULL, NULL, ts, base.range)))
      } else {
        stop("Go home and take your shitty input with you.")
      }
    }
  }
}

set.up.cluster <- function(parallel, type="SOCK", src) {
  ## Fire up the cluster...
  cluster <- NULL

  if(!is.logical(parallel)) {
    cat(paste("Creating cluster of", parallel, "nodes of type", type, "\n"))
    cat(paste("SRC:", src))
    cluster <- snow::makeCluster(parallel, type)
    snow::clusterCall(cluster, function() { source(src) })
    ##snow::clusterEvalQ(cluster, library(climdex.pcic.ncdf))
    ##snow::clusterEvalQ(cluster, try(getFromNamespace('nc_set_chunk_cache', 'ncdf4')(1024 * 2048, 1009), silent=TRUE))
  }
  cluster
}

#' Creates Climdex thresholds output file.
#'
#' Creates Climdex thresholds output file.
#'
#' This function creates a file suitable for outputting thresholds to, with all variables that can be created with the input data present in the file.
#'
#' @param thresholds.file The filename to be used for the thresholds file.
#' @param f The file(s) being used as sources for metadata.
#' @param ts The associated time data, as created by \code{nc.get.time.series}.
#' @param v.f.idx A mapping from variables to files, as created by \code{\link{get.var.file.idx}}.
#' @param variable.name.map A mapping from standardized names (tmax, tmin, prec) to NetCDF variable names.
#' @param base.range The base range; a vector of two numeric years.
#' @param dim.size Dimension sizes for the input.
#' @param dim.axes Dimension axes for the input.
#' @param threshold.dat Threshold metadata, as provided by \code{\link{get.thresholds.metadata}}.
#' @param author.data A vector containing named elements describing the author; see \code{\link{create.indices.from.files}}.
#' @return An object of class \code{ncdf4}.
#'
#' @examples
#' \donttest{
#' ## Establish basic inputs.
#' author.data <- list(institution="Looney Bin", institution_id="LBC")
#' input.files <- c("pr_NAM44_CanRCM4_ERAINT_r1i1p1_1989-2009.nc")
#'
#' ## Prepare derived inputs.
#' f <- lapply(input.files, ncdf4::nc_open)
#' variable.name.map <- c(tmax="tasmax", tmin="tasmin", prec="pr")
#' f.meta <- create.file.metadata(f, variable.name.map)
#' threshold.dat <- get.thresholds.metadata(names(f.meta$v.f.idx))
#'
#' ## Create output file
#' thresh.file <- create.thresholds.file("thresh.nc", f, f.meta$ts, f.meta$v.f.idx, variable.name.map,
#'                                       c(1981,1990), f.meta$dim.size, f.meta$dim.axes,
#'                                       threshold.dat, author.data)
#' }
#'
#' @export
create.thresholds.file <- function(thresholds.file, f, ts, v.f.idx, variable.name.map, base.range, dim.size, dim.axes, threshold.dat, author.data) {
  exemplar.file <- f[[v.f.idx[1]]]
  exemplar.var.name <- variable.name.map[names(v.f.idx)[1]]
  exemplar.var <- exemplar.file$var[[exemplar.var.name]]
  num.thresholds <- ifelse(is.null(attr(ts, "dpy")), 365, attr(ts, "dpy"))
  cal <- attr(ts, "cal")

  ## Get time metadata...
  old.time.dim <- exemplar.var$dim[[which(dim.axes == "T")]]
  time.units <- old.time.dim$units
  time.units.split <- strsplit(time.units, " ")[[1]]
  time.origin <- if(time.units.split[2] == "as") format(trunc(min(ts), units="days"), "%Y-%m-%d") else time.units.split[3]
  time.dim.name <- old.time.dim$name
  old.time.bnds.att <- ncdf4::ncatt_get(exemplar.file, time.dim.name, "bounds")
  time.bnds.name <- if(old.time.bnds.att$hasatt) old.time.bnds.att$value else "time_bnds"

  ## Set up time variables
  out.time <- as.numeric(julian(as.PCICt(paste(floor(mean(base.range)), 1:num.thresholds, sep="-"), attr(ts, "cal"), format="%Y-%j"), as.PCICt(time.origin, cal)), units="days")
  out.time.dim <- ncdf4::ncdim_def("time", paste("days since", time.origin), out.time, unlim=TRUE, calendar=cal, longname="time")

  ## Set up bounds
  input.bounds <- ncdf4.helpers::nc.get.dim.bounds.var.list(exemplar.file)
  input.bounds <- input.bounds[input.bounds != time.bnds.name]
  input.dim.names <- ncdf4.helpers::nc.get.dim.names(exemplar.file, exemplar.var.name)
  input.varname.list <- c(input.bounds, input.dim.names)

  bnds.dim <- ncdf4::ncdim_def("bnds", "", 1:2, create_dimvar=FALSE)
  if(length(input.bounds) > 0)
    bnds.dim <- exemplar.file$var[[input.bounds[1]]]$dim[[1]]
  out.time.bnds <- as.numeric(julian(as.PCICt(c(paste(base.range[1], 1:num.thresholds, sep="-"), paste(base.range[2], 1:num.thresholds, sep="-")), attr(ts, "cal"), format="%Y-%j"), as.PCICt(time.origin, cal)), units="days")
  dim(out.time.bnds) <- c(num.thresholds, 2)
  out.time.bnds <- t(out.time.bnds)
  out.time.bnds.var <- ncdf4::ncvar_def(time.bnds.name, '', list(bnds.dim, out.time.dim), longname='', prec="double")

  input.bounds.vars <- c(lapply(input.bounds, function(x) { exemplar.file$var[[x]] }), list(out.time.bnds.var))
  input.bounds.data <- c(lapply(input.bounds, function(x) { ncdf4::ncvar_get(exemplar.file, x) }), list(out.time.bnds))
  all.bounds <- c(input.bounds, time.bnds.name)
  names(input.bounds.data) <- names(input.bounds.vars) <- all.bounds

  ## Set up 2d and 3d dims
  out.dims.3d <- list(exemplar.var$dim[[which(dim.axes == 'X')]], exemplar.var$dim[[which(dim.axes == 'Y')]], out.time.dim)
  out.dims.2d <- list(exemplar.var$dim[[which(dim.axes == 'X')]], exemplar.var$dim[[which(dim.axes == 'Y')]])
  out.vars <- sapply(names(threshold.dat), function(n) {
    tinfo <- threshold.dat[[n]]
    if(tinfo$has.time)
      ncdf4::ncvar_def(n, tinfo$units, out.dims.3d, 1e20, tinfo$longname, prec="double")
    else
      ncdf4::ncvar_def(n, tinfo$units, out.dims.2d, 1e20, tinfo$longname, prec="double")
  }, simplify=FALSE)

  ## Tack bounds vars onto var list so they get created...
  all.vars <- c(input.bounds.vars, out.vars)

  ## Create file
  thresholds.netcdf <- ncdf4::nc_create(thresholds.file, all.vars, force_v4=TRUE)
  out.dim.axes <- c("X", "Y", "T")

  ## Copy attributes for all variables plus global attributes
  ncdf4::nc_redef(thresholds.netcdf)
  ncdf4.helpers::nc.copy.atts(exemplar.file, 0, thresholds.netcdf, 0, definemode=TRUE)
  for(v in input.varname.list) {
    ncdf4.helpers::nc.copy.atts(exemplar.file, v, thresholds.netcdf, v, definemode=TRUE)
  }

  put.ETCCDI.atts(thresholds.netcdf, "monClim", ncdf4::ncatt_get(exemplar.file, 0, "title")$value, author.data, definemode=TRUE)

  ## Attach history data to threshold data.
  lapply(out.vars, function(v) {
    put.history.att(thresholds.netcdf, v, definemode=TRUE)
    ncdf4::ncatt_put(thresholds.netcdf, v, "base_period", paste(base.range, collapse="-"), definemode=TRUE)
  })
  ncdf4::nc_enddef(thresholds.netcdf)

  ## Put bounds data.
  for(v in all.bounds) {
    ncdf4::ncvar_put(thresholds.netcdf, v, input.bounds.data[[v]])
  }

  return(thresholds.netcdf)
}

#' Create mapping from variables to files.
#'
#' Create mapping from variables to files.
#'
#' Given a variable name map and list of variables in each file, determine a mapping from variables to files.
#'
#' @param variable.name.map A mapping from standardized names (tmax, tmin, prec) to NetCDF variable names.
#' @param v.list A list containing a vector of variables in each file.
#' @return A vector mapping standardized variable names (tmax, tmin, prec) to indices in the file list.
#'
#' @examples
#' \dontrun{
#' ## Get mapping for a single file.
#' input.files <- c("pr_NAM44_CanRCM4_ERAINT_r1i1p1_1989-2009.nc")
#' f <- lapply(input.files, ncdf4::nc_open)
#' v.list <- lapply(f, ncdf4.helpers::nc.get.variable.list, min.dims=2)
#' v.f.idx <- get.var.file.idx(variable.name.map, v.list)
#' }
#'
#' @export
get.var.file.idx <- function(variable.name.map, v.list) {
  v.f.idx <- sapply(variable.name.map, function(v) { which(sapply(v.list, function(vl) { v %in% vl })) }, simplify=FALSE)
  v.f.idx <- unlist(v.f.idx[sapply(v.f.idx, length) > 0])
  return(v.f.idx)
}

#' Retrieve metadata about NetCDF-format files.
#'
#' Retrieve metadata about NetCDF-format files.
#'
#' Given a list of NetCDF files and a mapping from standard variable names (tmax, tmin, prec) to NetCDF variable names, retrieve a set of standardized metadata.
#'
#' @param f The list of NetCDF files.
#' @param variable.name.map A named character vector mapping standard variable names (tmax, tmin, prec) to NetCDF variable names.
#' @return A list containing time data (ts), dimension sizes (dim.size), dimension axes (dim.axes), source units (src.units), destination units (dest.units), a mapping from variables to files (v.f.idx), and a projection, if available.
#'
#' @examples
#' \dontrun{
#' ## Get metadata about a single input file.
#' input.files <- c("pr_NAM44_CanRCM4_ERAINT_r1i1p1_1989-2009.nc")
#' f <- lapply(input.files, ncdf4::nc_open)
#' f.meta <- create.file.metadata(f, variable.name.map)
#' }
#'
#' @export
create.file.metadata <- function(f, variable.name.map) {
  v.list <- lapply(f, ncdf4.helpers::nc.get.variable.list, min.dims=2)
  v.f.idx <- get.var.file.idx(variable.name.map, v.list)

  if(any(sapply(v.list, function(vl) { sum(variable.name.map %in% vl) }) == 0))
    stop("At least one input file doesn't contain any of the named variables.")
  if(anyDuplicated(unlist(names(v.f.idx))))
    stop("Variable(s) present in more than one input file.")

  ## Get units and specify destination units
  dest.units <- c(prec="kg m-2 d-1", tmax="degrees_C", tmin="degrees_C", tavg="degrees_C")
  dest.units <- dest.units[names(dest.units) %in% names(v.f.idx)]

  ## Get projection
  projection <- ncdf4.helpers::nc.get.proj4.string(f[[1]], v.list[[1]][1])
  stopifnot(!is.null(projection))
  if(projection == "")
    projection <- NULL

  return(list(ts=get.ts(f), dim.size=get.dim.size(f, v.f.idx, variable.name.map), dim.axes=get.dim.axes(f, v.f.idx, variable.name.map),
              src.units=sapply(names(v.f.idx), function(i) { f[[v.f.idx[i]]]$var[[variable.name.map[i]]]$units }),
              dest.units=dest.units, v.f.idx=v.f.idx, projection=projection))
}

#' Retrieve threshold metadata
#'
#' Retrieve threshold metadata
#'
#' Returns units, long names, locations within the climdexInput data structure, and whether time data should be included given the variable information available.
#'
#' @param var.names A vector containing names of available variables (tmax, tmin, prec).
#' @return A list containing metadata for each of the six thresholds.
#'
#' @examples
#' thresholds.meta <- get.thresholds.metadata("prec")
#'
#' @export
get.thresholds.metadata <- function(var.names) {
  threshold.dat <- list(tx10thresh=list(units="degrees_C", longname="10th_percentile_running_baseline_tasmax", has.time=TRUE, q.path=c("tmax", "outbase", "q10")),
                        tx90thresh=list(units="degrees_C", longname="90th_percentile_running_baseline_tasmax", has.time=TRUE, q.path=c("tmax", "outbase", "q90")),
                        tn10thresh=list(units="degrees_C", longname="10th_percentile_running_baseline_tasmin", has.time=TRUE, q.path=c("tmin", "outbase", "q10")),
                        tn90thresh=list(units="degrees_C", longname="90th_percentile_running_baseline_tasmin", has.time=TRUE, q.path=c("tmin", "outbase", "q90")),
                        r95thresh=list(units="kg m-2 d-1", longname="95th_percentile_baseline_wet_day_pr", has.time=FALSE, q.path=c("prec", "q95")),
                        r99thresh=list(units="kg m-2 d-1", longname="99th_percentile_baseline_wet_day_pr", has.time=FALSE, q.path=c("prec", "q99")))
  return(threshold.dat[sapply(threshold.dat, function(x) { x$q.path[1] %in% var.names })])
}

unsquash.dims <- function(dat.dim, subset, f, n) {
  dim.axes <- ncdf4.helpers::nc.get.dim.axes(f, n)
  return(sapply(dim.axes, function(x) { if(any(names(subset) == x)) length(subset[[x]]) else f$dim[[names(dim.axes)[dim.axes == x]]]$len }))
}

## Run Climdex to generate indices for computing Climdex on future data
#' Create Climdex thresholds used for computing threshold-based indices
#'
#' Create Climdex thresholds used for computing threshold-based indices
#'
#' For many applications, one may want to compute thresholds on one data set, then apply them to another. This is usually the case when comparing GCM (Global Climate Model) results for future time periods to either historical reanalysis data or historical / pre-industrial control runs from models. The purpose of this function is to compute these thresholds on the data supplied, saving them to the file specified. Then these thresholds can be used with \code{\link{create.indices.from.files}} to compute indices using the thresholds computed using this code.
#'
#' @param input.files A list of filenames of NetCDF files to be used as input. A NetCDF file may contain one or more variables.
#' @param output.file The name of the file to be created.
#' @param author.data A vector containing named elements describing the author; see \code{\link{create.indices.from.files}}.
#' @param variable.name.map A character vector mapping from standardized names (tmax, tmin, prec) to NetCDF variable names.
#' @param axis.to.split.on The axis to split up the data on for parallel / incremental processing.
#' @param fclimdex.compatible Whether the thresholds should be created to match fclimdex thresholds; affects padding at the ends of the base period.
#' @param base.range Vector of two numeric years specifying the start and end years.
#' @param parallel The number of parallel processing threads, or FALSE if no parallel processing is desired.
#' @param verbose Whether to be chatty.
 #' @param max.vals.millions The number of data values to process at one time (length of time dim * number of values * number of variables).
#' @param cluster.type The cluster type, as used by the \code{snow} library.
#'
#' @note NetCDF input files may contain one or more variables, named as per \code{variable.name.map}. The code will search the files for the named variables.
#'
#' @examples
#' \dontrun{
#' ## Prepare input data and calculate thresholds for file.
#' input.files <- c("pr_NAM44_CanRCM4_ERAINT_r1i1p1_1989-2009.nc")
#' author.data <- list(institution="Looney Bin", institution_id="LBC")
#' create.thresholds.from.file(input.files, "thresh.nc", author.data,
#'                             base.range=c(1991, 2000), parallel=FALSE)
#' }
#'
#' @export
create.thresholds.from.file <- function(input.files, output.file, author.data, variable.name.map=c(tmax="tasmax", tmin="tasmin", prec="pr", tavg="tas"), axis.to.split.on="Y", fclimdex.compatible=TRUE, base.range=c(1961, 1990), parallel=4, verbose=FALSE, max.vals.millions=10, cluster.type="SOCK", src="ncdf.R") {
  if(!(is.logical(parallel) || is.numeric(parallel)))
    stop("'parallel' option must be logical or numeric.")

  if(length(input.files) == 0)
    stop("Require at least one input file.")

  f <- lapply(input.files, ncdf4::nc_open)
  f.meta <- create.file.metadata(f, variable.name.map)

  ## Define what the threshold indices will look like...
  threshold.dat <- get.thresholds.metadata(names(f.meta$v.f.idx))

  ## Create the output file
  thresholds.netcdf <- create.thresholds.file(output.file, f, f.meta$ts, f.meta$v.f.idx, variable.name.map, base.range, f.meta$dim.size, f.meta$dim.axes, threshold.dat, author.data)

  cluster <- set.up.cluster(parallel, cluster.type, src)
  subsets <- ncdf4.helpers::get.cluster.worker.subsets(max.vals.millions * 1000000, f.meta$dim.size, f.meta$dim.axes, axis.to.split.on)

  write.thresholds.data <- function(out.list, out.sub) {
    lapply(names(threshold.dat), function(n) {
      d <- threshold.dat[[n]]
      if(d$has.time)
        dat <- t(sapply(out.list, function(y) { return(y[[d$q.path]]) }))
      else
        dat <- sapply(out.list, function(y) { return(y[[d$q.path[1]]][d$q.path[2]]) })
      dim(dat) <- unsquash.dims(dim(dat), out.sub, thresholds.netcdf, n)
      ncdf4.helpers::nc.put.var.subset.by.axes(thresholds.netcdf, n, dat, out.sub)
    })
    gc()
  }

  if(!is.null(cluster)) {
    lapply(f, ncdf4::nc_close)
    rm(f)

    snow::clusterExport(cluster, "input.files", environment())
    snow::clusterEvalQ(cluster, f <<- lapply(input.files, ncdf4::nc_open, readunlim=FALSE))

    ## Compute subsets and fire jobs off; collect and write out chunk-at-a-time
    parLapplyLBFiltered(cluster, subsets, get.quantiles.for.stripe, f.meta$ts, base.range, f.meta$dim.axes, f.meta$v.f.idx, variable.name.map, f.meta$src.units, f.meta$dest.units, local.filter.func=write.thresholds.data)

    snow::stopCluster(cluster)
  } else {
    ##try(getFromNamespace('nc_set_chunk_cache', 'ncdf4')(1024 * 2048, 1009), silent=TRUE)

    lapply(subsets, function(x) { write.thresholds.data(get.quantiles.for.stripe(x, f.meta$ts, base.range, f.meta$dim.axes, f.meta$v.f.idx, variable.name.map, f.meta$src.units, f.meta$dest.units, f), x) })

    lapply(f, ncdf4::nc_close)
  }

  ## Close all the files
  ncdf4::nc_close(thresholds.netcdf)

  cat("Finished computing thresholds\n")
  invisible(0)
}

#' Open thresholds file(s)
#'
#' Open thresholds file(s)
#'
#' This function opens one or more thresholds files and returns the \code{ncdf4} objects as a list.
#'
#' @param thresholds.files A character vector containing the names of thresholds files.
#' @return A list of objects of class \code{ncdf4}, or NULL if thresholds.files is NULL.
#'
#' @examples
#' \dontrun{
#' ## Open a single thresholds file
#' thresholds.files <- c("thresh.nc")
#' thresh <- thresholds.open(thresholds.files)
#' }
#'
#' @export
thresholds.open <- function(thresholds.files) {
  return(if(is.null(thresholds.files)) NULL else lapply(thresholds.files, ncdf4::nc_open))
}

#' Close thresholds file(s)
#'
#' Close thresholds file(s)
#'
#' This function closes one or more thresholds files.
#'
#' @param thresholds.nc A list of objects of class \code{ncdf4}, or NULL
#'
#' @examples
#' \dontrun{
#' ## Open a single thresholds file, then close it.
#' thresholds.files <- c("thresh.nc")
#' thresh <- thresholds.open(thresholds.files)
#' thresholds.close(thresh)
#' }
#'
#' @export
thresholds.close <- function(thresholds.nc) {
  if(!is.null(thresholds.nc)) lapply(thresholds.nc, ncdf4::nc_close)
  invisible(0)
}


get.time.origin <- function(f, dim.axes) {
  time.units <- f[[1]]$dim[[names(dim.axes)[which(dim.axes == "T")]]]$units
  time.units.split <- strsplit(gsub("[ ]+", " ", time.units), " ")[[1]]
  time.origin <- if(time.units.split[2] == "as") format(trunc(min(ts), units="days"), "%Y-%m-%d") else time.units.split[3]
  return(time.origin)
}

get.thresholds.f.idx <- function(thresholds.files, thresholds.name.map) {
  if(is.null(thresholds.files)) {
    return(NULL)
  } else {
    thresh <- thresholds.open(thresholds.files)
    t.f.idx <- get.var.file.idx(thresholds.name.map, lapply(thresh, ncdf4.helpers::nc.get.variable.list, min.dims=2))
    thresholds.close(thresh)
    return(t.f.idx)
  }
}

## Run Climdex and populate the output files
#' Create Climdex indices from NetCDF input files.
#'
#' Create Climdex indices from NetCDF input files.
#'
#' This function computes Climdex indices from NetCDF input files, writing out one file per variable named like the \code{template.filename}, which must follow the CMIP5 file naming conventions (this is a deficiency which will be corrected in later versions).
#'
#' The indices to be calculated can be specified; if not, they will be determined by data availability. Thresholds can be supplied (via \code{thresholds.files}) or, if there is data within the base period, calculated and used as part of the process. Note that in-base thresholds are separate from out-of-base thresholds; this is covered in more detail in the help for the \code{climdex.pcic} package.
#'
#' @param input.files A list of filenames of NetCDF files to be used as input. A NetCDF file may contain one or more variables.
#' @param out.dir The directory to put the output files in.
#' @param output.filename.template The output filename to be used as a template, which must follow the CMIP5 file naming conventions.
#' @param author.data Data describing the author; a character vector with 0 or more of the following named values:\describe{
#' \item{institution}{The institution generating the data.}
#' \item{institution_id}{An abbreviation for the institution generating the data.}
#' \item{indices_archive}{The URL the data is published at, if applicable.}
#' \item{contact}{The email address or contact info for the author.}
#' \item{references}{What to reference when citing this work.}
#' }
#' @param climdex.vars.subset A character vector of lower-case names of Climdex indices to calculate (eg: tr, fd, rx5day). See the list of 27 indices in the References section.
#' @param climdex.time.resolution The time resolution to compute indices at; one of "all" (both monthly and annual), "annual" (only annual), or "monthly" (only monthly).
#' @param variable.name.map A character vector mapping from standardized names (tmax, tmin, prec) to NetCDF variable names.
#' @param axis.to.split.on The axis to split up the data on for parallel / incremental processing.
#' @param fclimdex.compatible Whether the thresholds should be created to match fclimdex thresholds; affects padding at the ends of the base period.
#' @param base.range Vector of two numeric years specifying the start and end years.
#' @param parallel The number of parallel processing threads, or FALSE if no parallel processing is desired.
#' @param verbose Whether to be chatty.
#' @param thresholds.files A character vector of files containing thresholds to be used.
#' @param thresholds.name.map A mapping from threshold names to NetCDF variable names. The following names will be used: \describe{
#' \item{tx10thresh}{10th percentile for a 5-day running window of baseline daily maximum temperature.}
#' \item{tn10thresh}{10th percentile for a 5-day running window of baseline daily minimum temperature.}
#' \item{tx90thresh}{90th percentile for a 5-day running window of baseline daily maximum temperature.}
#' \item{tn90thresh}{90th percentile for a 5-day running window of baseline daily minimum temperature.}
#' \item{r95thresh}{95th percentile of daily precipitation in wet days (>=1 mm of rain).}
#' \item{r99thresh}{99th percentile of daily precipitation in wet days (>=1 mm of rain).}
#' }
#' @param max.vals.millions The number of data values to process at one time (length of time dim * number of values * number of variables).
#' @param cluster.type The cluster type, as used by the \code{snow} library.
#'
#' @note NetCDF input files may contain one or more variables, named as per \code{variable.name.map}. The code will search the files for the named variables. The same is true of thresholds files; one file may be supplied, or multiple files may be supplied, via the \code{thresholds.files} argument; and the name mapping may be supplied via the \code{thresholds.name.map} argument.
#'
#' @references \url{http://etccdi.pacificclimate.org/list_27_indices.shtml}
#' @examples
#' \dontrun{
#' ## Prepare input data and calculate indices for a single file
#' ## with a single thread (no parallelism).
#' input.files <- c("pr_NAM44_CanRCM4_ERAINT_r1i1p1_1989-2009.nc")
#' author.data <- list(institution="Looney Bin", institution_id="LBC")
#' create.indices.from.files(input.files, "out_dir/", input.files[1], author.data,
#'                           base.range=c(1991, 2000), parallel=FALSE)
#'
#' ## Prepare input data and calculate indices for two files
#' ## in parallel given thresholds.
#' input.files <- c("pr_NAM44_CanRCM4_ERAINT_r1i1p1_1989-2009.nc",
#'                  "tasmax_NAM44_CanRCM4_ERAINT_r1i1p1_1989-2009.nc")
#' author.data <- list(institution="Looney Bin", institution_id="LBC")
#' create.indices.from.files(input.files, "out_dir/", input.files[1], author.data,
#'                           base.range=c(1991, 2000), parallel=8, thresholds.files="thresh.nc")
#' }
#'
#' @export
create.indices.from.files <- function(input.files, out.dir, output.filename.template, author.data, climdex.vars.subset=NULL, climdex.time.resolution=c("all", "annual", "monthly"), variable.name.map=c(tmax="tasmax", tmin="tasmin", prec="pr", tavg="tas"), axis.to.split.on="Y", fclimdex.compatible=TRUE, base.range=c(1961, 1990), parallel=4, verbose=FALSE, thresholds.files=NULL, thresholds.name.map=c(tx10thresh="tx10thresh", tn10thresh="tn10thresh", tx90thresh="tx90thresh", tn90thresh="tn90thresh", r95thresh="r95thresh", r99thresh="r99thresh"), max.vals.millions=10, cluster.type="SOCK", src="ncdf.R") {
  if(!(is.logical(parallel) || is.numeric(parallel)))
    stop("'parallel' option must be logical or numeric.")

  if(length(input.files) == 0)
    stop("Require at least one input file.")

  ## Open files, determine mapping between files and variables.
  f <- lapply(input.files, ncdf4::nc_open)
  f.meta <- create.file.metadata(f, variable.name.map)

  ## Get thresholds variable-file mapping
  t.f.idx <- get.thresholds.f.idx(thresholds.files, thresholds.name.map)

  ## Get variable list, subset if necessary
  climdex.time.resolution <- match.arg(climdex.time.resolution)
  climdex.var.list <- get.climdex.variable.list(names(f.meta$v.f.idx), climdex.time.resolution, climdex.vars.subset)

  cdx.meta <- get.climdex.variable.metadata(climdex.var.list, output.filename.template)
  cdx.ncfile <- create.ncdf.output.files(cdx.meta, f, f.meta$v.f.idx, variable.name.map, f.meta$ts, get.time.origin(f, f.meta$dim.axes), base.range, out.dir, author.data)
  cdx.funcs <- get.climdex.functions(climdex.var.list)

  ## Compute indices, either single process or multi-process using 'parallel'
  subsets <- ncdf4.helpers::get.cluster.worker.subsets(max.vals.millions * 1000000, f.meta$dim.size, f.meta$dim.axes, axis.to.split.on)
  if(is.numeric(parallel)) {
    ## Setup...
    lapply(f, ncdf4::nc_close)
    rm(f)
    cluster <- set.up.cluster(parallel, cluster.type, src)
    snow::clusterExport(cluster, list("input.files", "thresholds.files"), environment())
    snow::clusterEvalQ(cluster, f <<- lapply(input.files, ncdf4::nc_open, readunlim=FALSE))
    snow::clusterEvalQ(cluster, thresholds.netcdf <<- thresholds.open(thresholds.files))

    ## Meat...
    parLapplyLBFiltered(cluster, subsets, compute.indices.for.stripe, cdx.funcs, f.meta$ts, base.range, f.meta$dim.axes, f.meta$v.f.idx, variable.name.map, f.meta$src.units, f.meta$dest.units, t.f.idx, thresholds.name.map, fclimdex.compatible, f.meta$projection, local.filter.func=function(x, x.sub) {
      write.climdex.results(x, x.sub, cdx.ncfile, f.meta$dim.size, cdx.meta$var.name)
    })

    ## Clean-up.
    snow::stopCluster(cluster)
  } else {
    ## Setup...
    thresholds.netcdf <- thresholds.open(thresholds.files)
    ##try(getFromNamespace('nc_set_chunk_cache', 'ncdf4')(1024 * 2048, 1009), silent=TRUE)

    ## Meat...
    lapply(subsets, function(x) { write.climdex.results(compute.indices.for.stripe(x, cdx.funcs, f.meta$ts, base.range, f.meta$dim.axes, f.meta$v.f.idx, variable.name.map, f.meta$src.units, f.meta$dest.units, t.f.idx, thresholds.name.map, fclimdex.compatible, f.meta$projection, f, thresholds.netcdf), x, cdx.ncfile, f.meta$dim.size, cdx.meta$var.name) })

    ## Clean-up.
    thresholds.close(thresholds.netcdf)
    lapply(f, ncdf4::nc_close)
  }

  ## Close all the output files
  lapply(cdx.ncfile, ncdf4::nc_close)

  cat("Finished computing indices\n")
  invisible(0)
}
# nolint end
