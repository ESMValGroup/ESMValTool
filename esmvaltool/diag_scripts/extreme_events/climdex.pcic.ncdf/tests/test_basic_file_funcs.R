#nolint start
author.data <- list(institution="Pacific Climate Impacts Consortium",
                    institution_id="PCIC",
                    indices_archive="Please check http://www.cccma.ec.gc.ca/data/climdex/climdex.shtml for errata or updates.",
                    contact="bronaugh@uvic.ca",
                    references="http://www.cccma.ec.gc.ca/data/climdex/"
                    )
x.subset <- 33:34
y.subset <- 17:18
correct.data.dir <- "correct_output/"
correct.thresh.file.6190 <- paste(correct.data.dir, "thresholds_monClim_CanESM2_historical_r1i1p1_1961-1990.nc", sep="")
thresh.omit.list <- c("tx10p", "tn10p", "tx10p", "tx90p", "wsdi", "csdi")

test.get.thresholds.chunk <- function() {
  ## Define mappings and filenames.
  thresholds.name.map <- c(tx10thresh="tx10thresh", tn10thresh="tn10thresh", tx90thresh="tx90thresh",
                           tn90thresh="tn90thresh", r95thresh="r95thresh", r99thresh="r99thresh")
  thresh.files <- "correct_output/thresholds_monClim_CanESM2_historical_r1i1p1_1961-1990.nc"

  if(all(file.exists(thresh.files))) {
    ## Open files, etc.
    cdx.funcs <- get.climdex.functions(get.climdex.variable.list("tmax"))
    thresholds.netcdf <- lapply(thresh.files, nc_open)
    t.f.idx <- get.var.file.idx(thresholds.name.map, lapply(thresholds.netcdf,
                                                            ncdf4.helpers::nc.get.variable.list, min.dims=2))

    ## Get thresholds chunk.
    dat <- get.thresholds.chunk(list(), cdx.funcs, thresholds.netcdf, t.f.idx, thresholds.name.map)
    checkEquals(thresholds.chunk.tmax.only, dat)

    lapply(thresholds.netcdf, nc_close)
  }
}

test.compute.indices.for.stripe <- function() {
  ## Define mappings and filenames.
  author.data <- list(institution="Looney Bin", institution_id="LBC")
  input.files <- list.files("test1/", full.names=TRUE)
  variable.name.map <- c(tmax="tasmax", tmin="tasmin", prec="pr")
  thresh.files <- "correct_output/thresholds_monClim_CanESM2_historical_r1i1p1_1961-1990.nc"

  if(length(input.files) > 0 && all(file.exists(input.files)) && all(file.exists(thresh.files))) {
    ## Open files, etc.
    cdx.funcs <- get.climdex.functions(get.climdex.variable.list("prec"))
    f <- lapply(input.files, ncdf4::nc_open)
    f.meta <- create.file.metadata(f, variable.name.map)
    climdex.var.list <- get.climdex.variable.list(names(f.meta$v.f.idx), "all", NULL)
    cdx.meta <- get.climdex.variable.metadata(climdex.var.list, input.files[1])

    thresholds.name.map <- c(tx10thresh="tx10thresh", tn10thresh="tn10thresh", tx90thresh="tx90thresh",
                             tn90thresh="tn90thresh", r95thresh="r95thresh", r99thresh="r99thresh")
    thresholds.netcdf <- lapply(thresh.files, nc_open)
    t.f.idx <- get.var.file.idx(thresholds.name.map, lapply(thresholds.netcdf,
                                                            ncdf4.helpers::nc.get.variable.list, min.dims=2))
    thresh.dat <- get.thresholds.chunk(list(), cdx.funcs, thresholds.netcdf, t.f.idx, thresholds.name.map)

    ## Compute indices for stripe
    cdx <- compute.indices.for.stripe(list(X=1:2, Y=1:2), cdx.funcs, f.meta$ts, c(1981, 1990), f.meta$dim.axes,
                                      f.meta$v.f.idx, variable.name.map, f.meta$src.units, f.meta$dest.units,
                                      t.f.idx, thresholds.name.map, f=f, thresholds.netcdf=thresholds.netcdf)

    lapply(thresholds.netcdf, nc_close)

    res <- lapply(names(cdx[[1]]), function(x) {
      fn <- list.files("correct_output/", pattern=paste("^", x, sep=""), full.names=TRUE)
      f.valid <- nc_open(fn, readunlim=FALSE)
      d.input <- ncvar_get(f.valid, strsplit(x, "_")[[1]][1])
      nc_close(f.valid)
      d.comparison <- t(do.call(cbind, lapply(cdx, function(cr) { cr[[x]] })))
      dim(d.comparison) <- dim(d.input)

      ## Apparently there are differences at the 3e-6 level between calculated and saved data... who knew?
      checkEquals(d.input, d.comparison, tolerance=1e-5)
      mean(abs(d.input - d.comparison))
    })

    lapply(f, nc_close)
  }
  invisible(0)
}

test.get.quantiles.for.stripe <- function() {
  historical.files <- list.files("historical/", full.names=TRUE)
  if(length(historical.files) > 0) {
    ## FIXME: This is untestable with the current input data.

    ## Establish basic inputs.
    author.data <- list(institution="Looney Bin", institution_id="LBC")
    input.files <- list.files("test1/", full.names=TRUE)

    ## Prepare derived inputs.
    f <- lapply(input.files, ncdf4::nc_open)
    variable.name.map <- c(tmax="tasmax", tmin="tasmin", prec="pr")
    f.meta <- create.file.metadata(f, variable.name.map)
    threshold.dat <- get.thresholds.metadata(names(f.meta$v.f.idx))

    ## Compute threshold quantiles for stripe
    q <- get.quantiles.for.stripe(list(Y=1), f.meta$ts, c(1981, 1990), f.meta$dim.axes,
                                  f.meta$v.f.idx, variable.name.map, f.meta$src.units,
                                  f.meta$dest.units, f)

    ## FIXME: Compare to valid data.

    lapply(f, nc_close)
  }
}

test.get.quantiles.object <- function() {
  ## Define mappings and filenames.
  thresholds.name.map <- c(tx10thresh="tx10thresh", tn10thresh="tn10thresh", tx90thresh="tx90thresh",
                           tn90thresh="tn90thresh", r95thresh="r95thresh", r99thresh="r99thresh")
  thresh.files <- "correct_output/thresholds_monClim_CanESM2_historical_r1i1p1_1961-1990.nc"

  if(all(file.exists(thresh.files))) {
    ## Open files, etc.
    cdx.funcs <- get.climdex.functions(get.climdex.variable.list("tmax"))
    thresholds.netcdf <- lapply(thresh.files, nc_open)
    t.f.idx <- get.var.file.idx(thresholds.name.map, lapply(thresholds.netcdf,
                                                            ncdf4.helpers::nc.get.variable.list, min.dims=2))
    ## Get thresholds chunk.
    dat <- get.thresholds.chunk(list(Y=1), cdx.funcs, thresholds.netcdf, t.f.idx, thresholds.name.map)

    ## Get quantiles object for index 2
    q <- get.quantiles.object(dat, 2)

    ## FIXME: Compare to a correct object.

    lapply(thresholds.netcdf, nc_close)
  }
}

test.get.northern.hemisphere.booleans <- function() {
  test.get.nh <- function(test.dir) {
    input.files <- list.files(test.dir, full.names=TRUE)
    f <- lapply(input.files, ncdf4::nc_open)
    f.v <- lapply(f, ncdf4.helpers::nc.get.variable.list, min.dims=2)
    bools <- get.northern.hemisphere.booleans(list(X=1:2, Y=1:2), f[[1]], f.v[[1]], NULL)
    lapply(f, ncdf4::nc_close)
    return(bools)
  }
  ## FIXME: Need test data.
  browser()
  if(file.exists("test3/"))
    checkEquals(test.get.nh("test3/"), rep(FALSE, 4))
  ## FIXME: Need test data, and a projected dataset.
  ##checkEquals(test.get.nh("test7/"), correct.data)
}

## FIXME: Needs proper test data. This is just a framework...
test.get.data <- function() {
  test.dir <- "test3/"
  variable.name.map <- c(tmax="tasmax", tmin="tasmin", prec="pr")
  if(file.exists(test.dir)) {
    input.files <- list.files(test.dir, full.names=TRUE)
    f <- lapply(input.files, ncdf4::nc_open)
    f.meta <- create.file.metadata(f, variable.name.map)
    d <- get.data(f[[f.meta$v.f.idx['prec']]], "pr", list(Y=2), "kg m-2 s-1", "kg m-2 s-1", c(lon="X",lat="Y",time="T"))
    lapply(f, ncdf4::nc_close)
  }
}

## FIXME: Needs proper test data. This is just a framework...
test.file.funcs <- function() {
  test.dir <- "test3/"
  if(file.exists(test.dir)) {
    input.files <- list.files(test.dir, full.names=TRUE)
    variable.name.map <- c(tmax="tasmax", tmin="tasmin", prec="pr")
    f <- lapply(input.files, ncdf4::nc_open)
    f.meta <- create.file.metadata(f, variable.name.map)
    lapply(f, ncdf4::nc_close)
  }
}

test.thresholds.create.and.indices <- function() {
  test.set <- paste("test", 1:6, "/", sep="")
  lapply(test.set[file.exists(test.set)], function(test) {
    input.file.list <- list.files(test, full.names=TRUE)
    print(file.exists(input.file.list))
    print(input.file.list)
    thresh.file <- tempfile()
    indices.dir.thresh <- tempdir()
    indices.dir.nothresh <- tempdir()
    create.thresholds.from.file(input.file.list, thresh.file, author.data, parallel=FALSE, base.range=c(2010, 2019))
    create.indices.from.files(input.file.list, indices.dir.thresh, input.file.list[1], author.data, parallel=FALSE, thresholds.files=correct.thresh.file.6190)

    ## Compare to base data.
    test.file.list <- list.files(indices.dir.thresh, pattern="ETCCDI")
    lapply(test.file.list, function(fn) {
      print(fn)
      f.test <- nc_open(paste(indices.dir.thresh, fn, sep="/"))
      f.correct <- nc_open(paste(correct.data.dir, fn, sep="/"))

      d.test <- ncvar_get(f.test, ncdf4.helpers::nc.get.variable.list(f.test)[1])
      d.correct <- ncvar_get(f.correct, ncdf4.helpers::nc.get.variable.list(f.correct)[1])

      checkEquals(d.test, d.correct)

      nc_close(f.test)
      nc_close(f.correct)
    })

    create.indices.from.files(input.file.list, indices.dir.thresh, input.file.list[1], author.data, parallel=FALSE, thresholds.files=thresh.file)
    create.indices.from.files(input.file.list, indices.dir.nothresh, input.file.list[1], author.data, parallel=FALSE, base.range=c(2010, 2019))

    unlink(paste(indices.dir.nothresh, "*", sep="/"))
    unlink(paste(indices.dir.thresh, "*", sep="/"))
    gc()
  })
}

parallel.thresholds.create.and.indices <- function() {
  test.set <- paste("test", 1:6, "/", sep="")
  lapply(test.set[file.exists(test.set)], function(test) {
    input.file.list <- list.files(test, full.names=TRUE)
    print(file.exists(input.file.list))
    print(input.file.list)
    thresh.file <- tempfile()
    indices.dir.thresh <- tempdir()
    indices.dir.nothresh <- tempdir()
    create.thresholds.from.file(input.file.list, thresh.file, author.data, parallel=4, base.range=c(2010, 2029))
    create.indices.from.files(input.file.list, indices.dir.thresh, input.file.list[1], author.data, parallel=4, thresholds.files=thresh.file)
    create.indices.from.files(input.file.list, indices.dir.nothresh, input.file.list[1], author.data, parallel=4, base.range=c(2010, 2029))
  })
}
#nolint end
