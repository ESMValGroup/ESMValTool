#nolint start
test.get.climdex.variable.list <- function() {
  checkEquals(climdex.var.list$tavg.all, get.climdex.variable.list(c("tavg")))
  checkEquals(climdex.var.list$tmax.all, get.climdex.variable.list(c("tmax")))
  checkEquals(climdex.var.list$tmax.mon, get.climdex.variable.list(c("tmax"), time.resolution="monthly"))
  checkEquals(climdex.var.list$tmax.yr, get.climdex.variable.list(c("tmax"), time.resolution="annual"))
  checkEquals(climdex.var.list$tmax.tmin.all, get.climdex.variable.list(c("tmax", "tmin")))
  checkEquals(climdex.var.list$tmax.prec.yr, get.climdex.variable.list(c("tmax", "prec"), time.resolution="annual"))
  checkEquals(climdex.var.list$prec.mon, get.climdex.variable.list(c("prec"), time.resolution="monthly"))
  checkEquals(climdex.var.list$prec.yr, get.climdex.variable.list(c("prec"), time.resolution="annual"))
  checkEquals(climdex.var.list$tmax.tmin.prec.all, get.climdex.variable.list(c("tmax", "tmin", "prec")))
  checkEquals(climdex.var.list$tmax.tmin.prec.sub, get.climdex.variable.list(c("tmax", "tmin", "prec"), climdex.vars.subset=c("su", "tr", "cdd", "gsl")))
  checkEquals(climdex.var.list$prec.sub, get.climdex.variable.list(c("prec"), climdex.vars.subset=c("su", "tr", "cdd", "gsl")))
  checkEquals(NULL, get.climdex.variable.list(c()))
}

test.get.climdex.variable.metadata <- function() {
  fn1 <- "tasmax_NAM44_CanRCM4_ERAINT_r1i1p1_1989-2009.nc"
  fn2 <- "pr_day_CanESM2_rcp85_r2i1p1_20060101-21001231.nc"
  checkEquals(climdex.var.meta$tavg.all.1, get.climdex.variable.metadata(climdex.var.list$tavg.all, fn1))
  checkEquals(climdex.var.meta$prec.yr.2, get.climdex.variable.metadata(climdex.var.list$prec.yr, fn2))
}

test.get.climdex.functions <- function() {
  checkEquals(climdex.functions$tmax.yr, get.climdex.functions(climdex.var.list$tmax.yr))
  checkEquals(climdex.functions$tmax.tmin.prec.all.fclimdex, get.climdex.functions(climdex.var.list$tmax.tmin.prec.all))
  checkEquals(climdex.functions$tmax.tmin.prec.all.notfclimdex, get.climdex.functions(climdex.var.list$tmax.tmin.prec.all, FALSE))
}

test.create.climdex.cmip5.filenames <- function() {
  fn.split <- c(model="CanESM2", emissions="rcp45", run="r1i1p1", tstart="20100101", tend="20991231")

  valid.tmax.mon.fn <- c("txxETCCDI_mon_CanESM2_rcp45_r1i1p1_201001-209912.nc", "txnETCCDI_mon_CanESM2_rcp45_r1i1p1_201001-209912.nc",
                         "tx10pETCCDI_mon_CanESM2_rcp45_r1i1p1_201001-209912.nc", "tx90pETCCDI_mon_CanESM2_rcp45_r1i1p1_201001-209912.nc")
  valid.tmax.all.fn <- c("suETCCDI_yr_CanESM2_rcp45_r1i1p1_2010-2099.nc", "idETCCDI_yr_CanESM2_rcp45_r1i1p1_2010-2099.nc",
                         "txxETCCDI_mon_CanESM2_rcp45_r1i1p1_201001-209912.nc", "txxETCCDI_yr_CanESM2_rcp45_r1i1p1_2010-2099.nc",
                         "txnETCCDI_mon_CanESM2_rcp45_r1i1p1_201001-209912.nc", "txnETCCDI_yr_CanESM2_rcp45_r1i1p1_2010-2099.nc",
                         "tx10pETCCDI_mon_CanESM2_rcp45_r1i1p1_201001-209912.nc", "tx10pETCCDI_yr_CanESM2_rcp45_r1i1p1_2010-2099.nc",
                         "tx90pETCCDI_mon_CanESM2_rcp45_r1i1p1_201001-209912.nc", "tx90pETCCDI_yr_CanESM2_rcp45_r1i1p1_2010-2099.nc",
                         "wsdiETCCDI_yr_CanESM2_rcp45_r1i1p1_2010-2099.nc", "altwsdiETCCDI_yr_CanESM2_rcp45_r1i1p1_2010-2099.nc")

  checkEquals(valid.tmax.mon.fn, create.climdex.cmip5.filenames(fn.split, climdex.var.list$tmax.mon))
  checkEquals(valid.tmax.all.fn, create.climdex.cmip5.filenames(fn.split, climdex.var.list$tmax.all))
}

test.flatten.dims <- function() {
  dat <- structure(1:8, .Dim=c(2, 2, 2))
  valid.flat <- structure(1:8, .Dim = c(2L, 4L))
  checkEquals(flatten.dims(dat, 2:3), valid.flat)
}
#nolint end
