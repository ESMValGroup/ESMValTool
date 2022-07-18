#nolint start
load("exemplar_data.rda")

if(require("RUnit", quietly=TRUE)) {
  ## Run all the tests
  library(climdex.pcic.ncdf)
  library(ncdf4)
  wd <- getwd()
  testsuite <- defineTestSuite("climdex.pcic.ncdf", dirs=wd, testFileRegexp = "^test_.+.R$", testFuncRegexp = "^test.+")
  climdex.pcic.ncdf.test.result <- runTestSuite(testsuite)
  printTextProtocol(climdex.pcic.ncdf.test.result)
  stopifnot(climdex.pcic.ncdf.test.result$climdex.pcic.ncdf$nFail == 0 && climdex.pcic.ncdf.test.result$climdex.pcic.ncdf$nErr == 0)
}
#nolint end
