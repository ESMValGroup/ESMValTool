log <- function(..., level = "INFO") {
  cat(format(Sys.time(), "%Y-%m-%d %X"), level, ":", ..., "\n")
}

# check for present library paths
RLIBPATH <- .libPaths()

# check if we can write in the present R libaries paths
if (any(file.access(RLIBPATH, 2) == 0)) {
  # if possible, use the standard one for following instalation
  RLIBLOC <- RLIBPATH[which(file.access(RLIBPATH, 2) == 0)[1]]
} else {
  # if not possible, create a local library in the home directory
  RLIBLOC <- Sys.getenv("R_LIBS_USER")
  dir.create(
    path = Sys.getenv("R_LIBS_USER"),
    showWarnings = FALSE,
    recursive = TRUE
  )
}

log("Installing packages to --> ", RLIBLOC)

# define the R mirror to download packages
pkg_mirror <- "https://cloud.r-project.org"
log("Using mirror: ", pkg_mirror)

# get the script path
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_name <- sub(
  file_arg_name, "",
  initial_options[grep(file_arg_name, initial_options)]
)
script_dirname <- dirname(script_name)

# read the dependencies
dependencies <- scan(paste(script_dirname, "r_requirements.txt", sep = "/"),
  what = "character"
)
# TODO: find a solution for script directory
inst_packages <- installed.packages()
package_list <-
  dependencies[!(dependencies %in% inst_packages[, "Package"])]

if (length(package_list) == 0) {
  log("All packages are already installed!")
} else {
  log("Number of packages to be installed: ", length(package_list))
}

Ncpus <- parallel::detectCores()
if (is.na(Ncpus)) {
  Ncpus <- 1
}

log("Installing packages:", package_list)
if (length(package_list) != 0) {
  install.packages(
    package_list,
    repos = pkg_mirror,
    Ncpus = Ncpus,
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}

failed <- list()
for (package_name in dependencies) {
  success <- library(package_name,
    character.only = TRUE,
    logical.return = TRUE
  )
  if (!success) {
    failed <- c(failed, package_name)
  }
}
if (length(failed) != 0) {
  log("Failed to install packages:", paste(failed, collapse = ", "))
  quit(status = 1, save = "no")
}

log("Successfully installed all packages")
