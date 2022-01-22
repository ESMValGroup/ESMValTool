log <- function(..., level = "INFO") {
  cat(format(Sys.time(), "%Y-%m-%d %X"), level, ":", ..., "\n")
}

# check for present library paths
rlibpath <- .libPaths()
Sys.setenv(R_INSTALL_STAGED = FALSE)

# check if we can write in the present R libraries paths
if (any(file.access(rlibpath, 2) == 0)) {
  # if possible, use the standard one for following installation
  rlibloc <- rlibpath[which(file.access(rlibpath, 2) == 0)[1]]
} else {
  # if not possible, create a local library in the home directory
  rlibloc <- Sys.getenv("R_LIBS_USER")
  dir.create(
    path = Sys.getenv("R_LIBS_USER"),
    showWarnings = FALSE,
    recursive = TRUE
  )
}

log("Installing packages to --> ", rlibloc)

# define the R mirror to download packages
pkg_mirror <- "https://cloud.r-project.org"
log("Using mirror: ", pkg_mirror)

# install the yaml package fist to handle the yaml-file
if (!require("yaml")) {
  install.packages("yaml", repos = pkg_mirror)
}

# get the script path
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_name <- sub(
  file_arg_name, "",
  initial_options[grep(file_arg_name, initial_options)]
)
script_dirname <- dirname(script_name)

# read the dependencies
dependencies <- yaml::read_yaml(
  paste(script_dirname, "r_requirements.yml", sep = "/"),
)
# TODO: find a solution for script directory

log(paste(unlist(names(dependencies)), collapse = ", "))
inst_packages <- installed.packages()
package_list <-
  names(dependencies)[!(names(dependencies) %in% inst_packages[, "Package"])]

if (length(package_list) == 0) {
  log("All packages are already installed!")
  quit()
} else {
  log("Number of packages to be installed: ", length(package_list))
}

n_cpus_override <- strtoi(Sys.getenv("R_INSTALL_N_CPUS", unset = "0"))
if (n_cpus_override) {
  n_cpus <- n_cpus_override
} else {
  n_cpus <- parallel::detectCores(all.tests = TRUE, logical = FALSE)
  if (is.na(n_cpus)) {
    n_cpus <- 1
  }
}
log("Using", n_cpus, "threads to compile packages")

# Install packages required for installation
install_dependencies <- list("remotes", "devtools")
install_dependencies <-
  install_dependencies[!install_dependencies %in% inst_packages[, "Package"]]
if (length(install_dependencies) > 0) {
  log(
    "Installing installation dependencies from CRAN:",
    paste(unlist(install_dependencies), collapse = ", ")
  )
  install.packages(
    unlist(install_dependencies),
    repos = pkg_mirror,
    Ncpus = n_cpus,
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}

# Special installations
install <- function(package, info) {
  if (!is.null(info)) {
    if (info$origin == "github") {
      log("Installing ", package, "from Github")
      remotes::install_github(info$path)
    }
    if (info$origin == "git") {
      log("Installing ", package, "from Git")
      remotes::install_git(info$url, ref = info$ref)
    }
    if (info$origin == "cran") {
      log("Installing version", info$version, "of", package, "from CRAN")
      devtools::install_version(
        package,
        version = info$version,
        repos = pkg_mirror
      )
    }
  }
}
for (package_name in names(dependencies)) {
  install(package_name, dependencies[[package_name]])
}

# Get missing packages from CRAN, last versions
inst_packages <- installed.packages()
package_list <-
  names(dependencies)[!(names(dependencies) %in% inst_packages[, "Package"])]
log(
  "Installing packages from CRAN:",
  paste(unlist(package_list), collapse = ", ")
)
if (length(package_list) != 0) {
  install.packages(
    package_list,
    repos = pkg_mirror,
    Ncpus = n_cpus,
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}

failed <- list()
for (package_name in names(dependencies)) {
  success <- library(
    package_name,
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
