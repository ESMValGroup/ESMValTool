# mini script to compile and install lintr (r-lintr) and yaml (r-yaml) from
# source without doing this from conda; conda installations introduce bugs in
# shared libs
log <- function(..., level = "INFO") {
  cat(format(Sys.time(), "%Y-%m-%d %X"), level, ":", ..., "\n")
}

# check for present library paths
rlibpath <- .libPaths()

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

pkg_mirror <- "https://cloud.r-project.org"
pkg_list <- c("devtools", "yaml", "lintr", "styler", "docopt")
install.packages(
  pkg_list,
  repos = pkg_mirror,
  Ncpus = 1,
  dependencies = c("Depends", "Imports", "LinkingTo")
)

failed <- list()
for (package_name in pkg_list) {
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
