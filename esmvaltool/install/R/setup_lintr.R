# mini script to compile and install lintr (r-lintr) from source
pkg_mirror <- "https://cloud.r-project.org"
install.packages(
    "lintr",
    repos = pkg_mirror,
    Ncpus = 1,
    dependencies = c("Depends", "Imports", "LinkingTo")
)
