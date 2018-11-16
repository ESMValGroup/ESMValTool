#check for present library paths
RLIBPATH <- .libPaths()

#check if we can write in the present R libaries paths
if (any(file.access(RLIBPATH, 2) == 0)) {
    #if possible, use the standard one for following instalation
    RLIBLOC <- RLIBPATH[which(file.access(RLIBPATH, 2) == 0)[1]]
} else {
    #if not possible, create a local library in the home directory
    RLIBLOC <- Sys.getenv("R_LIBS_USER")
    dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE,
               recursive = TRUE)
}

cat("\nINFO: Installing packages to --> ", RLIBLOC, "\n\n")

# define the R mirror to download packages
pkg_mirror <- "https://cloud.r-project.org"
print(paste("Using mirror: ", pkg_mirror))

# get the script path
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_name <- sub(file_arg_name, "",
                   initial_options[grep(file_arg_name, initial_options)])
script_dirname <- dirname(script_name)
print(script_dirname)

# read the dependencies
dependencies <- scan(
    paste(script.dirname, "r_requirements.txt", sep = "/"),
    what = "character"
)
# TODO: find a solution for script directory
inst_packages <- installed.packages()
package_list <- dependencies[!(dependencies %in% inst_packages[, "Package"])]

if (length(package_list) == 0) {
    print("All packages are installed!")
} else {
    print(paste("Number of packages to be installed: ", length(package_list)))
}

for (package_name in package_list) {
    print(paste("     Installing package --> ", package_name))
    install.packages(
        package_name,
        repos = pkg_mirror,
        dependencies = c("Depends", "Imports")
    )
}
