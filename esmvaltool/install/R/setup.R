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
dependencies <-scan(paste(script.dirname,"r_requirements.txt", sep="/"),what="character") # TODO: find a solution for script directory

packageList <- dependencies[!(dependencies %in% installed.packages()[,"Package"])]

if (length(packageList)==0) {
    print("All packages are installed!")
} else {
    print(paste("Number of packages to be installed: ", length(package_list)))
}

for (packageName in packageList) {
    print(paste("    Installing package --> ", packageName))
    install.packages(packageName, repos=pkgMirror, dependencies = c("Depends", "Imports"))
}

# print("List of installed packages:")
# installed.packages()
