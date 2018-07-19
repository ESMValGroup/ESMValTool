library(lintr)
args = commandArgs(trailingOnly=TRUE)
print(args)
check_paths <- list('esmvaltool', 'tests')

root_folder <- args(1)

has_errors <- FALSE

for(path in check_paths){
    for(file in list.files(file.path(root_folder, path))){
        print(file)
    }

    for (folder in list.dirs(file.path(root_folder, path))){
        print(folder)
            # if(os.path.splitext(filename)[1].lower() == '.r':
            #     errors <- lintr.lint(os.path.join(dirpath, filename),
            #                         linters)
            #     for error in errors:
            #         print(str(error)[0:-1])
            #         has_errors = True
    }
}
        
    