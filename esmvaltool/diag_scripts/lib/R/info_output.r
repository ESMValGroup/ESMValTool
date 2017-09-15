info_output <- function(output_string,
                        verbosity,
                        required_verbosity) {

    main_wd = Sys.getenv(c("ESMValTool_interface_data"))
    if (nchar(main_wd) == 0) {
        print(paste("info: ", output_string, sep = ""))
    }
    else {
        indent_step <- 3
        if (verbosity == 2) {
            padding_space <- 13
        }
        else {
            padding_space <- 1
        }

        if (verbosity >= required_verbosity) {
            entering_routine <- regexpr("<<<<<<<< Entering", output_string, fixed = T)[1]
            indent_file <- file.path(main_wd, "curr_trace_indent.txt")
            if (entering_routine != -1) {
                indent <- read_integer(indent_file)
                indent = indent + indent_step
                write_integer(indent, indent_file)
            }

            ## Add a prefix of ">>", "<<" or "  " to output_string
            indent_str = as.character(read_integer(indent_file))
            format <- paste("%0", indent_str, "d", sep = "")
            entering <- regexpr("<", output_string, fixed = T)[1]
            leaving <- regexpr(">", output_string, fixed = T)[1]
            if (entering != -1) {
                indent_str <- gsub("0", "<", sprintf(format, 0))
            }
            else if (leaving != -1) {
                indent_str <- gsub("0", ">", sprintf(format, 0))
            }
            else {
                indent_str <- sprintf(format, 0)
                indent_str <- sub("0", "", indent_str)
                indent_str <- gsub("0", " ", indent_str)
            }

            pasted_string <- paste("info: ", indent_str, output_string, sep = "")
            print(pasted_string)

            ## Decrease indentation if we're leaving an NCL routine
            leaving_routine <- regexpr(">>>>>>>> Leaving", output_string, fixed = T)[1]
            if (leaving_routine != -1) {
                indent <- read_integer(indent_file)
                indent = indent - indent_step
                write_integer(indent, indent_file)
            }
        }
    }
}

read_integer <- function (filename) 
{
    fileConn <- file(filename)
    indent <- as.integer(readLines(filename))
    close(fileConn)
    indent
}

write_integer <- function (indent, filename) 
{
    fileConn <- file(filename)
    writeLines(as.character(indent), fileConn)
    close(fileConn)
}
