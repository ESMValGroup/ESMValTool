# Wrappers to call external commands using system2
# Currently implemented: cdo, nco

cdo <-
  function(command,
             args = "",
             input = "",
             options = "",
             output = "",
             stdout = "",
             noout = F) {
    if (args != "") {
      args <- paste0(",", args)
    }
    if (stdout != "") {
      stdout <- paste0(" > '", stdout, "'")
      noout <- T
    }
    if (input[1] != "") {
      for (i in seq_along(input)) {
        input[i] <- paste0("'", input[i], "'")
      }
      input <- paste(input, collapse = " ")
    }
    output0 <- output
    if (output != "") {
      output <- paste0("'", output, "'")
    } else if (!noout) {
      output <- tempfile()
      output0 <- output
    }
    argstr <- paste0(
      options, " ", command, args, " ", input, " ", output,
      " ", stdout
    )
    print(paste("cdo", argstr))
    ret <- system2("cdo", args = argstr)
    if (ret != 0) {
      stop(paste("Failed (", ret, "): cdo", argstr))
    }
    return(output0)
  }

nco <- function(cmd, argstr) {
  ret <- system2(cmd, args = argstr)
  if (ret != 0) {
    stop(paste("Failed (", ret, "): ", cmd, " ", argstr))
  }
}
