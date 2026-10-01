.onAttach = function(libname, pkgname) {
    packageStartupMessage("**********\nThank you for using the simplePHENOTYPES\n",
    "For the reference publication, please run: citation(\"simplePHENOTYPES\")\n",
    "A Developmental version may be found at: https://github.com/samuelbfernandes/simplePHENOTYPES\n**********"
    )
}

# Variables injected into create_phenotypes() and qtn_from_user() by
# check_in() via assign(..., envir = parent.frame()).
utils::globalVariables(c(
  "add", "dom", "epi", "var",
  "add_QTN_num", "dom_QTN_num",
  "nonnumeric", "null_setting",
  "print1", "print2",
  "rep_by", "yes_no", "len_d",
  "mm", "tempdir", "path_out"
))

# data.table NSE in handle_finalreport().
utils::globalVariables(c(":=", "geno"))

#' Unwrap a `Result` returned by the Rust kernel
#'
#' The kernel entry points return `list(ok, err)` (extendr feature
#' `result_list`) instead of panicking, because a Rust panic aborts the R process
#' on toolchains whose unwinder cannot cross R's frames (a gcc-linked macOS
#' build). `err` becomes an ordinary R error here.
#' @param res the value returned by `.Call()`.
#' @param fn the kernel function name, for the message.
#' @keywords internal
#' @noRd
.unwrap_extendr <- function(res, fn) {
  if (inherits(res, "extendr_result")) {
    if (!is.null(res$err)) {
      stop(fn, "(): ", res$err, call. = FALSE)
    }
    return(res$ok)
  }
  res
}
