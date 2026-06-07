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
