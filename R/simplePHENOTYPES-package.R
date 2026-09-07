#' @useDynLib simplePHENOTYPES, .registration = TRUE
#' @keywords internal
"_PACKAGE"

# Package-level roxygen lives here, not in the generated R/extendr-wrappers.R.
# rextendr emits a `@docType package` block, which roxygen2 >= 7.3 deprecates and
# silently drops along with the @useDynLib directive it carries -- leaving a
# NAMESPACE with no useDynLib, so the compiled code never loads. Keeping the
# directive in a hand-maintained file makes it survive `rextendr::document()`.
