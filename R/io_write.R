#' Realized phenotypes in long format
#'
#' The canonical output of the v2 grammar: one row per individual x trait x rep
#' with columns `id`, `trait`, `rep`, `value`.
#'
#' @param sim a `phenotype_sim`.
#' @return a data frame in long format.
#' @export
phenotypes_long <- function(sim) {
  .check_sim(sim)
  sim$pheno
}

#' Realized phenotypes in wide format
#'
#' One row per individual x rep; one column per trait.
#'
#' @param sim a `phenotype_sim`.
#' @return a data frame in wide format.
#' @export
phenotypes_wide <- function(sim) {
  .check_sim(sim)
  long <- sim$pheno
  wide <- stats::reshape(
    long[, c("id", "rep", "trait", "value")],
    idvar = c("id", "rep"), timevar = "trait", direction = "wide"
  )
  names(wide) <- sub("^value\\.", "", names(wide))
  rownames(wide) <- NULL
  wide
}

#' Write realized phenotypes to disk
#'
#' Writes the long (default) or wide table as a delimited file. Specialized
#' exporters (gemma / plink / multi-file) are out of scope for the grammar core.
#'
#' @param sim a `phenotype_sim`.
#' @param file output path.
#' @param format "long" (default) or "wide".
#' @param sep field separator (default tab).
#' @return `file`, invisibly.
#' @export
write_phenotypes <- function(sim, file, format = c("long", "wide"),
                             sep = "\t") {
  .check_sim(sim)
  format <- match.arg(format)
  tab <- if (format == "long") phenotypes_long(sim) else phenotypes_wide(sim)
  data.table::fwrite(tab, file = file, sep = sep)
  invisible(file)
}
